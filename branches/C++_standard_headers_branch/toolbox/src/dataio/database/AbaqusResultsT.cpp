/* $Id: AbaqusResultsT.cpp,v 1.24 2003-11-21 22:41:46 paklein Exp $ */
/* created: S. Wimmer 9 Nov 2000 */

#include "AbaqusResultsT.h"
#include "fstreamT.h"
#include <time.h>

using namespace Tahoe;

AbaqusResultsT::AbaqusResultsT (ostream& message) :
  fMessage (message),
  fMarker("*"),
  fOutVersion("5.8-1"),
  fBinary(false),
  fBufferDone(0),
  fBufferSize(0),
  fCurrentLength(-1),
  fNumNodes(0),
  fNumElements (0),
  fNumNodeSets(0),
  fNumElementSets(0),
  fStartCount (0),
  fEndCount (0),
  fModalCount (0),
  fNumNodeVars (0),
  fNumElemVars (0),
  fNumQuadVars (0)
{
  SetVariableNames ();
}

bool AbaqusResultsT::Initialize (const char *filename)
{
  //TEMP - workaround for problem with CW7
  //fstreamT::FixPath(filename, fFileName); //not needed if fIn is ifstreamT
  fFileName = filename;
  fIn.open (fFileName);
  if (!fIn.is_open())
    {
      fMessage << "\n AbaqusResultsT::Initialize unable to open file " << filename << endl;
      return false;
    }
  
  int key = -1;
  ReadNextRecord (key);
  if (key != VERSION) 
    {
      fMessage << "\n AbaqusResultsT::Initialize not ASCII, trying Binary" << endl;
      fBinary = true;
      ResetFile ();

      if (ReadNextRecord (key) != OKAY)
	{
	  fMessage << "\n AbaqusResultsT::Initialize not Binary, unable to read next record" << endl;
	  cout << "key " << key << endl;
	  return false;
	}

      if (key != VERSION)
	{
	  fMessage << "\n AbaqusResultsT::Initialize not Binary" << endl;
	  return false;
	}
    }
  
  if (!ReadVersion ()) 
    {
      fMessage << "\n AbaqusResultsT::Initialize unreadable file" << endl;
      return false;
    }
    
    /* must be OK */
    return true;
}

void AbaqusResultsT::Create (const char* filename, bool binary, int numelems, int numnodes, double elemsize)
{
  fOut.open (filename);
  if (!fOut)
    {
      fMessage << "\n AbaqusResultsT::Create unable to open file " << filename << endl;
      throw ExceptionT::kDatabaseFail;
    }
  fFileName = filename;
  fBinary = binary;
  if (fBinary)
    fBufferSize = 512*kDoubleSize;
  else
    fBufferSize = 80;
  fBufferDone = 0;

  time_t now;
  time (&now);
  char date[40], time[20];
  strftime (time, 40, "%X", localtime (&now));
  strftime (date, 40, "%d-%b-%Y", localtime (&now));

  // write version record
  int length = 9;
  StringT temp;
  WriteASCII (fMarker);
  Write (length);
  Write (VERSION);
  Write (fOutVersion);
  temp = date;
  Write (temp, 2);
  temp = time;
  Write (temp);
  Write (numelems);
  Write (numnodes);
  Write (elemsize);
  fOut.flush();
}
 
void AbaqusResultsT::OpenWrite (const char *filename, bool binary, int bufferwritten)
{
  fOut.open (filename, ios::app);
  if (!fOut)
    {
      fMessage << "\n AbaqusResultsT::OpenWRite unable to open file " << filename << endl;
      throw ExceptionT::kDatabaseFail;
    }
  fFileName = filename;
  fBinary = binary;
  if (fBinary)
    fBufferSize = 512*kDoubleSize;
  else
    fBufferSize = 80;
  fBufferDone = bufferwritten;
}

int AbaqusResultsT::Close (void)
{
  int bufferdone = fBufferDone;
  fIn.close ();
  fFileName.Free ();
  fBufferDone = 0;
  fBufferSize = 0;
  fBuffer.Free ();
  fCurrentLength = 0;
  fNumNodes = 0;
  fNumElements = 0;
  fNumNodeSets = 0;
  fNumElementSets = 0;
  fStartCount = 0;
  fEndCount = 0;
  fModalCount = 0;
  for (int i=0; i < NVT; i++)
    fVariableTable[i].SetIOData (0, AbaqusVariablesT::kNotUsed, -1);
  fNumNodeVars = 0;
  fNumElemVars = 0;
  fNumQuadVars = 0;
  fElementNumber.Free ();
  fNodeNumber.Free ();
  fTimeIncs.Free ();
  fTimeSteps.Free ();
  fModeIncs.Free ();
  fModeSteps.Free ();
  fNodeSetNames.Free ();
  fElementSetNames.Free ();
  return bufferdone;
}

bool AbaqusResultsT::ScanFile (int &numelems, int &numnodes, int &numtimesteps, int &nummodes)
{
  int key = 0, oldkey = 0, error = OKAY, outputmode = -1, location = -1;
  fNumElements = 0;
  fNumNodes = 0;
  fModalCount = 0;
  fStartCount = 0;
  fEndCount = 0;
  error = ReadNextRecord (key);

  // kick out of loop if END or BAD error flag returned
  while (error == OKAY)
    {
      switch (key)
	{
	case ELEMENT: ScanElement (); break;
	case NODE: 
	  {
	    int number;
	    if (!Read (number)) return false;
	    fNodeNumber.Append (number);
	    fNumNodes++; 
	    break;
	  }
	case NODESET: 
	  {
	    StringT name;
	    if (!ReadSetName (name, 1)) return false;
	    fNodeSetNames.Append (name);
	    iArrayT tempset (0);
	    fNodeSets.Append (tempset);
	    StoreSet (fNodeSets[fNumNodeSets]);
	    fNumNodeSets++; 
	    break;
	  }
	case NODESETCONT:
	  {
	    StoreSet (fNodeSets[fNumNodeSets-1]);
	    break;
	  }
	case ELEMENTSET: 
	  {
	    StringT name;
	    if (!ReadSetName (name, 1)) return false;
	    fElementSetNames.Append (name);
	    iArrayT tempset (0);
	    fElementSets.Append (tempset);
	    StoreSet (fElementSets[fNumElementSets]);
	    fNumElementSets++; 
	    break;
	  }
	case ELEMSETCONT:
	  {
	    StoreSet (fElementSets[fNumElementSets-1]);
	    break;
	  }
	case MODAL: 
	  {
	    int number;
	    double mode;
	    if (!Read (number) || !Read (mode) ) return false;

	    fModeIncs.Append (number);
	    fModeSteps.Append (mode);
	    fModalCount++; 
	    break;
	  }
	case ENDINCREMENT: fEndCount++; break;
	case STARTINCREMENT: 
	  {
	    int procedure, step, number;
	    double steptime, creep, amp, time;
	    if (!Read (time) || !Read (steptime) || !Read (creep) || 
		!Read (amp) || !Read (procedure) || !Read (step) ||
		!Read (number) )
	      return false;

	    fTimeIncs.Append (number);
	    fTimeSteps.Append (time);
	    fStartCount++; 
	    break;
	  }
	case ELEMENTHEADER:
	  {
	    int objnum, intpt, sectpt;
	    if (fStartCount == 1)
	      ReadElementHeader (objnum, intpt, sectpt, location);
	    break;
	  }
	case OUTPUTDEFINE: 
	  {
	    StringT outsetname;
	    if (fStartCount == 1)
	      ReadOutputDefinitions (outputmode, outsetname);
	    break;
	  }
	}
      
      /* scan variable records, but only for the first time step */
      if (fStartCount == 1)
	{
	  /* if the record key is a variable, scan for variable within */
	  if (VariableKeyIndex (key) > 0)
	    ScanVariable (key, outputmode, location);
	}

      oldkey = key;
      error = ReadNextRecord (key);
      if (error == BAD)
	{
	  cout << "\nAbaqusResultsT::ScanFile unable to read next record \n";
	  cout << " or it was unsuccessful in identifying EOF.\n";
	  cout << "Old key " << oldkey << " Key " << key << endl;
	  return false;
	}
    }
  numelems = fNumElements;
  numnodes = fNumNodes;
  numtimesteps = fStartCount;
  nummodes = fModalCount;

  // error should equal END at this point, but add this statement for surety
  // if error == BAD, then a problem was encountered and false should be returned
  if (error == BAD)
    {
      cout << "\nAbaqusResultsT::ScanFile was unsuccessful in scanning to the end of the file.\n";
      return false;
    }

  // if error == END, then you scanned the file until EOF was encountered, 
  // and you have successfully scanned the file.
  return true;}

void AbaqusResultsT::ElementSetNames (ArrayT<StringT>& names) const
{
  for (int i=0; i < fNumElementSets; i++)
    names[i] = fElementSetNames[i];
}

void AbaqusResultsT::NodeSetNames (ArrayT<StringT>& names) const
{
  for (int i=0; i < fNumNodeSets; i++)
    names[i] = fNodeSetNames[i];
}

void AbaqusResultsT::NodeMap (iArrayT& n) const
{
  n.CopyPart (0, fNodeNumber, 0, fNumNodes);
}

int AbaqusResultsT::NumNodesInSet (const StringT& name) const
{
  int index = 0;
  bool found = false;
  for (int i = 0; i < fNumNodeSets && !found; i++)
    if (strncmp (fNodeSetNames[i].Pointer(), name.Pointer(), name.StringLength()) == 0)
      {
	found = true;
	index = i;
      }

  if (found)
    return fNodeSets[index].Length();
  else
    return 0;
}

void AbaqusResultsT::NodeSet (const StringT& name, iArrayT& nset) const
{
  int index = 0;
  bool found = false;
  for (int i = 0; i < fNumNodeSets && !found; i++)
    if (strncmp (fNodeSetNames[i].Pointer(), name.Pointer(), name.StringLength()) == 0)
      {
	found = true;
	index = i;
      }

  if (found)
    nset = fNodeSets[index];
  else
    nset.Free();
}

void AbaqusResultsT::ElementMap (iArrayT& e) const
{
  e.CopyPart (0, fElementNumber, 0, fNumElements);
}

int AbaqusResultsT::NumElements (const StringT& name) const
{
  int index = 0;
  bool found = false;
  for (int i = 0; i < fNumElementSets && !found; i++)
    if (strncmp (fElementSetNames[i].Pointer(), name.Pointer(), name.StringLength()) == 0)
      {
	found = true;
	index = i;
      }

  if (found)
    return fElementSets[index].Length();
  else
    return 0;
}

int AbaqusResultsT::NumElementNodes (const StringT& name)
{
  ResetFile ();
  StringT setname = "";

  while (strncmp (setname.Pointer(), name.Pointer(), name.StringLength()) != 0)
    {
      if (!AdvanceTo (ELEMENTSET))
	{
	  fMessage << "\nAbaqusResultsT::NumElementNodes, unable to advance to ELEMENTSET\n\n";
	  throw ExceptionT::kDatabaseFail;
	}
      if (!ReadSetName (setname, 1)) throw ExceptionT::kDatabaseFail;
    }

  int el;
  if (!Read (el)) throw ExceptionT::kDatabaseFail;

  ResetFile ();
  int cel=-1;
  while (cel != el)
    {
      if (!AdvanceTo (ELEMENT))
	{
	  fMessage << "\nAbaqusResultsT::NumElementNodes, unable to advance to ELEMENT\n\n";
	  throw ExceptionT::kDatabaseFail;
	}
      if (!Read (cel)) throw ExceptionT::kDatabaseFail;
    }

  return fCurrentLength - 1;
}

int AbaqusResultsT::NumElementQuadPoints (const StringT& name) const
{
  int index = 0;
  bool found = false;
  for (int i = 0; i < fNumElementSets && !found; i++)
    if (strncmp (fElementSetNames[i].Pointer(), name.Pointer(), name.StringLength()) == 0)
      {
	found = true;
	index = i;
      }

  if (found)
    {
      const iArrayT& set = fElementSets[index];
      for (int i=0; i < fElementNumber.Length(); i++)
	if (fElementNumber[i] == set[0])
	  return fNumElementQuadPoints [i];
    }
  return 0;
}

void AbaqusResultsT::ElementSet (const StringT& name, iArrayT& elset) const
{
  int index = 0;
  bool found = false;
  for (int i = 0; i < fNumElementSets && !found; i++)
    if (strncmp (fElementSetNames[i].Pointer(), name.Pointer(), name.StringLength()) == 0)
      {
	found = true;
	index = i;
      }

  if (found)
    elset = fElementSets[index];
  else
    elset.Free();
}

void AbaqusResultsT::GeometryCode (const StringT& name, GeometryT::CodeT& code)
{
  ResetFile ();
  StringT setname = "";

  while (strncmp (setname.Pointer(), name.Pointer(), name.StringLength()) != 0)
    {
      if (!AdvanceTo (ELEMENTSET))
	{
	  fMessage << "\nAbaqusResultsT::GeometryCode, unable to advance to ELEMENTSET\n\n";
	  throw ExceptionT::kDatabaseFail;
	}
      if (!ReadSetName (setname, 1)) throw ExceptionT::kDatabaseFail;
    }

  int el;
  if (!Read (el)) throw ExceptionT::kDatabaseFail;

  ResetFile ();
  int cel=-1;
  while (cel != el)
    {
      if (!AdvanceTo (ELEMENT))
	{
	  fMessage << "\nAbaqusResultsT::NumElementNodes, unable to advance to ELEMENT\n\n";
	  throw ExceptionT::kDatabaseFail;
	}
      if (!Read (cel)) throw ExceptionT::kDatabaseFail;
    }

  int numintpts;
  StringT elname;
  if (!Read (elname, 1)) throw ExceptionT::kDatabaseFail;
  if (TranslateElementName (elname.Pointer(), code, numintpts) == BAD)
    {
      fMessage << "\n AbaqusResultsT::GeometryCode Unable to translate element name" << endl;
      fMessage << "Name " << elname.Pointer() << endl;
      throw ExceptionT::kDatabaseFail;
    }
}

void AbaqusResultsT::ModeData (int index, int &number, double &mode) const
{
  if (index < 0 || index > fModalCount)
    throw ExceptionT::kOutOfRange;

  number = fModeIncs [index];
  mode = fModeSteps [index];
}

void AbaqusResultsT::TimeData (int index, int &number, double &time) const
{
  if (index < 0 || index > fStartCount)
    throw ExceptionT::kOutOfRange;

  number = fTimeIncs [index];
  time = fTimeSteps [index];
}

void AbaqusResultsT::NodeVariables (iArrayT& keys, iArrayT& dims) const
{
  int c = 0;
  for (int i=0; i < NVT; i++)
    if (fVariableTable[i].Type() == AbaqusVariablesT::kNode)
      {
	int dimension = fVariableTable[i].Dimension();
	for (int j=0; j < dimension; j++)
	  {
	    keys[c] = fVariableTable[i].Key();
	    dims[c] = dimension;
	    c++;
	  }
      }
}

void AbaqusResultsT::ElementVariables (iArrayT& keys, iArrayT& dims) const
{
  int c = 0;
  for (int i=0; i < NVT; i++)
    if (fVariableTable[i].Type() == AbaqusVariablesT::kElement)
      {
	int dimension = fVariableTable[i].Dimension();
	for (int j=0; j < dimension; j++)
	  {
	    keys[c] = fVariableTable[i].Key();
	    dims[c] = dimension;
	    c++;
	  }
      }
}

void AbaqusResultsT::QuadratureVariables (iArrayT& keys, iArrayT& dims) const
{
  int c = 0;
  for (int i=0; i < NVT; i++)
    if (fVariableTable[i].Type() == AbaqusVariablesT::kQuadrature)
      {
	int dimension = fVariableTable[i].Dimension();
	for (int j=0; j < dimension; j++)
	  {
	    keys[c] = fVariableTable[i].Key();
	    dims[c] = dimension;
	    c++;
	  }
      }
}

void AbaqusResultsT::VariablesUsed (const StringT& name, AbaqusVariablesT::TypeT vt, iArrayT& used)
{
  used = 0;
  iArrayT done (fVariableTable.Length());
  done = 0;

  /* examine the first time step/inc, variables should be same for all steps/incs */
  ResetFile ();
  int inc;
  double time;
  if (fModeIncs.Length() > 0)
    if (!NextMode (inc, time)) throw ExceptionT::kDatabaseFail;
  else
    if (!NextTimeSteps (inc, time)) throw ExceptionT::kDatabaseFail;
  
  int key;
  int ID, objnum, intpt, secpt, location, outputmode;
  StringT outsetname;
  while (ReadNextRecord (key) == OKAY)
    {
      switch (key)
	{
	case ENDINCREMENT:
	case MODAL:
	  return;
	case ELEMENTHEADER:
	  ReadElementHeader (objnum, intpt, secpt, location);
	  break;
	case OUTPUTDEFINE:
	  ReadOutputDefinitions (outputmode, outsetname);
	  break;
	default:
	  {
	    /* is this the setname we are interested in */
	    if (strncmp (outsetname.Pointer(), name.Pointer(), outsetname.StringLength()) == 0)
	      {
		/* make sure it is a variable */
		int index = VariableKeyIndex (key);
		if (index > 0 && index < done.Length() && done[index] < 1)
		  {
		    /* make sure the variable is of the type we are interested in */
		    if (CorrectType (outputmode, objnum, intpt, location, vt, ID))
		      {
			int offset = fVariableTable[index].IOIndex();
			int dim = fVariableTable[index].Dimension();
			for (int bj=0; bj < dim; bj++)
			  used [offset + bj] = fVariableTable[index].Dimension();
		      }
		  }
	      }
	  }
	}
    }
}

void AbaqusResultsT::ReadVariables (AbaqusVariablesT::TypeT vt, int step, dArray2DT& values, const StringT& name)
{
  /* read for all elements or nodes or just a set */
  iArrayT set;
//  bool subset = DataPoints (vt, name, set);

  /* number of quadrature points */
  int numquadpts = 0;
  if (vt == AbaqusVariablesT::kQuadrature)
    {
      numquadpts = NumElementQuadPoints (name);
      if (numquadpts < 1) throw ExceptionT::kDatabaseFail;
    }

  /* advance istream to time step */
  AdvanceToTimeIncrement (step);

  int key;
  int ID, objnum, intpt, secpt, location, outputmode, itemp;
  StringT outsetname, ctemp;
  while (ReadNextRecord (key) == OKAY)
    {
      switch (key)
	{
	case ENDINCREMENT:
	case MODAL:
	  return;
	case ELEMENTHEADER:
	  ReadElementHeader (objnum, intpt, secpt, location);
	  break;
	case OUTPUTDEFINE:
	  ReadOutputDefinitions (outputmode, outsetname);
	  break;
	default:
	  {
	    /* make sure it is a variable */
	    int index = VariableKeyIndex (key);
	    if (index > 0 && index < NVT)
	      {
		/* is the record found, one that you want to read */
		if (CorrectType (outputmode, objnum, intpt, location, vt, ID))
		  {
		    switch (fVariableTable[index].FirstAttribute())
		      {
		      case AbaqusVariablesT::kNodeNumber:
			if (!Read (ID)) throw ExceptionT::kDatabaseFail;
			break;
		      case AbaqusVariablesT::kIntType:
			if (!Read (itemp)) 
			  {
			    cout << fBuffer << endl;
			    cout << fBufferDone << endl;
			    throw ExceptionT::kDatabaseFail;
			  }
			break;
		      case AbaqusVariablesT::kCharType:
			if (!Read (ctemp, 1))
			  {
			    cout << fBuffer << endl;
			    cout << fBufferDone << endl;
			    throw ExceptionT::kDatabaseFail;
			  }
			break;
		      }
		    
		    bool save = true;
		    int row;
		    set.HasValue (ID, row);
		    if (row < 0 || row > set.Length())
		      save = false;
		    
		    // modify row by 
		    if (vt == AbaqusVariablesT::kQuadrature)
		      row = row*numquadpts + intpt - 1;
		    
		    if (save)
		      {
			int num = fCurrentLength;
			dArrayT v (num);
			for (int i=0; i < num; i++)
			  if (!Read (v[i])) throw ExceptionT::kDatabaseFail;
			
			int offset = fVariableTable[index].IOIndex();
			values.CopyPart (row*values.MinorDim() + offset, v, 0, num); 
		      }
		  }
	      }
	  }
	}
    }
}

const char *AbaqusResultsT::VariableName (int index) const
{
  if (index < 0 || index >= NVT)
    return "Unknown";
  else
    return fVariableTable[index].Name();
}

int AbaqusResultsT::VariableKey (const char* name) const
{
  for (int i=0; i < NVT; i++)
    {
      const StringT& n = fVariableTable[i].Name();
      if (strncmp (name, n.Pointer(), n.StringLength()) == 0)
	return fVariableTable[i].Key();
    }
  return -1;
}

int AbaqusResultsT::VariableKey (int index) const
{
  if (index < 0 || index >= NVT)
    return -1;
  else
    return fVariableTable[index].Key();
}

int AbaqusResultsT::VariableKeyIndex (int key) const
{
  for (int i=0; i < NVT; i++)
    if (key == fVariableTable[i].Key())
      return i;
  return -1;
}

bool AbaqusResultsT::NextCoordinate (int &number, dArrayT &nodes)
{
  if (!AdvanceTo (NODE))
    {
      fMessage << "\nAbaqusResultsT::NextCoordinate, unable to advance to NODE\n\n";
      throw ExceptionT::kDatabaseFail;
    }

  if (!Read (number))
    throw ExceptionT::kDatabaseFail;

  int dof = fCurrentLength;
  nodes.Allocate (dof);
  for (int ii=0; ii < dof; ii++)
    if (!Read (nodes[ii]))
      throw ExceptionT::kDatabaseFail;
  return true;
}

bool AbaqusResultsT::NextElement (int &number, GeometryT::CodeT &type, iArrayT &nodes)
{
  if (!AdvanceTo (ELEMENT))
    {
      fMessage << "\nAbaqusResultsT::NextElement, unable to advance to ELEMENT\n\n";
      throw ExceptionT::kDatabaseFail;
    }

  StringT name;
  if (!Read (number) || !Read (name, 1) )
    throw ExceptionT::kDatabaseFail;

  int numintpts;
  if (TranslateElementName (name.Pointer(), type, numintpts) == BAD) 
    {
      fMessage << "\n AbaqusResultsT::NextElement Unable to translate element name" << endl;
      fMessage << "Name " << name.Pointer() << endl;
      throw ExceptionT::kDatabaseFail;
    }

  int numnodes = fCurrentLength;
  nodes.Allocate (numnodes);
  nodes = 300;
  for (int i=0; i < numnodes; i++)
    {
      int temp;
      if (!Read (temp))
      throw ExceptionT::kDatabaseFail;
      nodes[i] = temp;
    }
  return true;
}

void AbaqusResultsT::WriteConnectivity (GeometryT::CodeT code, int startnumber, const iArray2DT& connects)
{
  StringT name;
  int numelemnodes;
  GetElementName (code, connects.MinorDim(), numelemnodes, name);

  int length = 4 + numelemnodes;
  for (int i=0; i < connects.MajorDim(); i++)
    {
      WriteASCII (fMarker);
      Write (length);
      Write (ELEMENT);
      Write (startnumber + i);
      Write (name);
      for (int j=0; j < numelemnodes; j++)
	Write (connects (i,j));
    }
}

void AbaqusResultsT::WriteCoordinates (const iArrayT& nodes_used, const dArray2DT& coords)
{
  int length = 3 + 6;
  const int *pn = nodes_used.Pointer();
  double zero = 0.;
  for (int i=0; i < nodes_used.Length(); i++)
    {
      WriteASCII (fMarker);
      Write (length);
      Write (NODE);
      Write (*pn);
		const double *pc = coords (*pn++ - 1);
      for (int j=0; j < coords.MinorDim(); j++)
	Write (*pc++);
      // fill to 6 degrees of freedom
      for (int k=coords.MinorDim(); k < 6; k++)
	Write (zero);
    }
}

void AbaqusResultsT::WriteElementSet (const StringT& name, const iArrayT& elms)
{
  AbaqusResultsT::GeneralKeys key = ELEMENTSET;
  int headerlength = 3;
  int num_vals_in_record = 80 - headerlength;
  const int *pe = elms.Pointer();
  for (int i=0; i < elms.Length(); i++)
    {
      if (i%num_vals_in_record == 0)
	{
	  int length = num_vals_in_record;
	  if (elms.Length() - i < num_vals_in_record)
	    length = elms.Length() - i;
	  if (i > 0)
	    {
	      key = ELEMSETCONT;
	      headerlength = 2;
	      num_vals_in_record = 80 - headerlength;
	    }
	  WriteASCII (fMarker);
	  Write (length + headerlength);
	  Write (key);
	  if (key == ELEMENTSET)
	    Write (name);
	}
      Write (*pe++);
    }
}

void AbaqusResultsT::WriteNodeSet (const StringT& name, const iArrayT& nodes)
{
  AbaqusResultsT::GeneralKeys key = NODESET;
  int headerlength = 3;
  int num_vals_in_record = 80 - headerlength;
  const int *pe = nodes.Pointer();
  for (int i=0; i < nodes.Length(); i++)
    {
      if (i%num_vals_in_record == 0)
	{
	  int length = num_vals_in_record;
	  if (nodes.Length() - i < num_vals_in_record)
	    length = nodes.Length() - i;
	  if (i > 0)
	    {
	      key = NODESETCONT;
	      headerlength = 2;
	      num_vals_in_record = 80 - headerlength;
	    }
	  WriteASCII (fMarker);
	  Write (length + headerlength);
	  Write (key);
	  if (key == NODESET)
	    Write (name);
	}
      Write (*pe++);
    }
}

void AbaqusResultsT::WriteActiveDOF (const iArrayT& active)
{
  int length = active.Length() + 2;
  WriteASCII (fMarker);
  Write (length);
  Write (ACTIVEDOF);
  for (int i=0; i < active.Length(); i++)
    Write (active[i]);
}

void AbaqusResultsT::WriteHeading (const StringT& heading)
{
  if (heading.Length() <= 1) return;
  int length = 10 + 2;
  WriteASCII (fMarker);
  Write (length);
  Write (HEADING);
  Write (heading, 10);
}

void AbaqusResultsT::WriteStartIncrement (int step, int inc, double totaltime, 
     double time, double timeincrement, AbaqusResultsT::AnalysisTypeT atype)
{
  double creep = 0, amplitude = 0, factor = 0, freq = 0;
  int perturb = 1, length = 23;
  WriteASCII (fMarker);
  Write (length);
  Write (STARTINCREMENT);
  Write (totaltime);
  Write (time);
  Write (creep);
  Write (amplitude);
  Write (atype);
  Write (step);
  Write (inc);
  Write (perturb);
  Write (factor);
  Write (freq);
  Write (timeincrement);
  StringT attribute ("default load case");
  Write (attribute, 10);
}

void AbaqusResultsT::WriteOutputDefinition (int key, const StringT& setname, GeometryT::CodeT code, int numelemnodes)
{
#pragma unused(setname)

  WriteASCII (fMarker);
  Write (5);
  Write (OUTPUTDEFINE);

  int index = VariableKeyIndex (key);
  if (fVariableTable[index].Point() != AbaqusVariablesT::kNodePoint)
    {
      StringT ename;
      int num_output_nodes;
      GetElementName (code, numelemnodes, num_output_nodes, ename);
      Write (ename);
    }
  else
    Write (0);
}

void AbaqusResultsT::WriteNodeVariables (int& i, const iArrayT& key, const dArray2DT& values, const iArrayT& nodes_used, int numdir, int numshear)
{
  // determine record length
  int count = 0;
  for (int j=i; j < key.Length(); j++)
    if (key[j] == key[i])
      count ++;
  int length = 2 + count;
  
  // account for node number
  if (VariableWrittenWithNodeNumber (key[i])) length++;

  // write data for this variable
  int index = VariableKeyIndex (key[i]);
  for (int n=0; n < nodes_used.Length(); n++)
    {
      // for element integration point data (assume node averaged)
      if (fVariableTable[index].Point() != AbaqusVariablesT::kNodePoint)
	WriteElementHeader (key[i], nodes_used[n], 0, 0, kElementNodeAveraged, numdir, numshear, 0, 0);
      
      // write record
      WriteASCII (fMarker);
      Write (length);
      if (VariableWrittenWithNodeNumber (key[i])) Write (nodes_used[n]);
      for (int m=i; m < i + count; m++)
	Write (values (n, m));
    }

  i += count - 1;
}

void AbaqusResultsT::WriteElementVariables (int& i, const iArrayT& key, const dArray2DT& values, const iArrayT& els_used, int numdir, int numshear)
{
  // determine record length
  int count = 0;
  for (int j=i; j < key.Length(); j++)
    if (key[j] == key[i])
      count ++;
  int length = 2 + count;
  
  // write data for this variable
  int index = VariableKeyIndex (key[i]);
  for (int n=0; n < els_used.Length(); n++)
    {
      // all element variables must have element header
      // no element data originates at node points
      if (fVariableTable[index].Point() != AbaqusVariablesT::kNodePoint)
	WriteElementHeader (key[i], els_used[n], 0, 0, kElementWhole, numdir, numshear, 0, 0);
      else
	throw ExceptionT::kDatabaseFail;
      
      // write record
      WriteASCII (fMarker);
      Write (length);
      Write (els_used[n]);
      for (int m=i; m < i + count; m++)
	Write (values (n, m));
    }

  i += count - 1;
}

void AbaqusResultsT::WriteEndIncrement (void)
{
  int length = 2;
  WriteASCII (fMarker);
  Write (length);
  Write (ENDINCREMENT);
  if (!fBinary)
    {
      iArrayT size (2);
      size [0] = fBufferSize - fBufferDone + 1;
      size [1] = fBufferSize + 1;
      for (int j=0; j < 2; j++)
	{
	  StringT space (size[j]);
	  for (int i=0; i < space.Length() - 1; i++)
	    space[i] = ' ';
	  space [space.Length() - 1] = '\0';
	  WriteASCII (space);
	}
    }
  fOut.flush();
}

void AbaqusResultsT::VersionNotes (ArrayT<StringT>& records)
{
  ResetFile ();
  if (!AdvanceTo (VERSION))
    {
      fMessage << "\nAbaqusResultsT::Version, unable to advance to VERSION\n\n";
      throw ExceptionT::kDatabaseFail;
    }
  records.Allocate (4);
  int numelems, numnodes;
  double elemleng;
  if (!Read (records[1], 1) || !Read (records[2], 2) || 
      !Read (records[3], 1) || !Read (numelems) || !Read (numnodes) || 
      !Read (elemleng) )
    throw ExceptionT::kDatabaseFail;
  records[0] = "ABAQUS";
}

void AbaqusResultsT::ResetFile (void)
{
  fIn.close ();
  fIn.open (fFileName);

  if (fBinary)
    fBufferSize = 512*kDoubleSize;
  else
    fBufferSize = 0;
  fBufferDone = 0;
  fCurrentLength = 0;
}

/******************* PRIVATE **********************/

bool AbaqusResultsT::ReadVersion (void)
{
  StringT version, date, time;
  int numelems, numnodes;
  double elemleng;
 
  if (!Read (version, 1) || !Read (date, 2) || !Read (time, 1) ||
      !Read (numelems) || !Read (numnodes) || !Read (elemleng) )
    {
      fMessage << "\n AbaqusResultsT::ScanFile Unable to read version record" << endl;
      return false;
    }
  //cout << version << "\n" << date << "\n" << time << "\n" << numelems
  //   << "\n" << numnodes << "\n" << elemleng << endl;

  if (strncmp (version.Pointer(), "5.8", 3) != 0 &&
      strncmp (version.Pointer(), "6.1", 3) != 0 &&
      strncmp (version.Pointer(), "6.2", 3) != 0 )
    {
      fMessage << "\n AbaqusResultsT::ScanFile Unrecognized version " << version << endl;
      return false;
    }
  return true;
}

bool AbaqusResultsT::NextMode (int &number, double &mode)
{
  if (!AdvanceTo (MODAL))
    return false;
  if (!Read (number) || !Read (mode) )
    return false;
  return true;
}

bool AbaqusResultsT::NextTimeSteps (int &number, double &time)
{
  int procedure, step;
  double steptime, creep, amp; 
  if (!AdvanceTo (STARTINCREMENT))
    return false;
  if (!Read (time) || !Read (steptime) || !Read (creep) || !Read (amp) ||
      !Read (procedure) || !Read (step) || !Read (number) )
    return false;
  return true;
}

void AbaqusResultsT::ScanElement (void)
{
  StringT name;
  GeometryT::CodeT type = GeometryT::kNone;
  int number, numintpts = 0;

  if (!Read (number) || !Read (name, 1)) 
    throw ExceptionT::kDatabaseFail;

  if (TranslateElementName (name.Pointer(), type, numintpts) == BAD) 
    {
      fMessage << "\n AbaqusResultsT::ScanElement Encountered unknown element type" << endl;
      fMessage << "Name " << name.Pointer() << endl;
      throw ExceptionT::kDatabaseFail;
    }
  fNumElements++;
  fElementNumber.Append (number);
  fNumElementQuadPoints.Append (numintpts);
}

void AbaqusResultsT::StoreSet (iArrayT& set)
{
  iAutoArrayT newset;
  newset.Append (set);

  int num = fCurrentLength;
  int temp;
  for (int j=0; j < num; j++)
    {
      if (!Read (temp)) throw ExceptionT::kDatabaseFail;
      newset.Append (temp);
    }

  set.Allocate (newset.Length());
  set.CopyPart (0, newset, 0, newset.Length());
}

void AbaqusResultsT::ReadOutputDefinitions (int &outputmode, StringT& setname)
{
  if (!Read (outputmode)|| !ReadSetName (setname, 1) )
    throw ExceptionT::kDatabaseFail;
}

void AbaqusResultsT::ReadElementHeader (int &objnum, int& intpt, int& secpt, int &location)
{
  if (!Read (objnum) || !Read (intpt) || !Read (secpt) || !Read (location) )
    throw ExceptionT::kDatabaseFail;
}

void AbaqusResultsT::ScanVariable (int key, int outputmode, int location)
{
  int index = VariableKeyIndex (key);
  if (index < 0)
    {
      fMessage << "\n AbaqusResultsT::ScanVariable Unable to map variable key: "
	       << key << endl;
      throw ExceptionT::kDatabaseFail;
    }

  /* quick exit, if we have already done this variable */
  if (fVariableTable[index].Dimension () != AbaqusVariablesT::kNotUsed &&
      fVariableTable[index].Type() != AbaqusVariablesT::kNotUsed &&
      fVariableTable[index].IOIndex() != AbaqusVariablesT::kNotUsed) return;

  /* number of components for this variable */
  int dim = fCurrentLength;
  switch (fVariableTable[index].FirstAttribute())
    {
    case AbaqusVariablesT::kNodeNumber:
    case AbaqusVariablesT::kIntType:
    case AbaqusVariablesT::kCharType:
      dim --;
      break;
    }

  int ioindex;
  AbaqusVariablesT::TypeT t;
  switch (outputmode)
    {
    case kElementOutput:
      {
	switch (location)
	  {
	  case kElementQuadrature:  /* quadrature point */
	    t = AbaqusVariablesT::kQuadrature;
	    ioindex = fNumQuadVars;
	    fNumQuadVars += dim;
	    break;
	  case kElementCentroidal: /* centroid of element */
	  case kElementWhole: /* whole element */
	    t = AbaqusVariablesT::kElement;
	    ioindex = fNumElemVars;
	    fNumElemVars += dim;
	    break;
	  case kElementNodal: /* node data */
	  case kElementNodeAveraged: /* nodal averaged */
	    t = AbaqusVariablesT::kNode;
	    ioindex = fNumNodeVars;
	    fNumNodeVars += dim;
	    break;
	  default:
	    {
	      fMessage << "\n AbaqusResultsT::Scan Variable, unrecognized location"
		       << location << endl;
	      throw ExceptionT::kDatabaseFail;
	    }
	  }
	break;
      }
    case kNodalOutput:
      {
	t = AbaqusVariablesT::kNode;
	ioindex = fNumNodeVars;
	fNumNodeVars += dim;
	break;
      }
    default:
      {
	fMessage << "\n AbaqusResultsT::Scan Variable, Unrecognize output mode."
		 << outputmode << endl;
	throw ExceptionT::kDatabaseFail;
      }
    }

  fVariableTable[index].SetIOData (dim, t, ioindex);
}

void AbaqusResultsT::WriteElementHeader (int key, int number, int intpt, int secpt, 
					 AbaqusResultsT::ElementVarType flag, int numdirect, 
					 int numshear, int numdir, int numsecforc)
{
#pragma unused(key)
//  int index = VariableKeyIndex (key);
  int rebarname = 0;
  WriteASCII ("*");
  Write (11);
  Write (ELEMENTHEADER);
  Write (number);
  Write (intpt);
  Write (secpt);
  Write (flag);
  Write (rebarname);
  Write (numdirect);
  Write (numshear);
  Write (numdir);
  Write (numsecforc);
}

/* this function tells the variable read function if there is an additional
   integer data, ususally a node number, in the variable record */
bool AbaqusResultsT::VariableWrittenWithNodeNumber (int key) const
{
  int index = VariableKeyIndex (key);
  if (index > 0 && index < NVT &&
      fVariableTable[index].FirstAttribute() == AbaqusVariablesT::kNodeNumber)
    return true;
  return false;
}

bool AbaqusResultsT::DataPoints (AbaqusVariablesT::TypeT vt, const StringT& name, iArrayT& set) const
{
  bool subset = true;
  if (name.Length() < 2) subset = false;
//  int numquadpts = 0;
  switch (vt)
    {
    case AbaqusVariablesT::kNode:
      {
	if (subset)
	  {
	    set.Allocate (NumNodesInSet (name));
	    NodeSet (name, set);
	  }
	else
	  {
	    set.Allocate (fNumNodes);
	    NodeMap (set);
	  }
	break;
      }
    case AbaqusVariablesT::kElement:
    case AbaqusVariablesT::kQuadrature:
      {
	if (subset)
	  {
	    set.Allocate (NumElements (name));
	    ElementSet (name, set);
//	    numquadpts = NumElementQuadPoints (name);
	  }
	else
	  {
	    set.Allocate (fNumElements);
	    ElementMap (set);
//	    numquadpts = NumElementQuadPoints (fElementSetNames[0]);
	  }
	break;
      }
    }
  return subset;
}

void AbaqusResultsT::AdvanceToTimeIncrement (int step)
{
  int number;
  double xtime;
  if (fModeIncs.Length() > 0)
    ModeData (step, number, xtime);
  else
    TimeData (step, number, xtime);

  int currentinc = -1;
  int numrewind = 0;
  double time;
  while (currentinc != number)
    {
      if (fModeIncs.Length() > 0)
	{
	  if (!NextMode (currentinc, time))
	    currentinc = -1;
	}
      else
	{
	  if (!NextTimeSteps (currentinc, time))
	    currentinc = -1;
	}

      /* starting in file after this time step */
      if ((currentinc == -1 || currentinc > number) && numrewind == 0)
	{
	  if (numrewind == 0)
	    {
	      ResetFile ();
	      numrewind ++;
	    }
	  else
	    {
	      fMessage << "\n AbaqusResultsT::ReadVariable: Searching for number = "
		       << number << ", " << xtime
		       << " currently at " << currentinc << ", " << time
		       << endl;
	      throw ExceptionT::kDatabaseFail;
	    }
	}
    }
}

bool AbaqusResultsT::CorrectType (int outputmode, int objnum, int intpt, int location, AbaqusVariablesT::TypeT vt, int& ID) const
{
  switch (outputmode)
    {
    case kElementOutput:
      switch (location)
	{
	case 0: /* quadrature point */
	  {
	    if (vt == AbaqusVariablesT::kQuadrature) 
	      {
		ID = objnum;
		return true;
	      }
	    break;
	  }
	case 1: /* centroid of element */
	case 5: /* whole element */
	  {
	    if (vt == AbaqusVariablesT::kElement) 
	      {
		ID = objnum;
		return true;
	      }
	    break;
	  }
	case 2: /* node data */
	  {
	    if (vt == AbaqusVariablesT::kNode) 
	      {
		ID = intpt;
		return true;
	      }
	    break;
	  }
	case 4: /* nodeal averaged */
	  {
	    if (vt == AbaqusVariablesT::kNode) 
	      {
		ID = objnum;
		return true;
	      }
	    break;
	  }
	}
    case kNodalOutput:
      {
	if (vt == AbaqusVariablesT::kNode) return true;
	break;
      }
    }
  return false;
}

int AbaqusResultsT::TranslateElementName (const char *name, GeometryT::CodeT &type, int &numintpts) const
{
  if (strncmp (name, "C", 1) == 0)
    return TranslateContinuum (name+1, type, numintpts);
  else if (strncmp (name, "DC", 2) == 0)
    return TranslateContinuum (name+2, type, numintpts);
  else if (strncmp (name, "AC", 2) == 0)
    return TranslateContinuum (name+2, type, numintpts);
  else if (strncmp (name, "DCC", 3) == 0)
    return TranslateContinuum (name+3, type, numintpts);
  else if (strncmp (name, "SPRING", 6) == 0)
    return TranslateSpring (name+6, type, numintpts);
  else if (strncmp (name, "S", 1) == 0) 
    return TranslateShell (name+1, type, numintpts);
  else if (strncmp (name, "R", 1) == 0)
    return TranslateRigid (name+1, type, numintpts);
  return BAD;
}

int AbaqusResultsT::TranslateContinuum (const char *name, GeometryT::CodeT &type, int &numintpts) const
{
  if (strncmp (name, "PE", 2) == 0)
    return Translate2D (name+2, type, numintpts);
  else if (strncmp (name, "PS", 2) == 0)
    return Translate2D (name+2, type, numintpts);
  else if (strncmp (name, "2D", 2) == 0)
    return Translate2D (name+2, type, numintpts);
  else if (strncmp (name, "GPE", 3) == 0)
    return Translate2D (name+3, type, numintpts);
  else if (strncmp (name, "3D", 2) == 0)
    return Translate3D (name+2, type, numintpts);
  return BAD;
}

int AbaqusResultsT::TranslateSpring (const char *name, GeometryT::CodeT &type, int &numintpts) const
{
  if (name[0] == '1')
    {
      type = GeometryT::kPoint;
      numintpts = 1;
      return OKAY;
    }
  else if (name[0] == '2' || name[0] == 'A')
    {
      type = GeometryT::kLine;
      numintpts = 2;
      return OKAY;
    }
  return BAD;
}

int AbaqusResultsT::Translate2D (const char *name, GeometryT::CodeT &type, int &numintpts) const
{
  if (strncmp (name, "3", 1) == 0)
    {
      type = GeometryT::kTriangle;
      numintpts = 1;
      return OKAY;
    }
  else if (strncmp (name, "4", 1) == 0)
    {
      type = GeometryT::kQuadrilateral;
      numintpts = 4;
      if (strchr (name, 'R') != NULL)
	numintpts = 1;
      return OKAY;
    }
  else if (strncmp (name, "6", 1) == 0)
    {
      type = GeometryT::kTriangle;
      numintpts = 3;
      return OKAY;
    }
  else if (strncmp (name, "8", 1) == 0)
    {
      type = GeometryT::kQuadrilateral;
      numintpts = 9;
      if (strchr (name, 'R') != NULL)
	numintpts = 4;
      return OKAY;
    }
  return BAD;
}

int AbaqusResultsT::Translate3D (const char *name, GeometryT::CodeT &type, int &numintpts) const
{
  if (strncmp (name, "4", 1) == 0)
    {
      type = GeometryT::kTetrahedron;
      numintpts = 1;
      return OKAY;
    }
  else if (strncmp (name, "6", 1) == 0)
    {
      type = GeometryT::kPentahedron;
      numintpts = 2;
      return OKAY;
    }
  else if (strncmp (name, "8", 1) == 0)
    {
      type = GeometryT::kHexahedron;
      numintpts = 8;
      if (strchr (name, 'R') != NULL)
	numintpts = 2;
      return OKAY;
    }
  else if (strncmp (name, "10", 2) == 0)
    {
      type = GeometryT::kTetrahedron;
      numintpts = 4;
      return OKAY;
    }
  else if (strncmp (name, "15", 2) == 0)
    {
      type = GeometryT::kPentahedron;
      numintpts = 6;
      return OKAY;
    }
  else if (strncmp (name, "20", 2) == 0)
    {
      type = GeometryT::kHexahedron;
      numintpts = 18;
      if (strchr (name, 'R') != NULL)
	numintpts = 8;
      return OKAY;
    }
  return BAD;
}

int AbaqusResultsT::TranslateShell (const char *name, GeometryT::CodeT &type, int &numintpts) const
{
  if (strncmp (name, "3", 1) == 0)
    {
      type = GeometryT::kTriangle;
      numintpts = 1;
      return OKAY;
    }
  else if (strncmp (name, "4", 1) == 0)
    {
      type = GeometryT::kQuadrilateral;
      numintpts = 4;
      if (strchr (name, 'R') != NULL)
	numintpts = 1;
      return OKAY;
    }
  else if (strncmp (name, "8R", 2) == 0)
    {
      type = GeometryT::kQuadrilateral;
      numintpts = 4;
      return OKAY;
    }
  return BAD;
}

int AbaqusResultsT::TranslateRigid (const char *name, GeometryT::CodeT &type, int &numintpts) const
{
  if (strncmp (name, "3D", 2) == 0) // send to 2D because it is rigid
    return Translate2D (name+2, type, numintpts);
  else if (strncmp (name, "2D", 2) == 0)
    {
      type = GeometryT::kPoint;
      numintpts = 1;
      return OKAY;
    }
  else if (strncmp (name, "AX", 2) == 0)
    {
      type = GeometryT::kPoint;
      numintpts = 1;
      return OKAY;
    }
  return BAD;
}

void AbaqusResultsT::GetElementName (GeometryT::CodeT geometry_code, int elemnodes, int& num_output_nodes, StringT& elem_name) const
{
  switch (geometry_code)
    {
    case GeometryT::kPoint:
      {
	elem_name = "MASS";	
	num_output_nodes = 1;
	break;
      }
    case GeometryT::kTriangle:
      {
	elem_name =  "CPE";
	num_output_nodes = (elemnodes < 6) ? 3 : 6;
	elem_name.Append (num_output_nodes);
	break;
      }
    case GeometryT::kQuadrilateral:
      {
	elem_name =  "CPE";
	num_output_nodes = (elemnodes < 8) ? 4 : 8;
	elem_name.Append (num_output_nodes);
	break;
      }
    case GeometryT::kHexahedron:
      {
	elem_name = "C3D";
	num_output_nodes = (elemnodes < 20) ? 8 : 20;
	elem_name.Append (num_output_nodes);
	break;
      }
    case GeometryT::kTetrahedron:
      {
	elem_name =  "C3D";
	num_output_nodes = (elemnodes < 10) ? 4 : 10;
	elem_name.Append (num_output_nodes);
	break;
      }
    case GeometryT::kPentahedron:
      {
	elem_name = "C3D";
	num_output_nodes = (elemnodes < 15) ? 6 : 15;
	elem_name.Append (num_output_nodes);
	break;
      }
    default:
      {
	fMessage << "\n AbaqusResultsT::GetElementName: cannot find name from geometry code: " << geometry_code << endl;
	throw ExceptionT::kDatabaseFail;
      }
    }
}

bool AbaqusResultsT::AdvanceTo (int target)
{
  int key = 0, error = OKAY;
  while (error == OKAY)
    {
      error = ReadNextRecord (key);
      if (key == target) return true;
    }
  //fMessage << "\n AbaqusResultsT::AdvanceTo: anable to advance to " << target << endl;
  return false;
}

bool AbaqusResultsT::SkipAttributes (void)
{
  if (fBinary)
    {
      if (fCurrentLength > 0)
	{
	  StringT temp;
	  if (!Read (temp, fCurrentLength))
	    return false;
	}
      fCurrentLength = 0;
      return true;
    }
  else
    {
      while (CheckBufferSize (fIn, 4))
	{
	  for (int i=fBufferDone; i < fBufferSize; i++)
	    {
	      if (fBuffer[fBufferDone++] == fMarker[0])
		{
		  fCurrentLength = 0;
		  return true;
		}
	    }
	}
      return false;
    }
}

int AbaqusResultsT::ReadNextRecord (int& key)
{
  if (!SkipAttributes ()) 
    return END;

  if (!fIn.good() || fIn.eof()) return END;

  int length = 0;
  if (!Read (length) || !Read (key)) 
    {
      if (!fIn.good() || fIn.eof()) return END;
      return BAD;
    }
  //cout << "length = " << length << " key=" << key << endl;
  fCurrentLength += length;
  return OKAY;
}

bool AbaqusResultsT::ReadSetName (StringT& s, int n)
{
  if (!Read (s, n)) return false;
  s.DropTrailingSpace ();
  s.DropLeadingSpace ();
  return true;
}

bool AbaqusResultsT::Read (StringT& s, int n)
{
  ArrayT<char> temp (n*kDoubleSize + 1);
  char *ps = temp.Pointer();
  for (int i=0; i < n; i++)
    {
      if (fBinary)
	{
	  CheckBufferSize (fIn);
	  fIn.read (ps, kDoubleSize);
	  if (fIn.eof ()) return false;
	  fBufferDone += kDoubleSize;
	  ps += kDoubleSize;
	}
      else
	{
	  if (!CheckBufferSize (fIn, 9)) 
	  	return false;
	  char c = fBuffer [fBufferDone++];
	  if (c != 'A') 
	  	return false;

	  for (int j=0; j < 8; j++)
	    *ps++ = fBuffer [fBufferDone++];
	}
    }
  *ps = '\0';
  s.Clear ();
  s.Append (temp.Pointer());
  fCurrentLength -= n;
  return true;
}

bool AbaqusResultsT::Read (int& i)
{
  if (fBinary)
    {
      CheckBufferSize (fIn);
      int temp;
      if (fIn.eof()) return false;
      fIn.read (reinterpret_cast<char *> (&temp), kDoubleSize);
      i = temp;
      fBufferDone += kDoubleSize;
    }
  else
    {
      // assume the maximum number of digits in a written integer is 10
      // the end increment (last entry) is always padded with spaces, 
      // so assumption should work
      if (!CheckBufferSize (fIn, 13)) return false;

      char c = fBuffer [fBufferDone++];
      if (c != 'I') return false;

      char w[3] = {'\0', '\0', '\0'};
      strncpy (&w[0], fBuffer.Pointer (fBufferDone), 2);
      fBufferDone += 2;
      int width = (int) atof (&w[0]);
 
      ArrayT<char> num (width+1);
      num.CopyPart (0, fBuffer, fBufferDone, width);
      fBufferDone += width;
      num [width] = '\0';
      i = (int) atof (num.Pointer());
    }
  fCurrentLength--;
  return true;
}

bool AbaqusResultsT::Read (double& d)
{
  if (fBinary)
    {
      CheckBufferSize (fIn);
      double temp;
      fIn.read (reinterpret_cast<char *> (&temp), kDoubleSize);
      d = temp;
      fBufferDone += kDoubleSize;
    }
  else
    {
      if (!CheckBufferSize (fIn, 23)) return false;

      char c = fBuffer [fBufferDone++];
      if (c != 'D') return false;

      ArrayT<char> num (19);
      num.CopyPart (0, fBuffer, fBufferDone, 18);
      fBufferDone += 18;
      num[18] = '\0';
      double base = atof (num.Pointer());
      
      c = fBuffer [fBufferDone++];
      if (c != 'D' && c != 'E') return false;

      char sign = fBuffer [fBufferDone++];

      ArrayT<char> expon (3);
      num.CopyPart (0, fBuffer, fBufferDone, 2);
      fBufferDone += 2;
      int exponent = (int) atof (num.Pointer());

      if (sign == '+' || sign == ' ')
	d = base * pow (10.0, exponent);
      else
	d = base * pow (10.0, -exponent);
    }
  fCurrentLength--;
  return true;
}

bool AbaqusResultsT::CheckBufferSize (istream& in, int numchars)
{
#pragma unused(in)
  
  if (fBinary) 
    return true;

  if (fBufferSize == 0 || (fBufferSize - fBufferDone) < numchars)
    {
      char temp [200];
      temp[0] = '\0';
      if (fBufferSize > 0) 
	strcpy (&temp[0], fBuffer.Pointer (fBufferDone));
      char nextline [90];

      if (!fIn.good() || fIn.eof()) return false;
      fIn.getline (&nextline[0], 89, '\n');
      
      strcat (temp, &nextline[0]);
      fBuffer = temp;
      fBufferSize = fBuffer.Length()-1;
      fBufferDone = 0;
      if (fBufferSize < numchars) return false;
    }
  return true;
}

void AbaqusResultsT::CheckBufferSize (istream& in)
{
  if (!fBinary) return;
  
  // FORTRAN footer
  if (fBufferDone == fBufferSize)
    {
      in.read (reinterpret_cast<char *> (&fBufferSize), kIntSize);
      fBufferDone = 0;
    }
  
  // FORTRAN header
  if (fBufferDone == 0)
    in.read (reinterpret_cast<char *> (&fBufferSize), kIntSize);
}

void AbaqusResultsT::Write (int i)
{
  if (fBinary)
    {
      CheckBufferSize (fOut);
      fOut.write (reinterpret_cast<char *> (&i), kDoubleSize);
      fBufferDone += kDoubleSize;
    }
  else
    {
      StringT s ("I ");
      StringT itext;
      itext.Append (i);
      s.Append (itext.Length() - 1);
      s.Append (itext);
      WriteASCII (s);
    }
}

void AbaqusResultsT::Write (double d)
{
  if (fBinary)
    {
      CheckBufferSize (fOut);
      fOut.write (reinterpret_cast<char *> (&d), kDoubleSize);
      fBufferDone += kDoubleSize;
    }
  else
    {
      double temp = d;
      StringT s ("D ");
      if (temp < 0)
	{
	  s = "D-";
	  temp = temp *-1;
	}

      int whole, exponent = 0;
      if (temp != 0)
	{
	  while (temp > 10)
	    {
	      temp = temp/10;
	      exponent++;
	    }
	  while (temp < 1)
	    {
	      temp = temp*10;
	      exponent--;
	    }
	}
      whole = (int) temp;
      temp = temp - whole;

      // append whole number
      s.Append (whole);

      // append fraction
      s.Append (".");
      for (int i=0; i < dprecision; i++)
	{
	  temp = temp *10;
	  int nextdigit = (int) temp;
	  s.Append (nextdigit);
	  temp = temp - nextdigit;
	}
      
      // append exponent
      if (exponent < 0)
	{
	  if (exponent < -99) exponent = -99;
	  exponent = exponent *-1;
	  s.Append ("D-");
	}
      else
	s.Append ("D+");
      s.Append (exponent, 2);

      WriteASCII (s);
    }
}

void AbaqusResultsT::Write (const StringT& s, int blocks)
{
  const char *ps = s.Pointer();
  if (fBinary)
    {
      CheckBufferSize (fOut);
      fOut.write (ps, kDoubleSize*blocks);
      fBufferDone += kDoubleSize*blocks;
    }
  else
    {
      for (int i=0, j=0; i < blocks; i++, j+= 8)
	{
	  StringT w = "A";
	  int copylength = 8;
	  if (j+8 > s.Length() - 1)
	    copylength = s.Length() - j - 1;
	  if (copylength < 0)
	    copylength = 0;
	  for (int k=j; k < j + copylength; k++)
	    w.Append (s[k]);
	  for (int m=copylength; m < 8; m++)
	    w.Append (' ');
	  WriteASCII (w);
	}
    }
}

void AbaqusResultsT::WriteASCII (const StringT& s)
{
  if (fBinary) return;
  
  if (fBufferDone + s.Length() < fBufferSize)
    {
      fOut << s;
      fBufferDone += s.Length() - 1;
    }
  else
    {
      for (int i=0; i < s.Length() - 1; i++)
	{
	  if (fBufferDone == fBufferSize)
	    {
	      fOut << '\n';
	      fBufferDone = 0;
	    }
	  fOut << s[i];
	  fBufferDone++;
	}
    }
}

void AbaqusResultsT::CheckBufferSize (ostream& out)
{
  if (!fBinary) return;

  // FORTRAN footer
  if (fBufferDone == fBufferSize)
    {
      out.write (reinterpret_cast<char *> (&fBufferSize), kIntSize);
      fBufferDone = 0;
    }

  // FORTRAN header
  if (fBufferDone == 0)
    out.write (reinterpret_cast<char *> (&fBufferSize), kIntSize);
}

void AbaqusResultsT::SetVariableNames (void)
{
  fVariableTable.Allocate (NVT);
  
  int i=0;
  /* Record Type, Record Key, First Attribute is a Node Number, Origin of Data */
  fVariableTable[i++].Set ("S", 11, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // stress components
  fVariableTable[i++].Set ("SP", 401, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // Principal stress components
  fVariableTable[i++].Set ("SINV", 12, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // mises, tresca, hydrostatic
  fVariableTable[i++].Set ("MISES", 75, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // 5.8 mises
  fVariableTable[i++].Set ("ALPHA", 86, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // back-stress tensor
  fVariableTable[i++].Set ("ALPHAP", 402, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // principal back-stress
  fVariableTable[i++].Set ("E", 21, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // strain components
  fVariableTable[i++].Set ("EP", 403, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // principle strain components
  fVariableTable[i++].Set ("LE", 89, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // logarithmic strains
  fVariableTable[i++].Set ("DG", 30, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // deformation gradient
  fVariableTable[i++].Set ("EE", 25, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // total elastic strains
  fVariableTable[i++].Set ("IE", 24, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // total inelastic strains
  fVariableTable[i++].Set ("PE", 22, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // plastic strains
  fVariableTable[i++].Set ("ENER", 14, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // energy densities
  fVariableTable[i++].Set ("SDV", 5, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // state dependent variables
  fVariableTable[i++].Set ("TEMP", 2, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // temperature
  fVariableTable[i++].Set ("FV", 9, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // field variables
  fVariableTable[i++].Set ("UVARM", 87, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // user defined
  fVariableTable[i++].Set ("LOCALDIR", 85, AbaqusVariablesT::kData, AbaqusVariablesT::kElementIntegration); // local coordinate directions
  	   
  fVariableTable[i++].Set ("LOADS", 3, AbaqusVariablesT::kCharType, AbaqusVariablesT::kElementWhole); // distributed loads
  	   
  fVariableTable[i++].Set ("U", 101, AbaqusVariablesT::kNodeNumber, AbaqusVariablesT::kNodePoint); // displacment
  fVariableTable[i++].Set ("V", 102, AbaqusVariablesT::kNodeNumber, AbaqusVariablesT::kNodePoint); // velocity
  fVariableTable[i++].Set ("A", 103, AbaqusVariablesT::kNodeNumber, AbaqusVariablesT::kNodePoint); // acceleration
  fVariableTable[i++].Set ("NT", 201, AbaqusVariablesT::kNodeNumber, AbaqusVariablesT::kNodePoint); // nodal temperature

  if (i != NVT)
    {
      fMessage << "\n AbaqusResultsT::SetVariableNames, incorrect allocation" << endl;
      throw ExceptionT::kDatabaseFail;
    }
}
