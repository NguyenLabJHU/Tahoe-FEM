/*
   CREATED: S. Wimmer 9 Nov 2000

*/

#include "AbaqusResultsT.h"

/* these variables are nodal and have a node number before the value list */
AbaqusResultsT::AbaqusResultsT (ostream& message) :
  fMessage (message),
  fNumElements (0),
  fNumNodes (0),
  fStartCount (0),
  fEndCount (0),
  fModalCount (0),
  fVarDimension (NVT),
  fVarType (NVT),
  fVarArrayColumn (NVT),
  fBinary (false),
  fBufferDone (0),
  fBufferSize (0),
  fCurrentLength (-1),
  fNumNodeSets (0),
  fNumElementSets (0),
  fNumNodeVars (0),
  fNumElemVars (0),
  fNumQuadVars (0)
{
  fVarDimension = BAD;
  fVarType = kNone;
  fVarArrayColumn = -1;
}

void AbaqusResultsT::Initialize (char *filename)
{
  fIn.open (filename);
  if (!fIn)
    {
      fMessage << "\n\nAbaqusResultsT::Initialize unable to open file " << filename << "\n\n";
      throw eDatabaseFail;
    }
  fFileName = filename;

  int key = -1;
  ReadNextRecord (key);
  if (key != VERSION) 
    {
      fMessage << "\n\nAbaqusResultsT::Initialize not ASCII, trying Binary.\n";
      fBinary = true;
      ResetFile ();
      ReadNextRecord (key);
      if (key != VERSION)
	{
	  fMessage << "\n\nAbaqusResultsT::Initialize not Binary.\n";
	  throw eDatabaseFail;
	}
    }

  if (!ReadVersion ()) 
    {
      fMessage << "\n\nAbaqusResultsT::Initialize unreadable file\n";
      throw eDatabaseFail;
    }
}

void AbaqusResultsT::Close (void)
{
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
  fVarDimension = BAD;
  fVarType = kNone;
  fVarArrayColumn = -1;
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
}

void AbaqusResultsT::ScanFile (int &numelems, int &numnodes, int &numtimesteps, int &nummodes)
{
  int key = 0, error = OKAY, outputmode = -1;
  int location = -1;
  fNumElements = 0;
  fNumNodes = 0;
  fModalCount = 0;
  fStartCount = 0;
  fEndCount = 0;
  int count = 0;
  error = ReadNextRecord (key);
  while (error == OKAY)
    {
      //cout << fCurrentLength << " " << key << endl;
      //if (count++ == 10) exit(0);
      switch (key)
	{
	case ELEMENT: ScanElement (); break;
	case NODE: 
	  {
	    int number;
	    if (!Read (number))
	      throw eDatabaseFail;
	    fNodeNumber.Append (number);
	    fNumNodes++; 
	    break;
	  }
	case NODESET: 
	  {
	    StringT name;
	    if (!Read (name, 1))
	      throw eDatabaseFail;
	    fNodeSetNames.Append (name);
	    fNumNodeSets++; 
	    break;
	  }
	case ELEMENTSET: 
	  {
	    StringT name;
	    if (!Read (name, 1))
	      throw eDatabaseFail;
	    fElementSetNames.Append (name);
	    fNumElementSets++; 
	    break;
	  }
	case MODAL: 
	  {
	    int number;
	    double mode;
	    if (!Read (number) || !Read (mode) )
	      throw eDatabaseFail;

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
	      throw eDatabaseFail;

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
	    if (fStartCount == 1)
	      ReadOutputDefinitions (outputmode);
	    break;
	  }
	}
      
      /* scan variable records */
      if (fStartCount == 1)
	{
	  VariableKeyT k = IntToVariableKey (key);
	  if (k != kNone)
	    ScanVariable (k, outputmode, location);
	}

      error = ReadNextRecord (key);
    }
  numelems = fNumElements;
  numnodes = fNumNodes;
  numtimesteps = fStartCount;
  nummodes = fModalCount;

  //cout << " NumNodeSets: " << fNumNodeSets 
  //<< "\n NumElemSets: " << fNumElementSets << "\n\n";
}

void AbaqusResultsT::ElementSetNames (ArrayT<StringT>& names) const
{
  for (int i=0; i < fNumElementSets; i++)
    names[i] = fElementSetNames[i];
}

void AbaqusResultsT::NodeSetNames (ArrayT<StringT>& names) const
{
  for (int i=0; i < fNumNodeSets; i++)
    {
      //cout << i << " " << fNumNodeSets << " " << fNodeSetNames[i] << endl;
      names[i] = fNodeSetNames[i];
    }
}

void AbaqusResultsT::NodeMap (iArrayT& n) const
{
  n.CopyPart (0, fNodeNumber, 0, fNumNodes);
}

int AbaqusResultsT::NumNodesInSet (StringT& name)
{
  ResetFile ();
  StringT setname = "";

  while (strncmp (setname.Pointer(), name.Pointer(), name.Length()) != 0)
    {
      //cout << setname << " " << name << endl;
      AdvanceTo (NODESET);
      if (!Read (setname, 1))
	throw eDatabaseFail;
    }

  int num = fCurrentLength;
  int key;
  ReadNextRecord (key);
  while (key == NODESETCONT)
    {
      num += fCurrentLength;
      ReadNextRecord (key);
    }
  
  return num;
}

void AbaqusResultsT::NodeSet (StringT& name, iArrayT& nset)
{
  ResetFile ();
  StringT setname = "";

  while (strncmp (setname.Pointer(), name.Pointer(), name.Length()) != 0)
    {
      AdvanceTo (NODESET);
      if (!Read (setname, 1))
	throw eDatabaseFail;
    }

  int e = 0;
  int num = fCurrentLength;
  for (int j=0; j < num; j++)
    if (!Read (nset[e++])) throw eDatabaseFail;

  int key;
  ReadNextRecord (key);
  while (key == NODESETCONT)
    {
      num = fCurrentLength;
      for (int i=0; i < num; i++)
	if (!Read (nset[e++])) throw eDatabaseFail;
      ReadNextRecord (key);
    }
}

void AbaqusResultsT::ElementMap (iArrayT& e) const
{
  e.CopyPart (0, fElementNumber, 0, fNumElements);
}

int AbaqusResultsT::NumElements (StringT& name) 
{
  ResetFile ();
  StringT setname = "";

  while (strncmp (setname.Pointer(), name.Pointer(), name.Length()) != 0)
    {
      AdvanceTo (ELEMENTSET);
      if (!Read (setname, 1))
	throw eDatabaseFail;
      //cout << setname << ". " << name << "." << endl;
    }

  int num = fCurrentLength;
  int key;
  ReadNextRecord (key);
  while (key == ELEMSETCONT)
    {
      num += fCurrentLength;
      ReadNextRecord (key);
    }
  
  return num;
}

int AbaqusResultsT::NumElementNodes (StringT& name)
{
  ResetFile ();
  StringT setname = "";

  while (strncmp (setname.Pointer(), name.Pointer(), name.Length()) != 0)
    {
      AdvanceTo (ELEMENTSET);
      if (!Read (setname, 1)) throw eDatabaseFail;
    }

  int el;
  if (!Read (el)) throw eDatabaseFail;

  ResetFile ();
  int cel=-1;
  while (cel != el)
    {
      AdvanceTo (ELEMENT);
      if (!Read (cel)) throw eDatabaseFail;
    }

  return fCurrentLength - 1;
}

int AbaqusResultsT::NumElementQuadPoints (StringT& name)
{
  ResetFile ();
  StringT setname = "";

  while (strncmp (setname.Pointer(), name.Pointer(), name.Length()) != 0)
    {
      AdvanceTo (ELEMENTSET);
      if (!Read (setname, 1)) throw eDatabaseFail;
    }

  int el;
  if (!Read (el)) throw eDatabaseFail;

  ResetFile ();
  int cel=-1;
  while (cel != el)
    {
      AdvanceTo (ELEMENT);
      if (!Read (cel)) throw eDatabaseFail;
    }

  int numintpts;
  GeometryT::CodeT code;
  StringT elname;
  if (!Read (elname, 1)) throw eDatabaseFail;
  if (TranslateElementName (elname.Pointer(), code, numintpts) == BAD)
    {
      fMessage << "\nAbaqusResultsT::GeometryCode Unable to translate element name.\n\n";
      throw eDatabaseFail;
    }

  return numintpts;
}

void AbaqusResultsT::ElementSet (StringT& name, iArrayT& elset) 
{
  ResetFile ();
  StringT setname = "";

  while (strncmp (setname.Pointer(), name.Pointer(), name.Length()) != 0)
    {
      AdvanceTo (ELEMENTSET);
      if (!Read (setname, 1))
	throw eDatabaseFail;
    }

  int e = 0;
  int num = fCurrentLength;
  for (int j=0; j < num; j++)
    if (!Read (elset[e++])) throw eDatabaseFail;

  int key;
  ReadNextRecord (key);
  while (key == ELEMSETCONT)
    {
      num = fCurrentLength;
      for (int i=0; i < num; i++)
	if (!Read (elset[e++])) throw eDatabaseFail;
      ReadNextRecord (key);
    }
}

void AbaqusResultsT::GeometryCode (StringT& name, GeometryT::CodeT& code)
{
  ResetFile ();
  StringT setname = "";

  while (strncmp (setname.Pointer(), name.Pointer(), name.Length()) != 0)
    {
      AdvanceTo (ELEMENTSET);
      if (!Read (setname, 1)) throw eDatabaseFail;
    }

  int el;
  if (!Read (el)) throw eDatabaseFail;

  ResetFile ();
  int cel=-1;
  while (cel != el)
    {
      AdvanceTo (ELEMENT);
      if (!Read (cel)) throw eDatabaseFail;
    }

  int numintpts;
  StringT elname;
  if (!Read (elname, 1)) throw eDatabaseFail;
  if (TranslateElementName (elname.Pointer(), code, numintpts) == BAD)
    {
      fMessage << "\nAbaqusResultsT::GeometryCode Unable to translate element name.\n\n";
      throw eDatabaseFail;
    }
}

void AbaqusResultsT::ModeData (int index, int &number, double &mode) const
{
  if (index < 0 || index > fModalCount)
    throw eOutOfRange;

  number = fModeIncs [index];
  mode = fModeSteps [index];
}

void AbaqusResultsT::TimeData (int index, int &number, double &time) const
{
  if (index < 0 || index > fStartCount)
    throw eOutOfRange;

  number = fTimeIncs [index];
  time = fTimeSteps [index];
}

void AbaqusResultsT::NodeVariables (ArrayT<VariableKeyT>& keys, iArrayT& dims) const
{
  int c = 0;
  for (int i=0; i < NVT; i++)
    if (fVarType[i] == kNodeVar)
      for (int j=0; j < fVarDimension[i]; j++)
	{
	  keys[c] = VariableKey (i);
	  dims[c] = fVarDimension[i];
	  c++;
	}
}

void AbaqusResultsT::ElementVariables (ArrayT<VariableKeyT>& keys, iArrayT& dims) const
{
  int c = 0;
  for (int i=0; i < NVT; i++)
    if (fVarType[i] == kElemVar)
      for (int j=0; j < fVarDimension[i]; j++)
	{
	  keys[c] = VariableKey (i);
	  dims[c] = fVarDimension[i];
	  c++;
	}
}

void AbaqusResultsT::QuadratureVariables (ArrayT<VariableKeyT>& keys, iArrayT& dims) const
{
  int c = 0;
  for (int i=0; i < NVT; i++)
    if (fVarType[i] == kQuadVar)
      for (int j=0; j < fVarDimension[i]; j++)
	{
	  keys[c] = VariableKey (i);
	  dims[c] = fVarDimension[i];
	  c++;
	}
}

void AbaqusResultsT::ReadVariables (VariableType vt, int step, dArray2DT& values, StringT& name)
{
  bool subset = true;
  if (name.Length() < 2) subset = false;
  iArrayT set;
  int numquadpts = 0;
  switch (vt)
    {
    case kNodeVar:
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
    case kElemVar:
    case kQuadVar:
      {
	if (subset)
	  {
	    set.Allocate (NumElements (name));
	    ElementSet (name, set);
	    numquadpts = NumElementQuadPoints (name);
	  }
	else
	  {
	    set.Allocate (fNumElements);
	    ElementMap (set);
	    numquadpts = NumElementQuadPoints (fElementSetNames[0]);
	  }
	break;
      }
    }

  // don't need to rewind after reading set
  int number;
  double time;
  if (fModeIncs.Length() > 0)
    ModeData (step, number, time);
  else
    TimeData (step, number, time);

  int currentinc = -1;
  while (currentinc != number)
    {
      if (fModeIncs.Length() > 0)
	NextMode (currentinc, time);
      else
	NextTimeSteps (currentinc, time);
    }

  int key;
  int ID, objnum, intpt, secpt, location, outputmode;
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
	  ReadOutputDefinitions (outputmode);
	  break;
	default:
	  {
	    if (CorrectType (outputmode, objnum, intpt, location, vt, ID))
	      {
		VariableKeyT k = IntToVariableKey (key);
		if (VariableWrittenWithNodeNumber (k)) 
		  if (!Read (ID)) throw eDatabaseFail;
		
		bool save = true;
		int row;
		set.HasValue (ID, row);
		if (row < 0 || row > set.Length())
		  save = false;

		// modify row by 
		if (vt == kQuadVar)
		  row = row*numquadpts + intpt - 1;

		if (save)
		  {
		    dArrayT v (fCurrentLength);
		    int num = fCurrentLength;
		    for (int i=0; i < num; i++)
		      if (!Read (v[i])) throw eDatabaseFail;

		    int index = VariableKeyIndex (k);
		    int offset = fVarArrayColumn [index];
		    //cout << ID << " " << intpt << " " << k << " ";
		    //cout << row << " " << values.MinorDim() << " " << offset << endl;
		    values.CopyPart (row*values.MinorDim() + offset, v, 0, v.Length()); 
		  }
	      }
	  }
	}
    }
}

const char *AbaqusResultsT::VariableName (int index) const
{
  AbaqusResultsT::VariableKeyT key = VariableKey (index);
  switch (key)
    {
    case kTEMP: return "TEMP";
    case kLOADS: return "LOADS";
    case kFLUXS: return "FLUXS";
    case kSDV: return "SDV";
    case kCOORD: return "COORD";
    case kS: return "S";
    case kSINV: return "SINV";
    case kE: return "E";
    case kPE: return "PE";
    case kCE: return "CE";
    case kIE: return "IE";
    case kEE: return "EE";
    case kU: return "U";
    case kV: return "V";
    case kA: return "A";
    case kNCOORD: return "COORD";
    case kSP: return "SP";
    case kEP: return "EP";
    case kNEP: return "NEP";
    case kLEP: return "LEP";
    case kERP: return "ERP";
    case kEEP: return "EEP";
    case kIEP: return "IEP";
    case kTHEP: return "THEP";
    case kPEP: return "PEP";
    case kCEP: return "CEP";
    }
  return "Unknown";
}

AbaqusResultsT::VariableKeyT AbaqusResultsT::VariableKey (int index) const
{
  switch (index)
    {
    case 0:  return kTEMP;	
    case 1:  return kLOADS;
    case 2:  return kFLUXS;
    case 3:  return kSDV;   
    case 4:  return kCOORD;
    case 5:  return kS;     
    case 6:  return kSINV;      
    case 7:  return kE;         
    case 8:  return kPE;        
    case 9:  return kCE;    
    case 10: return kIE;        
    case 11: return kEE;        
    case 12: return kU;     
    case 13: return kV;
    case 14: return kA;
    case 15: return kNCOORD;
    case 16: return kSP;    
    case 17: return kEP;    
    case 18: return kNEP;   
    case 19: return kLEP;   
    case 20: return kERP;   
    case 21: return kEEP;   
    case 22: return kIEP;   
    case 23: return kTHEP;  
    case 24: return kPEP;   
    case 25: return kCEP;   
    }
  return kNone;
}

AbaqusResultsT::VariableKeyT AbaqusResultsT::IntToVariableKey (int key) const
{
  switch (key)
    {
    case kTEMP: return kTEMP;	
    case kLOADS: return kLOADS;
    case kFLUXS: return kFLUXS;
    case kSDV: return kSDV;
    case kCOORD: return kCOORD;
    case kS: return kS;	
    case kSINV: return kSINV;
    case kE: return kE;	
    case kPE: return kPE;	
    case kCE: return kCE;	
    case kIE: return kIE;	
    case kEE: return kEE;	
    case kU: return kU;	
    case kV: return kV;
    case kA: return kA;
    case kNCOORD: return kNCOORD;
    case kSP: return kSP;	
    case kEP: return kEP;	
    case kNEP: return kNEP;
    case kLEP: return kLEP;
    case kERP: return kERP;
    case kEEP: return kEEP;
    case kIEP: return kIEP;
    case kTHEP: return kTHEP;
    case kPEP: return kPEP;
    case kCEP: return kCEP;
    }
  return kNone;
}

int AbaqusResultsT::VariableKeyIndex (AbaqusResultsT::VariableKeyT key) const
{
  switch (key)
    {
    case kTEMP:   return 0;
    case kLOADS:  return 1; 
    case kFLUXS:  return 2; 
    case kSDV:    return 3; 
    case kCOORD:  return 4; 
    case kS:      return 5; 
    case kSINV:   return 6; 
    case kE:      return 7; 
    case kPE:     return 8; 
    case kCE:     return 9; 
    case kIE:     return 10;
    case kEE:     return 11;
    case kU:      return 12;
    case kV:      return 13;
    case kA:	  return 14;
    case kNCOORD: return 15;
    case kSP:     return 16;
    case kEP:     return 17;
    case kNEP:    return 18;
    case kLEP:    return 19;
    case kERP:    return 20;
    case kEEP:    return 21;
    case kIEP:    return 22;
    case kTHEP:   return 23;
    case kPEP:    return 24;
    case kCEP:    return 25;
    }
  return -1;
}

bool AbaqusResultsT::NextCoordinate (int &number, dArrayT &nodes)
{
  AdvanceTo (NODE);

  if (!Read (number))
    throw eDatabaseFail;

  int dof = fCurrentLength;
  nodes.Allocate (dof);
  for (int ii=0; ii < dof; ii++)
    if (!Read (nodes[ii]))
      throw eDatabaseFail;
  return true;
}

bool AbaqusResultsT::NextElement (int &number, GeometryT::CodeT &type, iArrayT &nodes)
{
  AdvanceTo (ELEMENT);

  StringT name;
  if (!Read (number) || !Read (name, 1) )
    throw eDatabaseFail;

  //cout << number << " " << name << endl;

  int numintpts;
  if (TranslateElementName (name.Pointer(), type, numintpts) == BAD) 
    {
      fMessage << "\nAbaqusResultsT::NextElement Unable to translate element name.\n\n";
      throw eDatabaseFail;
    }
  //cout << type << " " << numintpts << endl;

  int numnodes = fCurrentLength;
  //cout << numnodes << endl;
  nodes.Allocate (numnodes);
  nodes = 300;
  //cout << nodes << endl;
  for (int i=0; i < numnodes; i++)
    {
      int temp;
      if (!Read (temp))
      throw eDatabaseFail;
      nodes[i] = temp;
    }
  return true;
}

void AbaqusResultsT::VersionNotes (ArrayT<StringT>& records)
{
  ResetFile ();
  AdvanceTo (VERSION);
  records.Allocate (4);
  int numelems, numnodes;
  double elemleng;
  if (!Read (records[1], 1) || !Read (records[2], 2) || 
      !Read (records[3], 1) || !Read (numelems) || !Read (numnodes) || 
      !Read (elemleng) )
    throw eDatabaseFail;
  records[0] = "ABAQUS";
}

void AbaqusResultsT::ResetFile (void)
{
  //printf ("\nRewinding File\n");
  fIn.close ();
  fIn.open (fFileName);

  if (fBinary)
    fBufferSize = 512 * sizeof (double);
  else
    fBufferSize = 0;
  fBufferDone = 0;
  fCurrentLength = 0;
}

/******************* PRIVATE **********************/

bool AbaqusResultsT::ReadVersion (void)
{
  //  cout << "ReadVersion: " << fCurrentLength << endl;
  StringT version, date, time;
  int numelems, numnodes;
  double elemleng;
  if (!Read (version, 1) || !Read (date, 2) || !Read (time, 1) ||
      !Read (numelems) || !Read (numnodes) || !Read (elemleng) )
    {
      fMessage << "\nAbaqusResultsT::ScanFile Unable to read version record\n\n";
      return false;
    }
  //cout << "ReadVersion: " << fCurrentLength << endl;

  /*cout << "\n Version: " << version << "\n    Date: "
    << date << "\n    Time: " << time << "\n NumElem: "
    << numelems << "\n NumNode: " << numnodes 
    << "\n ElemLen: " << elemleng << "\n\n";*/
  if (strncmp (version.Pointer(), "5.8", 3) != 0 &&
      strncmp (version.Pointer(), "6.1", 3) != 0 &&
      strncmp (version.Pointer(), "6.2", 3) != 0 )
    {
      fMessage << "\nAbaqusResultsT::ScanFile Unrecognized version " << version << "\n\n";
      return false;
    }
  return true;
}

void AbaqusResultsT::NextMode (int &number, double &mode)
{
  AdvanceTo (MODAL);
  if (!Read (number) || !Read (mode) )
    throw eDatabaseFail;
}

void AbaqusResultsT::NextTimeSteps (int &number, double &time)
{
  int procedure, step;
  double steptime, creep, amp; 
  AdvanceTo (STARTINCREMENT);
  if (!Read (time) || !Read (steptime) || !Read (creep) || !Read (amp) ||
      !Read (procedure) || !Read (step) || !Read (number) )
    throw eDatabaseFail;
}

void AbaqusResultsT::ScanElement (void)
{
  StringT name;
  GeometryT::CodeT type = GeometryT::kNone;
  int number, numintpts = 0;

  if (!Read (number) || !Read (name, 1)) 
    throw eDatabaseFail;

  if (TranslateElementName (name.Pointer(), type, numintpts) == BAD) 
    {
      fMessage << "\nAbaqusResultsT::ScanElement Encountered unknown element type\n\n";
      throw eDatabaseFail;
    }
  fNumElements++;
  fElementNumber.Append (number);
}

void AbaqusResultsT::ReadOutputDefinitions (int &outputmode)
{
  StringT setname, elemname;
  if (!Read (outputmode)|| !Read (setname, 1) )
    throw eDatabaseFail;
}

void AbaqusResultsT::ReadElementHeader (int &objnum, int& intpt, int& secpt, int &location)
{
  if (!Read (objnum) || !Read (intpt) || !Read (secpt) || !Read (location) )
    throw eDatabaseFail;
}

void AbaqusResultsT::ScanVariable (AbaqusResultsT::VariableKeyT key, int outputmode, int location)
{
  int index = VariableKeyIndex (key);
  if (index < 0)
    {
      fMessage << "\nAbaqusResultsT::ScanVariable Unable to map variable key\n\n";
      throw eDatabaseFail;
    }
  if (fVarDimension [index] != BAD) return;

  switch (outputmode)
    {
    case kElementOutput:
      {
	fVarDimension [index] = fCurrentLength;
	switch (location)
	  {
	  case 0:  /* quadrature point */
	    fVarType [index] = kQuadVar;
	    fVarArrayColumn [index] = fNumQuadVars;
	    fNumQuadVars += fVarDimension [index];
	    break;
	  case 1: /* centroid of element */
	  case 5: /* whole element */
	    fVarType [index] = kElemVar;
	    fVarArrayColumn [index] = fNumElemVars;
	    fNumElemVars += fVarDimension [index];
	    break;
	  case 2: /* node data */
	  case 4: /* nodal averaged */
	    fVarType [index] = kNodeVar;
	    fVarArrayColumn [index] = fNumNodeVars;
	    fNumNodeVars += fVarDimension [index];
	    break;
	  default:
	    {
	      fMessage << "\nAbaqusResultsT::Scan Variable, unrecognized location"
		       << location << "\n\n";
	      throw eDatabaseFail;
	    }
	  }
	/*printf ("  Variable: %i %s %i\n", key, 
	  VariableName (index), fVarDimension [index]);*/
	break;
      }
    case kNodalOutput:
      {
	fVarDimension [index] = fCurrentLength - 1;
	fVarType [index] = kNodeVar;
	fVarArrayColumn [index] = fNumNodeVars;
	fNumNodeVars += fVarDimension [index];
	/*printf ("  Variable: %i %s %i\n", key, 
	  VariableName (index), fVarDimension [index]);*/
	break;
      }
    default:
      {
	fMessage << "\nAbaqusResultsT::Scan Variable, Unrecognize output mode."
		 << outputmode << "\n\n";
	throw eDatabaseFail;
      }
    }
}

/* this function tells the variable read function if there is an additional
   integer data, ususally a node number, in the variable record */
bool AbaqusResultsT::VariableWrittenWithNodeNumber (AbaqusResultsT::VariableKeyT key) const
{
  switch (key)
    {
    case kLOADS: // load type not node number, but same idea
    case kFLUXS: // flux type not node number, but same idea
    case kU: 
    case kV:
    case kA:
    case kNCOORD:
      return true;
    }
  return false;
}

bool AbaqusResultsT::CorrectType (int outputmode, int objnum, int intpt, int location, VariableType vt, int& ID) const
{
  switch (outputmode)
    {
    case kElementOutput:
      switch (location)
	{
	case 0: /* quadrature point */
	  {
	    if (vt == kQuadVar) 
	      {
		ID = objnum;
		return true;
	      }
	    break;
	  }
	case 1: /* centroid of element */
	case 5: /* whole element */
	  {
	    if (vt == kElemVar) 
	      {
		ID = objnum;
		return true;
	      }
	    break;
	  }
	case 2: /* node data */
	  {
	    if (vt == kNodeVar) 
	      {
		ID = intpt;
		return true;
	      }
	    break;
	  }
	case 4: /* nodeal averaged */
	  {
	    if (vt == kNodeVar) 
	      {
		ID = objnum;
		return true;
	      }
	    break;
	  }
	}
    case kNodalOutput:
      {
	if (vt == kNodeVar) return true;
	break;
      }
    }
  return false;
}

int AbaqusResultsT::TranslateElementName (char *name, GeometryT::CodeT &type, int &numintpts)
{
  if (strncmp (name, "C", 1) == 0)
    return TranslateContinuum (name+1, type, numintpts);
  else if (strncmp (name, "DC", 2) == 0)
    return TranslateContinuum (name+2, type, numintpts);
  else if (strncmp (name, "AC", 2) == 0)
    return TranslateContinuum (name+2, type, numintpts);
  else if (strncmp (name, "DCC", 3) == 0)
    return TranslateContinuum (name+3, type, numintpts);
  else if (strncmp (name, "S", 1) == 0)
    return TranslateShell (name+1, type, numintpts);
  return BAD;
}

int AbaqusResultsT::TranslateContinuum (char *name, GeometryT::CodeT &type, int &numintpts)
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

int AbaqusResultsT::Translate2D (char *name, GeometryT::CodeT &type, int &numintpts)
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

int AbaqusResultsT::Translate3D (char *name, GeometryT::CodeT &type, int &numintpts)
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

int AbaqusResultsT::TranslateShell (char *name, GeometryT::CodeT &type, int &numintpts)
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

void AbaqusResultsT::AdvanceTo (int target)
{
  int key = 0, error = OKAY;
  while (error == OKAY)
    {
      error = ReadNextRecord (key);
      if (key == target) return;
    }
  fMessage << "Unable to advance to " << target << ".\n";
  throw eDatabaseFail;
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
	      if (fBuffer[fBufferDone++] == '*')
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

  int length;
  if (!Read (length) || !Read (key)) 
    {
      //fMessage << "\nAbaqusResultsT::ReadNextRecord Cannot Read Next Record.\n";
      return BAD;
    }
  fCurrentLength += length;
  return OKAY;
}

bool AbaqusResultsT::Read (StringT& s, int n)
{
  char c;
  ArrayT<char> temp (n * sizeof (double) + 1);
  char *ps = temp.Pointer();
  for (int i=0; i < n; i++)
    {
      if (fBinary)
	{
	  CheckBufferSize (fIn);
	  fIn.read (ps, sizeof (double));
	  if (fIn.eof ()) return false;
	  fBufferDone += sizeof (double);
	  ps += sizeof (double);
	}
      else
	{
	  if (!CheckBufferSize (fIn, 9)) return false;
	      
	  char c = fBuffer [fBufferDone++];
	  if (c != 'A') return false;

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
      fIn.read (reinterpret_cast<char *> (&temp), sizeof (double));
      i = temp;
      fBufferDone += sizeof (double);
    }
  else
    {
      // assume the maximum number of digits in a written integer is 10
      // the end increment (last entry) is always padded with spaces, 
      // so assumption should work
      if (!CheckBufferSize (fIn, 13)) return false;

      char c = fBuffer [fBufferDone++];
      if (c != 'I') return false;

      char w [3];
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
      fIn.read (reinterpret_cast<char *> (&temp), sizeof (double));
      d = temp;
      fBufferDone += sizeof (double);
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
	d = base * pow (10, exponent);
      else
	d = base * pow (10, -exponent);
    }
  fCurrentLength--;
  return true;
}

bool AbaqusResultsT::CheckBufferSize (istream& in, int numchars)
{
  if (fBinary) 
    return true;

  //cout << "Checking Buffer " << fBufferDone << " " << fBufferSize << " " << numchars << endl;
  if (fBufferSize == 0 || (fBufferSize - fBufferDone) < numchars)
    {
      //cout << "Updating Buffer" << endl;
      char temp [200];
      if (fBufferSize > 0) 
	strcpy (&temp[0], fBuffer.Pointer (fBufferDone));
      char nextline [90];
      if (!fIn.getline (&nextline[0], 89, '\n')) return false;
      
      strcat (temp, &nextline[0]);
      fBuffer = temp;
      //cout << "Buffer = " << fBuffer << ".end." << endl;
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
      in.read (reinterpret_cast<char *> (&fBufferSize), sizeof (int));
      fBufferDone = 0;
    }
  
  // FORTRAN header
  if (fBufferDone == 0)
    in.read (reinterpret_cast<char *> (&fBufferSize), sizeof (int));
}

