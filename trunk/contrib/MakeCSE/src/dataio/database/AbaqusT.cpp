// file: AbaqusT.cpp

// created:       SAW (05/30/00)

#include "AbaqusT.h"
#include "StringT.h"
#include "iosfwd.h"
#include <stdlib.h>
#include "iArrayT.h"
#include "GeometryT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "ArrayT.h"
#include "dArray2DT.h"
#include <time.h>

const char* kAbaqusVersion = "5.8-1";
const int kBufferSize = 512*sizeof (double);
const int kASCIIWrap = 80;

AbaqusT::AbaqusT (ostream& out, bool binary) :
  fBinary (binary),
  fOut (out),
  fBufferDone (0),
  fBufferSize (0),
  fModal (false)
{
  if (fBinary) fBufferSize = kBufferSize;
  else fBufferSize = kASCIIWrap;
}

void AbaqusT::ReadVersion (istream& in, StringT& version, StringT& date, StringT& time, int& elements, int& nodes)
{
  int length;
  double elemlength;
  if (!AdvanceStreamTo (in, aVersion, length))
    {
      cout << "\n\n AbaqusT unable to fine version record.\n\n";
      throw eGeneralFail;
    }

    Read (in, version);
    Read (in, date, 2);
    Read (in, time);
    Read (in, elements);
    Read (in, nodes);
    Read (in, elemlength);
}

void AbaqusT::IsModal (istream& in)
{
  int length;
  fModal = false;
  AbaqusT::KeyRecord key;
  bool search = true;
  if (AdvanceStreamTo (in, aModal, aEndIncrement, length, key, search)) 
    if (key == aModal) fModal = true;
}

void AbaqusT::ReadElementSets (istream& in, iArray2DT& data)
{
  int length, numsets = 0;
  AbaqusT::KeyRecord key;
  bool search = true;
  while (AdvanceStreamTo (in, aElementSet, aElemSetCont, length, key, search))
    {
      int num = length - 2;
      if (key == aElementSet)
	{
	  StringT name;
	  Read (in, name);
	  num--;
	  numsets++;
	}
      
      int element, dex;
      for (int i=0; i < num; i++)
	{
	  Read (in, element);
	  if (!data.ColumnHasValue (0, element, dex)) return;
	  data (dex, 2) = numsets;
	}
      search = false;
    }
}

void AbaqusT::ReadNodeSets (istream& in, iArray2DT& data)
{
  int length, numsets = 0;
  AbaqusT::KeyRecord key;
  bool search = true;
  while (AdvanceStreamTo (in, aNodeSet, aNodeSetCont, length, key, search))
    {
      int num = length - 2;
      if (key == aNodeSet)
	{
	  StringT name;
	  Read (in, name);
	  num--;
	  numsets++;
	}

      for (int i=0; i < num; i++)
	{
	  int node;
	  int dex = -1;
	  Read (in, node);
	  if (!data.ColumnHasValue (0, node, dex)) 
	    {
	      cout << "\nAbaqusT::ReadNodeSets unable to find node in fCoordinateData\n" << node << endl;
	      return;
	    }
	  if (dex > -1) data (dex, 1) = numsets;
	}
      search = false;
    }
}

void AbaqusT::ReadElement (istream& in, int& ID, GeometryT::GeometryCode& geocode, iArrayT& nodes)
{
  int length, numelemnodes;
  if (!AdvanceStreamTo (in, aElement, length))
    {
      cout << "\n\n AbaqusT unable to find element record.\n\n";
      return;
    }

  Read (in, ID);

  StringT name;
  Read (in, name);
  if (!TranslateElementName (name, geocode, numelemnodes))
    {
      cout << "\n\n AbaqusT::ReadElement unable to translate element type.\n";
      cout << name << endl;
      throw eGeneralFail;
    }

  nodes.Allocate (numelemnodes);
  for (int i=0; i < numelemnodes; i++)
    Read (in, nodes[i]);
}

void AbaqusT::ReadCoordinate (istream& in, int& ID, dArrayT& x)
{
  int length;
  if (!AdvanceStreamTo (in, aCoordinates, length))
    {
      cout << "\n\n AbaqusT::ReadCoordinate unable to fine coordinate record.\n\n";
      return;
    }

  Read (in, ID);

  x.Allocate (length - 3);
  for (int i=0; i < x.Length(); i++)
    Read (in, x[i]);
}

bool AbaqusT::NextTimeIncrement (istream& in, double& time)
{
  int length;
  AbaqusT::KeyRecord key;
  if (!AdvanceStreamTo (in, aIncrement, aModal, length, key, true))
    {
      cout << "\n\nAbaqusT::NextTimeIncrement, unable to find an increment or mode" << endl;
      return false;
    }

  if (key == aIncrement)
    {
      int clear = length - 2;
      double steptime, creep, amplitude, factor, freq, timeinc;
      int procedure, step, inc, perturb;
      StringT attribute;

      Read (in, time);
      Read (in, steptime);
      Read (in, creep);
      Read (in, amplitude);
      Read (in, procedure);
      Read (in, step);
      Read (in, inc);
      Read (in, perturb);
      Read (in, factor);
      Read (in, freq);
      Read (in, timeinc);
      clear -= 11;

      Read (in, attribute, clear);

      if (fModal)
	if (!AdvanceStreamTo (in, aModal, length)) 
	  {
	    cout << "AbaqusT::NextTimeIncrement, unable to find next mode" << endl;
	    return false;
	  }
	else key = aModal;
    }

  if (key == aModal)
    {
      cout << "\nFound Modal Data " << endl;
      int eigvalnumber;
      double dummy;
      Read (in, eigvalnumber);
      time = double (eigvalnumber);
      for (int i=0; i < length - 3; i++)
	Read (in, dummy);
    }

  return true;
}

void AbaqusT::ReadLabels (istream& in, ArrayT<AbaqusT::VariableKey>& nkeys, ArrayT<AbaqusT::VariableKey>& ekeys, bool& energydata)
{
  AbaqusT::KeyRecord key;
  AbaqusT::OutputFlag flag = aNoFlag;
  int length, clear = 0;
  while (NextRecord (in, key, length, clear))
    {
      clear = length - 2;
      switch (key)
	{
	case aOutputDef:
	  {
	    Read (in, flag);
	    clear--;
	    break;
	  }
	case aEndIncrement:
	case aModal:
	  return; // read only one time step
	case aElementHead:
	  {
	    int number, intpoint, section, location;
	    Read (in, number);
	    Read (in, intpoint);
	    Read (in, section);
	    Read (in, location);
	    clear -= 4;

	    switch (location)
	      {
	      case aQuadrature: // integration point
	      case aAtNodes: // data at element nodes
	      case aRebar: // rebar data
		break;
	      case aCentroidal:
	      case aWholeElement:
		flag = aElementFlag;
		break;
	      case aNodeAveraged:
		flag = aNodalFlag;
		break;
	      }
	    break;
	  }
	case aEnergy:
	  energydata = true;
	  break;
	default:
	  {
	    if (flag == aElementFlag || flag == aNodalFlag)
	      {
		int nodenumber;
		ReadPreliminary (in, key, nodenumber, clear);

		ArrayT<AbaqusT::VariableKey>* save;
		if (flag == aNodalFlag)
		  save = &nkeys;
		else
		  save = &ekeys;

		bool append = true;
		AbaqusT::VariableKey temp = Record2VariableKey (key);
		for (int i=0; i < save->Length() && append; i++)
		  if ((*save)[i] == temp) append = false;
		if (append)
		  {
		    int len = save->Length();
		    if (len == 0)
		      {
			save->Allocate (clear);
			(*save) = temp;
		      }
		    else
		      save->Resize (len + clear, temp, true, false);
		  }
	      }
	  }
	}
    }
}

void AbaqusT::ReadNodeVariables (istream& in, const ArrayT<AbaqusT::VariableKey>& nkeys, dArray2DT& values)
{
  int length, clear = 0;
  AbaqusT::VariableKey varkey = aUnknown;
  AbaqusT::KeyRecord nextkey;

  iArrayT row (nkeys.Length());
  row = 0;
  while (NextRecord (in, nextkey, length, clear))
    {
      clear = length - 2;
      varkey = aUnknown;
      
      if (nextkey == aModal || nextkey == aIncrement ||
	  nextkey == aEndIncrement) return;

      else
	{
	  varkey = Record2VariableKey (nextkey);
	  int column = -1;
	  for (int i=0; i < nkeys.Length() && column < 0; i++)
	    if (varkey == nkeys[i]) column = i;

	  if (column > -1)
	    {
	      int nodenumber;
	      ReadPreliminary (in, varkey, nodenumber, clear);
		  
	      double *pvalues = values(row[column]++) + column;
	      for (int k=0; k < clear; k++)
		Read (in, *pvalues++);
	      clear = 0;
	    }
	}
    }
}

void AbaqusT::ReadElementVariables (istream& in, const ArrayT<AbaqusT::VariableKey>& ekeys, dArray2DT& values, const iArrayT& elements)
{
  int length, clear = 0;
  AbaqusT::VariableKey varkey = aUnknown;
  AbaqusT::KeyRecord nextkey;

  int currentelement = -1;
  while (NextRecord (in, nextkey, length, clear))
    {
      clear = length - 2;
      varkey = aUnknown;

      if (nextkey == aModal || nextkey == aIncrement ||
	  nextkey == aEndIncrement) return;

      else if (nextkey == aElementHead)
	{
	  Read (in, currentelement);
	  clear--;
	}

      else
	{
	  int row;
	  if (elements.HasValue (currentelement, row))
	    {
	      varkey = Record2VariableKey (nextkey);
	      int column = -1;
	      for (int i=0; i < ekeys.Length() && column < 0; i++)
		if (varkey == ekeys[i]) column = i;
	      
	      if (column > -1)
		{
		  //int nodenumber; //only needed for NFLUX el var, I think
		  //ReadPreliminary (in, varkey, nodenumber, clear);
		  
		  double *pvalues = values (row) + column;
		  for (int k=0; k < clear; k++)
		    Read (in, *pvalues++);
		  clear = 0;
		}
	    }
	}
    }
}

void AbaqusT::Create (ostream& out, int numelems, int numnodes, double elemsize)
{
  time_t now;
  time(&now);
  char date[40], time[20];
  strftime(time, 40, "%X", localtime(&now));
  strftime(date, 20, "%d-%b-%Y", localtime(&now));
  
  StartRecord (out, 9, aVersion);
  StringT s = kAbaqusVersion;
  Write (out, s);
  s = date;
  Write (out, s, 2);
  s = time;
  Write (out, s);
  out.flush();
  Write (out, numelems);
  Write (out, numnodes);
  Write (out, elemsize);
}
void AbaqusT::WriteConnectivity (ostream& out, GeometryT::GeometryCode code, int startnumber, const iArray2DT& connects)
{
  StringT name;
  int numelemnodes;
  GetElementName (code, connects.MinorDim(), numelemnodes, name);

  int length = 4 + numelemnodes;
  for (int i=0; i < connects.MajorDim(); i++)
    {
      StartRecord (out, length, aElement);
      Write (out, i+1);
      Write (out, name);
      for (int j=0; j < numelemnodes; j++)
	Write (out, connects (i,j));
    }
}

void AbaqusT::WriteCoordinates (ostream& out, const iArrayT& nodes_used, const dArray2DT& coords)
{
  int dof = coords.MinorDim();
  int *pn = nodes_used.Pointer();
  for (int i=0; i < nodes_used.Length(); i++)
    {
      StartRecord (out, 3 + dof, aCoordinates);
      Write (out, *pn);
      int row = *pn++;
      double *pc = coords (row - 1);
      for (int j=0; j < dof; j++)
	Write (out, *pc++);
    }
}

void AbaqusT::WriteElementSet (ostream& out, const StringT& name, const iArrayT& elem_map)
{
  AbaqusT::KeyRecord key = aElementSet;
  int wrap = 80 - 3;
  int *pe = elem_map.Pointer();
  for (int i=0; i < elem_map.Length(); i++)
    {
      if (i%wrap == 0)
	{
	  int length = wrap;
	  if (elem_map.Length() - i < wrap) length = elem_map.Length() - i;
	  if (i > 0) key = aElemSetCont;
	  StartRecord (out, length + 3, key);
	  if (key == aElementSet) Write (out, name);
	}
      Write (out, *pe++);
    }
}

void AbaqusT::WriteNodeSet (ostream& out, const StringT& name, const iArrayT& nodes)
{
  AbaqusT::KeyRecord key = aNodeSet;
  int wrap = 80 - 3;
  int *pn = nodes.Pointer();
  for (int i=0; i < nodes.Length(); i++)
    {
      if (i%wrap == 0)
	{
	  int length = wrap;
	  if (nodes.Length() - i < wrap) length = nodes.Length() - i;
	  if (i > 0) key = aNodeSetCont;
	  StartRecord (out, length + 3, key);
	  if (key == aNodeSet) Write (out, name);
	}
      Write (out, *pn++);
    }
}

void AbaqusT::WriteActiveDOF (ostream& out, const iArrayT& active)
{
  int length = active.Length() + 2;
  StartRecord (out, length, aDOF);
  for (int i=0; i < active.Length(); i++)
    Write (out, active[i]);
}

void AbaqusT::WriteHeading (ostream& out, const StringT& heading)
{
  if (heading.Length() <= 1) return;
  StartRecord (out, 12, aHeading);
  Write (out, heading, 10);
}

void AbaqusT::WriteEndIncrement (ostream& out, bool endgeometry)
{
  StartRecord (out, 2, aEndIncrement);
  if (fBinary)
    {
      /*if (!endgeometry)
	{
	  StringT space (sizeof (double) + 1);
	  for (int i=0; i < space.Length(); i++)
	    space[i] = ' ';
	  space [space.Length() - 1] = '\0';
	  
	  int size = (fBufferSize - fBufferDone) / sizeof (double);
	  for (int j=0; j < size; j++)
	    Write (out, space);
	    }*/
    }
  else
    {
      iArrayT size (2);
      size[0] = fBufferSize - fBufferDone + 1;
      size[1] = fBufferSize + 1;
      for (int j=0; j < 2; j++)
	{
	  StringT space (size[j]);
	  for (int i=0; i < space.Length(); i++)
	    space[i] = ' ';
	  space [space.Length() - 1] = '\0';
	  WriteASCII (out, space);
	}
    }

  /*if (fBufferSize != fBufferDone && !endgeometry)
    {
      cout << "\n\nAbaqusT::EndIncrement: buffer size = "
	   << fBufferSize << " buffer written = " << fBufferDone << endl;
      throw eGeneralFail;
      }*/
}

void AbaqusT::WriteStartIncrement (ostream& out, int step, int inc, double totaltime, double time, double timeincrement, AbaqusT::Analysis analysistype)
{
  double creep  = 0, amplitude  = 0, factor  = 0, freq  = 0;
  int perturb  = 0;
  StartRecord (out, 23, aIncrement);
  Write (out, totaltime);
  Write (out, time);
  Write (out, creep);
  Write (out, amplitude);
  Write (out, analysistype);
  Write (out, step);
  Write (out, inc);
  Write (out, perturb);
  Write (out, factor);
  Write (out, freq);
  Write (out, timeincrement);
  StringT attribute ("default load case");
  Write (out, attribute, 10);
}

void AbaqusT::WriteNodalData (ostream& out, const ArrayT<AbaqusT::VariableKey>& key, const dArray2DT& values, const iArrayT& nodemap, GeometryT::GeometryCode code, int numnodes)
{
  StringT setname = "";
  int numdir = 0, numshear = 0;
  for (int i=0; i < key.Length(); i++)
    {
      // determine record length
      int length = 2;
      int count = 0;
      VariableRecordLength(key, i, length, count);

      // write output definition, account for element average nodal data
      WriteOutputDef (out, key[i], setname, code, numnodes, count, numdir, numshear);

      // write data records for that variable
      for (int k=0; k < nodemap.Length(); k++)
	{
	  // account for element averaged nodal data, determine component types
	  WriteElementHeader (out, key[i], nodemap[k], 0, 0, aNodeAveraged, numdir, numshear, 0, 0);
	    
	  StartRecord (out, length, key[i]);

	  // account for node number
	  WritePreliminary (out, key[i], nodemap[k]);

	  // write data
	  for (int m=i; m < i + count; m++)
	    if (key[m] == key[i])
	      Write (out, values (k, m));
	}
      
	  // adjust index
	  i += count - 1;
    }
}

void AbaqusT::WriteElementData (ostream& out, const ArrayT<AbaqusT::VariableKey>& key, const dArray2DT& values)
{
  cout << "element variables not supported" << endl;
  for (int i=0; i < key.Length(); i++)
    cout << key[i] << endl;
}

/*************************************************************************
 * Private
 *************************************************************************/

bool AbaqusT::AdvanceStreamTo (istream& in, AbaqusT::KeyRecord key, int& length)
{
  while (!in.eof())
    {
      if (!AdvanceToStar (in)) return false;
      
      // read
      int dkey, dlength;
      Read (in, dlength);
      Read (in, dkey);

      // compare
      if (dkey == key)
	{
	  length = dlength;
	  return true;
	}
      else
	ClearLine (in, dlength - 2);
    }
  return false;
}

bool AbaqusT::AdvanceStreamTo (istream& in, AbaqusT::KeyRecord key, AbaqusT::KeyRecord keycont, int& length, AbaqusT::KeyRecord& keyfound, bool search)
{
  while (!in.eof())
    {
      if (!AdvanceToStar (in)) return false;

      // read
      AbaqusT::KeyRecord dkey;
      int dlength;
      Read (in, dlength);
      Read (in, dkey);
      cout << "2 " << dlength << " " << dkey << endl;

      // compare
      if (dkey == key || dkey == keycont)
	{
	  length = dlength;
	  keyfound = dkey;
	  return true;
	}
      else if (search == false) return false; // shortcut
      else
	ClearLine (in, dlength - 2);
    }
  return false;
}

bool AbaqusT::NextRecord (istream& in, AbaqusT::KeyRecord& key, int& length, int clear)
{
  if (!AdvanceToStar (in)) return false;
  ClearLine (in, clear);
  
  // read
  Read (in, length);
  Read (in, key);
  return true;
}

bool AbaqusT::AdvanceToStar (istream& in) const
{
  if (fBinary) return true;

  // advance
  char c;
  in.get (c);
  while (c != '*') 
    {
      in.get (c);
      if (in.eof()) 
	{
	  cout << "\n\nAbaqusT::AdvanceToStar, unable to find start of next record" << endl;
	  return false;
	}
    }
  return true;
}

void AbaqusT::ClearLine (istream& in, int n)
{
  if (!fBinary) return;
  if (n < 1) return;
  StringT temp;
  Read (in, temp, n);
}

void AbaqusT::Read (istream& in, StringT& s, int n)
{
  char c;
  ArrayT<char> temp (n * sizeof (double) + 1);
  char *ps = temp.Pointer();
  for (int i=0; i < n; i++)
    {
      if (fBinary)
	{
	  CheckBufferSize (in);
	  in.read (ps, sizeof (double));
	  if (in.eof ()) return;
	  fBufferDone += sizeof (double);
	  ps += sizeof (double);
	}
      else
	{
	  in.get (c);
	  if (c == '\n') in.get(c);
	  if (c != 'A') return;
	  for (int j=0; j < sizeof (double); j++)
	    {
	      in.get (c);
	      if (c != '\n') *ps++ = c;
	      else j--;
	    }
	}
    }
  *ps = '\0';
  s.Append (temp.Pointer());
}

void AbaqusT::Read (istream& in, int& i)
{
  if (fBinary)
    {
      CheckBufferSize (in);
      int temp;
      if (in.eof()) return;
      in.read (reinterpret_cast<char *> (&temp), sizeof (double));
      i = temp;
      fBufferDone += sizeof (double);
    }
  else
    {
      char c;
      in.get(c);
      if (c == '\n') in.get(c);
      if (c != 'I') return;
      
      ArrayT<char> wtemp (3);
      char *pw = wtemp.Pointer();
      for (int w=0; w < 2; w++)
	{
	  in.get (c);
	  if (c != '\n') *pw++ = c;
	  else w--;
	}
      *pw = '\0';
      int width = atof (wtemp.Pointer());
      
      ArrayT<char> itemp (width+1);
      char *pi = itemp.Pointer();
      for (int j=0; j < width; j++)
	{
	  in.get (c);
	  if (c != '\n') *pi++ = c;
	  else j--;
	}
      *pi = '\0';
      i = atof (itemp.Pointer());
    }
}

void AbaqusT::Read (istream& in, AbaqusT::KeyRecord& k)
{
  int temp;
  Read (in, temp);
  k = (AbaqusT::KeyRecord) temp;
}

void AbaqusT::Read (istream& in, AbaqusT::OutputFlag& k)
{
  int temp;
  Read (in, temp);
  k = (AbaqusT::OutputFlag) temp;
}

void AbaqusT::Read (istream& in, double& d)
{
  if (fBinary)
    {
      CheckBufferSize (in);
      double temp;
      in.read (reinterpret_cast<char *> (&temp), sizeof (double));
      d = temp;
      fBufferDone += sizeof (double);
    }
  else
    {
      char c;
      in.get(c);
      if (c == '\n') in.get(c);
      if (c != 'D') return;
      
      ArrayT<char> temp (dprecision + 8);
      char *pt = temp.Pointer();
      for (int t=0; t < dprecision + 3; t++)
	{
	  in.get (c);
	  if (c != '\n') *pt++ = c;
	  else t--;
	}
      
      in.get (c);
      if (c == '\n') in.get(c);
      if (c != 'D' && c != 'E') return;
      *pt++ = 'E';
      
      for (int i=0; i < 3; i++)
	{
	  in.get (c);
	  if (c != '\n') *pt++ = c;
	  else i--;
	}
      *pt = '\0';
      d = atof (temp.Pointer());
    }
}

void AbaqusT::CheckBufferSize (istream& in)
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

bool AbaqusT::TranslateElementName (const StringT& name, GeometryT::GeometryCode& code, int& numelemnodes) const
{
  // continuum elements
  if (strncmp (name.Pointer(), "C", 1) == 0)
    if (TranslateContinuum (name.Pointer(1), code, numelemnodes)) return true;
    
  if (strncmp (name.Pointer(), "DCC", 3) == 0)
    if (TranslateContinuum (name.Pointer(3), code, numelemnodes)) return true;

  if (strncmp (name.Pointer(), "DC", 2) == 0 ||
      strncmp (name.Pointer(), "AC", 2) == 0)
    if (TranslateContinuum (name.Pointer(2), code, numelemnodes)) return true;

  // membrane elements
  if (strncmp (name.Pointer(), "M3D", 3) == 0)
    if (Translate2DName (name.Pointer(3), code, numelemnodes)) return true;

  // mass elements
  if (strncmp (name.Pointer(), "MASS", 4) == 0)
    {
      code = GeometryT::kPoint;
      numelemnodes = 1;
      return true;
    }
    

  cout << "\nAbaqusT does not support element type: " << name << endl;
  return false;
}

bool AbaqusT::TranslateContinuum (const char* name, GeometryT::GeometryCode& code, int& numelemnodes) const
{
  // 2D continuum elements
  if (strncmp (name, "PE", 2) == 0 ||
      strncmp (name, "PS", 2) == 0 ||
      strncmp (name, "2D", 2) == 0 ||
      strncmp (name, "AX", 2) == 0)
    if (Translate2DName (name += 2, code, numelemnodes)) return true;
    
  if (strncmp (name, "GPE", 3) == 0 ||
      strncmp (name, "GAX", 3) == 0 ||
      strncmp (name, "AXA", 3) == 0)
    if (Translate2DName (name += 3, code, numelemnodes)) return true;

  if (strncmp (name, "3D", 2) == 0)
    if (Translate3DName (name += 2, code, numelemnodes)) return true;
    
  return false;
}

bool AbaqusT::Translate2DName (const char* n, GeometryT::GeometryCode& code, int& numelemnodes) const
{
  numelemnodes = atoi (n);
  switch (numelemnodes)
    {
    case 3: 
    case 5:
    case 6:
      code = GeometryT::kTriangle;
      break;
    case 4: 
    case 8:
    case 9:
    case 10:
      code = GeometryT::kQuadrilateral;
      break;
    default:
      return false;
    }
  return true;
}

bool AbaqusT::Translate3DName (const char* n, GeometryT::GeometryCode& code, int& numelemnodes) const
{
  numelemnodes = atoi (n);
  switch (numelemnodes)
    {
    case 4: 
    case 10:
      code = GeometryT::kTetrahedron;
      break;
    case 8: 
    case 20:
    case 21:
    case 22:
    case 23:
    case 24:
    case 25:
    case 26:
    case 27:
      code = GeometryT::kHexahedron;
      break;
    case 6:
    case 15:
    case 16:
    case 17:
    case 18:
      code = GeometryT::kPentahedron;
    default:
      return false;
    }
  return true;
}

void AbaqusT::CheckBufferSize (ostream& out)
{
  if (!fBinary) return;
  
  if (fBufferDone == fBufferSize)
    {
      if (fBinary)
	// write buffer size
	out.write (reinterpret_cast<char *> (&fBufferSize), sizeof (int)); 
      fBufferDone = 0;
    }

  if (fBufferDone == 0)
    {
      if (fBinary)
	// write buffer size
	out.write (reinterpret_cast<char *> (&fBufferSize), sizeof (int)); 
    }
}

void AbaqusT::StartRecord (ostream& out, int length, int key)
{
  if (!fBinary)
    {
      StringT start ("*");
      WriteASCII (out, start);
    }
  Write (out, length);
  Write (out, key);
}

void AbaqusT::WriteASCII (ostream& out, const StringT& s)
{
  if (fBinary) return;

  if (fBufferDone + s.Length() < fBufferSize)
    {
      out << s;
      fBufferDone += s.Length() - 1;
    }
  else
    for (int i=0; i < s.Length() - 1; i++)
      {
	if (fBufferDone == fBufferSize) 
	  {
	    out << "\n";
	    fBufferDone = 0;
	  }
	out << s[i];
	fBufferDone++;
      }
}

void AbaqusT::Write (ostream& out, int i)
{
  if (fBinary)
    {
      CheckBufferSize (out);
      out.write (reinterpret_cast<char *> (&i), sizeof (double));
      fBufferDone += sizeof (double);
    }
  else
    {
      StringT s ("I ");
      StringT ints;
      ints.Append (i);
      s.Append (ints.Length()-1);
      s.Append (ints);
      WriteASCII (out, s);
    }
}

void AbaqusT::Write (ostream& out, double d)
{
  if (fBinary)
    {
      CheckBufferSize (out);
      out.write (reinterpret_cast<char *> (&d), sizeof (double));
      fBufferDone += sizeof (double);
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

      WriteASCII (out, s);
    }
}

void AbaqusT::Write (ostream& out, const StringT& s, int blocks)
{
  char *ps = s.Pointer();
  if (fBinary)
    {
      CheckBufferSize (out);
      out.write (ps, sizeof (double)*blocks);
      fBufferDone += sizeof (double)*blocks;
    }
  else
    {
      for (int i=0, j=0; i < blocks; i++, j+= 8)
	{
	  StringT w = "A";
	  int copylength = 8;
	  if (j+8 > s.Length() - 1) copylength = s.Length() - j - 1;
	  if (copylength < 0) copylength = 0;
	  for (int k=j; k < j+copylength; k++)
	    w.Append (s[k]);
	  for (int m=copylength; m < 8; m++)
	    w.Append (' ');
	  WriteASCII (out, w);
	}
    }
}

void AbaqusT::GetElementName (GeometryT::GeometryCode geometry_code, int elemnodes, int& num_output_nodes, StringT& elem_name) const
{
  switch (geometry_code)
    {
    case GeometryT::kPoint:
      elem_name = "MASS";	      
      num_output_nodes = 1;
      break;
      
    case GeometryT::kTriangle:
      elem_name =  "CPE";
      num_output_nodes = (elemnodes < 6) ? 3 : 6;
      elem_name.Append (num_output_nodes);
      break;
      
    case GeometryT::kQuadrilateral:
      elem_name =  "CPE";
      num_output_nodes = (elemnodes < 8) ? 4 : 8;
      elem_name.Append (num_output_nodes);
      break;
      
    case GeometryT::kHexahedron:
      elem_name = "C3D";
      num_output_nodes = (elemnodes < 20) ? 8 : 20;
      elem_name.Append (num_output_nodes);
      break;

    case GeometryT::kTetrahedron:
      elem_name =  "C3D";
      num_output_nodes = (elemnodes < 10) ? 4 : 10;
      elem_name.Append (num_output_nodes);
      break;

    case GeometryT::kPentahedron:
      elem_name = "C3D";
      num_output_nodes = (elemnodes < 15) ? 6 : 15;
      elem_name.Append (num_output_nodes);
      break;

    default:
      cout << "\nAbaqusT::GetElementName cannot find name\n\n";
      throw eGeneralFail;
    }
}

void AbaqusT::WriteOutputDef (ostream& out, int key, const StringT& setname, GeometryT::GeometryCode code, int numnodes, int count, int& numdir, int& numshear)
{
  AbaqusT::OutputFlag flag = aNodalFlag;
  GeometryT::GeometryCode geocode = GeometryT::kNone;
  int numelemnodes = 0;

  if (key == aSDV || key == aStress ||
      key == aTotalStrain || key == aPrinStress ||
      key == aUVARM)
    {
      flag = aElementFlag;
      geocode = code;
      numelemnodes = numnodes;
      if (key == aStress || key == aTotalStrain)
	switch (count)
	  {
	  case 3:
	    numdir = 2;
	    numshear = 1;
	    break;
	  case 6:
	    numdir = 3;
	    numshear = 3;
	    break;
	  default:
	    {
	      cout << "AbaqusT::WriteOutputDef, not programmed for " 
		   << count << " number of nodal averaged components\n";
	      throw eGeneralFail;
	    }
	  }
      else
	{
	  numdir = count;
	  numshear = 0;
	}
    }
      
  StartRecord (out, 5, aOutputDef);
  Write (out, flag);
  Write (out, setname);

  if (geocode != GeometryT::kNone)
    {
      StringT elemname;
      int num_output_nodes;
      GetElementName (geocode, numelemnodes, num_output_nodes, elemname);
      Write (out, elemname);
    }
  else
    {
      int none = 0;
      Write (out, none);
    }
}

void AbaqusT::WriteElementHeader (ostream& out, int key, int number, int integrationpoint, int sectionpoint, AbaqusT::ElementData flag, int numdirect, int numshear, int numdir, int numsecforc)
{
  if (key == aSDV || key == aStress ||
      key == aTotalStrain || key == aPrinStress ||
      key == aUVARM ||
      flag == aElementFlag)
    {
      int rebarname = 0;
      StartRecord (out, 11, aElementHead);
      Write (out, number);
      Write (out, integrationpoint);
      Write (out, sectionpoint);
      Write (out, flag);
      Write (out, rebarname);
      Write (out, numdirect);
      Write (out, numshear);
      Write (out, numdir);
      Write (out, numsecforc);
    }
}

void AbaqusT::ReadPreliminary (istream& in, int key, int& i, int& clear)
{
  if (key == aDisplacement || key == aVelocity ||
      key == aAcceleration || key == aCoordVariable)
    {
      Read (in, i);
      clear--;
    }
}

void AbaqusT::WritePreliminary (ostream& out, int key, int i)
{
  if (key == aDisplacement || key == aVelocity ||
      key == aAcceleration || key == aCoordVariable)
    Write (out, i);
}

void AbaqusT::VariableRecordLength (const ArrayT<AbaqusT::VariableKey>& key, int index, int& length, int& count)
{
      for (int j=index; j < key.Length(); j++)
	if (key[j] == key[index]) count++;
      length += count;

      // account for node number
      if (key[index] == aDisplacement || key[index] == aVelocity ||
	  key[index] == aAcceleration || key[index] == aCoordVariable)
	length++;
}

AbaqusT::VariableKey AbaqusT::Record2VariableKey (AbaqusT::KeyRecord key) const
{
  switch (key)
    {
    case aSDV: return aSDV;
    case aStress: return aStress;
    case aTotalStrain: return aTotalStrain;
    case aDisplacement: return aDisplacement;
    case aVelocity: return aVelocity;
    case aAcceleration: return aAcceleration;
    case aCoordVariable: return aCoordVariable;
    case aPrinStress: return aPrinStress;
    case aUVARM: return aUVARM;
    }
  return aUnknown;
}
