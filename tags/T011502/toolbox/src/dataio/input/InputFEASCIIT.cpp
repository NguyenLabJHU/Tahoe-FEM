/* $Id: InputFEASCIIT.cpp,v 1.7 2002-01-07 03:06:02 paklein Exp $ */
#include "InputFEASCIIT.h"
#include "ifstreamT.h"
#include "dArrayT.h"
#ifdef _MSC_VER
#include <strstrea.h>
#else
#include <strstream.h>
#endif

InputFEASCIIT::InputFEASCIIT (ostream& out) :
  InputBaseT (out),
  fFileRoot (""),
  fBlockID (0),
  fNumNodes (0),
  fNumElements (0),
  fNumDOF (0),
  fTimeSteps (0),
  fNodeVariable (0),
  fElementVariable (0)
{
}

bool InputFEASCIIT::Open (const StringT& filename)
{
  /* create file root */
  StringT suffix;
  suffix.Suffix (filename.Pointer());
  if (strncmp (suffix.Pointer(), ".geo", 4) == 0 ||
      strncmp (suffix.Pointer(), ".run", 4) == 0 ||
      strncmp (suffix.Pointer(), ".in", 3) == 0)
  fFileRoot.Root (filename);

	/* scan geometry file */
	ifstreamT geo;
	if (!OpenFile (geo, ".geo")) {
		cout << "\n InputFEASCIIT::Open: error opening geometry file: " << geo.filename() << endl;
		return false;
	}
	if (!ScanGeometryFile (geo)) {
		cout << "\n InputFEASCIIT::Open: error scanning geometry file: " << geo.filename() << endl;
		return false;
	}

	/* scan results file */
	ifstreamT run;
	if (!OpenFile (run, ".run")) {
		cout << "\n InputFEASCIIT::Open: error opening results file: " << run.filename() << endl;
		return false;
	}
    if (!ScanResultsFile (run)) {
		cout << "\n InputFEASCIIT::Open: error scanning results file: " << run.filename() << endl;
		return false;
    }
      
	/* must be OK */
	return true;
}

void InputFEASCIIT::Close (void)
{
  fFileRoot.Free();
  fBlockID.Free();
  fNumNodes = 0;
  fNumElements = 0;
  fNumDOF = 0;
  fTimeSteps.Free();
  fNodeVariable.Free();
  fElementVariable.Free();
}

void InputFEASCIIT::ElementGroupNames (ArrayT<StringT>& groupnames) const
{
  if (groupnames.Length() != fBlockID.Length()) throw eSizeMismatch;
  for (int i=0; i < groupnames.Length(); i++)
    {
      groupnames[i].Clear();
      groupnames[i].Append(fBlockID[i]);
    }
}

void InputFEASCIIT::ReadNodeMap (iArrayT& nodemap)
{
  nodemap.SetValueToPosition();
  nodemap++;
}

void InputFEASCIIT::ReadCoordinates (dArray2DT& coords)
{
  iArrayT nodemap (fNumNodes);
  ReadCoordinates (coords, nodemap);
}

void InputFEASCIIT::ReadCoordinates (dArray2DT& coords, iArrayT& nodemap)
{
  if (nodemap.Length() != fNumNodes ||
      coords.MajorDim() != fNumNodes ||
      coords.MinorDim() != fNumDOF ) throw eSizeMismatch;

  ifstreamT geo;
  OpenFile (geo, ".geo");
  
  StringT s;
  int count = 0;
  iArrayT used (fNodeVariable.Length());
  nodemap.SetValueToPosition();
  nodemap++;
  coords = 0;
  for (int i=0; i < fBlockID.Length(); i++)
    {
      if (!geo.FindString ("Nodal coordinates", s)) throw eDatabaseFail;

      dArray2DT c;
      iArrayT n;
      DataBlock (geo, used, n, c, true);
      n--;
 
      int *pn = n.Pointer();
      for (int i=0; i < c.MajorDim(); i++)
	coords.SetRow (*pn++, c(i));
    }
}

int InputFEASCIIT::NumElements (StringT& name)
{
  ifstreamT geo;
  OpenFile (geo, ".geo");

  if (!AdvanceToBlock (geo, name, "Connectivities")) return -1;

  int numels;
  StringT s;
  if (!geo.FindString ("Number of elements", s) ||
      !s.Tail ('=', numels)) throw eDatabaseFail;

  return numels;
}

int InputFEASCIIT::NumElementNodes (StringT& name)
{
  ifstreamT geo;
  OpenFile (geo, ".geo");

  if (!AdvanceToBlock (geo, name, "Connectivities")) return -1;

  int num;
  StringT s;
  if (!geo.FindString ("Number of element nodes", s) ||
      !s.Tail ('=', num)) throw eDatabaseFail;

  return num;
}

void InputFEASCIIT::ReadAllElementMap (iArrayT& elemmap)
{
  if (elemmap.Length() != fNumElements) throw eSizeMismatch;

  ifstreamT geo;
  OpenFile (geo, ".geo");

  StringT s;
  int numelms;
  int count = 0;
  char line [255];
  for (int i=0; i < fBlockID.Length(); i++)
    {
      if (!geo.FindString ("Connectivities", s) ||
	  !geo.FindString ("Number of elements", s) ||
	  !s.Tail ('=', numelms) ||
	  !geo.FindString ("element", s)) throw eDatabaseFail;

      for (int i=0; i < numelms; i++)
	{
	  geo >> elemmap[count++];
	  geo.getline (line, 254);
	}
    }
}

void InputFEASCIIT::ReadGlobalElementMap (StringT& name, iArrayT& elemmap)
{
  ifstreamT geo;
  OpenFile (geo, ".geo");

  StringT s;
  int numelms;
  char line [255];
  if (!AdvanceToBlock (geo, name, "Connectivities") ||
      !geo.FindString ("Number of elements", s) ||
      !s.Tail ('=', numelms)) throw eDatabaseFail;

  if (elemmap.Length() != numelms) throw eSizeMismatch;

  for (int i=0; i < numelms; i++)
    {
      geo >> elemmap[i];
      geo.getline (line, 254);
    }
}

void InputFEASCIIT::ReadGlobalElementSet (StringT& name, iArrayT& set)
{
  if (set.Length() != fNumElements) throw eSizeMismatch;

  ifstreamT geo;
  OpenFile (geo, ".geo");

  StringT s;
  int numelms;
  int count = 0;
  const int ID = atoi (name.Pointer());
  int found = -1;
  while (found != ID)
    {
      if (!geo.FindString ("Connectivities", s) ||
	  !geo.FindString ("Number of elements", s) ||
	  !s.Tail ('=', numelms) ||
	  !geo.FindString ("element", s)) throw eDatabaseFail;

      count += numelms;
    }
  
  if (set.Length() != numelms) throw eSizeMismatch;
  set.SetValueToPosition();
  set += count;
}

void InputFEASCIIT::ReadConnectivity (StringT& name, iArray2DT& connects)
{
  if (connects.Length() == 0) throw eSizeMismatch;

  ifstreamT geo;
  OpenFile (geo, ".geo");

  if (!AdvanceToBlock (geo, name, "Connectivities")) throw eDatabaseFail;

  StringT s;
  int numelms, numelnodes, elmid;
  iArrayT elms (connects.MinorDim());
  if (!geo.FindString ("Number of elements", s) ||
      !s.Tail ('=', numelms) ||
      !geo.FindString ("Number of element nodes", s) ||
      !s.Tail ('=', numelnodes) ||
      !geo.FindString ("element", s)) throw eDatabaseFail;

  if (numelnodes != connects.MinorDim()) throw eSizeMismatch;

  for (int i=0; i < numelms; i++)
    {
      geo >> elmid >> elms;
      connects.SetRow (i, elms);
    }
  
  connects--;
}

void InputFEASCIIT::ReadGeometryCode (StringT& name, GeometryT::CodeT& geocode)
{
  ifstreamT geo;
  OpenFile (geo, ".geo");

  if (!AdvanceToBlock (geo, name, "Connectivities")) throw eDatabaseFail;

  StringT s;
  if (!geo.FindString ("Geometry code", s)) throw eDatabaseFail;

  // either this or write an operator= for Geometry::CodeT
#ifdef _MSC_VER
  char *h = strstr ((char*) s, "=");
#else
  const char *h = strstr ((const char*) s, "=");
#endif
  istrstream istr (h+1);
  istr >> geocode;
}

void InputFEASCIIT::ReadTimeSteps (dArrayT& steps) 
{
  if (steps.Length() != fTimeSteps.Length()) throw eSizeMismatch;
  steps.Set (steps.Length(), fTimeSteps.Pointer());
}

void InputFEASCIIT::ReadNodeLabels (ArrayT<StringT>& nlabels) const
{
  if (nlabels.Length() != NumNodeVariables()) throw eSizeMismatch;
  for (int i=0; i < nlabels.Length(); i++)
    nlabels[i] = fNodeVariable[i];
}

void InputFEASCIIT::ReadElementLabels (ArrayT<StringT>& elabels) const
{
  if (elabels.Length() != NumElementVariables()) throw eSizeMismatch;
  for (int i=0; i < elabels.Length(); i++)
    elabels[i] = fElementVariable[i];
}

void InputFEASCIIT::NodeVariablesUsed (StringT& name, iArrayT& used)
{ 
  if (used.Length() != fNodeVariable.Length()) throw eSizeMismatch;
  used = 0;

  ifstreamT run;
  OpenFile (run, ".run");
  
  if (!AdvanceToBlock (run, name, "Nodal data")) throw eDatabaseFail;

  dArray2DT vals;
  iArrayT ids;
  DataBlock (run, used, ids, vals, true);
}

void InputFEASCIIT::ElementVariablesUsed (StringT& name, iArrayT& used)
{ 
  if (used.Length() != fElementVariable.Length()) throw eSizeMismatch;
  used = 0;

  ifstreamT run;
  OpenFile (run, ".run");
  
  if (!AdvanceToBlock (run, name, "Element data")) throw eDatabaseFail;

  dArray2DT vals;
  iArrayT ids;
  DataBlock (run, used, ids, vals, false);
}

void InputFEASCIIT::ReadAllNodeVariables (int step, dArray2DT& nvalues)
{
#pragma unused(step)

  if (nvalues.MajorDim() != fNumNodes) throw eSizeMismatch;

  ifstreamT run;
  OpenFile (run, ".run");
  
  StringT s;
  iArrayT used (fNodeVariable.Length()), ids;
  dArray2DT vals;
  nvalues = 0;
  for (int i=0; i < fBlockID.Length(); i++)
    {
      if (!run.FindString ("Nodal data", s)) throw eDatabaseFail;

      DataBlock (run, used, ids, vals, true);

      ids--;
      for (int i=0; i < ids.Length(); i++)
	for (int v=0, j=0; v < used.Length(); v++)
	  if (used[v] > 0)
	    nvalues (ids[i], v) = vals (i, j++);
    }
}

void InputFEASCIIT::ReadNodeVariables (int step, StringT& name, dArray2DT& nvalues)
{
#pragma unused(step)

  if (nvalues.Length() == 0) throw eSizeMismatch;

  ifstreamT run;
  OpenFile (run, ".run");

  if (!AdvanceToBlock (run, name, "Nodal data")) throw eDatabaseFail;

  iArrayT used (fNodeVariable.Length()), ids;
  dArray2DT vals;
  DataBlock (run, used, ids, vals, true);

  nvalues = 0;
  for (int i=0; i < ids.Length(); i++)
    for (int v=0, j=0; v < used.Length(); v++)
      if (used[v] > 0)
	nvalues (i, v) = vals (i, j++);
}

void InputFEASCIIT::ReadAllElementVariables (int step, dArray2DT& evalues)
{
#pragma unused(step)

  if (evalues.MajorDim() != fNumElements) throw eSizeMismatch;

  ifstreamT run;
  OpenFile (run, ".run");
  
  StringT s;
  iArrayT used (fElementVariable.Length()), ids;
  dArray2DT vals;
  evalues = 0;
  for (int i=0; i < fBlockID.Length(); i++)
    {
      if (!run.FindString ("Element data", s)) throw eDatabaseFail;

      DataBlock (run, used, ids, vals, true);

      ids--;
      for (int i=0; i < ids.Length(); i++)
	for (int v=0, j=0; v < used.Length(); v++)
	  if (used[v] > 0)
	    evalues (ids[i], v) = vals (i, j++);
    }
}

void InputFEASCIIT::ReadElementVariables (int step, StringT& name, dArray2DT& evalues)
{
#pragma unused(step)

  if (evalues.Length() == 0) throw eSizeMismatch;

  ifstreamT run;
  OpenFile (run, ".run");

  if (!AdvanceToBlock (run, name, "Element data")) throw eDatabaseFail;

  iArrayT used (fElementVariable.Length()), ids;
  dArray2DT vals;
  DataBlock (run, used, ids, vals, true);

  evalues = 0;
  for (int i=0; i < ids.Length(); i++)
    for (int v=0, j=0; v < used.Length(); v++)
      if (used[v] > 0)
	evalues (i, v) = vals (i, j++);
}

/******************** PRIVATE ***************************/

bool InputFEASCIIT::OpenFile (ifstreamT& in, const char* ext) const
{
  StringT file (fFileRoot);
  file.Append (ext);
  in.open (file);
  if (!in.is_open()) 
    {
      fout << "\nInputFEASCIIT::OpenFile unable to open " << file << "\n\n";
      return false;
    }
  return true;
}

bool InputFEASCIIT::ScanGeometryFile (ifstreamT& in)
{
  StringT s;
  if (!in.FindString ("G E O M E T R Y", s)) return false;

  int nid, eid, numelems;
  fNumNodes = 0;
  fNumElements = 0;
  fNumDOF = 0;
  iAutoArrayT nodes;
  iArrayT n, used;
  dArray2DT c;
  while (in.FindString ("Nodal coordinates", s))
    {
      if (!in.FindString ("Block number", s) ||
	  !s.Tail ('=', nid)) return false;
      fBlockID.Append (nid);

      DataBlock (in, used, n, c, true);
      int dof = c.MinorDim();
      int numnodes = c.MajorDim();
      if (fBlockID.Length() == 1) fNumDOF = dof;
      else if (fNumDOF != dof) return false;
      nodes.AppendUnique (n);

      if (!in.FindString ("Connectivities", s) ||
	  !in.FindString ("Block number", s) ||
	  !s.Tail ('=', eid)) return false;
      if (eid != nid) return false;

      if (!in.FindString ("Number of elements", s) ||
	  !s.Tail ('=', numelems)) return false;
      fNumElements += numelems;
    }

  fNumNodes = nodes.Length();
  if (fBlockID.Length() < 1) 
  	return false;
  else
  	return true;
}

bool InputFEASCIIT::ScanResultsFile (ifstreamT& in)
{
  StringT s;
  if (!in.FindString ("O U T P U T", s)) return false;

  int id, numb, vals;
  double t;
  fNodeVariable.Free();
  fElementVariable.Free();
  while (in.FindString ("Group number", s))
    {
      if (!in.FindString ("Time", s) ||
	  !s.Tail ('=', t)) return false;
      fTimeSteps.AppendUnique (t);

      if (!in.FindString ("Number of Blocks", s) ||
	  !s.Tail ('=', numb)) return false;

      for (int b=0; b < numb; b++)
	{
	  if (!in.FindString ("Nodal data", s) ||
	      !in.FindString ("Block number", s) ||
	      !s.Tail ('=', id) ||
	      !fBlockID.HasValue (id) ||
	      !in.FindString ("Number of values", s) ||
	      !s.Tail ('=', vals)) return false;
	  if (vals > 0 && fTimeSteps.Length() == 1)
	    for (int v=0; v < vals + 1; v++)
	      {
		in >> s;
		if (strncmp (s.Pointer(), "node", 4) != 0)
		  fNodeVariable.AppendUnique (s);
	      }

	  if (!in.FindString ("Element data", s) ||
	      !in.FindString ("Block number", s) ||
	      !s.Tail ('=', id) ||
	      !fBlockID.HasValue (id) ||
	      !in.FindString ("Number of values", s) ||
	      !s.Tail ('=', vals)) return false;
	  if (vals > 0 && fTimeSteps.Length() == 1)
	    for (int v=0; v < vals + 1; v++)
	      {
		in >> s;
		if (strncmp (s.Pointer(), "element", 7) != 0)
		  fElementVariable.AppendUnique (s);
	      }
	}
    }
  if (fTimeSteps.Length() < 1) 
  	return false;
  else
  	return true;
}

bool InputFEASCIIT::AdvanceToBlock (ifstreamT& in, const StringT& name, const char* tname) const
{
  const int ID = atoi (name.Pointer());
  int found = -1;
  StringT s;
  while (found != ID)
    {
      if (!in.FindString (tname, s) ||
	  !in.FindString ("Block number", s) ||
	  !s.Tail ('=', found)) return false;
    }
	return true;
}

void InputFEASCIIT::DataBlock (ifstreamT& in, iArrayT& used, iArrayT& ids, dArray2DT& vals, bool nodal) const
{
  StringT t;
  ArrayT<StringT> vars;
  if (nodal)
    {
      t = "Number of nodal points";
      vars.Allocate (fNodeVariable.Length());
      ReadNodeLabels (vars);
    }
  else
    {
      t = "Number of elements";
      vars.Allocate (fElementVariable.Length());
      ReadElementLabels (vars);
    }

  StringT s;
  int numvals;
  int num;
  if (!in.FindString (t.Pointer(), s) ||
      !s.Tail ('=', num) ||
      !in.FindString ("Number of values", s) ||
      !s.Tail ('=', numvals)) throw eDatabaseFail;

  /* read labels */
  used = 0;
  if (numvals > 0)
    {
      in >>s; // read "element" or "node"
      for (int v=0; v < numvals; v++)
	{
	  in >> s;
	  bool found = false;
	  for (int iv=0; iv < vars.Length() && !found; iv++)
	    {
	      int l = (vars[iv].Length() < s.Length()) ? vars[iv].Length() : s.Length();
	      if (strncmp (vars[iv].Pointer(), s.Pointer(), l-1) == 0)
		{
		  found = true;
		  used [iv] = 1;
		}
	    }
	}
    }

  vals.Allocate (num, numvals);
  ids.Allocate (num);
  int *pi = ids.Pointer();
  double *pv = vals.Pointer();
  for (int i=0; i < num; i++)
    {
      in >> *pi++;
      for (int j=0; j < numvals; j++)
	in >> *pv++;
    }
}
