/* $Id: InputFEASCIIT.cpp,v 1.12 2002-03-04 18:02:32 paklein Exp $ */

#include "InputFEASCIIT.h"

#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
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
	fBlockNumElem.Free();
	fBlockNumElemNode.Free();
	fBlockGeometry.Free();
  
  
  fNumNodes = 0;
  fNumElements = 0;
  fNumDOF = 0;
  fTimeSteps.Free();
  fNodeVariable.Free();
  fElementVariable.Free();
}

void InputFEASCIIT::QuadratureVariablesUsed (const StringT& name, iArrayT& used)
{
#pragma unused (name)
	used = 0;
}

void InputFEASCIIT::ReadNodeSet(const StringT& name, iArrayT& nodes)
{
#pragma unused (name)
	nodes.Free();
}

void InputFEASCIIT::ReadSideSetLocal(const StringT& setname, iArray2DT& sides) const
{
#pragma unused (setname)
	sides.Free ();
}

void InputFEASCIIT::ReadSideSetGlobal(const StringT& setname, iArray2DT& sides) const
{
#pragma unused (setname)
	sides.Free ();
}

void InputFEASCIIT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const
{
	qlabels.Free(); 
}

void InputFEASCIIT::ReadNodeSetVariables(int step, const StringT& name, dArray2DT& nvalues)
{
#pragma unused (step)
#pragma unused (name)
  nvalues.Free();
}

void InputFEASCIIT::ReadAllQuadratureVariable(int step, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (varindex)
	values.Free();
}

void InputFEASCIIT::ReadQuadratureVariable(int step, const StringT& name, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (name)
#pragma unused (varindex)
	values.Free();
}

void InputFEASCIIT::ReadAllQuadratureVariables(int step, dArray2DT& vals)
{
#pragma unused (step)
	vals.Free();
}

void InputFEASCIIT::ReadQuadratureVariables(int step, const StringT& name, dArray2DT& vals)
{
#pragma unused (step)
#pragma unused (name)
	vals.Free();
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

void InputFEASCIIT::ReadNodeID(iArrayT& node_id)
{
	dArray2DT coords(fNumNodes, fNumDOF);
	ReadCoordinates(coords, node_id);
}

void InputFEASCIIT::ReadCoordinates (dArray2DT& coords)
{
	iArrayT node_id(fNumNodes);
	ReadCoordinates(coords, node_id);
}

void InputFEASCIIT::ReadCoordinates (dArray2DT& coords, iArrayT& node_id)
{
	if (node_id.Length() != fNumNodes ||
       coords.MajorDim() != fNumNodes ||
       coords.MinorDim() != fNumDOF ) throw eSizeMismatch;

	ifstreamT geo;
	OpenFile(geo, ".geo");

	/* advance */
	StringT s;
	if (!geo.FindString ("G E O M E T R Y   D A T A", s)) throw eDatabaseFail;
	if (!geo.FindString ("Nodal coordinates", s)) throw eDatabaseFail;

	/* read data */
	iArrayT used;
	DataBlock(geo, used, node_id, coords, true);
	//NOTE: assumes you'll have at least as many nodal output variables
	//      as there are spatial dimensions
}

int InputFEASCIIT::NumElements(const StringT& name)
{
	int ID = atoi(name);
	int dex = fBlockID.PositionOf(ID);
	if (dex == -1) {
		cout << "\n InputFEASCIIT::NumElements: could not find block ID " << ID << endl;
		throw eDatabaseFail;
	}
	return fBlockNumElem[dex];
}

int InputFEASCIIT::NumElementNodes(const StringT& name)
{
	int ID = atoi(name);
	int dex = fBlockID.PositionOf(ID);
	if (dex == -1) {
		cout << "\n InputFEASCIIT::NumElementNodes: could not find block ID " << ID << endl;
		throw eDatabaseFail;
	}
	return fBlockNumElemNode[dex];
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

void InputFEASCIIT::ReadGlobalElementMap (const StringT& name, iArrayT& elemmap)
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

void InputFEASCIIT::ReadGlobalElementSet (const StringT& name, iArrayT& set)
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

void InputFEASCIIT::ReadConnectivity (const StringT& name, iArray2DT& connects)
{
  if (connects.Length() == 0) throw eSizeMismatch;

  ifstreamT geo;
  OpenFile (geo, ".geo");

  if (!AdvanceToBlock (geo, name, "Connectivities")) throw eDatabaseFail;

  StringT s;
  int numelms, numelnodes;
  iArrayT elms (connects.MinorDim());
  if (!geo.FindString ("Number of elements", s) ||
      !s.Tail ('=', numelms) ||
      !geo.FindString ("Number of element nodes", s) ||
      !s.Tail ('=', numelnodes) ||
      !geo.FindString ("element", s)) throw eDatabaseFail;

  if (numelnodes != connects.MinorDim()) throw eSizeMismatch;

	for (int i=0; i < numelms; i++)
	{
		int elm_dex, elm_id; 
		geo >> elm_dex>> elm_id >> elms;
		connects.SetRow (i, elms);
	}
	connects--;
}

void InputFEASCIIT::ReadGeometryCode (const StringT& name, GeometryT::CodeT& geocode)
{
	int ID = atoi(name);
	int dex = fBlockID.PositionOf(ID);
	if (dex == -1) {
		cout << "\n InputFEASCIIT::ReadGeometryCode: could not find block ID " << ID << endl;
		throw eDatabaseFail;
	}
	geocode = fBlockGeometry[dex];
}

void InputFEASCIIT::ReadTimeSteps (dArrayT& steps) 
{
	steps.Allocate(fTimeSteps.Length());
	fTimeSteps.CopyInto(steps);
}

void InputFEASCIIT::ReadNodeLabels (ArrayT<StringT>& nlabels) const
{
	/* allocate */
	nlabels.Allocate(NumNodeVariables());

	/* copy */
	for (int i=0; i < nlabels.Length(); i++)
		nlabels[i] = fNodeVariable[i];
}

void InputFEASCIIT::ReadElementLabels (ArrayT<StringT>& elabels) const
{
	/* allocate */
	elabels.Allocate(NumElementVariables());

	/* copy */
  	for (int i=0; i < elabels.Length(); i++)
		elabels[i] = fElementVariable[i];
}

void InputFEASCIIT::NodeVariablesUsed (const StringT& name, iArrayT& used)
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

void InputFEASCIIT::ElementVariablesUsed (const StringT& name, iArrayT& used)
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

void InputFEASCIIT::ReadAllNodeVariables(int step, dArray2DT& nvalues)
{
	if (step < 0) throw eDatabaseFail;
	if (nvalues.MajorDim() != fNumNodes) throw eSizeMismatch;
	ifstreamT run;
	OpenFile (run, ".run");

	StringT s;
	if (!run.FindString("O U T P U T", s)) throw eDatabaseFail;

	int count = 0;
	bool OK = run.FindString ("Group number", s);
	while (count < step && OK) {
		OK = run.FindString ("Group number", s);
		count++;
	}

	/* check */
	if (!OK) {
		cout << "\n InputFEASCIIT::ReadAllNodeVariables: could not find step index " 
		     << step << endl;
		throw eDatabaseFail;
	}
	
	/* advance to the edge of the nodal data block */
	if (!run.FindString ("Nodal data", s)) throw eDatabaseFail;

	/* read */
	iArrayT used (fNodeVariable.Length()), ids;
	DataBlock(run, used, ids, nvalues, true);
}

void InputFEASCIIT::ReadAllNodeVariable (int step, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (varindex)
  values.Free();
  cout << "InputFEASIIT::ReadAllNodeVariable not yet programmed\n\n";
  throw eGeneralFail;
}

void InputFEASCIIT::ReadNodeVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (name)
#pragma unused (varindex)
  values.Free();
  cout << "InputFEASIIT::ReadNodeVariable not yet programmed\n\n";
  throw eGeneralFail;
}

void InputFEASCIIT::ReadNodeVariables (int step, const StringT& name, dArray2DT& nvalues)
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

void InputFEASCIIT::ReadAllElementVariable (int step, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (varindex)
  values.Free();
  cout << "InputFEASIIT::ReadAllNodeVariable not yet programmed\n\n";
  throw eGeneralFail;
}

void InputFEASCIIT::ReadElementVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (name)
#pragma unused (varindex)
  values.Free();
  cout << "InputFEASIIT::ReadNodeVariable not yet programmed\n\n";
  throw eGeneralFail;
}

void InputFEASCIIT::ReadAllElementVariables (int step, dArray2DT& evalues)
{
	if (evalues.MajorDim() != fNumElements) throw eSizeMismatch;
	if (step < 0) throw eDatabaseFail;

	/* input stream */
	ifstreamT run;
	OpenFile(run, ".run");
	StringT s;
	if (!run.FindString("O U T P U T", s)) throw eDatabaseFail;

	/* advance to step */	
	int count = 0;
	bool OK = run.FindString ("Group number", s);
	while (count < step && OK) {
		OK = run.FindString ("Group number", s);
		count++;
	}

	/* check */
	if (!OK) {
		cout << "\n InputFEASCIIT::ReadAllElementVariables: could not find step index " 
		     << step << endl;
		throw eDatabaseFail;
	}

	/* advance to the edge of the nodal data block */
	if (!run.FindString ("Nodal data", s)) throw eDatabaseFail;

	iArrayT used (fElementVariable.Length()), ids;
	dArray2DT vals;
	evalues = 0;
	int dex = 0;
	for (int i=0; i < fBlockID.Length(); i++)
	{
		/* advance to start of block */
		if (!run.FindString ("Element data", s)) throw eDatabaseFail;

		/* read block */
		DataBlock (run, used, ids, vals, false);

		/* fill from top to bottom */
		for (int i = 0; i < ids.Length(); i++)
		{
			for (int v = 0, j = 0; v < used.Length(); v++)
	  			if (used[v] > 0) evalues(dex, v) = vals(i, j++);
	  		dex++;
	    }
    }
    
    /* check */
    if (evalues.MinorDim() > 0 && dex != fNumElements) {
    	cout << "\n InputFEASCIIT::ReadAllElementVariables: error joining values in blocks" << endl;
    	throw eDatabaseFail;
    }
}

void InputFEASCIIT::ReadElementVariables(int step, const StringT& name, dArray2DT& evalues)
{
	if (step < 0) throw eDatabaseFail;

	/* resolve block index */
	int ID = atoi(name);
	int dex = fBlockID.PositionOf(ID);
	if (dex == -1) {
		cout << "\n InputFEASCIIT::ReadElementVariables: could not find block ID " << ID << endl;
		throw eDatabaseFail;
	}

	/* input stream */
	ifstreamT run;
	OpenFile(run, ".run");
	StringT s;
	if (!run.FindString("O U T P U T", s)) throw eDatabaseFail;

	/* advance to step */	
	int count = 0;
	bool OK = run.FindString ("Group number", s);
	while (count < step && OK) {
		OK = run.FindString ("Group number", s);
		count++;
	}

	/* check */
	if (!OK) {
		cout << "\n InputFEASCIIT::ReadElementVariables: could not find step index " 
		     << step << endl;
		throw eDatabaseFail;
	}
	
	/* advance to the edge of the nodal data block */
	if (!run.FindString ("Nodal data", s)) throw eDatabaseFail;

	/* advance to block */
	for (int i = 0; i <= dex; i++)
		if (!run.FindString ("Element data", s)) 
			throw eDatabaseFail;

	/* verify block */
	int block_ID;
	if (!run.FindString ("Block ID", s) ||
        !s.Tail('=', block_ID)) throw eDatabaseFail;
	if (ID != block_ID) {
		cout << "\n InputFEASCIIT::ReadElementVariables: found block ID " << block_ID << '\n'
		     <<   "     at position " << dex << " instead of block ID " << ID << endl;
		throw eDatabaseFail;
	}

	/* read */
	iArrayT used (fElementVariable.Length()), ids;
	DataBlock(run, used, ids, evalues, false);	
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
  
	/* get number of element blocks */
	int num_blocks;
	if (!in.FindString ("Number of blocks", s) ||
	    !s.Tail ('=', num_blocks)) return false;

	/* coordinate data */
	if (!in.FindString ("Nodal coordinates", s)) return false;
	fNumNodes = 0;
	fNumDOF = 0;
	if (!in.FindString ("Number of nodal points", s) || !s.Tail ('=', fNumNodes)) 
		return false;
	if (!in.FindString ("Number of values", s) || !s.Tail ('=', fNumDOF)) 
		return false;
	
	/* scan block data */
	fNumElements = 0;
	if (!in.FindString ("Connectivities", s)) return false;
	for (int i = 0; i < num_blocks; i++)
	{
		int nid;
		if (!in.FindString ("Block ID", s) || !s.Tail ('=', nid)) 
			return false;
		fBlockID.Append(nid);
	
		int nel;
		if (!in.FindString ("Number of elements", s) || !s.Tail ('=', nel))
			return false;
		fBlockNumElem.Append(nel);
		fNumElements += nel;

		int nen;
		if (!in.FindString ("Number of element nodes", s) || !s.Tail ('=', nen))
			return false;
		fBlockNumElemNode.Append(nen);

		int icode;
		if (!in.FindString("Geometry code", s) || !s.Tail ('=', icode))
			return false;
		GeometryT::CodeT code = GeometryT::CodeT(icode);
		fBlockGeometry.Append(code);
	}

	/* return */
	if (fBlockID.Length() < 1) 
		return false;
	else
		return true;
}

bool InputFEASCIIT::ScanResultsFile (ifstreamT& in)
{
	StringT s;
	if (!in.FindString ("O U T P U T", s)) return false;

	fNodeVariable.Free();
	fElementVariable.Free();

	if (!in.FindString ("Group number", s)) return false;
	double t;
	if (!in.FindString ("Time", s) || !s.Tail ('=', t)) return false;
	fTimeSteps.Append(t);
	if (!in.FindString ("Number of blocks", s)) return false;

	/* scan nodal output labels */
	int vals;
	if (!in.FindString ("Nodal data", s) ||
        !in.FindString ("Number of values", s) ||
        !s.Tail ('=', vals)) {
	  cout << "\n InputFEASCIIT::ScanResultsFile: error scanning nodal values" << endl;
	  return false;
	}
	if (vals > 0) {
		in >> s >> s; // "index" and "node"
		for (int v = 0; v < vals; v++)
		{
			in >> s;
			fNodeVariable.Append(s);
		}
	}

	/* scan element output labels (from first block) */
	int id;
	if (!in.FindString ("Element data", s) ||
	    !in.FindString ("Block ID", s) ||
	    !s.Tail ('=', id) ||
	    !fBlockID.HasValue (id) ||
	    !in.FindString ("Number of values", s) ||
	    !s.Tail ('=', vals)) {
	  cout << "\n InputFEASCIIT::ScanResultsFile: error scanning element values" << endl;
	  return false;
	}
	if (vals > 0) {
		in >> s >> s; // "index" and "element"
		for (int v = 0; v < vals; v++)
		{
			in >> s;
			fElementVariable.Append(s);
		}
	}

	/* determine the remaining time steps */
	while (in.FindString ("Group number", s)) {
		double t;
		if (!in.FindString ("Time", s) || !s.Tail ('=', t)) return false;
		fTimeSteps.Append(t);
	}

	/* OK */
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
	  !in.FindString ("Block ID", s) ||
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
      
      /* bad patch - DataBlock used to read nodal output values as well
       * as nodal coordinates, which assume you have at least as many
       * output values as spatial dimensions */
      if (used.Length() < vars.Length()) used.Allocate(vars.Length());
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
		in >> s >> s; // read "index" and "element" | "node"
		for (int v = 0; v < numvals; v++)
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

	/* read values */
	vals.Allocate (num, numvals);
	if (numvals > 0)
	{
		ids.Allocate (num);
		int *pi = ids.Pointer();
		double *pv = vals.Pointer();
		for (int i=0; i < num; i++)
		{
			int index;
			in >> index >> *pi++;
			for (int j=0; j < numvals; j++)
			in >> *pv++;
		}
	}
	else ids.Allocate(0);
}
