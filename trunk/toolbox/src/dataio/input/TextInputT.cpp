/* $Id: TextInputT.cpp,v 1.3 2004-07-19 20:32:14 paklein Exp $ */
#include "TextInputT.h"

#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "ifstreamT.h"
#include "dArrayT.h"
#ifdef _MSC_VER
#include <strstrea.h>
#elif defined(__GCC_3__)
#include <strstream>
#else
#include <strstream.h>
#endif

using namespace Tahoe;

TextInputT::TextInputT (ostream& out) :
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

bool TextInputT::Open (const StringT& filename)
{
  /* create file root */
  StringT suffix;
  suffix.Suffix (filename.Pointer());
  if (strncmp (suffix.Pointer(), ".geo", 4) == 0 ||
      strncmp (suffix.Pointer(), ".run", 4) == 0 ||
      strncmp (suffix.Pointer(), ".in", 3) == 0)
  fFileRoot.Root (filename);
  fFilePath.FilePath(fFileRoot);

	/* scan geometry file */
	ifstreamT geo;
	if (!OpenFile (geo, ".geo")) {
		cout << "\n TextInputT::Open: error opening geometry file: " << geo.filename() << endl;
		return false;
	}
	if (!ScanGeometryFile (geo)) {
		cout << "\n TextInputT::Open: error scanning geometry file: " << geo.filename() << endl;
		return false;
	}

	/* scan results file */
	ifstreamT run;
	if (!OpenFile (run, ".run")) {
		cout << "\n TextInputT::Open: error opening results file: " << run.filename() << endl;
		return false;
	}
    if (!ScanResultsFile (run)) {
		cout << "\n TextInputT::Open: error scanning results file: " << run.filename() << endl;
		return false;
    }
      
	/* must be OK */
	return true;
}

void TextInputT::Close (void)
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

void TextInputT::QuadratureVariablesUsed (const StringT& name, iArrayT& used)
{
#ifdef __MWERKS__
#pragma unused (name)
#endif
	used = 0;
}

void TextInputT::ReadNodeSet(const StringT& name, iArrayT& nodes)
{
#ifdef __MWERKS__
#pragma unused (name)
#endif
	nodes.Free();
}

void TextInputT::ReadSideSetLocal(const StringT& setname, iArray2DT& sides) const
{
#ifdef __MWERKS__
#pragma unused (setname)
#endif
	sides.Free ();
}

void TextInputT::ReadSideSetGlobal(const StringT& setname, iArray2DT& sides) const
{
#ifdef __MWERKS__
#pragma unused (setname)
#endif
	sides.Free ();
}

void TextInputT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const
{
	qlabels.Free(); 
}

void TextInputT::ReadNodeSetVariables(int step, const StringT& name, dArray2DT& nvalues)
{
#ifdef __MWERKS__
#pragma unused (step)
#pragma unused (name)
#endif
  nvalues.Free();
}

void TextInputT::ReadAllQuadratureVariable(int step, int varindex, dArrayT& values)
{
#ifdef __MWERKS__
#pragma unused (step)
#pragma unused (varindex)
#endif
	values.Free();
}

void TextInputT::ReadQuadratureVariable(int step, const StringT& name, int varindex, dArrayT& values)
{
#ifdef __MWERKS__
#pragma unused (step)
#pragma unused (name)
#pragma unused (varindex)
#endif
	values.Free();
}

void TextInputT::ReadAllQuadratureVariables(int step, dArray2DT& vals)
{
#ifdef __MWERKS__
#pragma unused (step)
#endif
	vals.Free();
}

void TextInputT::ReadQuadratureVariables(int step, const StringT& name, dArray2DT& vals)
{
#ifdef __MWERKS__
#pragma unused (step)
#pragma unused (name)
#endif
	vals.Free();
}

void TextInputT::ElementGroupNames (ArrayT<StringT>& groupnames) const
{
  if (groupnames.Length() != fBlockID.Length()) throw ExceptionT::kSizeMismatch;
  for (int i=0; i < groupnames.Length(); i++)
    {
      groupnames[i] = fBlockID[i];
    }
}

void TextInputT::ReadNodeID(iArrayT& node_id)
{
	dArray2DT coords(fNumNodes, fNumDOF);
	ReadCoordinates(coords, node_id);
}

void TextInputT::ReadCoordinates (dArray2DT& coords)
{
	iArrayT node_id(fNumNodes);
	ReadCoordinates(coords, node_id);
}

void TextInputT::ReadCoordinates (dArray2DT& coords, iArrayT& node_id)
{
	if (node_id.Length() != fNumNodes ||
       coords.MajorDim() != fNumNodes ||
       coords.MinorDim() != fNumDOF ) throw ExceptionT::kSizeMismatch;

	ifstreamT geo;
	OpenFile(geo, ".geo");

	/* advance */
	StringT s;
	if (!geo.FindString ("G E O M E T R Y   D A T A", s)) throw ExceptionT::kDatabaseFail;
	if (!geo.FindString ("Nodal coordinates", s)) throw ExceptionT::kDatabaseFail;

	/* read data */
	iArrayT used;
	DataBlock(geo, used, node_id, coords, true);
	//NOTE: assumes you'll have at least as many nodal output variables
	//      as there are spatial dimensions
}

int TextInputT::NumElements(const StringT& name)
{
	int dex = fBlockID.PositionOf(name);
	if (dex == -1) {
		cout << "\n TextInputT::NumElements: could not find block ID " << name << endl;
		throw ExceptionT::kDatabaseFail;
	}
	return fBlockNumElem[dex];
}

int TextInputT::NumElementNodes(const StringT& name)
{
	int dex = fBlockID.PositionOf(name);
	if (dex == -1) {
		cout << "\n TextInputT::NumElementNodes: could not find block ID " << name << endl;
		throw ExceptionT::kDatabaseFail;
	}
	return fBlockNumElemNode[dex];
}

void TextInputT::ReadAllElementMap (iArrayT& elemmap)
{
  if (elemmap.Length() != fNumElements) throw ExceptionT::kSizeMismatch;

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
	  !geo.FindString ("element", s)) throw ExceptionT::kDatabaseFail;

      for (int i=0; i < numelms; i++)
	{
	  geo >> elemmap[count++];
	  geo.getline (line, 254);
	}
    }
}

void TextInputT::ReadGlobalElementMap (const StringT& name, iArrayT& elemmap)
{
  ifstreamT geo;
  OpenFile (geo, ".geo");

  StringT s;
  int numelms;
  char line [255];
  if (!AdvanceToBlock (geo, name, "Connectivities") ||
      !geo.FindString ("Number of elements", s) ||
      !s.Tail ('=', numelms)) throw ExceptionT::kDatabaseFail;
  if (elemmap.Length() != numelms) throw ExceptionT::kSizeMismatch;

	/* advance to the start of the connectivity block */
	if (!geo.FindString("index", s)) throw ExceptionT::kDatabaseFail;
	
	/* read map */
	elemmap = 0;
	for (int i=0; i < numelms; i++)
	{
		int index;
		geo >> index >> elemmap[i];
		geo.getline (line, 254);
    }
}

void TextInputT::ReadGlobalElementSet (const StringT& name, iArrayT& set)
{
  if (set.Length() != fNumElements) throw ExceptionT::kSizeMismatch;

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
	  !geo.FindString ("element", s)) throw ExceptionT::kDatabaseFail;

      count += numelms;
    }
  
  if (set.Length() != numelms) throw ExceptionT::kSizeMismatch;
  set.SetValueToPosition();
  set += count;
}

void TextInputT::ReadConnectivity (const StringT& name, iArray2DT& connects)
{
  ifstreamT geo;
  OpenFile (geo, ".geo");

  if (!AdvanceToBlock (geo, name, "Connectivities")) throw ExceptionT::kDatabaseFail;

  StringT s;
  int numelms, numelnodes;
  iArrayT elms (connects.MinorDim());
  if (!geo.FindString ("Number of elements", s) ||
      !s.Tail ('=', numelms) ||
      !geo.FindString ("Number of element nodes", s) ||
      !s.Tail ('=', numelnodes) ||
      !geo.FindString ("element", s)) throw ExceptionT::kDatabaseFail;

  if (numelms != connects.MajorDim() ||
      numelnodes != connects.MinorDim()) throw ExceptionT::kSizeMismatch;

	for (int i=0; i < numelms; i++)
	{
		int elm_dex, elm_id; 
		geo >> elm_dex>> elm_id >> elms;
		connects.SetRow (i, elms);
	}
	connects--;
}

void TextInputT::ReadGeometryCode (const StringT& name, GeometryT::CodeT& geocode)
{
	int dex = fBlockID.PositionOf(name);
	if (dex == -1) {
		cout << "\n TextInputT::ReadGeometryCode: could not find block ID " << name << endl;
		throw ExceptionT::kDatabaseFail;
	}
	geocode = fBlockGeometry[dex];
}

void TextInputT::ReadTimeSteps (dArrayT& steps) 
{
	steps.Dimension(fTimeSteps.Length());
	fTimeSteps.CopyInto(steps);
}

void TextInputT::ReadNodeLabels (ArrayT<StringT>& nlabels) const
{
	/* allocate */
	nlabels.Dimension(NumNodeVariables());

	/* copy */
	for (int i=0; i < nlabels.Length(); i++)
		nlabels[i] = fNodeVariable[i];
}

void TextInputT::ReadElementLabels (ArrayT<StringT>& elabels) const
{
	/* allocate */
	elabels.Dimension(NumElementVariables());

	/* copy */
  	for (int i=0; i < elabels.Length(); i++)
		elabels[i] = fElementVariable[i];
}

void TextInputT::NodeVariablesUsed (const StringT& name, iArrayT& used)
{ 
  if (used.Length() != fNodeVariable.Length()) throw ExceptionT::kSizeMismatch;
  used = 0;

  ifstreamT run;
  OpenFile (run, ".run");
  
  if (!AdvanceToBlock (run, name, "Nodal data")) throw ExceptionT::kDatabaseFail;

  dArray2DT vals;
  iArrayT ids;
  DataBlock (run, used, ids, vals, true);
}

void TextInputT::ElementVariablesUsed (const StringT& name, iArrayT& used)
{ 
  if (used.Length() != fElementVariable.Length()) throw ExceptionT::kSizeMismatch;
  used = 0;

  ifstreamT run;
  OpenFile (run, ".run");
  
  if (!AdvanceToBlock (run, name, "Element data")) throw ExceptionT::kDatabaseFail;

  dArray2DT vals;
  iArrayT ids;
  DataBlock (run, used, ids, vals, false);
}

void TextInputT::ReadAllNodeVariables(int step, dArray2DT& nvalues)
{
	const char caller[] = "TextInputT::ReadNodeVariables";

	StringT toc_file(fFileRoot);
	toc_file.Append(".run");
	if (is_old_format(toc_file)) {
		ReadAllNodeVariables_old(step, nvalues);
		return;
	}

	/* get file for specified step */
	StringT file;
	ResultsFile(toc_file, step, file);

	/* open results file */
	StringT results_file(fFilePath);
	results_file.Append(file);
  	ifstreamT run(results_file);
  	if (!run.is_open())
  		ExceptionT::GeneralFail(caller, "could not open file %s", results_file.Pointer());

	/* advance to the edge of the nodal data block */
	StringT s;
	if (!run.FindString ("Nodal data", s)) ExceptionT::DatabaseFail(caller);

	/* read */
	iArrayT used (fNodeVariable.Length()), ids;
	DataBlock(run, used, ids, nvalues, true);
}

void TextInputT::ReadAllNodeVariable (int step, int varindex, dArrayT& values)
{
#ifdef __MWERKS__
#pragma unused (step)
#pragma unused (varindex)
#endif
	values.Free();
	ExceptionT::GeneralFail("InputFEASIIT::ReadAllNodeVariable", "not yet programmed");
}

void TextInputT::ReadNodeVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
#ifdef __MWERKS__
#pragma unused (step)
#pragma unused (name)
#pragma unused (varindex)
#endif
	values.Free();
	ExceptionT::GeneralFail("InputFEASIIT::ReadNodeVariable", "not yet programmed");
}

void TextInputT::ReadNodeVariables(int step, const StringT& name, dArray2DT& nvalues)
{
	const char caller[] = "TextInputT::ReadNodeVariables";
		
	StringT toc_file(fFileRoot);
	toc_file.Append(".run");
	if (is_old_format(toc_file)) {
		ReadNodeVariables_old(step, name, nvalues);
		return;
	}
	
	if (nvalues.Length() == 0) ExceptionT::SizeMismatch(caller);

	/* get file for specified step */
	StringT file;
	ResultsFile(toc_file, step, file);

	/* open results file */
	StringT results_file(fFilePath);
	results_file.Append(file);
  	ifstreamT run(results_file);
  	if (!run.is_open())
  		ExceptionT::GeneralFail(caller, "could not open file %s", results_file.Pointer());

  	if (!AdvanceToBlock (run, name, "Nodal data")) ExceptionT::DatabaseFail(caller);

	iArrayT used (fNodeVariable.Length()), ids;
	dArray2DT vals;
	DataBlock (run, used, ids, vals, true);

	nvalues = 0;
	for (int i=0; i < ids.Length(); i++)
		for (int v=0, j=0; v < used.Length(); v++)
			if (used[v] > 0)
				nvalues (i, v) = vals (i, j++);
}

void TextInputT::ReadAllElementVariable (int step, int varindex, dArrayT& values)
{
#ifdef __MWERKS__
#pragma unused (step)
#pragma unused (varindex)
#endif
	values.Free();
	ExceptionT::GeneralFail("InputFEASIIT::ReadAllNodeVariable", "not yet programmed");
}

void TextInputT::ReadElementVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
#ifdef __MWERKS__
#pragma unused (step)
#pragma unused (name)
#pragma unused (varindex)
#endif
	values.Free();
	ExceptionT::GeneralFail("InputFEASIIT::ReadNodeVariable", "not yet programmed");
}

void TextInputT::ReadAllElementVariables(int step, dArray2DT& evalues)
{
	const char caller[] = "TextInputT::ReadAllElementVariables";
		
	StringT toc_file(fFileRoot);
	toc_file.Append(".run");
	if (is_old_format(toc_file)) {
		ReadAllElementVariables_old(step, evalues);
		return;
	}

	/* get file for specified step */
	StringT file;
	ResultsFile(toc_file, step, file);

	/* open results file */
	StringT results_file(fFilePath);
	results_file.Append(file);
  	ifstreamT run(results_file);
  	if (!run.is_open())
  		ExceptionT::GeneralFail(caller, "could not open file %s", results_file.Pointer());

	/* advance to the edge of the nodal data block */
	StringT s;
	if (!run.FindString ("Nodal data", s)) ExceptionT::DatabaseFail(caller);
	if (evalues.MajorDim() != fNumElements) ExceptionT::SizeMismatch(caller);

	iArrayT used (fElementVariable.Length()), ids;
	dArray2DT vals;
	evalues = 0;
	int dex = 0;
	for (int i=0; i < fBlockID.Length(); i++)
	{
		/* advance to start of block */
		if (!run.FindString ("Element data", s)) ExceptionT::DatabaseFail(caller);

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
    if (evalues.MinorDim() > 0 && dex != fNumElements)
    	ExceptionT::DatabaseFail(caller, "error joining values in blocks");
}

void TextInputT::ReadElementVariables(int step, const StringT& name, dArray2DT& evalues)
{
	const char caller[] = "TextInputT::ReadElementVariables";
		
	StringT toc_file(fFileRoot);
	toc_file.Append(".run");
	if (is_old_format(toc_file)) {
		ReadElementVariables_old(step, name, evalues);
		return;
	}

	/* resolve block index */
	int dex = fBlockID.PositionOf(name);
	if (dex == -1)
		ExceptionT::DatabaseFail(caller, "could not find block ID %s", name.Pointer());

	/* get file for specified step */
	StringT file;
	ResultsFile(toc_file, step, file);

	/* open results file */
	StringT results_file(fFilePath);
	results_file.Append(file);
  	ifstreamT run(results_file);
  	if (!run.is_open())
  		ExceptionT::GeneralFail(caller, "could not open file %s", results_file.Pointer());
	
	/* advance to the edge of the nodal data block */
	StringT s;
	if (!run.FindString ("Nodal data", s)) ExceptionT::DatabaseFail(caller);

	/* advance to block */
	for (int i = 0; i <= dex; i++)
		if (!run.FindString ("Element data", s)) 
			ExceptionT::DatabaseFail(caller);

	/* verify block */
	StringT block_ID;
	if (!run.FindString ("Block ID", s) ||
        !s.Tail('=', block_ID)) ExceptionT::DatabaseFail(caller);
	if (name != block_ID)
		ExceptionT::DatabaseFail(caller, "found block ID %s at position %d instead of ID %s",
			block_ID.Pointer(), dex, name.Pointer());

	/* read */
	iArrayT used (fElementVariable.Length()), ids;
	DataBlock(run, used, ids, evalues, false);	
}

/******************** PRIVATE ***************************/

bool TextInputT::OpenFile (ifstreamT& in, const char* ext) const
{
  StringT file (fFileRoot);
  file.Append (ext);
  in.open (file);
  if (!in.is_open()) 
    {
      fout << "\nTextInputT::OpenFile unable to open " << file << "\n\n";
      return false;
    }
  return true;
}

bool TextInputT::ScanGeometryFile (ifstreamT& in)
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
		StringT nid;
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

bool TextInputT::ScanResultsFile(ifstreamT& in)
{
	/* is old format */
	if (is_old_format(in.filename())) return ScanResultsFile_old(in);

	/* advance */
	StringT s;
	if (!in.FindString ("O U T P U T", s)) return false;

	/* first file name */
	StringT file;
	in >> file;
	if (!in.good()) return false;
	
	/* open first results file */
	StringT path, file_path;
	file_path.FilePath(in.filename());
	file_path.Append(file);
	ifstreamT dat(file_path);
	if (!dat.is_open()) {
		cout << "\n TextInputT::ScanResultsFile: error opening file: " << file_path << endl;
		return false;
	}

	/* get dimension from first file */
	if (!dat.FindString ("Group number", s)) return false;
	double t;
	if (!dat.FindString ("Time", s) || !s.Tail ('=', t)) return false;
	fTimeSteps.Append(t);
	if (!dat.FindString ("Number of blocks", s)) return false;

	/* scan nodal output labels */
	fNodeVariable.Free();
	int vals;
	if (!dat.FindString ("Nodal data", s) ||
        !dat.FindString ("Number of values", s) ||
        !s.Tail('=', vals)) {
		cout << "\n TextInputT::ScanResultsFile: error scanning nodal values" << endl;
		return false;
	}
	if (vals > 0) {
		dat >> s >> s; /* "index" and "node" */
		for (int v = 0; v < vals; v++) {
			dat >> s;
			fNodeVariable.Append(s);
		}
	}

	/* scan element output labels (from first block) */
	fElementVariable.Free();
	StringT id;
	if (!dat.FindString ("Element data", s) ||
	    !dat.FindString ("Block ID", s) ||
	    !s.Tail ('=', id) ||
	    !fBlockID.HasValue (id) ||
	    !dat.FindString ("Number of values", s) ||
	    !s.Tail ('=', vals)) {
		cout << "\n TextInputT::ScanResultsFile: error scanning element values" << endl;
		return false;
	}
	if (vals > 0) {
		dat >> s >> s; /* "index" and "element" */
		for (int v = 0; v < vals; v++) {
			dat >> s;
			fElementVariable.Append(s);
		}
	}
	dat.close();

	/* run through remaining file list */
	in >> file;
	while (in.good())
	{
		file_path.FilePath(in.filename());
		file_path.Append(file);
		dat.open(file_path);
		if (!dat.is_open()) {
			cout << "\n TextInputT::ScanResultsFile: error opening file: " << file_path << endl;
			return false;
		}
		double t;
		if (!dat.FindString ("Time", s) || !s.Tail ('=', t)) return false;
		fTimeSteps.Append(t);	
		dat.close();
		
		/* next file */
		in >> file;
	}
	
	return true;
}

bool TextInputT::ScanResultsFile_old(ifstreamT& in)
{
	/* advance */
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
	  cout << "\n TextInputT::ScanResultsFile: error scanning nodal values" << endl;
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
	StringT id;
	if (!in.FindString ("Element data", s) ||
	    !in.FindString ("Block ID", s) ||
	    !s.Tail ('=', id) ||
	    !fBlockID.HasValue (id) ||
	    !in.FindString ("Number of values", s) ||
	    !s.Tail ('=', vals)) {
	  cout << "\n TextInputT::ScanResultsFile: error scanning element values" << endl;
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

bool TextInputT::AdvanceToBlock (ifstreamT& in, const StringT& name, const char* tname) const
{
  //const int ID = atoi (name.Pointer());
  //int found = -1;
  bool found = false;
  StringT item;
  StringT s;
  while (!found && in.good())
    {
      if (!in.FindString (tname, s) ||
	  !in.FindString ("Block ID", s) ||
	  !s.Tail ('=', item)) 
	{
	  cout << "\n\nTextInputT::AdvanceToBlock unable to find:\n";
	  cout << "    Name = " << name << "\n    tname = " << tname;
	  cout << "\n    Last string read: " << s << endl;
	  return false;
	}
      if (strncmp (name.Pointer(), item.Pointer(), name.StringLength()) == 0)
	return true;
    }
  cout << "\n\nTextInputT::AdvanceToBlock end of file:\n";
  cout << "    Name = " << name << "\n    tname = " << tname;
  return false;
}

void TextInputT::DataBlock (ifstreamT& in, iArrayT& used, iArrayT& ids, dArray2DT& vals, bool nodal) const
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
      if (used.Length() < vars.Length()) used.Dimension(vars.Length());
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
      !s.Tail ('=', numvals)) throw ExceptionT::kDatabaseFail;

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
	else ids.Dimension(0);
}

/* return true if results file uses the pre-TOC format */
bool TextInputT::is_old_format(const StringT& file) const
{
	const char caller[] = "TextInputT::is_old_format";
	ifstreamT in(file);
	if (!in.is_open()) ExceptionT::DatabaseFail(caller, "could not open file: %s", file.Pointer());

	/* look for new format with table for contents */
	StringT s;
	if (!in.FindString ("O U T P U T", s)) ExceptionT::DatabaseFail(caller);
	StringT tail;
	s.Tail(':', tail);
	if (tail != "T O C")
		return true;
	else
		return false;
}

void TextInputT::ReadNodeVariables_old(int step, const StringT& name, dArray2DT& nvalues)
{
#ifdef __MWERKS__
#pragma unused(step)
#endif

	if (nvalues.Length() == 0) throw ExceptionT::kSizeMismatch;

	ifstreamT run;
	OpenFile (run, ".run");

	if (!AdvanceToBlock (run, name, "Nodal data")) throw ExceptionT::kDatabaseFail;

	iArrayT used (fNodeVariable.Length()), ids;
	dArray2DT vals;
	DataBlock (run, used, ids, vals, true);

	nvalues = 0;
	for (int i=0; i < ids.Length(); i++)
		for (int v=0, j=0; v < used.Length(); v++)
			if (used[v] > 0)
				nvalues (i, v) = vals (i, j++);
}

void TextInputT::ReadAllNodeVariables_old(int step, dArray2DT& nvalues)
{
	if (step < 0) throw ExceptionT::kDatabaseFail;
	if (nvalues.MajorDim() != fNumNodes) throw ExceptionT::kSizeMismatch;
	ifstreamT run;
	OpenFile (run, ".run");

	StringT s;
	if (!run.FindString("O U T P U T", s)) throw ExceptionT::kDatabaseFail;

	int count = 0;
	bool OK = run.FindString ("Group number", s);
	while (count < step && OK) {
		OK = run.FindString ("Group number", s);
		count++;
	}

	/* check */
	if (!OK) {
		cout << "\n TextInputT::ReadAllNodeVariables: could not find step index " 
		     << step << endl;
		throw ExceptionT::kDatabaseFail;
	}
	
	/* advance to the edge of the nodal data block */
	if (!run.FindString ("Nodal data", s)) throw ExceptionT::kDatabaseFail;

	/* read */
	iArrayT used (fNodeVariable.Length()), ids;
	DataBlock(run, used, ids, nvalues, true);
}

/* return the results file name for the given output step */
void TextInputT::ResultsFile(const StringT& toc_file, int step, StringT& file) const
{
	const char caller[] = "TextInputT::ResultsFile";
	if (is_old_format(toc_file))
		ExceptionT::DatabaseFail(caller, "data base contains no TOC: %s", toc_file.Pointer());
	
	ifstreamT in(toc_file);

	/* advance */
	StringT s;
	if (!in.FindString ("O U T P U T", s)) 
		ExceptionT::DatabaseFail(caller);

	/* get specified file name */
	int i = 0;
	for (int i = 0; in.good() && i <= step; i++)
		in >> file;
	
	/* check */
	if (!in.good())
		ExceptionT::DatabaseFail(caller, "step %d > %d: %s", step, i, toc_file.Pointer());
}

void TextInputT::ReadAllElementVariables_old(int step, dArray2DT& evalues)
{
	if (evalues.MajorDim() != fNumElements) throw ExceptionT::kSizeMismatch;
	if (step < 0) throw ExceptionT::kDatabaseFail;

	/* input stream */
	ifstreamT run;
	OpenFile(run, ".run");
	StringT s;
	if (!run.FindString("O U T P U T", s)) throw ExceptionT::kDatabaseFail;

	/* advance to step */	
	int count = 0;
	bool OK = run.FindString ("Group number", s);
	while (count < step && OK) {
		OK = run.FindString ("Group number", s);
		count++;
	}

	/* check */
	if (!OK) {
		cout << "\n TextInputT::ReadAllElementVariables: could not find step index " 
		     << step << endl;
		throw ExceptionT::kDatabaseFail;
	}

	/* advance to the edge of the nodal data block */
	if (!run.FindString ("Nodal data", s)) throw ExceptionT::kDatabaseFail;

	iArrayT used (fElementVariable.Length()), ids;
	dArray2DT vals;
	evalues = 0;
	int dex = 0;
	for (int i=0; i < fBlockID.Length(); i++)
	{
		/* advance to start of block */
		if (!run.FindString ("Element data", s)) throw ExceptionT::kDatabaseFail;

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
    	cout << "\n TextInputT::ReadAllElementVariables: error joining values in blocks" << endl;
    	throw ExceptionT::kDatabaseFail;
    }
}

void TextInputT::ReadElementVariables_old(int step, const StringT& name, dArray2DT& evalues)
{
	if (step < 0) throw ExceptionT::kDatabaseFail;

	/* resolve block index */
	int dex = fBlockID.PositionOf(name);
	if (dex == -1) {
		cout << "\n TextInputT::ReadElementVariables: could not find block ID " << name << endl;
		throw ExceptionT::kDatabaseFail;
	}

	/* input stream */
	ifstreamT run;
	OpenFile(run, ".run");
	StringT s;
	if (!run.FindString("O U T P U T", s)) throw ExceptionT::kDatabaseFail;

	/* advance to step */	
	int count = 0;
	bool OK = run.FindString ("Group number", s);
	while (count < step && OK) {
		OK = run.FindString ("Group number", s);
		count++;
	}

	/* check */
	if (!OK) {
		cout << "\n TextInputT::ReadElementVariables: could not find step index " 
		     << step << endl;
		throw ExceptionT::kDatabaseFail;
	}
	
	/* advance to the edge of the nodal data block */
	if (!run.FindString ("Nodal data", s)) throw ExceptionT::kDatabaseFail;

	/* advance to block */
	for (int i = 0; i <= dex; i++)
		if (!run.FindString ("Element data", s)) 
			throw ExceptionT::kDatabaseFail;

	/* verify block */
	StringT block_ID;
	if (!run.FindString ("Block ID", s) ||
        !s.Tail('=', block_ID)) throw ExceptionT::kDatabaseFail;
	if (name != block_ID) {
		cout << "\n TextInputT::ReadElementVariables: found block ID " << block_ID << '\n'
		     <<   "     at position " << dex << " instead of block ID " << name << endl;
		throw ExceptionT::kDatabaseFail;
	}

	/* read */
	iArrayT used (fElementVariable.Length()), ids;
	DataBlock(run, used, ids, evalues, false);	
}
