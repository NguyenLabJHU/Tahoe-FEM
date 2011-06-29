/* $Id: EnSightInputT.cpp,v 1.16 2004-06-17 06:38:24 paklein Exp $ */
/* created: sawimme (05/18/1998) */
#include "EnSightInputT.h"
#include "ifstreamT.h"
#include "dArray2DT.h"
#include "dArrayT.h"
#include "AutoArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"

using namespace Tahoe;

const int kDOF = 3; // all files are always written in 3D

EnSightInputT::EnSightInputT (ostream& out, bool binary) :
  InputBaseT (out),
  fData (out, binary, kDOF),
  fPartDimensions (0,0),
  fStartIncrement(-1),
  fIncrement(-1)
{

}

bool EnSightInputT::Open (const StringT& file)
{
	fCaseFile = file;
	fCaseFile.ToNativePathName();
	ifstreamT in ('#', fCaseFile);
	if (!in.is_open()) {
		cout << "\n EnSightInputT::Open: error opening case file: " << in.filename() << endl;	
		return false;
	}
	if (!fData.CaseFile (in, fGeometryFile)) {
		cout << "\n EnSightInputT::Open: case file format incorrect: " << in.filename() << endl;
		return false;
	}
	ScanGeometryFile ();
	
	/* must be OK */
	return true;
}

void EnSightInputT::ElementGroupNames (ArrayT<StringT>& groupnames) const
{
	ifstreamT in;
	OpenRelative(in, fGeometryFile, fCaseFile);

//  int id = -1;
  bool nodemap, elementmap;
  int numnodes, numelems, numelemnodes;
  fData.ReadGeometryHeader(in, nodemap, elementmap);
  
  iArrayT ids (fPartDimensions.MajorDim());
  for (int i=0; i < fPartDimensions.MajorDim(); i++)
    {
      fData.ReadPart (in, ids[i]);
      fData.SkipPart (in, nodemap, elementmap, numnodes, numelems, numelemnodes);
    }
  for (int j=0; j < ids.Length(); j++) {
  	groupnames[j].Clear();
    groupnames[j].Append(ids[j]);
    }
}

int EnSightInputT::NumElementGroups (void) const
{
	ifstreamT in;
	OpenRelative(in, fGeometryFile, fCaseFile);
  bool nodemap, elementmap;
  fData.ReadGeometryHeader(in, nodemap, elementmap);
  return fPartDimensions.MajorDim();
}

int EnSightInputT::NumNodes (void) const
{
  int num = 0;
  for (int i=0; i < fPartDimensions.MajorDim(); i++)
    num += fPartDimensions (i, kNumNodes);
  return num;
}

void EnSightInputT::ReadNodeID(iArrayT& node_id)
{
  if (node_id.Length() != NumNodes()) throw ExceptionT::kSizeMismatch;

  ScanGeometryFile();

  int id;
  bool nodemapgiven, elemmapgiven;

	ifstreamT in;
	OpenRelative(in, fGeometryFile, fCaseFile);
  fData.ReadGeometryHeader(in, nodemapgiven, elemmapgiven);

  if (nodemapgiven)
    {
      int num_elems, num_elem_nodes;
      dArray2DT tempcoords;
      iArrayT map;
      int num = 0;
      for (int i2=0; i2 < fPartDimensions.MajorDim(); i2++)
	{
	  fData.ReadPart (in, id);
	  
	  // read coordinate data
	  fData.ReadCoordinates (in, tempcoords, map, nodemapgiven);
	  node_id.CopyPart (num, map, 0, map.Length());
	  num += tempcoords.MajorDim();

	  // skip connectivity
	  fData.SkipConnectivity (in, num_elems, num_elem_nodes, elemmapgiven);
	}
    }
  else
    {
      node_id.SetValueToPosition ();
      node_id += 1;
    }
}

void EnSightInputT::ReadCoordinates (dArray2DT& coords)
{
  ScanGeometryFile();
  
  int id;
  bool nodemapgiven, elemmapgiven;
	ifstreamT in;
	OpenRelative(in, fGeometryFile, fCaseFile);
  fData.ReadGeometryHeader(in, nodemapgiven, elemmapgiven);
  
  if (coords.MajorDim() != NumNodes() ||
      coords.MinorDim() != 3) throw ExceptionT::kSizeMismatch;

  int num = 0;
  dArray2DT tempcoords;
  iArrayT map;
  int num_elems, num_elem_nodes;
  for (int i2=0; i2 < fPartDimensions.MajorDim(); i2++)
    {
      fData.ReadPart (in, id);
      
      // read coordinate data
      fData.ReadCoordinates (in, tempcoords, map, nodemapgiven);
      
      // store coordinate data
      coords.CopyPart (num*kDOF, tempcoords, 0, tempcoords.Length());
     num += tempcoords.MajorDim();
     
     // skip connectivity
     fData.SkipConnectivity (in, num_elems, num_elem_nodes, elemmapgiven);
    }
}

void EnSightInputT::ReadCoordinates (dArray2DT& coords, iArrayT& nodemap)
{
  ScanGeometryFile();
  
  int id;
  bool nodemapgiven, elemmapgiven;
	ifstreamT in;
	OpenRelative(in, fGeometryFile, fCaseFile);
  fData.ReadGeometryHeader(in, nodemapgiven, elemmapgiven);
  
  int num = 0;
  for (int i=0; i < fPartDimensions.MajorDim(); i++)
    num += fPartDimensions (i, kNumNodes);
  coords.Allocate (num, kDOF);
  if (nodemapgiven) nodemap.Allocate (num);
 
  dArray2DT tempcoords;
  iArrayT map;
  int num_elems, num_elem_nodes;
  num = 0;
  for (int i2=0; i2 < fPartDimensions.MajorDim(); i2++)
    {
      fData.ReadPart (in, id);
      
      // read coordinate data
      fData.ReadCoordinates (in, tempcoords, map, nodemapgiven);
      
      // store coordinate data
      coords.CopyPart (num*kDOF, tempcoords, 0, tempcoords.Length());
     if (nodemapgiven) nodemap.CopyPart (num, map, 0, map.Length());
     num += tempcoords.MajorDim();
     
     // skip connectivity
     fData.SkipConnectivity (in, num_elems, num_elem_nodes, elemmapgiven);
    }
}

int EnSightInputT::NumGlobalElements (void) const
{
  int num = 0;
  for (int i=0; i < fPartDimensions.MajorDim(); i++)
    num += fPartDimensions (i, kNumElements);
  return num;
}

int EnSightInputT::NumElements (const StringT& name)
{
  int part, ID = atoi (name.Pointer());
  //fPartDimensions.WriteNumbered (cout);
  fPartDimensions.ColumnHasValue (kPartID, ID, part); 
  return fPartDimensions (part, kNumElements);
}

int EnSightInputT::NumElementNodes (const StringT& name)
{
  int group = atoi (name.Pointer());

  ScanGeometryFile();
	ifstreamT in;
	OpenRelative(in, fGeometryFile, fCaseFile);
  int id = -1;

  // read header
  bool nodemapgiven, elemmapgiven;
  fData.ReadGeometryHeader(in, nodemapgiven, elemmapgiven);
  
  // forward down to part of interest
  int numnodes, numelems, numelemnodes, offset = 0;
  fData.ReadPart (in, id);
  while (id != group)
    {
      fData.SkipPart (in, nodemapgiven, elemmapgiven, numnodes, numelems, numelemnodes);
      offset += numelems;
      fData.ReadPart (in, id);
    }

  // skip coordinate data
  int num_nodes;
  fData.SkipCoordinates (in, num_nodes, nodemapgiven);
  
  // read connectivity
  iArray2DT connects;
  iArrayT elementmap;
  GeometryT::CodeT geocode;
  fData.ReadConnectivity (in, connects, elementmap, elemmapgiven, geocode);
  
  return connects.MinorDim();
}

void EnSightInputT::ReadAllElementMap (iArrayT& elemmap)
{
  if (elemmap.Length() != NumGlobalElements()) throw ExceptionT::kSizeMismatch;

  ScanGeometryFile();
	ifstreamT in;
	OpenRelative(in, fGeometryFile, fCaseFile);
  int id = -1;

  // read header
  bool nodemapgiven, elemmapgiven;
  fData.ReadGeometryHeader(in, nodemapgiven, elemmapgiven);
  
  if (elemmapgiven)
    {
      iArray2DT connects;
      iArrayT elementmap;
      GeometryT::CodeT geocode;
      int num=0;
      for (int i2=0; i2 < fPartDimensions.MajorDim(); i2++)
	{
	  fData.ReadPart (in, id);

	  // skip coordinate data
	  int num_nodes;
	  fData.SkipCoordinates (in, num_nodes, nodemapgiven);
  
	  // read connectivity
	  fData.ReadConnectivity (in, connects, elementmap, elemmapgiven, geocode);
	  elemmap.CopyPart (num, elementmap, 0, elementmap.Length());
	  num += elementmap.Length();
	}
    }
  else
    {
      elemmap.SetValueToPosition ();
      elemmap += 1;
    }
}

void EnSightInputT::ReadGlobalElementMap (const StringT& name, iArrayT& elementmap)
{
  int group = atoi (name.Pointer());
  ScanGeometryFile();
	ifstreamT in;
	OpenRelative(in, fGeometryFile, fCaseFile);
  int id = -1;
  
  // read header
  bool nodemapgiven, elemmapgiven;
  fData.ReadGeometryHeader(in, nodemapgiven, elemmapgiven);
  
  // forward down to part of interest
  int numnodes, numelems, numelemnodes, offset = 0;
  fData.ReadPart (in, id);
  while (id != group)
    {
      fData.SkipPart (in, nodemapgiven, elemmapgiven, numnodes, numelems, numelemnodes);
      offset += numelems;
      fData.ReadPart (in, id);
    }
  
  if (elemmapgiven)
    {
      // skip coordinate data
      int num_nodes;
      fData.SkipCoordinates (in, num_nodes, nodemapgiven);
      
      // read connectivity
      iArray2DT connects;
      GeometryT::CodeT geocode;
      fData.ReadConnectivity (in, connects, elementmap, elemmapgiven, geocode);
    }
  else
    {
      elementmap.SetValueToPosition ();
      elementmap += offset;
    }
}

void EnSightInputT::ReadGlobalElementSet (const StringT& name, iArrayT& set)
{
  ReadGlobalElementMap (name, set);
  set += -1;
}

void EnSightInputT::ReadConnectivity (const StringT& name, iArray2DT& connects)
{
  int group = atoi (name.Pointer());
  ScanGeometryFile();
	ifstreamT in;
	OpenRelative(in, fGeometryFile, fCaseFile);
  int id = -1;
  
  // read header
  bool nodemapgiven, elemmapgiven;
  fData.ReadGeometryHeader(in, nodemapgiven, elemmapgiven);
  
  // forward down to part of interest
  int numnodes, numelems, numelemnodes, offset = 0;
  fData.ReadPart (in, id);
  while (id != group)
    {
      fData.SkipPart (in, nodemapgiven, elemmapgiven, numnodes, numelems, numelemnodes);
      offset += numelems;
      fData.ReadPart (in, id);
    }
  
  // skip coordinate data
  int num_nodes;
  fData.SkipCoordinates (in, num_nodes, nodemapgiven);
  
  // read connectivity
  iArrayT elementmap;
  GeometryT::CodeT geocode;
  fData.ReadConnectivity (in, connects, elementmap, elemmapgiven, geocode);
  
  // renumber to global node numbering, offset to start at zero
  connects += offset - 1;
}

void EnSightInputT::ReadGeometryCode (const StringT& name, GeometryT::CodeT& geocode)
{
  int group = atoi (name.Pointer());
  ScanGeometryFile();
	ifstreamT in;
	OpenRelative(in, fGeometryFile, fCaseFile);
  int id = -1;
  
  // read header
  bool nodemapgiven, elemmapgiven;
  fData.ReadGeometryHeader(in, nodemapgiven, elemmapgiven);
  
  // forward down to part of interest
  int numnodes, numelems, numelemnodes, offset = 0;
  fData.ReadPart (in, id);
  while (id != group)
    {
      fData.SkipPart (in, nodemapgiven, elemmapgiven, numnodes, numelems, numelemnodes);
      offset += numelems;
      fData.ReadPart (in, id);
    }
  
  // skip coordinate data
  int num_nodes;
  fData.SkipCoordinates (in, num_nodes, nodemapgiven);
  
  // read connectivity
  iArray2DT connects;
  iArrayT elementmap;
  fData.ReadConnectivity (in, connects, elementmap, elemmapgiven, geocode);
}

void EnSightInputT::QARecords (ArrayT<StringT>& records)
{
  records.Allocate (4);
  records[0] = "EnSight";
  records[1] = "v6 Gold";
}

int  EnSightInputT::NumTimeSteps (void) const
{
  ifstreamT in ('#', fCaseFile);
  if (!AdvanceStream (in, "TIME")) return 0;
  return fData.NumTimeSteps (in);
}

void EnSightInputT::ReadTimeSteps (dArrayT& steps)
{
  ifstreamT in ('#', fCaseFile);
  steps.Allocate (0);
  if (!AdvanceStream (in, "TIME")) return;
  fData.ReadTimeSection (in, fStartIncrement, fIncrement, steps);
}

int  EnSightInputT::NumNodeVariables (void) const
{
  ifstreamT incase ('#', fCaseFile);
  if (!AdvanceStream (incase, "VARIABLE")) return 0;

  AutoArrayT<bool> nvector, evector;
  AutoArrayT<StringT> nl(20), el(20);
  fData.ReadVariableSection (incase, nl, el, nvector, evector, false);
  return nl.Length();
}

int  EnSightInputT::NumElementVariables (void) const
{
  ifstreamT incase ('#', fCaseFile);
  if (!AdvanceStream (incase, "VARIABLE")) return 0;

  AutoArrayT<bool> nvector, evector;
  AutoArrayT<StringT> nl(20), el(20);
  fData.ReadVariableSection (incase, nl, el, nvector, evector, false);
  return el.Length();
}

void EnSightInputT::ReadNodeLabels (ArrayT<StringT>& nlabels) const
{
  ifstreamT incase ('#', fCaseFile);
  if (!AdvanceStream (incase, "VARIABLE"))
    {
      nlabels.Allocate (0);
      return;
    }
  
  AutoArrayT<bool> nvector, evector;
  AutoArrayT<StringT> nl(20), el(20);
  fData.ReadVariableSection (incase, nl, el, nvector, evector, false);
  
  nlabels.Allocate (nl.Length());
  for (int j=0; j < nl.Length(); j++)
    nlabels[j] = nl[j];
}

void EnSightInputT::ReadElementLabels (ArrayT<StringT>& elabels) const
{
  ifstreamT incase ('#', fCaseFile);
  if (!AdvanceStream (incase, "VARIABLE"))
    {
      elabels.Allocate (0);
      return;
    }
  
  AutoArrayT<bool> nvector, evector;
  AutoArrayT<StringT> nl(20), el(20);
  fData.ReadVariableSection (incase, nl, el, nvector, evector, false);
  
  elabels.Allocate (el.Length());
  for (int j=0; j < el.Length(); j++)
    elabels[j] = el[j];
}

void EnSightInputT::NodeVariablesUsed (const StringT& name, iArrayT& used)
{ 
  used = 0;

  // gather variable labels and vector/scalar data
  ifstreamT incase ('#', fCaseFile);
  if (!AdvanceStream (incase, "VARIABLE")) return;
  
  AutoArrayT<bool> nvector, evector;
  AutoArrayT<StringT> nl(20), el(20);
  fData.ReadVariableSection (incase, nl, el, nvector, evector, false);
  if (used.Length() != nl.Length()) throw ExceptionT::kSizeMismatch;
  
  VariableUsed (name, used, nl, nvector, true);
}

void EnSightInputT::ElementVariablesUsed (const StringT& name, iArrayT& used)
{ 
  used = 0;

  // gather variable labels and vector/scalar data
  ifstreamT incase ('#', fCaseFile);
  if (!AdvanceStream (incase, "VARIABLE")) return;
  
  AutoArrayT<bool> nvector, evector;
  AutoArrayT<StringT> nl(20), el(20);
  fData.ReadVariableSection (incase, nl, el, nvector, evector, false);
  if (used.Length() != nl.Length()) throw ExceptionT::kSizeMismatch;
  
  VariableUsed (name, used, el, evector, false);
}

void EnSightInputT::ReadAllNodeVariable (int step, int varindex, dArrayT& values)
{
  int numg = NumElementGroups ();
  ArrayT<StringT> gnames (numg);
  ElementGroupNames (gnames);

  int offset = 0;
  dArrayT vals;
  for (int i=0; i < numg; i++)
    {
      ReadNodeVariable (step, gnames[i], varindex, vals);
      values.CopyPart (offset, vals, 0, vals.Length());
      offset += vals.Length();
    }
}

void EnSightInputT::ReadNodeVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
  int group_id, currentinc;
  VarPrelims_Geo (step, name, group_id, currentinc);

  // read variable filenames;
  AutoArrayT<bool> nvector;
  AutoArrayT<StringT> nl (20);
  VarPrelims_Case (nvector, nl, true);

  // if vector, determine component index
  int component = -1;
  if (nvector[varindex])
    component = ComponentIndex (varindex, nl);

  // read data
  ReadOneVariableData (component, nl[varindex], group_id, values, currentinc, true);
}

void EnSightInputT::ReadAllNodeVariables (int step, dArray2DT& nvalues)
{
  int numg = NumElementGroups ();
  ArrayT<StringT> gnames  (numg);
  ElementGroupNames (gnames);

  int offset = 0;
  dArray2DT vals;
  for (int i=0; i < numg; i++)
    {
      ReadNodeVariables (step, gnames[i], vals);

      if (vals.MinorDim() != nvalues.MinorDim()) throw ExceptionT::kSizeMismatch;

      nvalues.CopyPart (offset, vals, 0, vals.Length());
      offset += vals.Length();
    }
}

void EnSightInputT::ReadNodeVariables (int step, const StringT& name, dArray2DT& nvalues)
{
  int group_id, currentinc;
  VarPrelims_Geo (step, name, group_id, currentinc);

  // make sure variables exist and read filenames
  nvalues.Allocate (0,0);
  AutoArrayT<bool> nvector;
  AutoArrayT<StringT> nl (20);
  VarPrelims_Case (nvector, nl, true);
  
  // allocate space
  int dex;
  if (!fPartDimensions.ColumnHasValue (kPartID, group_id, dex)) return;
  nvalues.Allocate (fPartDimensions (dex, kNumNodes), nl.Length());
  nvalues = 0.0;
  
  // read data
  ReadVariableData (nvector, nl, group_id, nvalues, currentinc, true);
}

void EnSightInputT::ReadAllElementVariable (int step, int varindex, dArrayT& values)
{
  int numg = NumElementGroups ();
  ArrayT<StringT> gnames (numg);
  ElementGroupNames (gnames);

  int offset = 0;
  dArrayT vals;
  for (int i=0; i < numg; i++)
    {
      ReadElementVariable (step, gnames[i], varindex, vals);
      values.CopyPart (offset, vals, 0, vals.Length());
      offset += vals.Length();
    }
}

void EnSightInputT::ReadElementVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
  int group_id, currentinc;
  VarPrelims_Geo (step, name, group_id, currentinc);

  // read variable filenames
  AutoArrayT<bool> evector;
  AutoArrayT<StringT> el (20);
  VarPrelims_Case (evector, el, false);

  // if vector, determine component index
  int component = -1;
  if (evector [varindex])
    component = ComponentIndex (varindex, el);

  //read data
  ReadOneVariableData (component, el[varindex], group_id, values, currentinc, false);
}

void EnSightInputT::ReadAllElementVariables (int step, dArray2DT& evalues)
{
  int numg = NumElementGroups ();
  ArrayT<StringT> gnames  (numg);
  ElementGroupNames (gnames);

  int offset = 0;
  dArray2DT vals;
  for (int i=0; i < numg; i++)
    {
      ReadElementVariables (step, gnames[i], vals);

      if (vals.MinorDim() != evalues.MinorDim()) throw ExceptionT::kSizeMismatch;

      evalues.CopyPart (offset, vals, 0, vals.Length());
      offset += vals.Length();
    }
}

void EnSightInputT::ReadElementVariables (int step, const StringT& name, dArray2DT& evalues)
{
  int group_id, currentinc;
  VarPrelims_Geo (step, name, group_id, currentinc);
  
  // make sure variables exist
  evalues.Allocate (0,0);
  AutoArrayT<bool> evector;
  AutoArrayT<StringT> el (20);
  VarPrelims_Case (evector, el, false);
  
  // allocate space
  int dex;
  if (!fPartDimensions.ColumnHasValue (kPartID, group_id, dex)) return;
  evalues.Allocate (fPartDimensions (dex, kNumElements), el.Length());
  evalues = 0.0;
  
  // read data
  ReadVariableData (evector, el, group_id, evalues, currentinc, false);
}

/*************************************************************************
*
* Private
*
*************************************************************************/

bool EnSightInputT::AdvanceStream (istream& in, const char* key) const
{
  StringT word (81);
  StringT line (256);
  while (in.good())
    {
      in >> word;
      in.getline (line.Pointer(), 254); // clear line
      if (word == key) return true;
    }
  return false;
}

void EnSightInputT::ScanGeometryFile (void)
{
  bool nm, em;
  int partID;
  int numnodes, numelems, numelemnodes;
  if (fPartDimensions.Length() == 0)
    {
      ifstreamT in;
      OpenRelative(in, fGeometryFile, fCaseFile);
      fData.ReadGeometryHeader (in, nm, em);
      while (fData.ReadPart (in, partID))
	{
	  fData.SkipPart (in, nm, em, numnodes, numelems, numelemnodes);

	  int length = fPartDimensions.MajorDim();
	  if (length == 0)
	    fPartDimensions.Allocate (1,kPartDims);
	  else
	    fPartDimensions.Resize(length + 1, 0);

	  fPartDimensions (length, kNumNodes) = numnodes;
	  fPartDimensions (length, kNumElements) = numelems;
	  fPartDimensions (length, kPartID) = partID;
	}
    }
}

StringT EnSightInputT::CreateVariableFile (const StringT& old, int inc) const
{
  // determine number of wild cards;
  int numwild = 0;
  int index = 0;
  for (int i=0; i < old.Length(); i++)
    {
      if (old[i] == '*') numwild++;
      if (numwild == 0) index++;
    }
  
  StringT fileinc;
  fileinc.Append (inc, numwild);
  StringT filename = old;
  filename.CopyPart (index, fileinc, 0, numwild);

  return filename;
}

int EnSightInputT::ComponentIndex (int varindex, ArrayT<StringT>& labels) const
{
  for (int i=varindex+1; i < varindex+kDOF-1; i++)
    {
      if (strncmp (labels[varindex].Pointer(),labels[i].Pointer(), labels[i].StringLength()) != 0) 
	return i-varindex;
    }
  return 0;
}

void EnSightInputT::VarPrelims_Geo (int step, const StringT& name, int& group_id, int& currentinc)
{
  group_id = atoi (name.Pointer());
  // set fStartIncrement
  if (fStartIncrement < 0)
    {
      dArrayT temp;
      ReadTimeSteps (temp);
      if (temp.Length() == 0) {
      	currentinc = 0;
      	return;
      }
    }
  currentinc = fStartIncrement + fIncrement*step;

  // set fPartDimensions
  ScanGeometryFile();  
}

void EnSightInputT::VarPrelims_Case (AutoArrayT<bool>& vector, AutoArrayT<StringT>& labels, bool node)
{
  // make sure variable exists
  ifstreamT incase ('#', fCaseFile);
  if (!AdvanceStream (incase, "VARIABLE")) return;

  // read variable filenames;
  AutoArrayT<bool> tempvector;
  AutoArrayT<StringT> templabel (20);
  if (node)
    fData.ReadVariableSection (incase, labels, templabel, vector, tempvector, true);
  else
    fData.ReadVariableSection (incase, templabel, labels, tempvector, vector, true);
}

void EnSightInputT::ReadOneVariableData (int component, const StringT& label, int group_id, dArrayT& values, int currentinc, bool nodal) const
{
	StringT filename = CreateVariableFile (label, currentinc);
	ifstreamT in;
	OpenRelative(in, filename, fCaseFile);

  // read variable file header
  StringT header;
  fData.ReadVariableHeader (in, header);

  // search file for part 
  bool found = false;
  int id, dex, num;
  dArray2DT temp;
  while (fData.ReadPart (in, id) && !found)
    {
      if (!fPartDimensions.ColumnHasValue (kPartID, id, dex)) throw ExceptionT::kGeneralFail;

      if (nodal)
	num = fPartDimensions (dex, kNumNodes);
      else
	num = fPartDimensions (dex, kNumElements);

      if (component > -1 || component < kDOF)
	temp.Allocate (num, kDOF);
      else
	temp.Allocate (num, 1);
      temp = 0.0;

      fData.ReadVariable (in, temp);

      // found part
      if (id == group_id)
	{
	  // figure out which component to copy
	  if (component > -1 || component < kDOF)
	    {
	      dArray2DT temp2 (temp.MinorDim(), temp.MajorDim());
	      temp2.Transpose (temp);
	      values.CopyPart (0, temp2, temp.MinorDim()*component, temp.MinorDim());
	    }
	  else
	    values.CopyPart (0, temp, 0, temp.Length());
	  found = true;
	}
    }
}

void EnSightInputT::ReadVariableData (ArrayT<bool>& vector, ArrayT<StringT>& labels, int group_id, dArray2DT& values, int currentinc, bool nodal) const
{
  for (int i=0; i < labels.Length(); i++)
    {
		StringT filename = CreateVariableFile (labels[i], currentinc);
		ifstreamT in;
		OpenRelative(in, filename, fCaseFile);

      // read variable file header
      StringT header;
      fData.ReadVariableHeader (in, header);
      
      // search file for part
      bool found = false;
      int id, dex, num;
      dArray2DT temp;
      while (fData.ReadPart (in, id) && !found)
	{
	  if (!fPartDimensions.ColumnHasValue (kPartID, id, dex)) throw ExceptionT::kGeneralFail;
	  
	  if (nodal)
	    num = fPartDimensions (dex, kNumNodes);
	  else
	    num = fPartDimensions (dex, kNumElements);

	  if (vector[i])
	    temp.Allocate (num, kDOF);
	  else
	    temp.Allocate (num, 1);
	  temp = 0.0;

	  fData.ReadVariable (in, temp);

	  // found part
	  if (id == group_id)
	    {
	      values.BlockColumnCopyAt (temp, i);
	      found = true;
	    }
	}

      if (vector[i])
	i += kDOF - 1;
    }
}

void EnSightInputT::VariableUsed (const StringT& name, iArrayT& used, ArrayT<StringT>& labels, ArrayT<bool>& vector, bool nodal) const
{
  // examine each variable file for entries with this part
  int group_id = atoi (name.Pointer());
  for (int i=0; i < labels.Length(); i++)
    {
		StringT filename = CreateVariableFile (labels[i], 0);
		ifstreamT in;
		OpenRelative(in, filename, fCaseFile);

      // read variable file header
      StringT header;
      fData.ReadVariableHeader (in, header);

      // search file for part
      bool found = false;
      int id, dex, num;
      dArray2DT temp;
      while (fData.ReadPart (in, id) && !found)
	{
	  if (!fPartDimensions.ColumnHasValue (kPartID, id, dex)) throw ExceptionT::kGeneralFail;

	  if (nodal)
	    num = fPartDimensions (dex, kNumNodes);
	  else
	    num = fPartDimensions (dex, kNumElements);

	  if (vector[i])
	    temp.Allocate (num, kDOF);
	  else
	    temp.Allocate (num, 1);
	  temp = 0.0;

	  fData.ReadVariable (in, temp);

	  // found part
	  if (id == group_id)
	    {
	      if (vector[i])
		for (int hg=0; hg < kDOF; hg++)
		  used[i+hg] = kDOF;
	      else
		used[i] = 1;
	    }
	}

      if (vector[i])
	i += kDOF - 1;
    }
}
