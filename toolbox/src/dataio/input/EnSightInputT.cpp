/* $Id: EnSightInputT.cpp,v 1.4 2001-08-07 23:11:54 paklein Exp $ */
/* created: sawimme (05/18/1998)                                          */

#include "EnSightInputT.h"
#include "ios_fwd_decl.h"
#include "fstreamT.h"
#include "dArray2DT.h"
#include "dArrayT.h"
#include "AutoArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"

const int kDOF = 3; // all files are always written in 3D

EnSightInputT::EnSightInputT (ostream& out, bool binary) :
  InputBaseT (out),
  fData (out, binary, kDOF),
  fPartDimensions (0,0)
{
}

void EnSightInputT::Open (const StringT& file)
{
  fCaseFile = file;
  ifstreamT in ('#', fCaseFile);
  if (!fData.CaseFile (in, fGeometryFile))
    {
      cout << "\n\nEnSightInputT, case file format incorrect.\n ";
      throw eGeneralFail;
    }
}

int EnSightInputT::NumElementGroups (void)
{
  ifstream in (fGeometryFile);
  bool nodemap, elementmap;
  fData.ReadGeometryHeader(in, nodemap, elementmap);
  return fPartDimensions.MajorDim();
}

int EnSightInputT::NumNodes (void)
{
  int num = 0;
  for (int i=0; i < fPartDimensions.MajorDim(); i++)
    num += fPartDimensions (i, 0);
  return num;
}

void EnSightInputT::ReadNodeMap (iArrayT& nodemap)
{
  if (nodemap.Length() != NumNodes()) throw eSizeMismatch;

  ScanGeometryFile();

  int id;
  bool nodemapgiven, elemmapgiven;
  ifstream in (fGeometryFile);
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
	  nodemap.CopyPart (num, map, 0, map.Length());
	  num += tempcoords.MajorDim();

	  // skip connectivity
	  fData.SkipConnectivity (in, num_elems, num_elem_nodes, elemmapgiven);
	}
    }
  else
    {
      nodemap.SetValueToPosition ();
      nodemap += 1;
    }
}

void EnSightInputT::ReadCoordinates (dArray2DT& coords)
{
  ScanGeometryFile();
  
  int id;
  bool nodemapgiven, elemmapgiven;
  ifstream in (fGeometryFile);
  fData.ReadGeometryHeader(in, nodemapgiven, elemmapgiven);
  
  if (coords.MajorDim() != NumNodes() ||
      coords.MinorDim() != 3) throw eSizeMismatch;

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
  ifstream in (fGeometryFile);
  fData.ReadGeometryHeader(in, nodemapgiven, elemmapgiven);
  
  int num = 0;
  for (int i=0; i < fPartDimensions.MajorDim(); i++)
    num += fPartDimensions (i, 0);
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

int EnSightInputT::NumGlobalElements (void)
{
  int num = 0;
  for (int i=0; i < fPartDimensions.MajorDim(); i++)
    num += fPartDimensions (i, 1);
  return num;
}

void EnSightInputT::ReadAllElementMap (iArrayT& elemmap)
{
  if (elemmap.Length() != NumGlobalElements()) throw eSizeMismatch;

  ScanGeometryFile();
  ifstream in (fGeometryFile);
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

int  EnSightInputT::NumTimeSteps (void)
{
  dArrayT steps;
  ReadTimeSteps (steps);
  return steps.Length();
}

void EnSightInputT::ReadTimeSteps (dArrayT& steps)
{
  ifstreamT in ('#', fCaseFile);
  steps.Allocate (0);
  if (!AdvanceStream (in, "TIME")) return;
  fData.ReadTimeSection (in, fStartIncrement, fIncrement, steps);
}

int  EnSightInputT::NumNodeVariables (void)
{
  ifstreamT incase ('#', fCaseFile);
  if (!AdvanceStream (incase, "VARIABLE")) return 0;

  AutoArrayT<bool> nvector, evector;
  AutoArrayT<StringT> nl(20), el(20);
  fData.ReadVariableSection (incase, nl, el, nvector, evector, false);
  return nl.Length();
}

int  EnSightInputT::NumElementVariables (void)
{
  ifstreamT incase ('#', fCaseFile);
  if (!AdvanceStream (incase, "VARIABLE")) return 0;

  AutoArrayT<bool> nvector, evector;
  AutoArrayT<StringT> nl(20), el(20);
  fData.ReadVariableSection (incase, nl, el, nvector, evector, false);
  return el.Length();
}

/*************************************************************************
*
* Protected
*
*************************************************************************/

void EnSightInputT::ElementGroupIDs (iArrayT& groupnums)
{
  ifstream in (fGeometryFile);
  int id = -1;
  bool nodemap, elementmap;
  int numnodes, numelems, numelemnodes;
  fData.ReadGeometryHeader(in, nodemap, elementmap);
  
  iArrayT ids (fPartDimensions.MajorDim());
  for (int i=0; i < fPartDimensions.MajorDim(); i++)
    {
      fData.ReadPart (in, ids[i]);
      fData.SkipPart (in, nodemap, elementmap, numnodes, numelems, numelemnodes);
    }
  groupnums.Allocate (ids.Length());
  groupnums.CopyPart (0, ids, 0, ids.Length());
}

int EnSightInputT::NumElements_ID (int ID)
{
  int part;
  fPartDimensions.ColumnHasValue (2, ID, part); 
  return fPartDimensions (1, part);
}

int EnSightInputT::NumElementNodes_ID (int group)
{
  ScanGeometryFile();
  ifstream in (fGeometryFile);
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

void EnSightInputT::ReadGlobalElementMap_ID (int group, iArrayT& elementmap)
{
  ScanGeometryFile();
  ifstream in (fGeometryFile);
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

void EnSightInputT::ReadConnectivity_ID (int group, iArray2DT& connects)
{
  ScanGeometryFile();
  ifstream in (fGeometryFile);
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

void EnSightInputT::ReadGeometryCode_ID (int group, GeometryT::CodeT& geocode)
{
  ScanGeometryFile();
  ifstream in (fGeometryFile);
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

void EnSightInputT::ReadElementLabels_ID (int ID, ArrayT<StringT>& elabels)
{
#pragma unused (ID)
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

void EnSightInputT::ReadElementVariables_ID (int step, int group_id, dArray2DT& evalues)
{
  // set fStartIncrement
  if (fStartIncrement < 0)
    {
      dArrayT temp;
      ReadTimeSteps (temp);
      if (temp.Length() == 0) return;
    }
  int currentinc = fStartIncrement + fIncrement*step;
  
  // set fPartDimensions
  ScanGeometryFile();
  
  // make sure variables exist
  evalues.Allocate (0,0);
  ifstreamT incase ('#', fCaseFile);
  if (!AdvanceStream (incase, "VARIABLE"))  return;
  
  // read variable filenames;
  AutoArrayT<bool> nvector, evector;
  AutoArrayT<StringT> nl (20), el (20);
  fData.ReadVariableSection (incase, nl, el, nvector, evector, true);
  
  // allocate space
  int dex;
  if (!fPartDimensions.ColumnHasValue (2, group_id, dex)) return;
  evalues.Allocate (fPartDimensions (dex, 1), el.Length());
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
      ifstream in (fGeometryFile);
      fData.ReadGeometryHeader (in, nm, em);
      while (fData.ReadPart (in, partID))
	{
	  fData.SkipPart (in, nm, em, numnodes, numelems, numelemnodes);

	  int length = fPartDimensions.MajorDim();
	  if (length == 0)
	    fPartDimensions.Allocate (1,3);
	  else
	    fPartDimensions.Resize(length + 1, 0);

	  fPartDimensions (length, 0) = numnodes;
	  fPartDimensions (length, 1) = numelems;
	  fPartDimensions (length, 2) = partID;
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

void EnSightInputT::ReadVariableData (ArrayT<bool>& vector, ArrayT<StringT>& labels, int group_id, dArray2DT& values, int currentinc, bool nodal) const
{
  for (int i=0; i < labels.Length(); i++)
    {
      StringT filename = CreateVariableFile (labels[i], currentinc);
      ifstream in (filename);
      
      if (!in)
	{
	  cout << "\n\nEnSightInputT::ReadVariableData Unable to open: "
	       << filename << endl;
	}

      // read variable file header
      StringT header;
      fData.ReadVariableHeader (in, header);
      
      // search file for part
      bool found = false;
      int id, dex, num;
      dArray2DT temp;
      while (fData.ReadPart (in, id) && !found)
	{
	  if (!fPartDimensions.ColumnHasValue (2, id, dex)) throw eGeneralFail;
	  
	  if (nodal)
	    num = fPartDimensions (dex, 0);
	  else
	    num = fPartDimensions (dex, 1);

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
