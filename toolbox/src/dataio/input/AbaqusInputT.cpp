/* $Id: AbaqusInputT.cpp,v 1.3 2001-08-08 11:47:12 sawimme Exp $ */
/* created: sawimme (05/18/1998)                                          */

#include "AbaqusInputT.h"
#include "ios_fwd_decl.h"
#include <fstream.h>
#include "dArray2DT.h"
#include "dArrayT.h"
#include "iAutoArrayT.h"

AbaqusInputT::AbaqusInputT (ostream& out, bool binary) :
InputBaseT (out),
	fData (out, binary),
	fNumElements (0),
	fNumNodes (0),
	fDimensions (0),
	fEnergyData (false)
{
}

void AbaqusInputT::Open (const StringT& file)
{
  fResultFile = file;
  fNumElements = 0;
  fNumNodes = 0;
  Initialize ();
  if (fNumElements > 0) ReadElementSets ();
  if (fNumNodes > 0) ReadNodeSets ();
}

void AbaqusInputT::Close (void)
{
  fNumElements = 0;
  fNumNodes = 0;
  fDimensions = 0;

  fElementData.Free ();
  fCoordinateData.Free();
  fElementSetNames.Free ();
  fNodeSetNames.Free ();
}

void AbaqusInputT::ElementGroupNames (ArrayT<StringT>& groupnames)
{
  if (groupnames.Length() != fElementSetNames.Length()) throw eSizeMismatch;
  for (int k=0; k < groupnames.Length(); k++)
    groupnames[k] = fElementSetNames [k];
}

void AbaqusInputT::NodeSetNames (ArrayT<StringT>& nodenames)
{
  if (nodenames.Length() != fNodeSetNames.Length()) throw eSizeMismatch;
  for (int k=0; k < nodenames.Length(); k++)
    nodenames[k] = fNodeSetNames[k];
}

void AbaqusInputT::ReadNodeMap (iArrayT& nodemap)
{
  if (nodemap.Length() != fNumNodes) throw eSizeMismatch;
  for (int n=0; n < nodemap.Length(); n++)
    nodemap[n] = fCoordinateData (n, 0);
}

void AbaqusInputT::ReadCoordinates (dArray2DT& coords)
{
  if (coords.MajorDim() != fNumNodes || coords.MinorDim() != fDimensions) throw eSizeMismatch;
  dArrayT x;
  int ID;
  ifstream in;
  if (!OpenFile (in)) throw eGeneralFail;
  for (int i=0; i < fNumNodes; i++)
    {
      fData.ReadCoordinate (in, ID, x);
      coords.SetRow (i, x);
    }
}

void AbaqusInputT::ReadCoordinates (dArray2DT& coords, iArrayT& nodemap)
{
  if (coords.MajorDim() != fNumNodes || coords.MinorDim() != fDimensions ||
      nodemap.Length() != fNumNodes ) throw eSizeMismatch;

  dArrayT x;
  int ID;
  ifstream in;
  if (!OpenFile (in)) throw eGeneralFail;
  for (int i=0; i < fNumNodes; i++)
    {
      fData.ReadCoordinate (in, ID, x);
      coords.SetRow (i, x);
      nodemap [i] = ID;
    }
}

int AbaqusInputT::NumElements (StringT& name)
{
  int group = 0;
  for (int k=0; k < fElementSetNames.Length() && group == 0; k++)
    if (fElementSetNames [k] == name)
      group = k+1;

  int count = 0;
  int numelemnodes;
  for (int i=0, j=0; i < fNumElements; i++)
    if (fElementData (i, 2) == group)
      {
	count++;
	if (j++ == 0) numelemnodes = fElementData (i, 3);
      }
  return count;
}

int AbaqusInputT::NumElementNodes (StringT& name)
{
  int group = 0;
  for (int k=0; k < fElementSetNames.Length() && group == 0; k++)
    if (fElementSetNames [k] == name)
      group = k+1;

  for (int i=0; i < fNumElements; i++)
    if (fElementData (i, 2) == group)
      return fElementData (i, 3);
  return -1;
}

void AbaqusInputT::ReadAllElementMap (iArrayT& elemmap)
{
  if (elemmap.Length() != fNumElements) throw eSizeMismatch;
  for (int i=0; i < fNumElements; i++)
    elemmap[i] = fElementData (i, 0);
}

void AbaqusInputT::ReadGlobalElementMap (StringT& name, iArrayT& elemmap)
{
  if (elemmap.Length() != NumElements (name)) throw eSizeMismatch;

  int group = 0;
  for (int k=0; k < fElementSetNames.Length() && group == 0; k++)
    if (fElementSetNames [k] == name)
      group = k+1;

  int count = 0;
  for (int i=0; i < fNumElements; i++)
    if (fElementData (i, 2) == group)
      elemmap [count++] = fElementData (i, 0);
}

void AbaqusInputT::ReadConnectivity (StringT& name, iArray2DT& connects)
{
  int group = 0;
  for (int k=0; k < fElementSetNames.Length() && group == 0; k++)
    if (fElementSetNames [k] == name)
      group = k+1;

  int count = 0;
  int numelemnodes;
  for (int i=0, j=0; i < fNumElements; i++)
    if (fElementData (i, 2) == group)
      {
	count++;
	if (j++ == 0) numelemnodes = fElementData (i, 3);
      }
  
  int ID;
  GeometryT::CodeT geocode, code;
  iArrayT nodes, internal;
  ifstream in;
  if (!OpenFile (in)) throw eGeneralFail;
  for (int k=0, m=0; k < fNumElements; k++)
    if (fElementData (k,2) == group)
      {
	fData.ReadElement (in, ID, code, nodes);
	
	// change to consecutive numbering, offset from zero
	internal.Allocate (nodes.Length());
	int dex;
	for (int n=0; n < nodes.Length(); n++)
	  {
	    if (!fCoordinateData.ColumnHasValue (0, nodes[n], dex))
	      throw eGeneralFail;
	    internal[n] = dex;
	  }
	if (m == 0) geocode = code;
	else if (geocode != code)
	  {
	    cout << "\n\nAbaqusInputT::ReadConnectivity, mismatched geo. code\n";
	    throw eGeneralFail;
	  }

	connects.SetRow (m, internal);
	m++;
      }
}

void AbaqusInputT::ReadGeometryCode (StringT& name, GeometryT::CodeT& geocode)
{
  geocode = GeometryT::kNone;
  int group = 0;
  for (int k=0; k < fElementSetNames.Length() && group == 0; k++)
    if (fElementSetNames [k] == name)
      group = k+1;

  for (int i=0; i < fNumElements; i++)
    if (fElementData (i, 2) == group)
      {
	switch (fElementData (i,1))
	  {
	  case GeometryT::kPoint:
	    geocode = GeometryT::kPoint;
	    break;
	  case GeometryT::kLine:
	    geocode = GeometryT::kLine;
	    break;
	  case GeometryT::kQuadrilateral:
	    geocode = GeometryT::kQuadrilateral;
	    break;
	  case GeometryT::kTriangle:
	    geocode = GeometryT::kTriangle;
	    break;
	  case GeometryT::kHexahedron:
	    geocode = GeometryT::kHexahedron;
	    break;
	  case GeometryT::kTetrahedron:
	    geocode = GeometryT::kTetrahedron;
	    break;
	  case GeometryT::kPentahedron:
	    geocode = GeometryT::kPentahedron;
	    break;
	  }
	return;
      }
}

int AbaqusInputT::NumNodesInSet (StringT& name)
{
  int set_num = 0;
  for (int k=0; k < fNodeSetNames.Length() && set_num == 0; k++)
    if (fNodeSetNames [k] == name)
      set_num = k+1;

  int count = 0;
  for (int i=0; i < fNumNodes; i++)
    if (fCoordinateData (i,1) == set_num) count++;
  
  return count;
}

void AbaqusInputT::ReadNodeSet (StringT& name, iArrayT& nodes)
{
  int group = 0;
  for (int k=0; k < fNodeSetNames.Length() && group == 0; k++)
    if (fNodeSetNames [k] == name)
      group = k+1;

  if (nodes.Length() != NumNodesInSet (name)) throw eSizeMismatch;
  
  int *pn = nodes.Pointer();
  for (int j=0; j < fNumNodes; j++)
    if (fCoordinateData (j,1) == group)
      *pn++ = j; //change to global consecutive numbering, offset to zero
}

void AbaqusInputT::QARecords (ArrayT<StringT>& records)
{
  records.Allocate (4);
  records[0] = "Abaqus";
  records[1] = fVersion;
  records[2] = fDate;
  records[3] = fTime;
}

int AbaqusInputT::NumTimeSteps (void)
{
  ifstream in;
  if (!OpenFile (in)) throw eGeneralFail;

  int count = 0;
  double time;
  while (fData.NextTimeIncrement (in, time))
    count ++;

  return count;
}

void AbaqusInputT::ReadTimeSteps (dArrayT& steps)
{
  ifstream in;
  if (!OpenFile (in)) throw eGeneralFail;

  double time;
  AutoArrayT<double> temp;
  while (fData.NextTimeIncrement (in, time))
    temp.Append (time);
  
  steps.Allocate (temp.Length());
  steps.CopyPart (0, temp, 0, temp.Length());
  
  // check time steps to see if duplicate time value is used
  for (int i=1; i < steps.Length(); i++)
    if (steps[i] == steps[i-1])
      {
	double shift;
	if (i+1 != steps.Length())
	  shift = (steps[i+1]-steps[i])/1.e-5;
	else
	  shift = steps[i]/1.e-5;
	steps[i] += shift;
      }
}

int AbaqusInputT::NumNodeVariables (void)
{
  ifstream in;
  if (!OpenFile (in)) throw eGeneralFail;
  double time;
  if (!fData.NextTimeIncrement (in, time)) throw eGeneralFail;
  ArrayT<AbaqusT::VariableKeyT> nkeys, ekeys;
  fData.ReadLabels (in, nkeys, ekeys, fEnergyData);
  return nkeys.Length();
}

int AbaqusInputT::NumElementVariables (void)
{
  ifstream in;
  if (!OpenFile (in)) throw eGeneralFail;
  double time;
  if (!fData.NextTimeIncrement (in, time)) throw eGeneralFail;
  ArrayT<AbaqusT::VariableKeyT> nkeys, ekeys;
  fData.ReadLabels (in, nkeys, ekeys, fEnergyData);
  return ekeys.Length();
}

void AbaqusInputT::ReadElementLabels (StringT& name, ArrayT<StringT>& elabels)
{
#pragma unused (name)
  ifstream in;
  if (!OpenFile (in)) throw eGeneralFail;
  double time;
  if (!fData.NextTimeIncrement (in, time)) throw eGeneralFail;
  ArrayT<AbaqusT::VariableKeyT> nkeys, ekeys;
  fData.ReadLabels (in, nkeys, ekeys, fEnergyData);
  elabels.Allocate (ekeys.Length());
  
  StringT newlabel;
  char incrementor = '1';
  int unknown = 0;
  for (int e=0; e < ekeys.Length(); e++)
    {
      if (e == 0 || ekeys[e] != ekeys[e-1]) incrementor = '1';
      SetLabelName (ekeys[e], unknown, elabels[e], incrementor);
      incrementor++;
    }
}

void AbaqusInputT::ReadElementVariables (int step, StringT& name, dArray2DT& evalues)
{
  int group = 0;
  for (int k=0; k < fElementSetNames.Length() && group == 0; k++)
    if (fElementSetNames [k] == name)
      group = k+1;

  // determine number of variables
  ifstream in;
  if (!OpenFile (in)) throw eGeneralFail;
  ArrayT<AbaqusT::VariableKeyT> nkeys, ekeys;
  double time;
  for (int i=0; i < step + 1; i++)
    if (!fData.NextTimeIncrement (in, time)) throw eGeneralFail;
  fData.ReadLabels (in, nkeys, ekeys, fEnergyData);
  
  // read elemental variables
  if (ekeys.Length() > 0)
    {
      if (!OpenFile (in)) throw eGeneralFail; // rewind file
      for (int i2=0; i2 < step + 1; i2++)
	if (!fData.NextTimeIncrement (in, time)) throw eGeneralFail;
      
      iAutoArrayT elems;
      for (int i=0, j=0; i < fNumElements; i++)
	if (fElementData (i, 2) == group)
	  elems.Append (fElementData (i, 0));
      
      evalues.Allocate (elems.Length(), ekeys.Length());
      iArrayT els (elems.Length());
      els.Set (elems.Length(), elems.Pointer());
      fData.ReadElementVariables (in, ekeys, evalues, els);
    }
}

/*************************************************************************
*
* Private
*
*************************************************************************/

bool AbaqusInputT::OpenFile (ifstream& in)
{
  if (in) in.close();
  fData.ResetBufferSize ();
  in.open (fResultFile);
  if (in) return true;
  else
    {
      cout << "\n\nAbaqusInputT::OpenFile unable to open file: "
	   << fResultFile << endl;
      return false;
    }
}

void AbaqusInputT::Initialize (void)
{
  ifstream in;
  if (!OpenFile (in))
    {
      cout << "\n\nAbaqusInputT, unable to open file: " << fResultFile << endl;
      throw eGeneralFail;
    }
  
  // read num elements, num dimensions, etc.
  fData.ReadVersion (in, fVersion, fDate, fTime, fNumElements, fNumNodes);
  
  // set up element dat
  if (!OpenFile (in)) throw eGeneralFail; // rewind file
  int ID;
  GeometryT::CodeT geocode;
  iArrayT nodes;
  fElementData.Allocate (fNumElements, 4);
  fElementData = -1;
  int *p = fElementData.Pointer();
  for (int e=0; e < fNumElements; e++)
    {
      fData.ReadElement (in, ID, geocode, nodes);
      *p++ = ID;
      *p++ = geocode;
      *p++;
      *p++ = nodes.Length();
    }
  
  // set up coord data
  if (!OpenFile (in)) throw eGeneralFail; // rewind file
  fCoordinateData.Allocate (fNumNodes, 2);
  fCoordinateData = -1;
  dArrayT x;
  p = fCoordinateData.Pointer();
  for (int n=0; n < fNumNodes; n++)
    {
      fData.ReadCoordinate (in, ID, x);
      *p++ = ID;
      *p++;
    }
  fDimensions = x.Length();
  
  // determine if modal analysis
  if (!OpenFile (in)) throw eGeneralFail; // rewind file
  fData.IsModal (in);
}

void AbaqusInputT::ReadElementSets (void)
{
  ifstream (in);
  if (!OpenFile (in)) throw eGeneralFail;
  fData.ReadElementSets (in, fElementData, fElementSetNames);
}

void AbaqusInputT::ReadNodeSets (void)
{
  ifstream (in);
  if (!OpenFile (in)) throw eGeneralFail;
  fData.ReadNodeSets (in, fCoordinateData, fNodeSetNames);
}

void AbaqusInputT::SetLabelName (AbaqusT::VariableKeyT key, int& unknown, StringT& name, char incrementor) const
{
  switch (key)
    {
    case AbaqusT::aSDV:           name = "SDV_"; break;
    case AbaqusT::aStress:        name = "S_";   break;
    case AbaqusT::aTotalStrain:   name = "E_";   break;
    case AbaqusT::aDisplacement:  name = "D_";   break;
    case AbaqusT::aVelocity:      name = "V_";   break;
    case AbaqusT::aAcceleration:  name = "A_";   break;
    case AbaqusT::aCoordVariable: name = "X_";   break;
    case AbaqusT::aPrinStress:    name = "PS_";  break;
    default:
      name = "unk";
      name.Append (unknown++);
      name.Append ("_");
    }
  name.Append (incrementor);
}

