/* $Id: TahoeInputT.cpp,v 1.2 2001-08-07 23:11:55 paklein Exp $ */
/* created: sawimme July 2001 */

#include "TahoeInputT.h"

TahoeInputT::TahoeInputT (ostream& out) :
  InputBaseT (out),
  fModel ()
{
}

void TahoeInputT::Open (const StringT& file)
{
  fModel.OpenRead (file);
}

void TahoeInputT::Close (void)
{
  fModel.Close ();
}

int TahoeInputT::NumElementGroups (void)
{
  iArrayT ids;
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw eGeneralFail;
  return ids.Length();
}

int TahoeInputT::NumSideSets (void)
{
  iArrayT ids;
  if (fModel.GetSideSetID (ids) == ModelFileT::kFail)
    throw eGeneralFail;
  return ids.Length();
}

int TahoeInputT::NumNodeSets (void)
{
  iArrayT ids;
  if (fModel.GetNodeSetID (ids) == ModelFileT::kFail)
    throw eGeneralFail;
  return ids.Length();
}

int TahoeInputT::NumNodes (void)
{
  int numnodes, dims;
  if (fModel.GetDimensions (numnodes, dims) == ModelFileT::kFail)
    throw eGeneralFail;
  return numnodes;
}

int TahoeInputT::NumDimensions (void)
{
  int numnodes, dims;
  if (fModel.GetDimensions (numnodes, dims) == ModelFileT::kFail)
    throw eGeneralFail;
  return dims;
}

void TahoeInputT::ReadNodeMap (iArrayT& nodemap)
{
  if (nodemap.Length() != NumNodes()) throw eSizeMismatch;
  nodemap.SetValueToPosition ();
}

void TahoeInputT::ReadCoordinates (dArray2DT& coords)
{
  if (coords.MajorDim() != NumNodes() ||
      coords.MinorDim() != NumDimensions()) throw eSizeMismatch;
  if (fModel.GetCoordinates (coords) == ModelFileT::kFail) 
    throw eGeneralFail;
}

void TahoeInputT::ReadCoordinates (dArray2DT& coords, iArrayT& nodemap)
{
  ReadCoordinates (coords);
  ReadNodeMap (nodemap);
}

int TahoeInputT::NumGlobalElements (void)
{
  int numelems = 0;
  iArrayT ids;
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw eGeneralFail;
  for (int i=0; i < ids.Length(); i++)
    numelems += NumElements_ID (ids[i]);
  return numelems;
}

void TahoeInputT::ReadAllElementMap (iArrayT& elemmap)
{
  if (elemmap.Length() != NumGlobalElements()) throw eSizeMismatch;
  elemmap.SetValueToPosition ();
  elemmap += 1;
}

/******************* PROTECTED ********************/

void TahoeInputT::ElementGroupIDs (iArrayT& groupnums)
{
  if (groupnums.Length() != NumElementGroups()) throw eSizeMismatch;
  if (fModel.GetElementSetID (groupnums) == ModelFileT::kFail)
    throw eGeneralFail;
}

void TahoeInputT::SideSetIDs (iArrayT& nums)
{
  if (nums.Length() != NumSideSets()) throw eSizeMismatch;
  if (fModel.GetSideSetID (nums) == ModelFileT::kFail)
    throw eGeneralFail;
}

void TahoeInputT::NodeSetIDs (iArrayT& nums)
{
  if (nums.Length() != NumNodeSets()) throw eSizeMismatch;
  if (fModel.GetNodeSetID (nums) == ModelFileT::kFail)
    throw eGeneralFail;
}

int TahoeInputT::NumElements_ID (int ID)
{
  int numelems, dims;
  if (fModel.GetElementSetDimensions (ID, numelems, dims) == ModelFileT::kFail)
    throw eGeneralFail;
  return numelems;
}

int TahoeInputT::NumElementNodes_ID (int ID)
{
  int numelems, dims;
  if (fModel.GetElementSetDimensions (ID, numelems, dims) == ModelFileT::kFail)
    throw eGeneralFail;
  return dims;
}

void TahoeInputT::ReadGlobalElementMap_ID (int ID, iArrayT& elemmap)
{
  if (elemmap.Length() != NumElements_ID (ID)) throw eSizeMismatch;

  int numelems = 0;
  iArrayT ids;
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw eGeneralFail;
  for (int i=0; i < ids.Length(); i++)
    {
      numelems += NumElements_ID (ids[i]);
      if (ids[i] == ID) break;
    }

  elemmap.SetValueToPosition ();
  elemmap += 1 + numelems;
}

void TahoeInputT::ReadConnectivity_ID (int ID, iArray2DT& connects)
{
  if (fModel.GetElementSet (ID, connects) == ModelFileT::kFail) 
    throw eGeneralFail;

  connects += -1;
}

void TahoeInputT::ReadGeometryCode_ID (int ID, GeometryT::CodeT& code)
{
  int length, numelemnodes;
  int numnodes, dims;
  if (fModel.GetDimensions (numnodes, dims) == ModelFileT::kFail)
    throw eGeneralFail;
  if (fModel.GetElementSetDimensions (ID, length, numelemnodes) == ModelFileT::kFail) 
    throw eGeneralFail;
  SetCode (numelemnodes, dims, code);
}

int TahoeInputT::NumNodesInSet_ID (int id)
{
  int num;
  if (fModel.GetNodeSetDimensions (id, num) == ModelFileT::kFail)
    throw eGeneralFail;
  return num;
}

void TahoeInputT::ReadNodeSet_ID (int id, iArrayT& nodes)
{
  if (fModel.GetNodeSet (id, nodes) == ModelFileT::kFail) 
    throw eGeneralFail;

  nodes += -1;
}

int TahoeInputT::NumSidesInSet_ID (int id)
{
  int num;
  if (fModel.GetSideSetDimensions (id, num) == ModelFileT::kFail)
    throw eGeneralFail;
  return num;
}

int TahoeInputT::SideSetGroupIndex_ID (int id)
{
  int elsetid;
  iArray2DT sides (NumSidesInSet_ID(id), 2);
  if (fModel.GetSideSet (id, elsetid, sides) == ModelFileT::kFail)
    throw eGeneralFail;
  return elsetid;
}

void TahoeInputT::ReadSideSetLocal_ID (int id, iArray2DT& sides)
{
  int elsetid;
  if (fModel.GetSideSet (id, elsetid, sides) == ModelFileT::kFail)
    throw eGeneralFail;
}

void TahoeInputT::ReadSideSetGlobal_ID (int id, iArray2DT& sides)
{
  int elsetid;
  if (fModel.GetSideSet (id, elsetid, sides) == ModelFileT::kFail)
    throw eGeneralFail;

  iArrayT ids;
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw eGeneralFail;

  int num_elems, dim, offset = 0;
  for (int i=0; i < ids.Length(); i++)
    {
      if (ids[i] == elsetid) break;
      fModel.GetElementSetDimensions (ids[i], num_elems, dim);
      offset += num_elems;
    }

  int *pelem = sides.Pointer();
  for (int j=0; j < sides.MajorDim(); j++, pelem += 2)
    *pelem += offset;

  sides += -1;
}

/******************* PRIVATE ********************/

void TahoeInputT::SetCode (int numelemnodes, int dof, GeometryT::CodeT& code)
{
  code = GeometryT::kNone;
  if (dof == 2)
    switch (numelemnodes)
      {
      case 6: case 3: code = GeometryT::kTriangle; break;
      case 8: case 4: code = GeometryT::kQuadrilateral; break;
      }
  else if (dof == 3)
    switch (numelemnodes)
      {
      case 4: case 10: code = GeometryT::kTetrahedron; break;
      case 8: case 20: code = GeometryT::kHexahedron; break;
      case 6: case 15: code = GeometryT::kPentahedron; break;
      }
}
