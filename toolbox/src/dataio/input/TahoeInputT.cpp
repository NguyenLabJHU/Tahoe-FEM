/* $Id: TahoeInputT.cpp,v 1.5 2001-12-16 23:53:45 paklein Exp $ */
/* created: sawimme July 2001 */

#include "TahoeInputT.h"

TahoeInputT::TahoeInputT (ostream& out) :
  InputBaseT (out),
  fModel ()
{
}

void TahoeInputT::Open (const StringT& file)
{
  if (fModel.OpenRead (file) == ModelFileT::kFail)
    {
      fout << "\n\nTahoeInputT::Open unable to open file: ";
      fout << file << "\n\n";
      cout << "\n\nTahoeInputT::Open unable to open file: ";
      cout << file << "\n\n";
      throw eDatabaseFail;
    }
}

void TahoeInputT::Close (void)
{
  fModel.Close ();
}

void TahoeInputT::ElementGroupNames (ArrayT<StringT>& groupnames) const
{
  if (groupnames.Length() != NumElementGroups()) throw eSizeMismatch;
  iArrayT ids (groupnames.Length());
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw eDatabaseFail;
  for (int i=0; i < ids.Length(); i++)
    groupnames[i].Append (ids[i]);
}

void TahoeInputT::SideSetNames (ArrayT<StringT>& sidenames) const
{
  if (sidenames.Length() != NumSideSets()) throw eSizeMismatch;
  iArrayT nums (sidenames.Length());
  if (fModel.GetSideSetID (nums) == ModelFileT::kFail)
    throw eDatabaseFail;
  for (int i=0; i < nums.Length(); i++)
    sidenames[i].Append (nums[i]);
}

void TahoeInputT::NodeSetNames (ArrayT<StringT>& nodenames) const
{
  if (nodenames.Length() != NumNodeSets()) throw eSizeMismatch;
  iArrayT nums (nodenames.Length());
  if (fModel.GetNodeSetID (nums) == ModelFileT::kFail)
    throw eDatabaseFail;
  for (int i=0; i < nums.Length(); i++)
    nodenames[i].Append (nums[i]);
}

int TahoeInputT::NumElementGroups (void) const
{
  iArrayT ids;
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw eDatabaseFail;
  return ids.Length();
}

int TahoeInputT::NumSideSets (void) const
{
  iArrayT ids;
  if (fModel.GetSideSetID (ids) == ModelFileT::kFail)
    throw eDatabaseFail;
  return ids.Length();
}

int TahoeInputT::NumNodeSets (void) const
{
  iArrayT ids;
  if (fModel.GetNodeSetID (ids) == ModelFileT::kFail)
    throw eDatabaseFail;
  return ids.Length();
}

int TahoeInputT::NumNodes (void) const
{
  int numnodes, dims;
  if (fModel.GetDimensions (numnodes, dims) == ModelFileT::kFail)
    throw eDatabaseFail;
  return numnodes;
}

int TahoeInputT::NumDimensions (void) const
{
  int numnodes, dims;
  if (fModel.GetDimensions (numnodes, dims) == ModelFileT::kFail)
    throw eDatabaseFail;
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
    throw eDatabaseFail;
}

void TahoeInputT::ReadCoordinates (dArray2DT& coords, iArrayT& nodemap)
{
  ReadCoordinates (coords);
  ReadNodeMap (nodemap);
}

int TahoeInputT::NumGlobalElements (void) const
{
  int numelems = 0;
  iArrayT ids;
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw eDatabaseFail;
  int nume, dims;
  for (int i=0; i < ids.Length(); i++)
    {
      if (fModel.GetElementSetDimensions (ids[i], nume, dims) == ModelFileT::kFail)
	throw eDatabaseFail;
      numelems += nume;
    }
  return numelems;
}

int TahoeInputT::NumElements (StringT& name)
{
  int ID = atoi (name.Pointer());
  int numelems, dims;
  if (fModel.GetElementSetDimensions (ID, numelems, dims) == ModelFileT::kFail)
    throw eDatabaseFail;
  return numelems;
}

int TahoeInputT::NumElementNodes (StringT& name)
{
  int ID = atoi (name.Pointer());
  int numelems, dims;
  if (fModel.GetElementSetDimensions (ID, numelems, dims) == ModelFileT::kFail)
    throw eDatabaseFail;
  return dims;
}

void TahoeInputT::ReadAllElementMap (iArrayT& elemmap)
{
  if (elemmap.Length() != NumGlobalElements()) throw eSizeMismatch;
  elemmap.SetValueToPosition ();
  elemmap += 1;
}

void TahoeInputT::ReadGlobalElementMap (StringT& name, iArrayT& elemmap)
{
  int ID = atoi (name.Pointer());
  if (elemmap.Length() != NumElements (name)) throw eSizeMismatch;

  int numelems = 0;
  iArrayT ids;
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw eDatabaseFail;
  int nume, dims;
  for (int i=0; i < ids.Length(); i++)
    {
      if (fModel.GetElementSetDimensions (ids[i], nume, dims) == ModelFileT::kFail)
	throw eDatabaseFail;
      numelems += nume;
      if (ids[i] == ID) break;
    }

  elemmap.SetValueToPosition ();
  elemmap += 1 + numelems;
}

void TahoeInputT::ReadGlobalElementSet (StringT& name, iArrayT& set)
{
  ReadGlobalElementMap (name, set);
  set += -1;
}

void TahoeInputT::ReadConnectivity (StringT& name, iArray2DT& connects)
{
  int ID = atoi (name.Pointer());
  if (fModel.GetElementSet (ID, connects) == ModelFileT::kFail) 
    throw eDatabaseFail;

  connects += -1;
}

void TahoeInputT::ReadGeometryCode (StringT& name, GeometryT::CodeT& code)
{
  int ID = atoi (name.Pointer());
  int length, numelemnodes;
  int numnodes, dims;
  if (fModel.GetDimensions (numnodes, dims) == ModelFileT::kFail)
    throw eDatabaseFail;
  if (fModel.GetElementSetDimensions (ID, length, numelemnodes) == ModelFileT::kFail) 
    throw eDatabaseFail;
  SetCode (numelemnodes, dims, code);
}

int TahoeInputT::NumNodesInSet (StringT& name)
{
  int id = atoi (name.Pointer());
  int num;
  if (fModel.GetNodeSetDimensions (id, num) == ModelFileT::kFail)
    throw eDatabaseFail;
  return num;
}

void TahoeInputT::ReadNodeSet (StringT& name, iArrayT& nodes)
{
  int id = atoi (name.Pointer());
  if (fModel.GetNodeSet (id, nodes) == ModelFileT::kFail) 
    throw eDatabaseFail;

  nodes += -1;
}

int TahoeInputT::NumSidesInSet (StringT& name) const
{
  int id = atoi (name.Pointer());
  int num;
  if (fModel.GetSideSetDimensions (id, num) == ModelFileT::kFail)
    throw eDatabaseFail;
  return num;
}

StringT TahoeInputT::SideSetGroupName (StringT& name) const
{
  int id = atoi (name.Pointer());
  int elsetid;
  iArray2DT sides (NumSidesInSet(name), 2);
  if (fModel.GetSideSet (id, elsetid, sides) == ModelFileT::kFail)
    throw eDatabaseFail;

  StringT elname;
  elname.Append (elsetid);
  return elname;
}

void TahoeInputT::ReadSideSetLocal (StringT& name, iArray2DT& sides) const
{
  int id = atoi (name.Pointer());
  int elsetid;
  if (fModel.GetSideSet (id, elsetid, sides) == ModelFileT::kFail)
    throw eDatabaseFail;

  sides += -1;
}

void TahoeInputT::ReadSideSetGlobal (StringT& name, iArray2DT& sides) const
{
  int id = atoi (name.Pointer());
  int elsetid;
  if (fModel.GetSideSet (id, elsetid, sides) == ModelFileT::kFail)
    throw eDatabaseFail;

  iArrayT ids;
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw eDatabaseFail;

  int num_elems, dim, offset = 0;
  for (int i=0; i < ids.Length(); i++)
    {
      if (ids[i] == elsetid) break;
      if (fModel.GetElementSetDimensions (ids[i], num_elems, dim) == ModelFileT::kFail)
	throw eDatabaseFail;
      offset += num_elems;
    }

  int *pelem = sides.Pointer();
  for (int j=0; j < sides.MajorDim(); j++, pelem += 2)
    *pelem += offset;

  sides += -1;
}

/******************* PRIVATE ********************/

void TahoeInputT::SetCode (int numelemnodes, int dof, GeometryT::CodeT& code) const
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
