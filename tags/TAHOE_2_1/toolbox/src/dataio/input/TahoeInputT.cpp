/* $Id: TahoeInputT.cpp,v 1.14 2004-04-15 23:55:41 paklein Exp $ */
/* created: sawimme July 2001 */
#include "TahoeInputT.h"

using namespace Tahoe;

TahoeInputT::TahoeInputT (ostream& out) :
  InputBaseT (out),
  fModel ()
{
}

bool TahoeInputT::Open (const StringT& file)
{
  if (fModel.OpenRead (file) == ModelFileT::kFail)
    {
      fout << "\n\nTahoeInputT::Open unable to open file: ";
      fout << file << "\n\n";
      cout << "\n\nTahoeInputT::Open unable to open file: ";
      cout << file << "\n\n";
      return false;
    } else return true;
}

void TahoeInputT::Close (void)
{
  fModel.Close ();
}

void TahoeInputT::ElementGroupNames (ArrayT<StringT>& groupnames) const
{
	const char caller[] = "TahoeInputT::ElementGroupNames";
  if (groupnames.Length() != NumElementGroups()) ExceptionT::SizeMismatch(caller);
  iArrayT ids (groupnames.Length());
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail) ExceptionT::DatabaseFail(caller);
  for (int i=0; i < ids.Length(); i++) {
  	groupnames[i].Clear();
    groupnames[i].Append(ids[i]);
    }
}

void TahoeInputT::SideSetNames (ArrayT<StringT>& sidenames) const
{
	const char caller[] = "TahoeInputT::SideSetNames";
  if (sidenames.Length() != NumSideSets()) ExceptionT::SizeMismatch(caller);
  iArrayT nums (sidenames.Length());
  if (fModel.GetSideSetID (nums) == ModelFileT::kFail) ExceptionT::DatabaseFail(caller);
  for (int i=0; i < nums.Length(); i++) {
  	sidenames[i].Clear();
    sidenames[i].Append(nums[i]);
    }
}

void TahoeInputT::NodeSetNames (ArrayT<StringT>& nodenames) const
{
	const char caller[] = "TahoeInputT::NodeSetNames";
  if (nodenames.Length() != NumNodeSets()) ExceptionT::SizeMismatch(caller);
  iArrayT nums (nodenames.Length());
  if (fModel.GetNodeSetID (nums) == ModelFileT::kFail)
    ExceptionT::DatabaseFail(caller);
  for (int i=0; i < nums.Length(); i++) {
  	nodenames[i].Clear();
    nodenames[i].Append(nums[i]);
    }
}

int TahoeInputT::NumElementGroups (void) const
{
  iArrayT ids;
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;
  return ids.Length();
}

int TahoeInputT::NumSideSets (void) const
{
  iArrayT ids;
  if (fModel.GetSideSetID (ids) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;
  return ids.Length();
}

int TahoeInputT::NumNodeSets (void) const
{
  iArrayT ids;
  if (fModel.GetNodeSetID (ids) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;
  return ids.Length();
}

int TahoeInputT::NumNodes (void) const
{
  int numnodes, dims;
  if (fModel.GetDimensions (numnodes, dims) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;
  return numnodes;
}

int TahoeInputT::NumDimensions (void) const
{
  int numnodes, dims;
  if (fModel.GetDimensions (numnodes, dims) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;
  return dims;
}

void TahoeInputT::ReadNodeID(iArrayT& node_id)
{
	if (node_id.Length() != NumNodes()) throw ExceptionT::kSizeMismatch;
	node_id.SetValueToPosition ();
	node_id++;
}

void TahoeInputT::ReadCoordinates (dArray2DT& coords)
{
	const char caller[] = "TahoeInputT::ReadCoordinates";
	if (coords.MajorDim() != NumNodes() ||
        coords.MinorDim() != NumDimensions()) 
        ExceptionT::SizeMismatch(caller);

	if (fModel.GetCoordinates (coords) == ModelFileT::kFail) 
		ExceptionT::DatabaseFail(caller);
}

void TahoeInputT::ReadCoordinates (dArray2DT& coords, iArrayT& node_id)
{
	ReadCoordinates(coords);
	ReadNodeID(node_id);
}

int TahoeInputT::NumGlobalElements (void) const
{
  int numelems = 0;
  iArrayT ids;
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;
  int nume, dims;
  for (int i=0; i < ids.Length(); i++)
    {
      if (fModel.GetElementSetDimensions (ids[i], nume, dims) == ModelFileT::kFail)
	throw ExceptionT::kDatabaseFail;
      numelems += nume;
    }
  return numelems;
}

int TahoeInputT::NumElements (const StringT& name)
{
  int ID = atoi (name.Pointer());
  int numelems, dims;
  if (fModel.GetElementSetDimensions (ID, numelems, dims) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;
  return numelems;
}

int TahoeInputT::NumElementNodes (const StringT& name)
{
  int ID = atoi (name.Pointer());
  int numelems, dims;
  if (fModel.GetElementSetDimensions (ID, numelems, dims) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;
  return dims;
}

void TahoeInputT::ReadAllElementMap (iArrayT& elemmap)
{
  if (elemmap.Length() != NumGlobalElements()) throw ExceptionT::kSizeMismatch;
  elemmap.SetValueToPosition ();
  elemmap += 1;
}

void TahoeInputT::ReadGlobalElementMap (const StringT& name, iArrayT& elemmap)
{
  int ID = atoi (name.Pointer());
  if (elemmap.Length() != NumElements (name)) throw ExceptionT::kSizeMismatch;

  int numelems = 0;
  iArrayT ids;
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;
  int nume = 0, dims;
  for (int i=0; i < ids.Length(); i++)
    {
      numelems += nume;
      if (fModel.GetElementSetDimensions (ids[i], nume, dims) == ModelFileT::kFail) throw ExceptionT::kDatabaseFail;
      if (ids[i] == ID) break;
    }

  elemmap.SetValueToPosition ();
  elemmap += 1 + numelems;
}

void TahoeInputT::ReadGlobalElementSet (const StringT& name, iArrayT& set)
{
  ReadGlobalElementMap (name, set);
  set += -1;
}

void TahoeInputT::ReadConnectivity (const StringT& name, iArray2DT& connects)
{
  int ID = atoi (name.Pointer());
  if (fModel.GetElementSet (ID, connects) == ModelFileT::kFail) 
    throw ExceptionT::kDatabaseFail;

  connects += -1;
}

void TahoeInputT::ReadGeometryCode (const StringT& name, GeometryT::CodeT& code)
{
  int ID = atoi (name.Pointer());
  int length, numelemnodes;
  int numnodes, dims;
  if (fModel.GetDimensions (numnodes, dims) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;
  if (fModel.GetElementSetDimensions (ID, length, numelemnodes) == ModelFileT::kFail) 
    throw ExceptionT::kDatabaseFail;
  SetCode (numelemnodes, dims, code);
}

int TahoeInputT::NumNodesInSet (const StringT& name)
{
  int id = atoi (name.Pointer());
  int num;
  if (fModel.GetNodeSetDimensions (id, num) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;
  return num;
}

void TahoeInputT::ReadNodeSet (const StringT& name, iArrayT& nodes)
{
  int id = atoi (name.Pointer());
  if (fModel.GetNodeSet (id, nodes) == ModelFileT::kFail) 
    throw ExceptionT::kDatabaseFail;

  nodes += -1;
}

int TahoeInputT::NumSidesInSet (const StringT& name) const
{
  int id = atoi (name.Pointer());
  int num;
  if (fModel.GetSideSetDimensions (id, num) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;
  return num;
}

StringT TahoeInputT::SideSetGroupName (const StringT& name) const
{
  int id = atoi (name.Pointer());
  int elsetid;
  iArray2DT sides (NumSidesInSet(name), 2);
  if (fModel.GetSideSet (id, elsetid, sides) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;

  StringT elname;
  elname.Append (elsetid);
  return elname;
}

void TahoeInputT::ReadSideSetLocal (const StringT& name, iArray2DT& sides) const
{
  int id = atoi (name.Pointer());
  int elsetid;
  if (fModel.GetSideSet (id, elsetid, sides) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;

  sides += -1;
}

void TahoeInputT::ReadSideSetGlobal (const StringT& name, iArray2DT& sides) const
{
  int id = atoi (name.Pointer());
  int elsetid;
  if (fModel.GetSideSet (id, elsetid, sides) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;

  iArrayT ids;
  if (fModel.GetElementSetID (ids) == ModelFileT::kFail)
    throw ExceptionT::kDatabaseFail;

  int num_elems, dim, offset = 0;
  for (int i=0; i < ids.Length(); i++)
    {
      if (ids[i] == elsetid) break;
      if (fModel.GetElementSetDimensions (ids[i], num_elems, dim) == ModelFileT::kFail)
	throw ExceptionT::kDatabaseFail;
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
	if (dof == 1 && numelemnodes == 1)
		code = GeometryT::kPoint;
	else if (dof == 1 && numelemnodes == 2)
		code = GeometryT::kLine;	
	else if (dof == 2)
    switch (numelemnodes)
      {
      case 1: code = GeometryT::kPoint; break;
      case 6: case 3: code = GeometryT::kTriangle; break;
      case 8: case 4: code = GeometryT::kQuadrilateral; break;
      }
  else if (dof == 3)
    switch (numelemnodes)
      {
      case 1: code = GeometryT::kPoint; break;
      case 4: case 10: code = GeometryT::kTetrahedron; break;
      case 8: case 20: code = GeometryT::kHexahedron; break;
      case 6: case 15: code = GeometryT::kPentahedron; break;
      }
}
