/* $Id: ExodusInputT.cpp,v 1.2 2001-08-03 19:16:43 sawimme Exp $ */
/* created: sawimme (12/04/1998)                                          */

#include "ExodusInputT.h"
#include "InputBaseT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "dArrayT.h"

ExodusInputT::ExodusInputT (ostream& out) :
  InputBaseT (out),
  fData (out)
{
}

void ExodusInputT::Open (const StringT& filename)
{
  fData.OpenRead (filename);
}

void ExodusInputT::ReadNodeMap (iArrayT& nodemap)
{
  if (nodemap.Length() != NumNodes()) throw eSizeMismatch;
  fData.ReadNodeMap (nodemap);
}

void ExodusInputT::ReadCoordinates (dArray2DT& coords)
{
  if (coords.MajorDim() != NumNodes() ||
      coords.MinorDim() != fData.NumDimensions()) throw eSizeMismatch;
  fData.ReadCoordinates (coords);
}

void ExodusInputT::ReadCoordinates (dArray2DT& coords, iArrayT& nodemap)
{
  ReadNodeMap (nodemap);
  ReadCoordinates (coords);
}

int ExodusInputT::NumGlobalElements (void)
{
  iArrayT ids (NumElementGroups());
  fData.ElementBlockID (ids);
  int count = 0;
  int numelems, numelemnodes;
  for (int i=0; i < ids.Length(); i++)
    {
      fData.ReadElementBlockDims (ids[i], numelems, numelemnodes);
      count += numelems;
    }
  return count;
}

void ExodusInputT::ReadAllElementMap (iArrayT& elemmap)
{
  if (elemmap.Length() != NumGlobalElements()) throw eSizeMismatch;
  elemmap.SetValueToPosition ();
  elemmap += 1;
}

void ExodusInputT::ReadTimeSteps (dArrayT& steps)
{
  if (steps.Length() != NumTimeSteps()) throw eSizeMismatch;
  for (int i=0; i < steps.Length(); i++)
    fData.ReadTime (i+1, steps[i]);
}

/*************************************************************************
* Protected
*************************************************************************/

void ExodusInputT::ElementGroupIDs (iArrayT& groupnums)
{
  if (groupnums.Length() != NumElementGroups()) throw eSizeMismatch;
  fData.ElementBlockID (groupnums);
}

void ExodusInputT::SideSetIDs (iArrayT& sidenums)
{
  if (sidenums.Length() != NumSideSets ()) throw eSizeMismatch;
  fData.SideSetID (sidenums);
}

void ExodusInputT::NodeSetIDs (iArrayT& nodenums)
{
  if (nodenums.Length() != NumNodeSets ()) throw eSizeMismatch;
  fData.NodeSetID (nodenums);
}

int ExodusInputT::NumElements_ID (int ID)
{
  int numelems, numelemnodes;
  fData.ReadElementBlockDims (ID, numelems, numelemnodes);
  return numelems;
}

int ExodusInputT::NumElementNodes_ID (int ID)
{
  int numelems, numelemnodes;
  fData.ReadElementBlockDims (ID, numelems, numelemnodes);
  return numelemnodes;
}

void ExodusInputT::ReadGlobalElementMap_ID (int id, iArrayT& elemmap)
{
  if (elemmap.Length() != NumElements_ID (id)) throw eSizeMismatch;

  int offset = 0;
  iArrayT eid (NumElementGroups());
  ElementGroupIDs (eid);
  for (int i=0; i < eid.Length(); i++)
    {
      if (eid[i] == id) break;
      offset += NumElements_ID (eid[i]);
    }

  elemmap.SetValueToPosition ();
  elemmap += offset + 1;
}

void ExodusInputT::ReadConnectivity_ID (int group, iArray2DT& connects)
{
  int numelems, numelemnodes;
  GeometryT::CodeT geocode;
  fData.ReadElementBlockDims (group, numelems, numelemnodes);

  if (connects.MajorDim() != numelems ||
      connects.MinorDim() != numelemnodes) throw eSizeMismatch;
  fData.ReadConnectivities (group, geocode, connects);

  connects += -1;
}

void ExodusInputT::ReadGeometryCode_ID (int group, GeometryT::CodeT& code)
{
  int numelems, numelemnodes;
  fData.ReadElementBlockDims (group, numelems, numelemnodes);
  iArray2DT connects (numelems, numelemnodes);
  fData.ReadConnectivities (group, code, connects);
}

void ExodusInputT::ReadNodeSet_ID (int set_num, iArrayT& nodes)
{
  if (nodes.Length() != fData.NumNodesInSet(set_num)) throw eSizeMismatch;
  fData.ReadNodeSet (set_num, nodes);
  
  nodes += -1;
}

int ExodusInputT::SideSetGroupIndex_ID (int setnum)
{
  int block_ID;
  iArray2DT sides (fData.NumSidesInSet (setnum), 2);
  fData.ReadSideSet (setnum, block_ID, sides);
  return block_ID;
}

void ExodusInputT::ReadSideSetLocal_ID (int set_num, iArray2DT& sides)
{
  if (sides.MajorDim() != fData.NumSidesInSet (set_num) ||
      sides.MinorDim() != 2) throw eSizeMismatch;
  int block_ID;
  fData.ReadSideSet (set_num, block_ID, sides);
  
  sides += -1;
}

void ExodusInputT::ReadSideSetGlobal_ID (int set_num, iArray2DT& sides)
{
  if (sides.MajorDim() != fData.NumSidesInSet (set_num) ||
      sides.MinorDim() != 2) throw eSizeMismatch;
  int block_ID;
  fData.ReadSideSet (set_num, block_ID, sides);
  
  iArrayT ids;
  ElementGroupIDs (ids);
  int num_elems, dim, offset = 0;
  for (int i=0; i < ids.Length(); i++)
    {
      if (ids[i] == block_ID) break;
      fData.ReadElementBlockDims (ids[i], num_elems, dim);
      offset += num_elems;
    }
  
  int *pelem = sides.Pointer();
  for (int j=0; j < sides.MajorDim(); j++, pelem += 2)
    *pelem += offset;
  
  sides += -1;
}

void ExodusInputT::ReadElementLabels_ID (int group_id, ArrayT<StringT>& elabels)
{
#pragma unused (group_id)
  fData.ReadElementLabels (elabels);
}

void ExodusInputT::ReadElementVariables_ID (int step, int group_id, dArray2DT& evalues)
{
  int num = NumElementVariables ();

  int numelems, dim;
  fData.ReadElementBlockDims (group_id, numelems, dim);
  
  // read element variables
  if (evalues.MajorDim() != numelems ||
      evalues.MinorDim() != num) throw eSizeMismatch;

  dArrayT temp (numelems);
  for (int e=0; e < num; e++)
    {
      fData.ReadElementVariable (step+1, group_id, e+1, temp);
      evalues.SetColumn (e, temp);
    }
}

/*************************************************************************
* Private
*************************************************************************/

void ExodusInputT::NodesUsed(const iArray2DT& connects, iArrayT& nodesused) const
{
	/* quick exit */
	if (connects.Length() == 0) return;

	/* compressed number range */
	int min, max;
	connects.MinMax(min, max);
	int range = max - min + 1;

	/* local map */
	iArrayT node_map(range);

	/* determine used nodes */
	node_map = 0;
	for (int i = 0; i < connects.Length(); i++)
		node_map[connects[i] - min] = 1;

	/* collect list */
	nodesused.Allocate(node_map.Count(1));
	int dex = 0;
	int*  p = node_map.Pointer();
	for (int j = 0; j < node_map.Length(); j++)
		if (*p++ == 1) nodesused[dex++] = j + min;
}
