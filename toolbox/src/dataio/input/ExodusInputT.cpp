/* $Id: ExodusInputT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (12/04/1998)                                          */

#include "ExodusInputT.h"
#include "InputBaseT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "dArrayT.h"

ExodusInputT::ExodusInputT (ostream& out, const char* filename) :
InputBaseT (out),
	fData (out)
{
fData.OpenRead (filename);
}

void ExodusInputT::GroupNumbers (iArrayT& groupnums) const
{
groupnums.Allocate (NumElementGroups ());
fData.ElementBlockID (groupnums);
}

void ExodusInputT::SideSetNumbers (iArrayT& sidenums) const
{
sidenums.Allocate (NumSideSets ());
fData.SideSetID (sidenums);
}

void ExodusInputT::NodeSetNumbers (iArrayT& nodenums) const
{
nodenums.Allocate (NumNodeSets ());
fData.NodeSetID (nodenums);
}

void ExodusInputT::ReadCoordinates (dArray2DT& coords, iArrayT& nodemap)
{
int numnodes = fData.NumNodes();
coords.Allocate (numnodes, fData.NumDimensions());
fData.ReadCoordinates (coords);
nodemap.Allocate (numnodes);
fData.ReadNodeMap (nodemap);
}

void ExodusInputT::ReadConnectivity (int group, GeometryT::CodeT& geocode, iArray2DT& connects, iArrayT& elementmap)
{
#pragma unused (elementmap)
int numelems, numelemnodes;
fData.ReadElementBlockDims (group, numelems, numelemnodes);
connects.Allocate (numelems, numelemnodes);
fData.ReadConnectivities (group, geocode, connects);

connects += -1;
}

void ExodusInputT::ReadNodeSet (int set_num, iArrayT& nodes) const
{
nodes.Allocate (fData.NumNodesInSet (set_num));
fData.ReadNodeSet (set_num, nodes);

nodes += -1;
}

void ExodusInputT::ReadSideSet (int set_num, iArray2DT& sides) const
{
sides.Allocate (fData.NumSidesInSet (set_num), 2);
int block_ID;
fData.ReadSideSet (set_num, block_ID, sides);

sides += -1;
}

void ExodusInputT::ReadSideSetGlobal (int set_num, iArray2DT& sides) const
{
sides.Allocate (fData.NumSidesInSet (set_num), 2);
int block_ID;
fData.ReadSideSet (set_num, block_ID, sides);

iArrayT ids;
GroupNumbers (ids);
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

void ExodusInputT::Close (void)
{
fData.Close ();
}

void ExodusInputT::QARecords (ArrayT<StringT>& records) const
{
fData.ReadQA (records);
}

void ExodusInputT::ReadTimeSteps (dArrayT& steps)
{
int numsteps = fData.NumTimeSteps ();
steps.Allocate (numsteps);
for (int i=0; i < steps.Length(); i++)
fData.ReadTime (i+1, steps[i]);
}

void ExodusInputT::ReadLabels (ArrayT<StringT>& nlabels, ArrayT<StringT>& elabels, int group_id)
{
#pragma unused (group_id)
fData.ReadNodeLabels (nlabels);
fData.ReadElementLabels (elabels);
}

void ExodusInputT::ReadVariables (int step, int group_id, dArray2DT& nvalues, dArray2DT& evalues)
{
ArrayT<StringT> nlabel, elabel;
ReadLabels (nlabel, elabel, group_id);

// read nodal variables
dArray2DT ntemp (fData.NumNodes(), nlabel.Length());
dArrayT temp (nvalues.MajorDim());
for (int n=0; n < nlabel.Length(); n++)
{
fData.ReadNodalVariable (step+1, n+1, temp);
ntemp.SetColumn (n, temp);
}

// only want to return nodal variable values for nodes used
int numelems, numelemnodes;
iArray2DT connects;
GeometryT::CodeT geocode;
fData.ReadElementBlockDims (group_id, numelems, numelemnodes);
connects.Allocate (numelems, numelemnodes);
fData.ReadConnectivities (group_id, geocode, connects);
connects += -1;
iArrayT nodesused;
NodesUsed (connects, nodesused);
nvalues.Allocate (nodesused.Length(), nlabel.Length());
nvalues.RowCollect (nodesused, ntemp);

// read element variables
evalues.Allocate (numelems, elabel.Length());
temp.Allocate (numelems);
for (int e=0; e < elabel.Length(); e++)
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
