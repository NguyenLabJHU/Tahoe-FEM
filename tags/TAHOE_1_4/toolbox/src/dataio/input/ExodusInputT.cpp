/* $Id: ExodusInputT.cpp,v 1.19 2004-05-07 22:34:12 paklein Exp $ */
/* created: sawimme (12/04/1998) */

#include "ExodusInputT.h"
#include "InputBaseT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "dArrayT.h"


using namespace Tahoe;

ExodusInputT::ExodusInputT (ostream& out) :
  InputBaseT (out),
  fData (out)
{
}

bool ExodusInputT::Open (const StringT& filename)
{
	if (!fData.OpenRead (filename)) {
		cout << "\n ExodusInputT::Open: error opening file: " << filename << endl;
		return false;
	}
	else return true;
}

void ExodusInputT::ElementGroupNames (ArrayT<StringT>& groupnames) const
{
  if (groupnames.Length() != NumElementGroups()) throw ExceptionT::kSizeMismatch;
  iArrayT ids (groupnames.Length());
  fData.ElementBlockID (ids);
  for (int i=0; i < groupnames.Length(); i++) {
  	groupnames[i].Clear();
    groupnames[i].Append(ids[i]);
    }
}

void ExodusInputT::SideSetNames (ArrayT<StringT>& sidenames) const
{
  if (sidenames.Length() != NumSideSets ()) throw ExceptionT::kSizeMismatch;
  iArrayT sidenums (sidenames.Length());
  fData.SideSetID (sidenums);
  for (int i=0; i < sidenames.Length(); i++) {
  	sidenames[i].Clear();
    sidenames[i].Append(sidenums[i]);
    }
}

void ExodusInputT::NodeSetNames (ArrayT<StringT>& nodenames) const
{
  if (nodenames.Length() != NumNodeSets ()) throw ExceptionT::kSizeMismatch;
  iArrayT nodenums (nodenames.Length());
  fData.NodeSetID (nodenums);
  for (int i=0; i < nodenames.Length(); i++) {
  	nodenames[i].Clear();
    nodenames[i].Append(nodenums[i]);
    }
}

void ExodusInputT::ReadNodeID(iArrayT& node_id)
{
	if (node_id.Length() != NumNodes()) throw ExceptionT::kSizeMismatch;
	fData.ReadNodeMap(node_id);
}

void ExodusInputT::ReadCoordinates (dArray2DT& coords)
{
  if (coords.MajorDim() != NumNodes() ||
      coords.MinorDim() != fData.NumDimensions()) throw ExceptionT::kSizeMismatch;
  fData.ReadCoordinates (coords);
}

void ExodusInputT::ReadCoordinates (dArray2DT& coords, iArrayT& node_id)
{
  ReadNodeID(node_id);
  ReadCoordinates (coords);
}

int ExodusInputT::NumGlobalElements (void) const
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

int ExodusInputT::NumElements (const StringT& name)
{
  int ID = atoi (name.Pointer());
  int numelems, numelemnodes;
  fData.ReadElementBlockDims (ID, numelems, numelemnodes);
  return numelems;
}

int ExodusInputT::NumElementNodes (const StringT& name)
{
  int ID = atoi (name.Pointer());
  int numelems, numelemnodes;
  fData.ReadElementBlockDims (ID, numelems, numelemnodes);
  return numelemnodes;
}

void ExodusInputT::ReadAllElementMap (iArrayT& elemmap)
{
  if (elemmap.Length() != NumGlobalElements()) throw ExceptionT::kSizeMismatch;
  elemmap.SetValueToPosition ();
  elemmap += 1;
}

void ExodusInputT::ReadGlobalElementMap (const StringT& name, iArrayT& elemmap)
{
  if (elemmap.Length() != NumElements (name)) throw ExceptionT::kSizeMismatch;
//  int id = atoi (name.Pointer());
  int offset = 0;
  ArrayT<StringT> eid (NumElementGroups());
  ElementGroupNames (eid);
  for (int i=0; i < eid.Length(); i++)
    {
      if (eid[i] == name) break;
      offset += NumElements (eid[i]);
    }

  elemmap.SetValueToPosition ();
  elemmap += offset + 1;
}

void ExodusInputT::ReadGlobalElementSet (const StringT& name, iArrayT& set)
{
  ReadGlobalElementMap (name, set);
  set += -1;
}

void ExodusInputT::ReadConnectivity (const StringT& name, iArray2DT& connects)
{
  int group = atoi (name.Pointer());
  int numelems, numelemnodes;
  GeometryT::CodeT geocode;
  fData.ReadElementBlockDims (group, numelems, numelemnodes);

  if (connects.MajorDim() != numelems ||
      connects.MinorDim() != numelemnodes) throw ExceptionT::kSizeMismatch;
  fData.ReadConnectivities (group, geocode, connects);

  connects += -1;
}

void ExodusInputT::ReadGeometryCode (const StringT& name, GeometryT::CodeT& code)
{
  int group = atoi (name.Pointer());
  int numelems, numelemnodes;
  fData.ReadElementBlockDims (group, numelems, numelemnodes);
  iArray2DT connects (numelems, numelemnodes);
  fData.ReadConnectivities (group, code, connects);
}

void ExodusInputT::ReadNodeSet (const StringT& name, iArrayT& nodes)
{
  int set_num = atoi (name.Pointer());
  if (nodes.Length() != fData.NumNodesInSet(set_num)) throw ExceptionT::kSizeMismatch;
  fData.ReadNodeSet (set_num, nodes);
  
  nodes += -1;
}

StringT ExodusInputT::SideSetGroupName (const StringT& name) const
{
  int setnum = atoi (name.Pointer());
  int block_ID;
  iArray2DT sides (fData.NumSidesInSet (setnum), 2);
  fData.ReadSideSet (setnum, block_ID, sides);

  StringT elname;
  elname.Append (block_ID);
  return elname;
}

void ExodusInputT::ReadSideSetLocal (const StringT& name, iArray2DT& sides) const
{
  int set_num = atoi (name.Pointer());
  if (sides.MajorDim() != fData.NumSidesInSet (set_num) ||
      sides.MinorDim() != 2) throw ExceptionT::kSizeMismatch;
  int block_ID;
  fData.ReadSideSet (set_num, block_ID, sides);
  
  sides += -1;
}

void ExodusInputT::ReadSideSetGlobal (const StringT& name, iArray2DT& sides) const
{
  int set_num = atoi (name.Pointer());
  if (sides.MajorDim() != fData.NumSidesInSet (set_num) ||
      sides.MinorDim() != 2) throw ExceptionT::kSizeMismatch;
  int block_ID;
  fData.ReadSideSet (set_num, block_ID, sides);
  
  iArrayT ids (fData.NumElementBlocks());
  fData.ElementBlockID (ids);
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

void ExodusInputT::ReadTimeSteps (dArrayT& steps)
{
  steps.Dimension(NumTimeSteps());
  for (int i=0; i < steps.Length(); i++)
    fData.ReadTime (i+1, steps[i]);
}

void ExodusInputT::NodeVariablesUsed (const StringT& name, iArrayT& used)
{
#pragma unused(name)
#pragma unused(used)

  // TEMP
  used = 1;
}

void ExodusInputT::ElementVariablesUsed (const StringT& name, iArrayT& used)
{ 
#pragma unused(name)
#pragma unused(used)

  // TEMP: I think there is a variable table that could be used ???
  used = 1;
}

void ExodusInputT::ReadAllNodeVariable (int step, int varindex, dArrayT& values)
{
  values.Allocate (NumNodes());
  fData.ReadNodalVariable (step+1, varindex+1, values);
}

void ExodusInputT::ReadNodeVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
  iArray2DT connects (NumElements (name), NumElementNodes (name));
  ReadConnectivity (name, connects);

  iArrayT nodesused;
  NodesUsed (connects, nodesused);

  values.Allocate (nodesused.Length());
  dArrayT temp (NumNodes());
  fData.ReadNodalVariable (step+1, varindex+1, temp);
  values.Collect (nodesused, temp);
}

void ExodusInputT::ReadAllNodeVariables (int step, dArray2DT& values)
{
	if (values.MajorDim() != NumNodes() ||
     	values.MinorDim() != NumNodeVariables()) {
      	cout << "\n ExodusInputT::ReadAllNodeVariables: dimension of array {"
      	     << values.MajorDim() << "," << values.MinorDim() << "}\n"
      	     <<   "     does not match the file dimensions {"
      	     << NumNodes() << "," << NumNodeVariables() << "}." << endl;
     	throw ExceptionT::kSizeMismatch;
	}

  dArrayT temp (NumNodes());
  for (int n=0; n < values.MinorDim(); n++)
    {
      fData.ReadNodalVariable (step+1, n+1, temp);
      values.SetColumn (n, temp);
    }
}


void ExodusInputT::ReadNodeVariables (int step, const StringT& name, dArray2DT& values)
{
  iArray2DT connects (NumElements (name), NumElementNodes (name));
  ReadConnectivity (name, connects);

  iArrayT nodesused;
  NodesUsed (connects, nodesused);

  values.Allocate (nodesused.Length(), NumNodeVariables ());
  //if (values.MajorDim() != nodesused ||
  //  values.MinorDim() != NumNodeVariables ()) throw ExceptionT::kSizeMismatch;

  dArrayT temp (NumNodes());
  dArray2DT temp2 (NumNodes(), NumNodeVariables ());
  for (int n=0; n < values.MinorDim(); n++)
    {
      fData.ReadNodalVariable (step+1, n+1, temp);
      temp2.SetColumn (n, temp);
    }

  values.RowCollect (nodesused, temp2);
}

void ExodusInputT::ReadNodeSetVariables (int step, const StringT& nsetname, dArray2DT& values)
{
  iArrayT ns (NumNodesInSet (nsetname));
  ReadNodeSet (nsetname, ns);

  iArrayT nodesused;
  NodesUsed (ns, nodesused);

  values.Allocate (nodesused.Length(), NumNodeVariables ());
  //if (values.MajorDim() != nodesused ||
  //  values.MinorDim() != NumNodeVariables ()) throw ExceptionT::kSizeMismatch;

  dArrayT temp (NumNodes());
  dArray2DT temp2 (NumNodes(), NumNodeVariables ());
  for (int n=0; n < values.MinorDim(); n++)
    {
      fData.ReadNodalVariable (step+1, n+1, temp);
      temp2.SetColumn (n, temp);
    }

  values.RowCollect (nodesused, temp2);
}

void ExodusInputT::ReadAllElementVariable (int step, int varindex, dArrayT& values)
{
  if (values.Length() != NumGlobalElements ()) throw ExceptionT::kSizeMismatch;

  int ng = NumElementGroups ();
  iArrayT ids (ng);
  fData.ElementBlockID (ids);

  for (int i=0, offset=0; i < ng; i++)
    {
      int numelems, dim;
      fData.ReadElementBlockDims (ids[i], numelems, dim);
      dArrayT temp (numelems);
      fData.ReadElementVariable (step+1, ids[i], varindex+1, temp);
      values.CopyPart (offset, temp, 0, numelems);
      offset += numelems;
    }
}

void ExodusInputT::ReadElementVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
  int group_id = atoi (name.Pointer());
  int numelems, dim;
  fData.ReadElementBlockDims (group_id, numelems, dim);
  if (values.Length() != numelems) throw ExceptionT::kSizeMismatch;
  fData.ReadElementVariable (step+1, group_id, varindex+1, values);
}

void ExodusInputT::ReadAllElementVariables(int step, dArray2DT& values)
{
	int nvar = values.MinorDim();
	if (values.MajorDim() != NumGlobalElements() || nvar != NumElementVariables())
		ExceptionT::SizeMismatch("ExodusInputT::ReadAllElementVariables");

	/* quick exit */
	if (nvar == 0) return;

	/* collect all element block id's */
	int ng = NumElementGroups();
	iArrayT ids(ng);
	fData.ElementBlockID(ids);
  
	dArray2DT vt(nvar, values.MajorDim());
	int offset = 0;
	for (int i = 0; i < ng; i++)
	{
		int numelems, dim;
      	fData.ReadElementBlockDims (ids[i], numelems, dim);

		dArrayT var(numelems);
		for (int e = 0; e < nvar; e++) {
			fData.ReadElementVariable(step+1, ids[i], e+1, var);
			vt.CopyPart(e*numelems + offset, var, 0, numelems);
		}

		offset += numelems;
    }
    
    /* copy/transpose to return */
  	values.Transpose(vt);
}

void ExodusInputT::ReadElementVariables (int step, const StringT& name, dArray2DT& evalues)
{
  int group_id = atoi (name.Pointer());
  int num = NumElementVariables ();

  if (num ==0)
    {
      evalues = 0;
      return;
    }

  int numelems, dim;
  fData.ReadElementBlockDims (group_id, numelems, dim);
  
  // read element variables
  if (evalues.MajorDim() != numelems ||
      evalues.MinorDim() != num) throw ExceptionT::kSizeMismatch;

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

void ExodusInputT::NodesUsed(const nArrayT<int>& connects, iArrayT& nodesused) const
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
	nodesused.Dimension(node_map.Count(1));
	int dex = 0;
	int*  p = node_map.Pointer();
	for (int j = 0; j < node_map.Length(); j++)
		if (*p++ == 1) nodesused[dex++] = j + min;
}
