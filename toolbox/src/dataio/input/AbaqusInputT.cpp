/* $Id: AbaqusInputT.cpp,v 1.7 2001-12-16 23:53:44 paklein Exp $ */
/* created: sawimme (05/18/1998)                                          */

#include "AbaqusInputT.h"
#include "ios_fwd_decl.h"
#include <fstream.h>
#include "dArray2DT.h"
#include "dArrayT.h"
#include "iAutoArrayT.h"

AbaqusInputT::AbaqusInputT (ostream& out) :
  InputBaseT (out),
  fData (out)
{
}

void AbaqusInputT::Open (const StringT& file)
{
  fData.Initialize (file);
  fData.ScanFile (fNumElements, fNumNodes, fNumTimeSteps, fNumModes);
}

void AbaqusInputT::Close (void) 
{ 
  fData.Close (); 
  fNumElements = 0;
  fNumNodes = 0;
  fNumTimeSteps = 0;
  fNumModes = 0;
}

void AbaqusInputT::ElementGroupNames (ArrayT<StringT>& groupnames) const
{
  if (groupnames.Length() != fData.NumElementSets ()) throw eSizeMismatch;
  fData.ElementSetNames (groupnames);
}

void AbaqusInputT::NodeSetNames (ArrayT<StringT>& nodenames) const
{
  if (nodenames.Length() != fData.NumNodeSets ()) throw eSizeMismatch;
  fData.NodeSetNames (nodenames);
}

void AbaqusInputT::ReadNodeMap (iArrayT& nodemap)
{
  if (nodemap.Length() != fNumNodes) throw eSizeMismatch;
  fData.NodeMap (nodemap);
}

void AbaqusInputT::ReadCoordinates (dArray2DT& coords)
{
  if (coords.MajorDim() != fNumNodes) throw eSizeMismatch;
  fData.ResetFile ();
  
  dArrayT c;
  int n;
  int cm = coords.MinorDim();
  for (int i=0, j=0; i < fNumNodes; i++, j+= cm)
    {
      fData.NextCoordinate (n, c);
      coords.CopyPart (j, c, 0, cm);
    }
}

void AbaqusInputT::ReadNodeSet (StringT& name, iArrayT& nodes)
{
  if (nodes.Length() != NumNodesInSet (name)) 
    {
      fout << "\nAbaqusInputT::ReadNodeSet, array size mismatch\n";
      throw eSizeMismatch;
    }
  fData.NodeSet (name, nodes);

  // offset and map to start numbering at zero
  // account for discontinuous numbering
  iArrayT map (fNumNodes);
  ReadNodeMap (map);
  MapOffset (nodes, map);
}

void AbaqusInputT::ReadCoordinates (dArray2DT& coords, iArrayT& nodemap)
{
  ReadNodeMap (nodemap);
  ReadCoordinates (coords);
}

void AbaqusInputT::ReadAllElementMap (iArrayT& elemmap)
{
  if (elemmap.Length() != fNumElements) throw eSizeMismatch;
  fData.ElementMap (elemmap);
}

void AbaqusInputT::ReadGlobalElementMap (StringT& name, iArrayT& elemmap)
{
  if (elemmap.Length() != NumElements (name)) throw eSizeMismatch;
  fData.ElementSet (name, elemmap);
}

void AbaqusInputT::ReadGlobalElementSet (StringT& name, iArrayT& set)
{
  ReadGlobalElementMap (name, set);

  // offset and map to start numbering at zero
  // account for discontinuous numbering
  iArrayT map (fNumElements);
  ReadAllElementMap (map);
  MapOffset (set, map);
}

void AbaqusInputT::ReadConnectivity (StringT& name, iArray2DT& connects)
{
  iArrayT ellist (connects.MajorDim());
  ReadGlobalElementMap (name, ellist);
  fData.ResetFile ();
  
  iArrayT map (fNumNodes);
  ReadNodeMap (map);
  //cout << "map length " << map.Length() << endl;

  iArrayT n;
  int el;
  int cm = connects.MinorDim();
  GeometryT::CodeT code;
  GeometryT::CodeT firstcode;
  for (int i=0, j=0; i < fNumElements; i++)
    {
      fData.NextElement (el, code, n);

      /* check element type */
      if (i==0) 
	firstcode = code;
      else if (code != firstcode) 
	{
	  fout << "AbaqusInputT::ReadConnectivity, geo code does not match\n";
	  fout << "Group " << name << " firstcode = " << firstcode << " code= " << code;
	  fout << "\nElement " << el << "\n\n";
	  throw eDatabaseFail;
	}

      //cout << n << endl;
      int kdex;
      ellist.HasValue (el, kdex);
      if (kdex > -1 && kdex < connects.MajorDim())
	{
	  // offset and map to start numbering at zero
	  // account for discontinuous numbering
	  MapOffset (n, map);
	  connects.CopyPart (j, n, 0, cm);
	  j += cm;

	  // quick escape
	  if (j == connects.Length()) return;
	}
    }
}

void AbaqusInputT::ReadTimeSteps (dArrayT& steps)
{
  int num = NumTimeSteps ();
  if (steps.Length() != num) throw eSizeMismatch;
  
  int number;
  if (fNumModes > 0)
    for (int i=0; i < fNumModes; i++)
      fData.ModeData (i, number, steps[i]);
  else
    for (int i=0; i < fNumTimeSteps; i++)
      fData.TimeData (i, number, steps[i]);
}

void AbaqusInputT::ReadNodeLabels (ArrayT<StringT>& nlabels) const
{
  if (nlabels.Length() != NumNodeVariables()) throw eSizeMismatch;

  iArrayT keys (nlabels.Length());
  iArrayT dims (nlabels.Length());
  fData.NodeVariables (keys, dims);
  SetLabelName (keys, dims, nlabels);
}

void AbaqusInputT::ReadElementLabels (ArrayT<StringT>& elabels) const
{
  if (elabels.Length() != NumElementVariables()) throw eSizeMismatch;

  iArrayT keys (elabels.Length());
  iArrayT dims (elabels.Length());
  fData.ElementVariables (keys, dims);
  SetLabelName (keys, dims, elabels);
}

void AbaqusInputT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const
{
  if (qlabels.Length() != NumQuadratureVariables()) throw eSizeMismatch;

  iArrayT keys (qlabels.Length());
  iArrayT dims (qlabels.Length());
  fData.QuadratureVariables (keys, dims);
  SetLabelName (keys, dims, qlabels);
}

void AbaqusInputT::NodeVariablesUsed (StringT& name, iArrayT& used)
{
  if (used.Length() != NumNodeVariables()) throw eSizeMismatch;
  fData.VariablesUsed (name, AbaqusVariablesT::kNode, used);
}

void AbaqusInputT::ElementVariablesUsed (StringT& name, iArrayT& used)
{ 
  if (used.Length() != NumElementVariables()) throw eSizeMismatch;
  fData.VariablesUsed (name, AbaqusVariablesT::kElement, used);
}

void AbaqusInputT::QuadratureVariablesUsed (StringT& name, iArrayT& used)
{ 
  if (used.Length() != NumQuadratureVariables()) throw eSizeMismatch;
  fData.VariablesUsed (name, AbaqusVariablesT::kQuadrature, used);
}

void AbaqusInputT::ReadAllNodeVariables (int step, dArray2DT& values)
{
  StringT name ("\0");
  int numv = NumNodeVariables ();
  if (values.MinorDim() != numv) throw eSizeMismatch;
  fData.ReadVariables (AbaqusVariablesT::kNode, step, values, name);
}

void AbaqusInputT::ReadNodeVariables (int step, StringT& elsetname, dArray2DT& values)
{
  iArray2DT connects (NumElements (elsetname), NumElementNodes (elsetname));
  ReadConnectivity (elsetname, connects);

  iArrayT nodesused;
  NodesUsed (connects, nodesused);

  // read all values
  dArray2DT temp (NumNodes(), values.MinorDim());
  ReadAllNodeVariables (step, temp);

  values.Allocate (nodesused.Length(), NumNodeVariables());
  values.RowCollect (nodesused, temp);
}

void AbaqusInputT::ReadNodeSetVariables (int step, StringT& nsetname, dArray2DT& values)
{
  int numv = NumNodeVariables ();
  if (values.MinorDim() != numv) throw eSizeMismatch;
  fData.ReadVariables (AbaqusVariablesT::kNode, step, values, nsetname);
}

void AbaqusInputT::ReadAllElementVariables (int step, dArray2DT& values)
{
  StringT name ("\0");
  int numv = NumNodeVariables ();
  if (values.MinorDim() != numv) throw eSizeMismatch;
  fData.ReadVariables (AbaqusVariablesT::kElement, step, values, name);
}

void AbaqusInputT::ReadElementVariables (int step, StringT& name, dArray2DT& evalues)
{
  int numv = NumElementVariables ();
  if (evalues.MinorDim() != numv) throw eSizeMismatch;
  fData.ReadVariables (AbaqusVariablesT::kElement, step, evalues, name);
}

void AbaqusInputT::ReadAllQuadratureVariables (int step, dArray2DT& values)
{
  int numv = NumQuadratureVariables ();
  if (values.MinorDim() != numv) throw eSizeMismatch;
  StringT name ("\0");
  fData.ReadVariables (AbaqusVariablesT::kQuadrature, step, values, name);
}

void AbaqusInputT::ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues)
{
  int numv = NumQuadratureVariables ();
  if (qvalues.MinorDim() != numv) throw eSizeMismatch;
  fData.ReadVariables (AbaqusVariablesT::kQuadrature, step, qvalues, name);
}

/*************************************************************************
*
* Private
*
*************************************************************************/

void AbaqusInputT::SetLabelName (const iArrayT& key, const iArrayT& dims, ArrayT<StringT>& name) const
{
  int u=1;
  for (int i=0; i < name.Length();)
    {
      int d = dims[i];
      for (int j=0; j < d; j++, i++)
	{
	  int index = fData.VariableKeyIndex (key[i]);
	  name[i] = fData.VariableName (index);
	  name[i].Append ("_");
	  if (index > 0)
	    name[i].Append (j+1);
	else
	  name[i].Append (u++);
	}
    }
}

void AbaqusInputT::MapOffset (ArrayT<int>& set, const iArrayT& map) const
{
  int index;
  for (int n=0; n < set.Length(); n++)
    {
      map.HasValue (set[n], index);
      if (index < 0 || index >= map.Length()) throw eOutOfRange;
      set[n] = index;
    }
}

void AbaqusInputT::NodesUsed (const nArrayT<int>& connects, iArrayT& nodesused) const
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
