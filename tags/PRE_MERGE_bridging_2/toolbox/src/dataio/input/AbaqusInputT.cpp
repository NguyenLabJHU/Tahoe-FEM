/* $Id: AbaqusInputT.cpp,v 1.14 2002-10-20 22:36:54 paklein Exp $ */
/* created: sawimme (05/18/1998) */

#include "AbaqusInputT.h"
#include "ios_fwd_decl.h"
#include <fstream.h>
#include "dArray2DT.h"
#include "dArrayT.h"
#include "iAutoArrayT.h"


using namespace Tahoe;

AbaqusInputT::AbaqusInputT (ostream& out) :
  InputBaseT (out),
  fData (out)
{
}

bool AbaqusInputT::Open (const StringT& file)
{
	if (!fData.Initialize (file))
	{
		cout << "\n AbaqusInputT::Open: error initializing file: " << file << endl;
		return false;
	}
	if (!fData.ScanFile (fNumElements, fNumNodes, fNumTimeSteps, fNumModes))
	{
		cout << "\n AbaqusInputT::Open: error scanning file: " << file << endl;  
		return false;  	
	} return true;
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
  if (groupnames.Length() != fData.NumElementSets ()) throw ExceptionT::kSizeMismatch;
  fData.ElementSetNames (groupnames);
}

void AbaqusInputT::NodeSetNames (ArrayT<StringT>& nodenames) const
{
  if (nodenames.Length() != fData.NumNodeSets ()) throw ExceptionT::kSizeMismatch;
  fData.NodeSetNames (nodenames);
}

void AbaqusInputT::ReadNodeID (iArrayT& node_id)
{
	if (node_id.Length() != fNumNodes) throw ExceptionT::kSizeMismatch;
	fData.NodeMap(node_id);
}

void AbaqusInputT::ReadCoordinates (dArray2DT& coords)
{
  if (coords.MajorDim() != fNumNodes) throw ExceptionT::kSizeMismatch;
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

void AbaqusInputT::ReadNodeSet (const StringT& name, iArrayT& nodes)
{
  if (nodes.Length() != NumNodesInSet (name)) 
    {
      fout << "\nAbaqusInputT::ReadNodeSet, array size mismatch\n";
      throw ExceptionT::kSizeMismatch;
    }
  fData.NodeSet (name, nodes);

  // offset and map to start numbering at zero
  // account for discontinuous numbering
  iArrayT map(fNumNodes);
  ReadNodeID(map);
  MapOffset(nodes, map);
}

void AbaqusInputT::ReadCoordinates (dArray2DT& coords, iArrayT& node_id)
{
	ReadNodeID(node_id);
	ReadCoordinates(coords);
}

void AbaqusInputT::ReadAllElementMap (iArrayT& elemmap)
{
  if (elemmap.Length() != fNumElements) throw ExceptionT::kSizeMismatch;
  fData.ElementMap (elemmap);
}

void AbaqusInputT::ReadGlobalElementMap (const StringT& name, iArrayT& elemmap)
{
  if (elemmap.Length() != NumElements (name)) throw ExceptionT::kSizeMismatch;
  fData.ElementSet (name, elemmap);
}

void AbaqusInputT::ReadGlobalElementSet (const StringT& name, iArrayT& set)
{
  ReadGlobalElementMap (name, set);

  // offset and map to start numbering at zero
  // account for discontinuous numbering
  iArrayT map (fNumElements);
  ReadAllElementMap (map);
  MapOffset (set, map);
}

void AbaqusInputT::ReadConnectivity (const StringT& name, iArray2DT& connects)
{
  iArrayT ellist (connects.MajorDim());
  ReadGlobalElementMap (name, ellist);
  fData.ResetFile ();
  
  iArrayT map (fNumNodes);
  ReadNodeID (map);
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
	  throw ExceptionT::kDatabaseFail;
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
  if (steps.Length() != num) throw ExceptionT::kSizeMismatch;
  
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
  if (nlabels.Length() != NumNodeVariables()) throw ExceptionT::kSizeMismatch;

  iArrayT keys (nlabels.Length());
  iArrayT dims (nlabels.Length());
  fData.NodeVariables (keys, dims);
  SetLabelName (keys, dims, nlabels);
}

void AbaqusInputT::ReadElementLabels (ArrayT<StringT>& elabels) const
{
  if (elabels.Length() != NumElementVariables()) throw ExceptionT::kSizeMismatch;

  iArrayT keys (elabels.Length());
  iArrayT dims (elabels.Length());
  fData.ElementVariables (keys, dims);
  SetLabelName (keys, dims, elabels);
}

void AbaqusInputT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const
{
  if (qlabels.Length() != NumQuadratureVariables()) throw ExceptionT::kSizeMismatch;

  iArrayT keys (qlabels.Length());
  iArrayT dims (qlabels.Length());
  fData.QuadratureVariables (keys, dims);
  SetLabelName (keys, dims, qlabels);
}

void AbaqusInputT::NodeVariablesUsed (const StringT& name, iArrayT& used)
{
  if (used.Length() != NumNodeVariables()) throw ExceptionT::kSizeMismatch;
  fData.VariablesUsed (name, AbaqusVariablesT::kNode, used);
}

void AbaqusInputT::ElementVariablesUsed (const StringT& name, iArrayT& used)
{ 
  if (used.Length() != NumElementVariables()) throw ExceptionT::kSizeMismatch;
  fData.VariablesUsed (name, AbaqusVariablesT::kElement, used);
}

void AbaqusInputT::QuadratureVariablesUsed (const StringT& name, iArrayT& used)
{ 
  if (used.Length() != NumQuadratureVariables()) throw ExceptionT::kSizeMismatch;
  fData.VariablesUsed (name, AbaqusVariablesT::kQuadrature, used);
}

void AbaqusInputT::ReadAllNodeVariable (int step, int varindex, dArrayT& value)
{
  int n = NumNodeVariables ();
  dArray2DT  vals (value.Length(), n);
  ReadAllNodeVariables (step, vals);

  double *v = vals.Pointer(varindex);
  double *t = value.Pointer();
  for (int i=0; i < value.Length(); i++)
    {
      *t++ = *v;
      v += n;
    }
}

void AbaqusInputT::ReadNodeVariable (int step, const StringT& elsetname, int varindex, dArrayT& value)
{
  int n = NumNodeVariables ();
  dArray2DT  vals (value.Length(), n);
  ReadNodeVariables (step, elsetname, vals);

  double *v = vals.Pointer(varindex);
  double *t = value.Pointer();
  for (int i=0; i < value.Length(); i++)
    {
      *t++ = *v;
      v += n;
    }
}

void AbaqusInputT::ReadAllNodeVariables (int step, dArray2DT& values)
{
  StringT name ("\0");
  int numv = NumNodeVariables ();
  if (values.MinorDim() != numv) throw ExceptionT::kSizeMismatch;
  fData.ReadVariables (AbaqusVariablesT::kNode, step, values, name);
}

void AbaqusInputT::ReadNodeVariables (int step, const StringT& elsetname, dArray2DT& values)
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

void AbaqusInputT::ReadNodeSetVariables (int step, const StringT& nsetname, dArray2DT& values)
{
  int numv = NumNodeVariables ();
  if (values.MinorDim() != numv) throw ExceptionT::kSizeMismatch;
  fData.ReadVariables (AbaqusVariablesT::kNode, step, values, nsetname);
}

void AbaqusInputT::ReadAllElementVariable (int step, int varindex, dArrayT& value)
{
  int n = NumElementVariables ();
  dArray2DT  vals (value.Length(), n);
  ReadAllElementVariables (step, vals);

  double *v = vals.Pointer(varindex);
  double *t = value.Pointer();
  for (int i=0; i < value.Length(); i++)
    {
      *t++ = *v;
      v += n;
    }
}

void AbaqusInputT::ReadElementVariable (int step, const StringT& elsetname, int varindex, dArrayT& value)
{
  int n = NumElementVariables ();
  dArray2DT  vals (value.Length(), n);
  ReadElementVariables (step, elsetname, vals);

  double *v = vals.Pointer(varindex);
  double *t = value.Pointer();
  for (int i=0; i < value.Length(); i++)
    {
      *t++ = *v;
      v += n;
    }
}

void AbaqusInputT::ReadAllElementVariables (int step, dArray2DT& values)
{
  StringT name ("\0");
  int numv = NumElementVariables ();
  if (values.MinorDim() != numv) throw ExceptionT::kSizeMismatch;
  fData.ReadVariables (AbaqusVariablesT::kElement, step, values, name);
}

void AbaqusInputT::ReadElementVariables (int step, const StringT& name, dArray2DT& evalues)
{
  int numv = NumElementVariables ();
  if (evalues.MinorDim() != numv) throw ExceptionT::kSizeMismatch;
  fData.ReadVariables (AbaqusVariablesT::kElement, step, evalues, name);
}

void AbaqusInputT::ReadAllQuadratureVariable (int step, int varindex, dArrayT& value)
{
  int n = NumQuadratureVariables ();
  dArray2DT  vals (value.Length(), n);
  ReadAllQuadratureVariables (step, vals);

  double *v = vals.Pointer(varindex);
  double *t = value.Pointer();
  for (int i=0; i < value.Length(); i++)
    {
      *t++ = *v;
      v += n;
    }
}

void AbaqusInputT::ReadQuadratureVariable (int step, const StringT& elsetname, int varindex, dArrayT& value)
{
  int n = NumQuadratureVariables ();
  dArray2DT vals (value.Length(), n);
  ReadQuadratureVariables (step, elsetname, vals);

  double *v = vals.Pointer(varindex);
  double *t = value.Pointer();
  for (int i=0; i < value.Length(); i++)
    {
      *t++ = *v;
      v += n;
    }
}

void AbaqusInputT::ReadAllQuadratureVariables (int step, dArray2DT& values)
{
  int numv = NumQuadratureVariables ();
  if (values.MinorDim() != numv) throw ExceptionT::kSizeMismatch;
  StringT name ("\0");
  fData.ReadVariables (AbaqusVariablesT::kQuadrature, step, values, name);
}

void AbaqusInputT::ReadQuadratureVariables (int step, const StringT& name, dArray2DT& qvalues)
{
  int numv = NumQuadratureVariables ();
  if (qvalues.MinorDim() != numv) throw ExceptionT::kSizeMismatch;
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
      if (index < 0 || index >= map.Length()) throw ExceptionT::kOutOfRange;
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
	nodesused.Dimension(node_map.Count(1));
	int dex = 0;
	int*  p = node_map.Pointer();
	for (int j = 0; j < node_map.Length(); j++)
		if (*p++ == 1) nodesused[dex++] = j + min;
}
