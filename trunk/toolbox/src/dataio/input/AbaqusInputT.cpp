/* $Id: AbaqusInputT.cpp,v 1.4 2001-09-04 14:46:37 sawimme Exp $ */
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
  iArrayT map;
  ReadNodeMap (map);
  for (int n=0; n < nodes.Length(); n++)
    {
      int index;
      map.HasValue (nodes[n], index);
      if (index < 0 || index >= nodes.Length()) throw eOutOfRange;
      nodes[n] = index;
    }
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
  iArrayT map;
  ReadAllElementMap (map);
  for (int n=0; n < set.Length(); n++)
    {
      int index;
      map.HasValue (set[n], index);
      if (index < 0 || index >= set.Length()) throw eOutOfRange;
      set[n] = index;
    }
}

void AbaqusInputT::ReadConnectivity (StringT& name, iArray2DT& connects)
{
  iArrayT ellist (connects.MajorDim());
  ReadGlobalElementMap (name, ellist);
  fData.ResetFile ();
  
  iArrayT map;
  ReadNodeMap (map);

  iArrayT n;
  int el;
  int cm = connects.MinorDim();
  GeometryT::CodeT code;
  for (int i=0, j=0; i < fNumElements; i++)
    {
      fData.NextElement (el, code, n);
      if (ellist.HasValue (el))
	{
	  // offset and map to start numbering at zero
	  // account for discontinuous numbering
	  for (int k=0; k < n.Length(); k++)
	    {
	      int index;
	      map.HasValue (n[k], index);
	      if (index < 0 || index >= n.Length()) throw eOutOfRange;
	      n[k] = index;
	    }
	  connects.CopyPart (j, n, 0, cm);
	  j += cm;
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

  ArrayT<AbaqusResultsT::VariableKeyT> keys (nlabels.Length());
  iArrayT dims (nlabels.Length());
  fData.NodeVariables (keys, dims);
  SetLabelName (keys, dims, nlabels);
}

void AbaqusInputT::ReadElementLabels (ArrayT<StringT>& elabels) const
{
  if (elabels.Length() != NumElementVariables()) throw eSizeMismatch;

  ArrayT<AbaqusResultsT::VariableKeyT> keys (elabels.Length());
  iArrayT dims (elabels.Length());
  fData.ElementVariables (keys, dims);
  SetLabelName (keys, dims, elabels);
}

void AbaqusInputT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const
{
  if (qlabels.Length() != NumQuadratureVariables()) throw eSizeMismatch;

  ArrayT<AbaqusResultsT::VariableKeyT> keys (qlabels.Length());
  iArrayT dims (qlabels.Length());
  fData.QuadratureVariables (keys, dims);
  SetLabelName (keys, dims, qlabels);
}

void AbaqusInputT::ReadAllNodeVariables (int step, dArray2DT& values)
{
  StringT name ("\0");
  int numv = NumNodeVariables ();
  if (values.MinorDim() != numv) throw eSizeMismatch;
  fData.ReadVariables (AbaqusResultsT::kNodeVar, step, values, name);
}

void AbaqusInputT::ReadNodeVariables (int step, StringT& elsetname, dArray2DT& values)
{
  fout << "AbaqusInputT::ReadNodeVariables not yet programmed\n";
  throw eDatabaseFail;
}

void AbaqusInputT::ReadNodeSetVariables (int step, StringT& nsetname, dArray2DT& values)
{
  int numv = NumNodeVariables ();
  if (values.MinorDim() != numv) throw eSizeMismatch;
  fData.ReadVariables (AbaqusResultsT::kNodeVar, step, values, nsetname);
}

void AbaqusInputT::ReadAllElementVariables (int step, dArray2DT& values)
{
  StringT name ("\0");
  int numv = NumNodeVariables ();
  if (values.MinorDim() != numv) throw eSizeMismatch;
  fData.ReadVariables (AbaqusResultsT::kElemVar, step, values, name);
}

void AbaqusInputT::ReadElementVariables (int step, StringT& name, dArray2DT& evalues)
{
  int numv = NumElementVariables ();
  if (evalues.MinorDim() != numv) throw eSizeMismatch;
  fData.ReadVariables (AbaqusResultsT::kElemVar, step, evalues, name);
}

void AbaqusInputT::ReadAllQuadratureVariables (int step, dArray2DT& values)
{
  //StringT name ("\0");
  //fData.ReadVariables (AbaqusResultsT::kQuadVar, step, values, name);
  fout << "AbaqusInputT::ReadAllQuadratureVariables not yet programmed\n";
  throw eDatabaseFail;
}

void AbaqusInputT::ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues)
{
  int numv = NumQuadratureVariables ();
  if (qvalues.MinorDim() != numv) throw eSizeMismatch;
  fData.ReadVariables (AbaqusResultsT::kQuadVar, step, qvalues, name);
}

/*************************************************************************
*
* Private
*
*************************************************************************/

void AbaqusInputT::SetLabelName (const ArrayT<AbaqusResultsT::VariableKeyT>& key, const iArrayT& dims, ArrayT<StringT>& name) const
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

