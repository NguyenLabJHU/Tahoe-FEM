/* $Id: AVSOutputT.cpp,v 1.7 2004-06-17 06:41:07 paklein Exp $ */
/* created: sawimme (05/10/2001) */
#include "AVSOutputT.h"

#include "AVST.h"
#include "OutputSetT.h"
#include "dArray2DT.h"
#include "ifstreamT.h"
#include "iArray2DT.h"

using namespace Tahoe;

AVSOutputT::AVSOutputT (ostream& out, const ArrayT<StringT>& out_strings, bool binary) :
  OutputBaseT (out, out_strings),
  fBinary (binary)
{
}

void AVSOutputT::WriteGeometry (void)
{
  AVST avs (fout, fBinary);

  int num_sets = fElementSets.Length();
  CreateElementBlockIDs ();
  for (int e=0; e < num_sets; e++)
    if (fElementSets[e]->NumNodes() > 0)
      {
	StringT avsfile = CreateFileName (e);
	ofstream avsout (avsfile);

	int num_elems = fElementSets[e]->NumElements();
	int num_nodes = fElementSets[e]->NumNodes();
	int num_nvars = 0;
	int num_evars = 0;
	int num_gvars = 0;

	avs.WriteHeader (avsout, num_nodes, num_elems, num_nvars, num_evars, num_gvars);

	iArrayT nodes_used;
	nodes_used.Alias (fElementSets[e]->NodesUsed());

	WriteCoordinates (avsout, avs, e, nodes_used);
	WriteConnectivity (avsout, avs, e, nodes_used);
      }
}

void AVSOutputT::WriteOutput (double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values)
{
  OutputBaseT::WriteOutput (time, ID, n_values, e_values);

  AVST avs (fout, fBinary);

  if (fElementSets[ID]->NumNodes() > 0)
    {
      StringT avsfile = CreateFileName (ID);
      ofstream avsout (avsfile);
      const ArrayT<StringT>& node_labels = fElementSets[ID]->NodeOutputLabels();
      const ArrayT<StringT>& elem_labels = fElementSets[ID]->ElementOutputLabels();

      int num_elems = fElementSets[ID]->NumElements();
      int num_nodes = fElementSets[ID]->NumNodes();
      int num_nvars = 0;
      int num_evars = 0;
      int num_gvars = 0;
      CountVariables (num_nvars, node_labels);
      CountVariables (num_evars, elem_labels);
      avs.WriteHeader (avsout, num_nodes, num_elems, num_nvars, num_evars, num_gvars);

      iArrayT nodes_used;
      nodes_used.Alias (fElementSets[ID]->NodesUsed());
      WriteCoordinates (avsout, avs, ID, nodes_used);
      WriteConnectivity (avsout, avs, ID, nodes_used);

      WriteVariable (avsout, avs, node_labels, n_values, num_nvars);
      WriteVariable (avsout, avs, elem_labels, e_values, num_evars);
    }
}

// *************** PRIVATE ********************

StringT AVSOutputT::CreateFileName (int index) const
{
  StringT var (fOutroot);

  /* tack on sequence number */
  if (fSequence > 0) var.Append (".sq", fSequence + 1);
  
  /* tack on output set */
  var.Append (".io", fElementBlockIDs[index][0]);

  /* tack on print increment */
  var.Append (".ps", fElementSets[index]->PrintStep() + 1);

  /* tack on extension */
  var.Append (".inp");

  return var;
}

void AVSOutputT::CountVariables (int &num, const ArrayT<StringT>& labels) const
{
      ArrayT<StringT> extension (3);
      num = labels.Length();
      int dof = fCoordinates->MinorDim();

      /* account for filling vector to 3 dimensions */
      if (dof < 3)
	for (int i=0; i < labels.Length(); i++)
	  if (IsVector (labels, i, extension, dof)) 
	    num += 3 - dof;
}

void AVSOutputT::WriteCoordinates (ostream &avsout, AVST &avs, int index, iArrayT &nodes_used) const
{
#pragma unused(index)

  dArray2DT local (nodes_used.Length(), fCoordinates->MinorDim());
  for (int i=0; i < nodes_used.Length(); i++)
    local.SetRow (i, (*fCoordinates)(nodes_used[i]));

  int firstnodeID = 1;
  avs.WriteCoordinates (avsout, local, firstnodeID);
}

void AVSOutputT::WriteConnectivity (ostream &avsout, AVST &avs, int index, iArrayT &nodes_used) const
{
  int num_blocks = fElementSets[index]->NumBlocks();
  int firstelemID = 1;
  for (int i=0; i < num_blocks; i++)
    {
      const iArray2DT* connects = fElementSets[index]->Connectivities(fElementSets[index]->BlockID(i));
      iArray2DT localconn (connects->MajorDim(), connects->MinorDim());
      LocalConnectivity (nodes_used, *connects, localconn);
      localconn++;

      avs.WriteCells (avsout, fElementSets[index]->Geometry(), localconn, fElementBlockIDs[index][i], firstelemID);
      firstelemID += connects->MajorDim();
    }
}

void AVSOutputT::WriteVariable (ostream &avsout, AVST &avs, const ArrayT<StringT>& labels, const dArray2DT& values, int num_vars) const
{
  int dof = fCoordinates->MinorDim();
  int firstnodeID = 1;

  if (labels.Length () <= 0) return;

  /* if there are no vector variables, or vectors are 3D */
  if (num_vars == labels.Length())
    {
      avs.WriteDataHeader (avsout, labels);
      avs.WriteData (avsout, values, firstnodeID);
      return;
    }

  /* else if vectors, create labels */
  int count = 0;
  AutoArrayT<StringT> actual_labels (20);
  AutoArrayT<int> var_dim;
  ArrayT<StringT> extension (3);
  for (int i=0; i < labels.Length(); i++)
    if (IsVector (labels, i, extension, dof))
      {
	actual_labels.Append (extension);
	var_dim.Append (3);
	count += 3;
	i += dof - 1;
      }
    else
      {
	actual_labels.Append (labels[i]);
	var_dim.Append (1);
	count ++;
      }
  avs.WriteDataHeader (avsout, actual_labels);

  /* fill vector slots to 3D */
  int length = values.MajorDim();
  int wv = values.MinorDim();
  dArray2DT av (num_vars, length);
  av = 0.;
  dArray2DT temp (wv, length);
  temp.Transpose (values);
  for (int k=0, j=0, l=0; k < wv; k++)
    {
      av.CopyPart (j, temp, l, length);
      l += length;
      j += length;

      if (var_dim[k] == 3)
	{
	  for (int m=1; m < dof && k < wv; m++, k++)
	    {
	      av.CopyPart (j, temp, l, length);
	      l += length;
	      j += length;
	    }
	  for (int n=dof; n < 3; n++)
	    j += length;
	}
    }
  temp.Free ();
  temp.Allocate (length, num_vars);
  temp.Transpose (av);
  avs.WriteData (avsout, temp, firstnodeID);
}

// is the variable a vector, if so, provide new labels for a 3D vector.
bool AVSOutputT::IsVector (const ArrayT<StringT>& inlabels, int index, ArrayT<StringT>& extension, int dof) const
{
  if ( (strstr ((const char*) inlabels[index], "_x")) == NULL &&
       (strstr ((const char*) inlabels[index], "_X")) == NULL) 
    return false;
  
  if (dof >= 2)
    {
      if (inlabels.Length() < index+1)
	return false;
      if ( (strstr ((const char*) inlabels[index+1], "_y")) == NULL &&
	   (strstr ((const char*) inlabels[index+1], "_Y")) == NULL)
	return false;
    }

  if (dof == 3)
    {
      if (inlabels.Length() < index+2)
	return false;
      if ( (strstr ((const char*) inlabels[index+2], "_z")) == NULL &&
	   (strstr ((const char*) inlabels[index+2], "_Z")) == NULL)
	return false;
    }

  // create vector extensions, must fill to 3D
  StringT temp (inlabels[index]);
  temp.Drop (-2);
  extension[0] = temp;
  extension[0].Append ("_X");
  extension[1] = temp;
  extension[1].Append ("_Y");
  extension[2] = temp;
  extension[2].Append ("_Z");

  return true;
}
