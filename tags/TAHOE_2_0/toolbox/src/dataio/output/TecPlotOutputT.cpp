/* $Id: TecPlotOutputT.cpp,v 1.4 2002-07-02 19:57:07 cjkimme Exp $ */
/* created: sawimme (06/06/2000)                                          */

#include "TecPlotOutputT.h"
#include "TecPlotT.h"
#include "OutputSetT.h"
#include "ios_fwd_decl.h"
#include <fstream.h>
#include "iArray2DT.h"
#include "dArray2DT.h"


using namespace Tahoe;

TecPlotOutputT::TecPlotOutputT(ostream& out, const ArrayT<StringT>& out_strings, int digits) :
  OutputBaseT(out, out_strings),
  fNumDigits (digits)
{

}

/* print geometry from multiple element groups to one file */
void TecPlotOutputT::WriteGeometry(void)
{
  // create file name
  StringT filename = fOutroot;
  FileName (0, filename, -1);
  ofstream out (filename);
  TecPlotT tec (fout, false);
  
  // write header data
  int dof = fCoordinates->MinorDim();
  ArrayT<StringT> vars (dof);
  char x= 'X';
  for (int j=0; j < dof; j++)
    vars[j].Append (x++);
  tec.WriteHeader (out, fTitle, vars);
  
  // write element sets
  for (int e=0; e < fElementSets.Length(); e++)
    if (fElementSets[e]->NumNodes() > 0)
      {
	const ArrayT<StringT>& blockIDs = fElementSets[e]->BlockID();
	for (int b=0; b < fElementSets[e]->NumBlocks(); b++)
	  {
	    // nodes used by this block
	    const iArrayT& nodes_used = fElementSets[e]->BlockNodesUsed(blockIDs[b]);

	    // write zone header
	    StringT zonetitle = "Grp ";
	    zonetitle.Append (fElementSets[e]->ID());
	    zonetitle.Append (".", blockIDs[b]);
	    int numnodes = nodes_used.Length();
	    int numelems = fElementSets[e]->NumBlockElements(blockIDs[b]);
	    tec.WriteFEZone (out, zonetitle, numnodes, numelems, fElementSets[e]->Geometry(), true);
	
	    // write only the nodes used by that connectivity block
	    dArray2DT local_coords(nodes_used.Length(), fCoordinates->MinorDim());
	    local_coords.RowCollect(nodes_used, *fCoordinates);
	    tec.WriteData (out, local_coords);
	
	    // write the connectivity block
	    const iArray2DT* c = fElementSets[e]->Connectivities(blockIDs[b]);
	    iArray2DT local_connects(c->MajorDim(), c->MinorDim());
	    LocalConnectivity(nodes_used, *c, local_connects);
	    local_connects++;
	    tec.WriteConnectivity (out, fElementSets[e]->Geometry(), local_connects);
	    local_connects--;
	  }
      }
}

void TecPlotOutputT::WriteOutput(double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values)
{
  OutputBaseT::WriteOutput (time, ID, n_values, e_values);
  if (fElementSets[ID]->NumNodes() == 0) return;
  
  // open file
  StringT filename;
  FileName (ID, filename, fElementSets[ID]->PrintStep());
  ofstream out (filename);
  TecPlotT tec (fout, false);
  
  // write header
  int dof = fCoordinates->MinorDim();
  const ArrayT<StringT>& node_labels = fElementSets[ID]->NodeOutputLabels();
  ArrayT<StringT> vars (dof + node_labels.Length());
  char x= 'X';
  for (int j=0; j < dof; j++)
    vars[j].Append (x++);
  for (int i=0; i < node_labels.Length(); i++)
    vars[i+dof] = node_labels[i];
  tec.WriteHeader (out, fTitle, vars);
  
  // one zone per block
  const ArrayT<StringT>& blockIDs = fElementSets[ID]->BlockID();
  for (int b=0; b < fElementSets[ID]->NumBlocks(); b++)
    {
      // nodes used by this block
      const iArrayT& nodes_used = fElementSets[ID]->BlockNodesUsed(blockIDs[b]);
      const iArray2DT* c = fElementSets[ID]->Connectivities(blockIDs[b]);

      // write zone header
      StringT zonetitle = "Grp ";
      zonetitle.Append (fElementSets[ID]->ID());
      zonetitle.Append (".", blockIDs[b]);
      int numnodes = nodes_used.Length();
      int numelems = fElementSets[ID]->NumBlockElements (blockIDs[b]);
      tec.WriteFEZone (out, zonetitle, numnodes, numelems, fElementSets[ID]->Geometry(), true);
  
      // write coordinates
      dArray2DT local_coords (nodes_used.Length(), fCoordinates->MinorDim());
      local_coords.RowCollect (nodes_used, *fCoordinates);
      tec.WriteData (out, local_coords);
  
      // write variable data, since we are using BLOCK format,
      // can write separately from coordinate list
      if (n_values.MajorDim() > 0)
	{
	  iArrayT block_nodes;
	  dArray2DT local_vars;
	  NodalBlockValues (ID, b, n_values, local_vars, block_nodes);
	  tec.WriteData (out, local_vars);
	}

      // write the connectivity block
      iArray2DT local_connects(c->MajorDim(), c->MinorDim());
      LocalConnectivity(nodes_used, *c, local_connects);
      local_connects++;
      tec.WriteConnectivity (out, fElementSets[ID]->Geometry(), local_connects);
      local_connects--;
    }
  
  out.flush ();
}

/*************************************************************************
* Protected
*************************************************************************/

/*************************************************************************
* Private
*************************************************************************/

/* generate database file name for the given ID */
void TecPlotOutputT::FileName(int ID, StringT& filename, int printstep) const
{
	/* root */
	filename = fOutroot;

	/* tack on sequence number */
	if (fSequence > 0) filename.Append(".sq", fSequence + 1);

	/* I/O ID */
	filename.Append(".io", ID);

	/* print step */
	if (printstep > -1) filename.Append(".ps", printstep, fNumDigits);

	/* extension */
	filename.Append(".dat");
}

