/* $Id: ExodusOutputT.cpp,v 1.15 2003-08-08 00:27:55 paklein Exp $ */
/* created: sawimme (05/18/1999) */

#include "ExodusOutputT.h"
#include "ExodusT.h"
#include "OutputSetT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include <time.h>


using namespace Tahoe;

ExodusOutputT::ExodusOutputT(ostream& out, const ArrayT<StringT>& out_strings):
OutputBaseT(out, out_strings)
{

}

/* print geometry from multiple element groups to one file */
void ExodusOutputT::WriteGeometry(void)
{
	ExodusT exo(cout);
	CreateElementBlockIDs ();
	CreateGeometryFile (exo);
	exo.Close ();
}

void ExodusOutputT::WriteOutput(double time, int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
	/* inherited */
	OutputBaseT::WriteOutput(time, ID, n_values, e_values);

	/* ExodusII does not like empty files with no nodes */
	if (fElementSets[ID]->NumNodes() == 0) return;

	/* ExodusII database */
	ExodusT exo(cout);
	if (fElementSets[ID]->PrintStep() == 0 ||
	    fElementSets[ID]->Changing())
		/* create new file */
	        CreateResultsFile(ID, exo);
	else
	{
		/* database file name */
		StringT filename;
		FileName(ID, filename);
	
		/* append output to existing results */
		if (!exo.OpenWrite(filename))
		  {
		    cout << "\n\nExodusOutputT::WriteOutput call to ExodusT::OpenWrite failed\n";
		    cout << "filename = " << filename << "\n";
		    cout << "time = " << time << "\n";
		    cout << "output ID = " << ID << "\n\n";
		    throw ExceptionT::kDatabaseFail;
		  }
	}
	
	/* print step - changing implies 1 result per file */
	int print_step = (fElementSets[ID]->Changing()) ? 1 : fElementSets[ID]->PrintStep() + 1;

	/* write time */
	exo.WriteTime(print_step, time);

	/* write nodal data */
	if (n_values.Length() > 0)
	{
		/* separate values by variable */
		dArrayT values(n_values.MajorDim());
		for (int i = 0; i < n_values.MinorDim(); i++)
		{
			n_values.ColumnCopy(i, values);
			exo.WriteNodalVariable(print_step, i + 1, values);
		}
	}

	/* write element data */
	const ArrayT<StringT>& blockIDs = fElementSets[ID]->BlockID ();
	if (e_values.Length() > 0)
	{
		/* separate values by block */
		for (int b=0; b < blockIDs.Length(); b++)
	    {
			dArray2DT e_block (fElementSets[ID]->NumBlockElements(blockIDs[b]), 
					   e_values.MinorDim());
			ElementBlockValues (ID, b, e_values, e_block);
	      
			/* separate values by variable */
			dArrayT values(e_block.MajorDim());
			for (int i = 0; i < e_block.MinorDim(); i++)
			{
				e_block.ColumnCopy(i, values);
				exo.WriteElementVariable(print_step, fElementBlockIDs[ID][b], i + 1, values);
			}
	    }
	}
}

/*************************************************************************
* Protected
*************************************************************************/

/*************************************************************************
* Private
*************************************************************************/

/* generate database file name for the given ID */
void ExodusOutputT::FileName(int ID, StringT& filename) const
{
	/* root */
	filename = fOutroot;

	/* tack on sequence number */
	if (fSequence > 0) filename.Append(".sq", fSequence + 1);

	/* group number */
	//filename.Append(".gp", fElementSets[ID]->ID());
	//NOTE: it's hard to resolve the number of element groups from the
	//      number of io groups based solely on what's in the database,
	//      so skip it for now.	

	/* I/O ID */
	filename.Append(".io", ID);

	/* changing geometry */
	if (fElementSets[ID]->Changing())
		filename.Append(".ps", fElementSets[ID]->PrintStep() + 1);

	/* extension */
	filename.Append(".exo");
}

/* create results file */
void ExodusOutputT::CreateResultsFile(int ID, ExodusT& exo)
{
	/* database file name */
	StringT filename;
	FileName(ID, filename);

	/* set initialization parameters */
	int dim = fCoordinates->MinorDim();
	int num_nodes = fElementSets[ID]->NumNodes();
	int num_elem = fElementSets[ID]->NumElements();
	int num_blks = fElementSets[ID]->NumBlocks();
	int num_node_sets = 0;
	int num_side_sets = 0;
	
	/* create new file */
	ArrayT<StringT> info, qa;
	AssembleQA (qa);
	if (!exo.Create(filename, fTitle, info, qa, dim, num_nodes,
			num_elem, num_blks, num_node_sets, num_side_sets))
	  {
	    cout << "\n\nExodusOutputT::WriteOutput call to ExodusT::Create failed\n";
	    cout << "output ID = " << ID << "\n\n";
	    throw ExceptionT::kDatabaseFail;
	  }

	/* write geometry */
	iArrayT nodes_used;
	nodes_used.Alias(fElementSets[ID]->NodesUsed());
	WriteCoordinates (exo, nodes_used);
	WriteConnectivity (ID, exo, nodes_used);

	/* write nodal variable labels */
	const ArrayT<StringT>& node_labels = fElementSets[ID]->NodeOutputLabels();
	exo.WriteLabels(node_labels, ExodusT::kNode);

	/* write element variable labels */
	const ArrayT<StringT>& elem_labels = fElementSets[ID]->ElementOutputLabels();
	exo.WriteLabels(elem_labels, ExodusT::kElement);
}

void ExodusOutputT::CreateGeometryFile(ExodusT& exo)
{
  /* create integer ID values from string values */
  String2IntIDs (fNodeSetNames, fNodeSetIntIDs);
  String2IntIDs (fSideSetNames, fSideSetIntIDs);

  StringT filename = fOutroot;

  /* changing geometry */
  bool change = false;
  for (int j=0; j < fElementSets.Length() && !change; j++)
    if (fElementSets[j]->Changing()) change = true;
  if (change)
    filename.Append(".ps", fElementSets[0]->PrintStep() + 1);
  filename.Append(".exo");
  
  int dim = fCoordinates->MinorDim();
  int num_nodes = fCoordinates->MajorDim();
  int num_node_sets = fNodeSets.Length();
  int num_side_sets = fSideSets.Length();
  
  int num_elem = 0, num_blks = 0;
  for (int e=0; e < fElementSets.Length(); e++)
      {
	num_blks += fElementSets[e]->NumBlocks();
	num_elem += fElementSets[e]->NumElements();
      }
  
  ArrayT<StringT> info, qa;
  AssembleQA (qa);
  exo.Create (filename, fTitle, info, qa, dim, num_nodes,
	      num_elem, num_blks, num_node_sets, num_side_sets);
  
  // write coordinates
  iArrayT nodes_used (num_nodes);
  nodes_used.SetValueToPosition();
  WriteCoordinates (exo, nodes_used);
  
  // write connectivities
  for (int i=0; i < fElementSets.Length(); i++)
      WriteConnectivity (i, exo, nodes_used);
  
  // write node sets
  for (int n=0; n < fNodeSets.Length(); n++)
    {
      iArrayT& set = *((iArrayT*) fNodeSets[n]);
      set++;
      // exodus does not support string labels, use index instead
      exo.WriteNodeSet (fNodeSetIntIDs[n], set);
      set--;
    }
  
  // write side sets, send local element numbering
  // send element block ID, not group index
  for (int s=0; s < fSideSets.Length(); s++)
    {
      /* search for group name */
      StringT& gname = fSSGroupNames [s];
      int gindex, bindex;
      ElementGroupBlockIndex (gname, gindex, bindex);

      int block_ID = fElementBlockIDs[gindex][bindex];
      iArray2DT& set = *((iArray2DT*) fSideSets[s]);
      set++;
      // exodus does not support string labels, use index instead
      exo.WriteSideSet (fSideSetIntIDs[s], block_ID, set);
      set--;
    }
}

void ExodusOutputT::AssembleQA (ArrayT<StringT>& qa) const
{
	time_t now;
	time(&now);
	char date[40], time[20];
	strftime(date, 40, "%x", localtime(&now));
	strftime(time, 20, "%X", localtime(&now));

	qa.Allocate (4);
	qa[0] = fCodeName;
	qa[1] = fVersion;
	qa[2] = date;
	qa[3] = time;
}

void ExodusOutputT::WriteCoordinates (ExodusT& exo, const iArrayT& nodes_used)
{
	dArray2DT local_coords(nodes_used.Length(), fCoordinates->MinorDim());
	local_coords.RowCollect(nodes_used, *fCoordinates);

	/* collect set ids */
	iArrayT node_id;
	if (!fNodeID) {
		node_id = nodes_used;
		node_id++;
	}
	else {
		node_id.Dimension(nodes_used.Length());
		node_id.Collect(nodes_used, *fNodeID);
	}

	/* write */
	exo.WriteCoordinates(local_coords, &node_id);
}

void ExodusOutputT::WriteConnectivity (int ID, ExodusT& exo, const iArrayT& nodes_used)
{
	iArray2DT local_connects;
	const ArrayT<StringT>& blockIDs = fElementSets[ID]->BlockID();
	for (int i = 0; i < fElementSets[ID]->NumBlocks(); i++)
    {
		const iArray2DT& connects = *(fElementSets[ID]->Connectivities(blockIDs[i]));

		/* generate connectivities in block-local numbering */
		local_connects.Dimension(connects);
		LocalConnectivity(nodes_used, connects, local_connects);

		local_connects++;
		exo.WriteConnectivities(fElementBlockIDs[ID][i], fElementSets[ID]->Geometry(), local_connects);
    }
}
