/* $Id: ExodusOutputT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (05/18/1999)                                          */

#include "ExodusOutputT.h"
#include "ExodusT.h"
#include "OutputSetT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include <time.h>

ExodusOutputT::ExodusOutputT(ostream& out, const ArrayT<StringT>& out_strings):
OutputBaseT(out, out_strings)
{

}

/* print geometry from multiple element groups to one file */
void ExodusOutputT::WriteGeometry(void)
{
ExodusT exo (fout);
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
	ExodusT exo(fout);
	if (fElementSets[ID]->PrintStep() == 0)
		/* create new file */
		CreateResultsFile(ID, exo);
	else
	{
		/* database file name */
		StringT filename;
		FileName(ID, filename);
	
		/* append output to existing results */
		exo.OpenWrite(filename);
	}

	/* write time */
	exo.WriteTime(fElementSets[ID]->PrintStep() + 1, time);

	/* write nodal data */
	if (n_values.Length() > 0)
	{
		dArrayT values(n_values.MajorDim());
		for (int i = 0; i < n_values.MinorDim(); i++)
		{
			n_values.ColumnCopy(i, values);
			exo.WriteNodalVariable(fElementSets[ID]->PrintStep() + 1,
				i + 1, values);
		}
	}

	/* write element data */
	if (e_values.Length() > 0)
	{
		dArrayT values(e_values.MajorDim());
		for (int i = 0; i < e_values.MinorDim(); i++)
		{
			e_values.ColumnCopy(i, values);
			exo.WriteElementVariable(fElementSets[ID]->PrintStep() + 1,
				fElementSets[ID]->ID(), i + 1, values);
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
	int num_blks = 1;
	int num_node_sets = 0;
	int num_side_sets = 0;
	
	/* create new file */
	ArrayT<StringT> info, qa;
	AssembleQA (qa);
	exo.Create(filename, fTitle, info, qa, dim, num_nodes,
		num_elem, num_blks, num_node_sets, num_side_sets);

	/* write geometry */
	iArrayT nodes_used;
	nodes_used.Alias(fElementSets[ID]->NodesUsed());
	WriteCoordinates (exo, nodes_used);
	WriteConnectivity (ID, exo, nodes_used);

	/* write nodal variable labels */
	const ArrayT<StringT>& node_labels = fElementSets[ID]->NodeOutputLabels();
	exo.WriteNodeLabels(node_labels);

	/* write element variable labels */
	const ArrayT<StringT>& elem_labels = fElementSets[ID]->ElementOutputLabels();
	exo.WriteElementLabels(elem_labels);
}

void ExodusOutputT::CreateGeometryFile(ExodusT& exo)
{
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
if (fElementSets[e]->NumNodes() > 0)
{
	num_blks++;
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
if (fElementSets[i]->NumNodes() > 0)
WriteConnectivity (i, exo, nodes_used);

// write node sets
for (int n=0; n < fNodeSets.Length(); n++)
{
iArrayT& set = *((iArrayT*) fNodeSets[n]);
set++;
exo.WriteNodeSet (fNodeSetIDs[n], set);
set--;
}

// write side sets, send local element numbering
// send element block ID, not group index
for (int s=0; s < fSideSets.Length(); s++)
{
iArray2DT& set = *((iArray2DT*) fSideSets[s]);
set++;
int block_ID = fElementSets[fSSGroupID[s]]->ID();
exo.WriteSideSet (fSideSetIDs[s], block_ID, set);
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

void ExodusOutputT::WriteCoordinates (ExodusT& exo, iArrayT& nodes_used)
{
	dArray2DT local_coords(nodes_used.Length(), fCoordinates->MinorDim());
	local_coords.RowCollect(nodes_used, *fCoordinates);
	nodes_used++;
	exo.WriteCoordinates(local_coords, &nodes_used);
	nodes_used--;
}

void ExodusOutputT::WriteConnectivity (int ID, ExodusT& exo, const iArrayT& nodes_used)
{
	const iArray2DT& connects = fElementSets[ID]->Connectivities();
	iArray2DT local_connects(connects.MajorDim(), connects.MinorDim());
	LocalConnectivity(nodes_used, connects, local_connects);
	local_connects++;
	exo.WriteConnectivities(fElementSets[ID]->ID(), fElementSets[ID]->Geometry(),
		local_connects);
	local_connects--;
}
