/* $Id: TecPlotOutputT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (06/06/2000)                                          */

#include "TecPlotOutputT.h"
#include "TecPlotT.h"
#include "OutputSetT.h"
#include "ios_fwd_decl.h"
#include <fstream.h>
#include "iArray2DT.h"
#include "dArray2DT.h"

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
FileName (0, filename);
ofstream out (filename);

// write header data
TecPlotT tec (fout, false);
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
	StringT title = "Grp ";
	title.Append (fElementSets[e]->ID());
	tec.WriteFEZone (out, title, fElementSets[e]->NumNodes(), fElementSets[e]->NumElements(), fElementSets[e]->Geometry(), true);
	
	iArrayT nodes_used;
	nodes_used.Alias(fElementSets[e]->NodesUsed());
	dArray2DT local_coords(nodes_used.Length(), fCoordinates->MinorDim());
	local_coords.RowCollect(nodes_used, *fCoordinates);
	tec.WriteData (out, local_coords);
	
	const iArray2DT& connects = fElementSets[e]->Connectivities();
	iArray2DT local_connects(connects.MajorDim(), connects.MinorDim());
	LocalConnectivity(nodes_used, connects, local_connects);
	local_connects++;
	tec.WriteConnectivity (out, fElementSets[e]->Geometry(), local_connects);
	local_connects--;
}
}

void TecPlotOutputT::WriteOutput(double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values)
{
OutputBaseT::WriteOutput (time, ID, n_values, e_values);
if (fElementSets[ID]->NumNodes() == 0) return;

// open file
StringT filename;
FileName (ID, filename);
ofstream out;
TecPlotT tec (fout, false);
if (fElementSets[ID]->PrintStep() == 0)
out.open (filename);
else
out.open (filename, ios::app);

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

// writing one zone to the file per print increment
StringT title = "Grp ";
title.Append (fElementSets[ID]->ID());
if (fElementSets[ID]->PrintStep() == 0 || fElementSets[ID]->Changing())
tec.WriteFEZone (out, title, fElementSets[ID]->NumNodes(), fElementSets[ID]->NumElements(), fElementSets[ID]->Geometry(), true);
else
tec.WriteFEZone (out, title, fElementSets[ID]->NumNodes(), fElementSets[ID]->NumElements(), fElementSets[ID]->Geometry(), false);

// write coordinates
iArrayT nodes_used;
nodes_used.Alias (fElementSets[ID]->NodesUsed());
dArray2DT local_coords (nodes_used.Length(), fCoordinates->MinorDim());
local_coords.RowCollect (nodes_used, *fCoordinates);
tec.WriteData (out, local_coords);

// write variable data, since we are using BLOCK format,
// can write separately from coordinate list
tec.WriteData (out, n_values);

// write connectivity
if (fElementSets[ID]->PrintStep() == 0 || fElementSets[ID]->Changing())
{
const iArray2DT& connects = fElementSets[ID]->Connectivities();
iArray2DT local_connects(connects.MajorDim(), connects.MinorDim());
LocalConnectivity(nodes_used, connects, local_connects);
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
void TecPlotOutputT::FileName(int ID, StringT& filename) const
{
	/* root */
	filename = fOutroot;

	/* tack on sequence number */
	if (fSequence > 0) filename.Append(".sq", fSequence + 1);

	/* I/O ID */
	filename.Append(".io", ID);

	/* print step */
	//filename.Append(".ps", fElementSets[ID]->PrintStep() + 1, fNumDigits);

	/* extension */
	filename.Append(".dat");
}

