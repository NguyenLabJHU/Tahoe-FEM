/* $Id: OutputSetT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (03/07/2000)                                          */

#include "OutputSetT.h"
#include "iArrayT.h"
#include "iArray2DT.h"

/* array behavior */
const bool ArrayT<OutputSetT*>::fByteCopy = true;

/* constructor */
OutputSetT::OutputSetT(int ID, GeometryT::CodeT geometry_code,
	const iArray2DT& connectivities, const ArrayT<StringT>& n_labels,
	const ArrayT<StringT>& e_labels, bool changing):
	fPrintStep(-1),
	fID(ID),
	fChanging(changing),
	fGeometry(geometry_code),
	fConnectivities(connectivities)
{
	fNodeOutputLabels.Allocate(n_labels.Length());
	for (int i = 0; i < fNodeOutputLabels.Length(); i++)
		fNodeOutputLabels[i] = n_labels[i];

	fElementOutputLabels.Allocate(e_labels.Length());
	for (int j = 0; j < fElementOutputLabels.Length(); j++)
		fElementOutputLabels[j] = e_labels[j];
		
	/* set the nodes used array */
	SetNodesUsed();
}

OutputSetT::OutputSetT(const OutputSetT& source):
	fPrintStep(-1),
	fID(source.fID),
	fChanging(source.fChanging),
	fGeometry(source.fGeometry),
	fConnectivities(source.fConnectivities),
	fNodesUsed(source.fNodesUsed)
{
	fNodeOutputLabels.Allocate(source.fNodeOutputLabels.Length());
	for (int i = 0; i < fNodeOutputLabels.Length(); i++)
		fNodeOutputLabels[i] = source.fNodeOutputLabels[i];

	fElementOutputLabels.Allocate(source.fElementOutputLabels.Length());
	for (int j = 0; j < fElementOutputLabels.Length(); j++)
		fElementOutputLabels[j] = source.fElementOutputLabels[j];
}

/* dimensions */
int OutputSetT::NumNodes(void) const { return fNodesUsed.Length(); }
int OutputSetT::NumElements(void) const { return fConnectivities.MajorDim(); }

/* set nodes used */
void OutputSetT::SetNodesUsed(void)
{
	/* quick exit */
	if (fConnectivities.Length() == 0) return;

	/* compressed number range */
	int min, max;
	fConnectivities.MinMax(min, max);
	int range = max - min + 1;

	/* local map */
	iArrayT node_map(range);

	/* determine used nodes */
	node_map = 0;
	for (int i = 0; i < fConnectivities.Length(); i++)
		node_map[fConnectivities[i] - min] = 1;

	/* collect list */
	fNodesUsed.Allocate(node_map.Count(1));
	int dex = 0;
	int*  p = node_map.Pointer();
	for (int j = 0; j < node_map.Length(); j++)
		if (*p++ == 1) fNodesUsed[dex++] = j + min;
}
