/* $Id: PointInCellDataT.cpp,v 1.3 2003-05-05 00:58:26 paklein Exp $ */
#include "PointInCellDataT.h"
#include "ContinuumElementT.h"
#include "InverseMapT.h"

using namespace Tahoe;

/* collect a list of the nodes used in cells containing a non-zero number
 * of points */
void PointInCellDataT::GenerateCellConnectivities(void)
{
	const char caller[] = "PointInCellDataT::GenerateCellConnectivities";
	if (!fContinuumElement) ExceptionT::GeneralFail(caller, "element pointer not set");

	const ElementSupportT& elem_support = fContinuumElement->ElementSupport();

	/* mark nodes used in non-empty cells */
	iArrayT nodes_used(elem_support.NumNodes());
	nodes_used = 0;
	int cell_count = 0;
	for (int i = 0; i < fPointInCell.MajorDim(); i++) 
	{
		if (fPointInCell.MinorDim(i) > 0) 
		{
			cell_count++;
			const iArrayT& nodes = fContinuumElement->ElementCard(i).NodesX();
			for (int j = 0; j < nodes.Length(); j++)
				nodes_used[nodes[j]] = 1;
		}
	}
	
	/* collect nodes */
	fCellNodes.Dimension(nodes_used.Count(1));
	int dex = 0;
	for (int i = 0; i < nodes_used.Length(); i++)
		if (nodes_used[i] == 1)
			fCellNodes[dex++] = i;

	/* dimension local connectivities */
	fCellConnectivities.Dimension(cell_count, fContinuumElement->NumElementNodes());

	/* inverse numbering map */
	InverseMapT global_to_local;
	global_to_local.SetMap(fCellNodes);

	/* create connectivities in local numbering */
	dex = 0;
	for (int i = 0; i < fPointInCell.MajorDim(); i++) 
	{
		if (fPointInCell.MinorDim(i) > 0) 
		{
			int* local = fCellConnectivities(dex++);
			const iArrayT& nodes = fContinuumElement->ElementCard(i).NodesX();
			for (int j = 0; j < nodes.Length(); j++)
				local[j] = global_to_local.Map(nodes[j]);
		}
	}
}
