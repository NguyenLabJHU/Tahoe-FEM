/* $Id: PointInCellDataT.cpp,v 1.1.2.2 2003-02-10 09:25:37 paklein Exp $ */
#include "PointInCellDataT.h"

#include "ContinuumElementT.h"

/* collect a list of the nodes used in cells containing a non-zero number
 * of points */
void PointInCellDataT::CollectCellNodes(iArrayT& cell_nodes) const
{
	const char caller[] = "PointInCellDataT::CollectCellNodes";
	if (!fContinuumElement) ExceptionT::GeneralFail(caller, "element pointer not set");

	const ElementSupportT& elem_support = fContinuumElement->ElementSupport();

	/* mark nodes used in non-empty cells */
	iArrayT nodes_used(elem_support.NumNodes());
	nodes_used = 0;
	for (int i = 0; i < fContinuumElement->NumElements(); i++) 
	{
		if (fPointInCell.MinorDim(i) > 0) 
		{
			const iArrayT& nodes = fContinuumElement->ElementCard(i).NodesX();
			for (int j = 0; j < nodes.Length(); j++) {
				nodes_used[nodes[j]] = 1;
			}
		}
	}
	
	/* collect nodes */
	cell_nodes.Dimension(nodes_used.Count(1));
	int dex = 0;
	for (int i = 0; i < nodes_used.Length(); i++)
		if (nodes_used[i] == 1)
			cell_nodes[dex++] = i;
}
