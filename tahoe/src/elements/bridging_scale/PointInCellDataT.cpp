/* $Id: PointInCellDataT.cpp,v 1.6 2004-06-26 05:53:20 paklein Exp $ */
#include "PointInCellDataT.h"
#include "ContinuumElementT.h"
#include "InverseMapT.h"

using namespace Tahoe;

/* collect the list of nodes in cells containing points */
int PointInCellDataT::CollectCellNodes(void)
{
	const char caller[] = "PointInCellDataT::CollectCellNodes";
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

	return cell_count;
}

/* collect a list of the nodes used in cells containing a non-zero number
 * of points */
void PointInCellDataT::GenerateCellConnectivities(void)
{
	const char caller[] = "PointInCellDataT::GenerateCellConnectivities";
	if (!fContinuumElement) ExceptionT::GeneralFail(caller, "element pointer not set");

	/* collect nodes */
	int cell_count = CollectCellNodes();

	/* dimension local connectivities */
	fCellConnectivities.Dimension(cell_count, fContinuumElement->NumElementNodes());

	/* inverse numbering map */
	InverseMapT global_to_local;
	global_to_local.SetMap(fCellNodes);

	/* create connectivities in local numbering */
	int dex = 0;
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

/* translate interpolation data to sparse matrix data */
void PointInCellDataT::InterpolationDataToMatrix(iArrayT& r, iArrayT& c, dArrayT& v) const
{
	int nv = fInterpolationWeights.Length();
	r.Dimension(nv);
	c.Dimension(nv);
	v.Dimension(nv);

	int index = 0;
	iArrayT fwd(fInterpolationWeights.MajorDim());
	fGlobalToLocal.Forward(fwd);
	for (int i = 0; i < fInterpolatingCell.Length(); i++) {
		const iArrayT& nodes = fContinuumElement->ElementCard(i).NodesU();
		const double* w = fInterpolationWeights(i);
		for (int j = 0; j < nodes.Length(); j++) {
			r[index] = fwd[i];
			c[index] = nodes[j];
			v[index] = w[j];
			index++;
		}
	}
}
