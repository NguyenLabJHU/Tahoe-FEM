/* $Id: PointInCellDataT.h,v 1.7 2005-03-11 20:35:21 paklein Exp $ */
#ifndef _POINT_IN_CELL_DATA_T_H_
#define _POINT_IN_CELL_DATA_T_H_

/* direct members */
#include "RaggedArray2DT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "InverseMapT.h"
#include "InterpolationDataT.h"

namespace Tahoe {

/* forward declarations */
class ContinuumElementT;

/** data structures needed for evaluating field data at arbitrary
 * points in a mesh. */
class PointInCellDataT
{
public:

	/** constructor */
	PointInCellDataT(void);
	
	/** \name associated element group */
	/*@{*/
	void SetContinuumElement(const ContinuumElementT& element) { fContinuumElement = &element; };
	const ContinuumElementT* ContinuumElement(void) const { return fContinuumElement; };
	/*@}*/

	/** \name particle in cell data */
	/*@{*/
	RaggedArray2DT<int>& PointInCell(void) { return fPointInCell; };
	const RaggedArray2DT<int>& PointInCell(void) const { return fPointInCell; };
	RaggedArray2DT<double>& PointInCellCoords(void) { return fPointInCellCoords; };
	const RaggedArray2DT<double>& PointInCellCoords(void) const { return fPointInCellCoords; };
	/*@}*/
	
	/** \name interpolation data
	 * The major dimension of these arrays is the number of interpolating points. */
	/*@{*/
	/** interpolation data */
	dArray2DT& InterpolationWeights(void) { return fInterpolationWeights; };

	/** const access to interpolation data */
	const dArray2DT& InterpolationWeights(void) const { return fInterpolationWeights; };

	/** cell containing each point of interpolation with points referred to in
	 * local numbering */
	iArrayT& InterpolatingCell(void) { return fInterpolatingCell; };

	/** const access to the cells containing each point of interpolation */
	const iArrayT& InterpolatingCell(void) const { return fInterpolatingCell; };

	/** global to local map */
	InverseMapT& GlobalToLocal(void) { return fGlobalToLocal; };

	/** const access to global to local map */
	const InverseMapT& GlobalToLocal(void) const { return fGlobalToLocal; };
	
	/** translate interpolation data to sparse matrix data */
	void InterpolationDataToMatrix(iArrayT& r, iArrayT& c, dArrayT& v) const;
	/*@}*/

	/** \name interpolation to nodes from points in filled cells */
	/*@{*/
	InterpolationDataT& PointToNode(void) { return fPointToNode; };
	const InterpolationDataT& PointToNode(void) const { return fPointToNode; };
	/*@}*/

	/** \name projection from points to points */
	/*@{*/
	InterpolationDataT& PointToPoint(void) { return fPointToPoint; };
	const InterpolationDataT& PointToPoint(void) const { return fPointToPoint; };
	/*@}*/

	/** collect the list of nodes in cells containing points. Returns the number of non-empty
	 * cells. The list of nodes is accessed with PointInCellDataT::CellNodes */
	int CollectCellNodes(void);

	/** generate non-empty cell connectivities in local numbering. Generates the data
	 * accessed with PointInCellDataT::CellNodes and PointInCellDataT::CellConnectivities.
	 * The numbering of nodes corresponds to the index of the real node number in the 
	 * PointInCellDataT::CellNodes array. Also calls PointInCellDataT::CollectCellNodes,
	 * to collect the nodes in non-empty cells. */
	void GenerateCellConnectivities(void);

	/** nodes in cells containing points. The list corresponds to the latest
	 * call to PointInCellDataT::GenerateCellConnectivities */
	const iArrayT& CellNodes(void) const { return fCellNodes; };

	/** cell connectivities in local ordering. The connectivities correspond to the latest
	 * call to PointInCellDataT::GenerateCellConnectivities */
	const iArray2DT& CellConnectivities(void) const { return fCellConnectivities; };

private:

	/** associated continuum group */
	const ContinuumElementT* fContinuumElement;

	/** list of points per element: [n_cell] x [n_part_i] */
	RaggedArray2DT<int> fPointInCell;
	
	/** take fPointInCell, now have list of inverse mappings per element:
	 *  [n_cell] x [n_inversemap_i] */
	RaggedArray2DT<double> fPointInCellCoords;

	/** nodes in cells containing points */
	iArrayT fCellNodes;

	/** connectivities of cells containing points with local numbering.The numbering of nodes
	 * corresponds to the index of the real node number in the PointInCellDataT::fCellNodes
	 * array */
	iArray2DT fCellConnectivities;
	 	
	/** \name interpolation data. Information needed to interpolate data
	 * from the elements in PointInCellDataT::fContinuumElement onto an arbitrary
	 * set of points contained in PointInCellDataT::fGlobalToLocal. */
	/*@{*/
	/** map from global id of interpolation point to the index in the
	 * interpolation data */
	InverseMapT fGlobalToLocal;
	
	/** cell containing each point of interpolation with points referred to in
	 * local numbering */
	iArrayT fInterpolatingCell;
	
	/** interpolation weights. Dimension [np] x [nen]. Assumes all cells have the same 
	 * number of nodes, and therefore the same number of weights. */
	dArray2DT fInterpolationWeights;	
	/*@}*/

	/** interpolation to nodes from points in filled cells */
	InterpolationDataT fPointToNode;

	/** projection from points to points */
	InterpolationDataT fPointToPoint;
};

} /* namespace Tahoe */

#endif /* _POINT_IN_CELL_DATA_T_H_ */
