/* $Id: PointInCellDataT.h,v 1.1.2.5 2003-02-23 02:39:28 paklein Exp $ */
#ifndef _POINT_IN_CELL_DATA_T_H_
#define _POINT_IN_CELL_DATA_T_H_

/* direct members */
#include "RaggedArray2DT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "InverseMapT.h"

namespace Tahoe {

/* forward declarations */
class ContinuumElementT;

/** data structures needed for evaluating field data at arbitrary
 * points in a mesh. */
class PointInCellDataT
{
public:

	/** constructor */
	PointInCellDataT(void) { };
	
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
	/*@}*/
	
	/** \name interpolation data
	 * The major dimension of these arrays is the number of interpolating points. */
	/*@{*/
	/** interpolation data */
	dArray2DT& InterpolationWeights(void) { return fInterpolationWeights; };

	/** const access to interpolation data */
	const dArray2DT& InterpolationWeights(void) const { return fInterpolationWeights; };

	/** cell containing each point of interpolation */
	iArrayT& InterpolatingCell(void) { return fInterpolatingCell; };

	/** const access to the cells containing each point of interpolation */
	const iArrayT& InterpolatingCell(void) const { return fInterpolatingCell; };

	/** global to local map */
	InverseMapT& GlobalToLocal(void) { return fGlobalToLocal; };

	/** const access to global to local map */
	const InverseMapT& GlobalToLocal(void) const { return fGlobalToLocal; };
	/*@}*/

	/** generate non-empty cell connectivities in local numbering. Generates the data
	 * accessed with PointInCellDataT::CellNodes and PointInCellDataT::CellConnectivities.
	 * The numbering of nodes corresponds to the index of the real node number in the 
	 * PointInCellDataT::CellNodes array. */
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
	 	
	/** \name interpolation data */
	/*@{*/
	/** map from global id of interpolating point to the index in the
	 * interpolation data */
	InverseMapT fGlobalToLocal;
	
	/** cell containing each point of interpolation */
	iArrayT fInterpolatingCell;
	
	/** interpolation weights. Dimension [np] x [nen]. Assumes all cells have the same 
	 * number of nodes, and therefore the same number of weights. */
	dArray2DT fInterpolationWeights;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _POINT_IN_CELL_DATA_T_H_ */
