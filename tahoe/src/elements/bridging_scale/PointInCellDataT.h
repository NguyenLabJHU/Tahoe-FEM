/* $Id: PointInCellDataT.h,v 1.1.2.1 2003-02-10 02:16:42 paklein Exp $ */
#ifndef _POINT_IN_CELL_DATA_T_H_
#define _POINT_IN_CELL_DATA_T_H_

/* direct members */
#include "RaggedArray2DT.h"
#include "iArray2DT.h"

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
	const ContinuumElementT* ContinuumElement(void) { return fContinuumElement; };
	/*@}*/

	/** \name particle in cell data */
	/*@{*/
	RaggedArray2DT<int>& PointInCell(void) { return fPointInCell; };
	RaggedArray2DT<double>& PointInCellCoords(void) { return fPointInCellCoords; };
	/*@}*/

	/** collect a list of the nodes used in cells containing a non-zero number
	 * of points */
	void CollectCellNodes(iArrayT& cell_nodes) const;

private:

	/** associated continuum group */
	const ContinuumElementT* fContinuumElement;

//	iArray2DT fElementNodesUsed;

	/** list of points per element: [n_cell] x [n_part_i] */
	RaggedArray2DT<int> fPointInCell;
	
	/** take fPointInCell, now have list of inverse mappings per element:
	 *  [n_cell] x [n_inversemap_i] */
	RaggedArray2DT<double> fPointInCellCoords;
};

} /* namespace Tahoe */

#endif /* _POINT_IN_CELL_DATA_T_H_ */
