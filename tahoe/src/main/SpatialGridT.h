/* $Id: SpatialGridT.h,v 1.2 2004-10-06 21:28:57 paklein Exp $ */
#ifndef _SPATIAL_GRID_T_H_
#define _SPATIAL_GRID_T_H_

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class CommunicatorT;

/** class to sort points into bins for 1D/2D/3D */
class SpatialGridT
{
public:

	/** define handling of points lying outside the grid */
	enum GridBoundT {
		kCutOff, /**< ignore point beyond grid */
		kExtended, /**< extend cells along edge of grid */
		kError /**< throw exception if point lies beyond the grid */
	};

	/** constructor */
	SpatialGridT(GridBoundT grid_bound = kExtended);

	/** set number of grid cells in each dimension */
	void Dimension(const iArrayT& num_cells);

	/** \name set grid bounds */
	/*@{*/
	/** set bounds */
	void SetBounds(const dArray2DT& min_max);

	/** set bounds based on points.
	 * \param comm multiprocess communicator
	 * \param points coordinates of the points
	 * \param points_used array of row indicies of the points to consider or NULL is all
	 *        points should be considered */
	void SetBounds(CommunicatorT& comm, const dArray2DT& points, const ArrayT<int>* points_used = NULL);
	/*@}*/

	/** sort coordinates into bins.
	 * \param points coordinates of the points
	 * \param bin returns with the bin number associated with every point
	 * \param bin_counts returns with the number of points in each bin
	 * \param points_used array of row indicies of the points to consider or NULL is all
	 *        points should be considered */
	void Bin(const dArray2DT& points, iArrayT& bin, iArrayT& bin_counts, 
		const ArrayT<int>* points_used = NULL);

private:

	/** \name binning methods */
	/*@{*/
	void Bin2D(const dArray2DT& points, iArrayT& bin, iArrayT& bin_counts,  const ArrayT<int>* points_used);
	void Bin3D(const dArray2DT& points, iArrayT& bin, iArrayT& bin_counts,  const ArrayT<int>* points_used);
	/*@}*/

protected:

	/** handling of outliers */
	GridBoundT fGridBound;

	/** grid dimensions */
	iArrayT fNx;

	/** index offsets in flattened numbering of cells */
	iArrayT fN_offset;

	/** grid bounds */
	dArray2DT fMinMax;
	
	/** grid sizes */
	dArrayT fdx;
};

} /* namespace Tahoe */

#endif /* _SPATIAL_GRID_T_H_ */
