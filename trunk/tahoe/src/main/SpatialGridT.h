/* $Id: SpatialGridT.h,v 1.4 2004-11-18 00:17:25 paklein Exp $ */
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

	/** return the bounds of the given grid. The bounds are associated with the last call
	 * to SpatialGridT::SetBounds and the grid dimensions set with SpatialGridT::Dimension. */
	void GridBounds(const ArrayT<int>& grid_position, dArray2DT& bounds) const;

	/** sort coordinates into bins.
	 * \param points coordinates of the points
	 * \param bin returns with the bin number associated with every point
	 * \param bin_counts returns with the number of points in each bin
	 * \param points_used array of row indicies of the points to consider or NULL is all
	 *        points should be considered */
	void Bin(const dArray2DT& points, iArrayT& bin, iArrayT& bin_counts, 
		const ArrayT<int>* points_used = NULL);

	/** \name grid position to processor maps
	 * Inverse methods to CommManagerT::Processor2Grid. Methods return -1 if
	 * the grid indicies are off the grid. */
	/*@{*/
	int Grid2Processor(int i) const;
	int Grid2Processor(int i, int j) const;
	int Grid2Processor(int i, int j, int k) const;
	
	int Grid2Processor(const ArrayT<int>& grid_position) const;
	/*@}*/

	/** \name processor to grid position */
	/*@{*/
	void Processor2Grid(int p, int& i) const;
	void Processor2Grid(int p, int& i, int& j) const;
	void Processor2Grid(int p, int& i, int& j, int& k) const;

	void Processor2Grid(int p, ArrayT<int>& grid_position) const;
	/*@}*/

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
	
	/** total number of processes */
	int fSize;

	/** index offsets in flattened numbering of cells */
//	iArrayT fN_offset;

	/** grid bounds */
	dArray2DT fMinMax;
	
	/** grid sizes */
	dArrayT fdx;
};

inline int SpatialGridT::Grid2Processor(const ArrayT<int>& grid_position) const {
	int nsd = grid_position.Length();
	const int* gp = grid_position.Pointer();
	if (nsd == 1)
		return Grid2Processor(gp[0]);
	else if (nsd == 2)
		return Grid2Processor(gp[0], gp[1]);
	else if (nsd == 3)
		return Grid2Processor(gp[0], gp[1], gp[2]);
	else
		ExceptionT::GeneralFail("SpatialGridT::Grid2Processor", "unsupported dimension %d", nsd);
	return -1;	
}

inline void SpatialGridT::Processor2Grid(int p, ArrayT<int>& grid_position) const {
	int nsd = grid_position.Length();
	int* gp = grid_position.Pointer();
	if (nsd == 1)
		Processor2Grid(p, gp[0]);
	else if (nsd == 2)
		Processor2Grid(p, gp[0], gp[1]);
	else if (nsd == 3)
		Processor2Grid(p, gp[0], gp[1], gp[2]);
	else
		ExceptionT::GeneralFail("SpatialGridT::Processor2Grid", "unsupported dimension %d", nsd);
}

} /* namespace Tahoe */

#endif /* _SPATIAL_GRID_T_H_ */
