/* $Id: CartesianGridT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (11/10/2000)                                          */

#ifndef _CARTESIAN_GRID_T_H_
#define _CARTESIAN_GRID_T_H_

/* direct members class */
#include "iArrayT.h"
#include "RaggedArray2DT.h"
#include "iArray2DT.h"

class CartesianGridT
{
public:

	/* grid boundary conditions */
	enum BoundaryConditionT {kFree = 0,
	                     kPeriodic = 1};

	/* constructor */
	CartesianGridT(void);

	/* (re-) set dimensions */
	void SetDimensions(const iArrayT& dimensions,
		const ArrayT<BoundaryConditionT>& bc);

	/* accessor */
	int Dimensionality(void) const;
	int Dimension(int axis) const;
	BoundaryConditionT BoundaryCondition(int axis) const;
	const iArray2DT& NeighborList(void) const;
	const iArrayT Partition(void) const;

	/* partition the grid (num_parts <= num_cells) */
	void PartitionGrid(int num_parts, const iArrayT& cell_weight);

	/* indexed distances */
	int Shift(int axis) const;
	int Width(int axis) const;

private:

	/* set neighbor lists */
	void SetNeighborLists(void);

private:

	/* dimensions along each axis */
	iArrayT fDimensions; // {n_1, n_2,..., n_i}
	
	/* boundary conditions for each axis */
	ArrayT<BoundaryConditionT> fBC;

	/* grid cell map */
	iArrayT fCellMap; // [n_1] x [n_2] x ... x [n_i] x 1, in generalized
	                  // column major ordering from n_i -> n_1, i.e.,
	                  // n_i is sequential
	
	/* grid cell neighbor lists */
	iArray2DT fNeighborList; // [num_cells] x [num_neighbors + 1]
		// [num_neighbors] = [{{n_1-, n_1+}, {n_2-, n_2+},..., {n_i-, n_i+}}]
		//                 = 2^[num_dimensions]
		// self tacked on end to make "connectivity"

	/* indexed distances */
	iArrayT fShift; // indexed distances between adjacent cells
	iArrayT fWidth; // indexed distances across the grid
};

/* inlines */
inline int CartesianGridT::Dimensionality(void) const
{
	return fDimensions.Length();
}

inline int CartesianGridT::Dimension(int axis) const
{
	return fDimensions[axis];
}

inline int CartesianGridT::Shift(int axis) const { return fShift[axis]; }
inline int CartesianGridT::Width(int axis) const { return fWidth[axis]; }


inline CartesianGridT::BoundaryConditionT
	CartesianGridT::BoundaryCondition(int axis) const
{
	return fBC[axis];
}

inline const iArray2DT& CartesianGridT::NeighborList(void) const
{
	return fNeighborList;
}

inline const iArrayT CartesianGridT::Partition(void) const
{
	return fCellMap;
}

#endif /* _CARTESIAN_GRID_T_H_ */
