/* $Id: iGridManagerT.h,v 1.11 2004-03-18 01:15:39 paklein Exp $ */
/* created: paklein (09/13/1998) */
#ifndef _I_GRIDMANAGER_T_H_
#define _I_GRIDMANAGER_T_H_

/* environment */
#include "Environment.h"
#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class iArrayT;
template <class TYPE> class ArrayT;
class dArray2DT;
class iNodeT;
class iGridManager1DT;
class iGridManager2DT;
class iGridManager3DT;
template <class TYPE> class AutoArrayT;

/** iNodeT grid with unified interface for 1D/2D/3D and lightweight
 * file dependencies */
class iGridManagerT
{
public:

	/** \name constructors. The grid isn't actually filled with data
	 * by the constructor. Data is (re)-filled with a call to
	 * iGridManagerT::Reset. */
	/*@{*/
	/** construct grid with specified grid layout.
	 * \param n_grid number of grid cells in each coordinate direction
	 * \param coords coordinates of search objects
	 * \param nodes_used optional list of subset of coords to include
	 *        in the search grid */
	iGridManagerT(const iArrayT& n_grid, const dArray2DT& coords,
		const ArrayT<int>* nodes_used);

	/** construct grid without explicit grid layout 
	 * \param avg_cell_nodes approximate average number of points per
	 *        grid cell assuming a even distribution in space.
	 * \param max_cells maximum allowable number of cells. -1 for
	 *        no maximum.
	 * \param coords coordinates of search objects
	 * \param nodes_used optional list of subset of coords to include
	 *        in the search grid */
	iGridManagerT(int avg_cell_nodes, int max_cells, const dArray2DT& coords,
		const ArrayT<int>* nodes_used);
	/*@}*/
	
	/* destructor */
	~iGridManagerT(void);	 	
	
	/** reconfigure grid with stored coordinate data */
	void Reset(void);
	//void Reset(const dArray2DT& coords, const iArrayT* nodes_used);
	//TODO - add version that allows you to change the coordinate data
	//       but use the same grid allocation?

	/** return the coordinate array */
	const dArray2DT& Coordinates(void) const;

	/** \name neighbors
	 * returns neighbors coords(n) (SELF not included) */
	/*@{*/
	void Neighbors(int n, double tol, AutoArrayT<int>& neighbors);
	void Neighbors(int n, const ArrayT<double>& tol_xyz, AutoArrayT<int>& neighbors);
	/*@}*/

	/** \name hits in the neighborhood
	 * return list of data falling within the defined region */
	/*@{*/
	const AutoArrayT<iNodeT>& HitsInRegion(const double* coords, double distance);
	const AutoArrayT<iNodeT>& HitsInRegion(const double* coords, int cell_span);
	const AutoArrayT<iNodeT>& HitsInRegion(const double* coords, const ArrayT<double>& tol_xyz);
	/*@}*/

	/** the distance covered by the given cell span */
	double CellSpan(int cell_span) const;

	/** write grid statistics to the output stream */
	void WriteStatistics(ostream& out) const;

	/** dump grid contents to output stream */
	void DumpGrid(ostream& out) const;	

private:

	/** \name 1D/2D/3D search grids */
	/*@{*/
	iGridManager1DT* fGrid1D;
	iGridManager2DT* fGrid2D;
	iGridManager3DT* fGrid3D;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _I_GRIDMANAGER_T_H_ */
