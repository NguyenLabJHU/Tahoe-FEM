/* $Id: iGridManagerT.h,v 1.3 2001-06-19 00:52:19 paklein Exp $ */
/* created: paklein (09/13/1998)                                          */
/* iNodeT grid with unified interface for 2D/3D and lightweight           */
/* file dependencies                                                      */

#ifndef _I_GRIDMANAGER_T_H_
#define _I_GRIDMANAGER_T_H_

/* environment */
#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class iArrayT;
template <class TYPE> class ArrayT;
class dArray2DT;
class iNodeT;
class iGridManager2DT;
class iGridManager3DT;
template <class TYPE> class AutoArrayT;

class iGridManagerT
{
public:

	/* constructor */
	iGridManagerT(const iArrayT& n_grid, const dArray2DT& coords,
		const iArrayT* nodes_used);
	
	/* destructor */
	~iGridManagerT(void);	 	
	
	/* reconfigure grid with stored coordinate data */
	void Reset(void);
	//void Reset(const dArray2DT& coords, const iArrayT* nodes_used);
	//TODO - add version that allows you to change the coordinate data
	//       but use the same grid allocation?

	/* neighbors - returns neighbors coords(n) (SELF not included) */
	void Neighbors(int n, double tol, AutoArrayT<int>& neighbors);
	void Neighbors(int n, const ArrayT<double>& tol_xyz, AutoArrayT<int>& neighbors);

	/* return list of data falling within the defined region */
	const AutoArrayT<iNodeT>& HitsInRegion(double* coords, double distance);
	const AutoArrayT<iNodeT>& HitsInRegion(double* coords, int cell_span);
	const AutoArrayT<iNodeT>& HitsInRegion(double* coords, const ArrayT<double>& tol_xyz);

	/* the distance covered by the given cell span */
	double CellSpan(int cell_span) const;

	/* grid statistics */
	void WriteStatistics(ostream& out) const;

private:

	/* 2D/3D search grids */
	iGridManager2DT* fGrid2D;
	iGridManager3DT* fGrid3D;
};

#endif /* _I_GRIDMANAGER_T_H_ */
