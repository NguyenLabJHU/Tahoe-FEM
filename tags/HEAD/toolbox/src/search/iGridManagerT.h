/* $Id: iGridManagerT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
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

	/* neighbors - returns neighbors coords(n) (SELF not included) */
	void Neighbors(int n, double tol, AutoArrayT<int>& neighbors);

	/* return list of data falling within the defined region */
	const AutoArrayT<iNodeT>& HitsInRegion(double* coords, double distance);
	const AutoArrayT<iNodeT>& HitsInRegion(double* coords, int cell_span);

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
