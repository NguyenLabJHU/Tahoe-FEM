/* $Id: iPeriodicGrid3DT.h,v 1.1.1.1.6.1 2002-06-27 18:01:16 cjkimme Exp $ */
/* created: paklein (12/18/1997)                                          */

#ifndef _I_PER_GRID3D_T_H_
#define _I_PER_GRID3D_T_H_

/* base class */
#include "iGridManager3DT.h"


namespace Tahoe {

class iPeriodicGrid3DT: public iGridManager3DT
{
public:

	/* constructor */
	iPeriodicGrid3DT(int nx, int ny, int nz, const dArray2DT& coords,
		const iArrayT* nodes_used, const dArrayT& periodicity);
	
	/* neighbors - returns the number of neighbors */
	int PeriodicNeighbors(int n, double tol, iArrayT& neighbors);

protected:

	/* periodicity dimensions */
	const dArrayT& fPeriodicity; //pass negative for NOT periodic
	
	/* sorted hit list */
	AutoArrayT<int> fSortedHits;
};

} // namespace Tahoe 
#endif /* _I_PER_GRID3D_T_H_ */
