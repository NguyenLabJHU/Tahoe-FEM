/* $Id: iGridManager1DT.h,v 1.3 2002-11-22 01:53:12 paklein Exp $ */
#ifndef _I_GRIDMANAGER1D_T_H_
#define _I_GRIDMANAGER1D_T_H_

/* base class */
#include "GridManager1DT.h"

/* direct members */
#include "iNodeT.h"

namespace Tahoe {

class iGridManager1DT: public GridManager1DT<iNodeT>
{
public:

	/* constructor */
	iGridManager1DT(int nx, const dArray2DT& coords, const iArrayT* nodes_used);
	
	/* neighbors - returns neighbors coords(n) (SELF not included) */
	void Neighbors(int n, double tol, AutoArrayT<int>& neighbors);
	void Neighbors(int n, const ArrayT<double>& tol_xy, AutoArrayT<int>& neighbors);
	
	/* reconfigure grid with stored coordinate data */
	void Reset(void);
	
protected:

	const dArray2DT& fCoords;
	const iArrayT* fNodesUsed;
};

} // namespace Tahoe

#endif /* _I_GRIDMANAGER1D_T_H_ */
