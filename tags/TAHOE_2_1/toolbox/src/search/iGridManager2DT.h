/* $Id: iGridManager2DT.h,v 1.5 2003-01-27 06:42:48 paklein Exp $ */
/* created: paklein (12/09/1997) */
#ifndef _I_GRIDMANAGER2D_T_H_
#define _I_GRIDMANAGER2D_T_H_

/* base class */
#include "GridManager2DT.h"

/* direct members */
#include "iNodeT.h"

namespace Tahoe {

/** 2D iNodeT grid */
class iGridManager2DT: public GridManager2DT<iNodeT>
{
public:

	/** constructor */
	iGridManager2DT(int nx, int ny, const dArray2DT& coords, const ArrayT<int>* nodes_used);
	
	/** neighbors - returns neighbors coords(n) (SELF not included) */
	void Neighbors(int n, double tol, AutoArrayT<int>& neighbors);
	void Neighbors(int n, const ArrayT<double>& tol_xy, AutoArrayT<int>& neighbors);
	
	/** reconfigure grid with stored coordinate data */
	void Reset(void);

	/** return the coordinate array */
	const dArray2DT& Coordinates(void) const { return fCoords; };

protected:

	const dArray2DT& fCoords;
	const ArrayT<int>* fNodesUsed;
};

} // namespace Tahoe 
#endif /* _I_GRIDMANAGER2D_T_H_ */
