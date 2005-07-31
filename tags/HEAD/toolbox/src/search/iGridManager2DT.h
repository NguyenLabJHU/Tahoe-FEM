/* $Id: iGridManager2DT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (12/09/1997)                                          */
/* iNodeT grid. offset specifies the starting number for the iNodeT       */
/* tags.                                                                  */

#ifndef _I_GRIDMANAGER2D_T_H_
#define _I_GRIDMANAGER2D_T_H_

/* base class */
#include "GridManager2DT.h"

/* direct members */
#include "iNodeT.h"

class iGridManager2DT: public GridManager2DT<iNodeT>
{
public:

	/* constructor */
	iGridManager2DT(int nx, int ny, const dArray2DT& coords, const iArrayT* nodes_used);
	
	/* neighbors - returns neighbors coords(n) (SELF not included) */
	void Neighbors(int n, double tol, AutoArrayT<int>& neighbors);
	
	/* reconfigure grid with stored coordinate data */
	void Reset(void);
	
protected:

	const dArray2DT& fCoords;
	const iArrayT* fNodesUsed;
};

#endif /* _I_GRIDMANAGER2D_T_H_ */
