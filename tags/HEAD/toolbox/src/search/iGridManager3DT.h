/* $Id: iGridManager3DT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (12/09/1997)                                          */
/* iNodeT grid. offset specifies the starting number for the iNodeT       */
/* tags.                                                                  */

#ifndef _I_GRIDMANAGER3D_T_H_
#define _I_GRIDMANAGER3D_T_H_

/* base class */
#include "GridManager3DT.h"

/* direct members */
#include "iNodeT.h"

/* forward declarations */
class iArrayT;

class iGridManager3DT: public GridManager3DT<iNodeT>
{
public:

	/* constructor */
	iGridManager3DT(int nx, int ny, int nz, const dArray2DT& coords,
		const iArrayT* nodes_used);
	
	/* neighbors - returns neighbors coords(n) (SELF not included) */
	void Neighbors(int n, double tol, AutoArrayT<int>& neighbors);

	/* reconfigure grid with stored coordinate data */
	void Reset(void);

protected:

	/* process hit list - returns all, and repeated, tags != skiptag */
	void ProcessHits(double* target, double tol, int skiptag, int& count,
		iArrayT& neighbors);

	/* process hit list - only returns unique tags > skiptag */
	void ProcessHitsSorted(double* target, double tol,
		int skiptag, AutoArrayT<int>& neighbors);
	
protected:

	const dArray2DT& fCoords;
	const iArrayT* fNodesUsed;
};

#endif /* _I_GRIDMANAGER3D_T_H_ */
