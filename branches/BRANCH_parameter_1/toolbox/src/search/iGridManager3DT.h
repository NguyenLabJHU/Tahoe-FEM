/* $Id: iGridManager3DT.h,v 1.6 2003-01-27 06:42:48 paklein Exp $ */
/* created: paklein (12/09/1997) */
#ifndef _I_GRIDMANAGER3D_T_H_
#define _I_GRIDMANAGER3D_T_H_

/* base class */
#include "GridManager3DT.h"

/* direct members */
#include "iNodeT.h"

namespace Tahoe {

/* forward declarations */
class iArrayT;

/** 3D iNodeT grid */
class iGridManager3DT: public GridManager3DT<iNodeT>
{
public:

	/** constructor */
	iGridManager3DT(int nx, int ny, int nz, const dArray2DT& coords,
		const ArrayT<int>* nodes_used);
	
	/** neighbors - returns neighbors coords(n) (SELF not included) */
	void Neighbors(int n, double tol, AutoArrayT<int>& neighbors);
	void Neighbors(int n, const ArrayT<double>& tol_xyz, AutoArrayT<int>& neighbors);

	/** reconfigure grid with stored coordinate data */
	void Reset(void);

	/** return the coordinate array */
	const dArray2DT& Coordinates(void) const { return fCoords; };

protected:

	/* process hit list - returns all, and repeated, tags != skiptag */
	void ProcessHits(double* target, double tol, int skiptag, int& count,
		iArrayT& neighbors);

	/* process hit list - only returns unique tags > skiptag */
	void ProcessHitsSorted(double* target, double tol,
		int skiptag, AutoArrayT<int>& neighbors);
	
protected:

	const dArray2DT& fCoords;
	const ArrayT<int>* fNodesUsed;
};

} // namespace Tahoe 
#endif /* _I_GRIDMANAGER3D_T_H_ */
