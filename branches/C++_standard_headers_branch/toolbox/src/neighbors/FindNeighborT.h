/* $Id: FindNeighborT.h,v 1.5 2002-10-20 22:39:06 paklein Exp $ */
/* created: paklein (03/21/1997)                                          */
/* FindNeighborT.h                                                        */

#ifndef _FIND_NEIGHBORT_H_
#define _FIND_NEIGHBORT_H_

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class iGridManager1DT;
class iGridManager2DT;
class iGridManager3DT;

class FindNeighborT
{
public:

	/* Constructor */
	FindNeighborT(const dArray2DT& coords, int maxneighbors);
	FindNeighborT(const iArrayT& nodesused, const dArray2DT& coords,
		int maxneighbors);

	/* Destructor */
	virtual ~FindNeighborT(void);

	/* Print neighbors to output stream */
	virtual void OutputNeighors(ostream& out, double tolerance);
	void GetNeighors(iArray2DT& edges, double tolerance);

	/* Print coordinate data to the output stream */
	void PrintCoords(ostream& out) const;

private:

	/* Determine neighbors */
	void FindNeighors1D(double tolerance);
	void FindNeighors2D(double tolerance);
	void FindNeighors3D(double tolerance);
	
	/* allocate memory */
	void Dimension(int numpts, int nsd);

	/* Determine number of 2 body interactions */
	int Count2Body(void) const;
	int Count2BodyMapped(void) const;

	/* copy unique 2 body into edges */
	void Set2Body(iArray2DT& edges) const;
	void Set2BodyMapped(iArray2DT& edges) const;

protected:

	int fNumPts;
	int fnsd;
	int fMaxNeighbors;

	const iArrayT* const fNodeMap; //also used as mode switch

	iArray2DT	fNeighbors; //neighbor list for every node
	iArrayT		fCount;		//neighbor count

private:

	const dArray2DT& fglCoords;	//all nodal coordinates
	dArray2DT fCoords;			//searched nodal coordinates
	
	iGridManager1DT*        fGrid1D;        //search grid
	iGridManager2DT*	fGrid2D;	//search grid
	iGridManager3DT*	fGrid3D;	//search grid
};

} // namespace Tahoe 
#endif /* _FIND_NEIGHBORT_H_ */
