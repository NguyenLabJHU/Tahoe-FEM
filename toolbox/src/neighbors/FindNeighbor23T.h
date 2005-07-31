/* $Id: FindNeighbor23T.h,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (03/21/1997)                                          */
/* FindNeighbor23T.h                                                      */
/* Interface for finding 2 and 3 body neighbors                           */

#ifndef _FIND_NEIGHBOR23T_H_
#define _FIND_NEIGHBOR23T_H_

/* base class */
#include "FindNeighborT.h"

class FindNeighbor23T: public FindNeighborT
{
public:

	/* Constructor */
	FindNeighbor23T(const dArray2DT& coords, int maxneighbors);
	FindNeighbor23T(const iArrayT& nodesused, const dArray2DT& coords,
		int maxneighbors);

	/* Print neighbors to output stream */
	virtual void OutputNeighors(ostream& out, double tolerance);
	void GetNeighors(iArray2DT& edges, iArray2DT& angles, double tolerance);

private:

	/* Determine number of 3 body interactions */
	int Count3Body(void) const;

	/* Determine 3 body interactions */
	void Set3Body(iArray2DT& angles) const;
	void Set3BodyMapped(iArray2DT& angles) const;
};

#endif /* _FIND_NEIGHBOR23T_H_ */
