/* $Id: FindNeighbor23T.h,v 1.2 2002-07-02 19:57:20 cjkimme Exp $ */
/* created: paklein (03/21/1997)                                          */
/* FindNeighbor23T.h                                                      */
/* Interface for finding 2 and 3 body neighbors                           */

#ifndef _FIND_NEIGHBOR23T_H_
#define _FIND_NEIGHBOR23T_H_

/* base class */
#include "FindNeighborT.h"


namespace Tahoe {

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

} // namespace Tahoe 
#endif /* _FIND_NEIGHBOR23T_H_ */
