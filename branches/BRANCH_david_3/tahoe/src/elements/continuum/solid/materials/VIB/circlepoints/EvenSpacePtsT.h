/* $Id: EvenSpacePtsT.h,v 1.4 2004-07-15 08:28:03 paklein Exp $ */
/* created: paklein (11/02/1997) */
#ifndef _EVENSPACE_PTS_T_H_
#define _EVENSPACE_PTS_T_H_

/* base class */
#include "CirclePointsT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

class EvenSpacePtsT: public CirclePointsT
{
public:

	/* constructor */
	EvenSpacePtsT(int n);

	/* generate points with the given orientation angle theta */
	virtual const dArray2DT& CirclePoints(double theta);

private:

	/* parameters */
	int fNtheta;
			
};

} // namespace Tahoe 
#endif /* _EVENSPACE_PTS_T_H_ */
