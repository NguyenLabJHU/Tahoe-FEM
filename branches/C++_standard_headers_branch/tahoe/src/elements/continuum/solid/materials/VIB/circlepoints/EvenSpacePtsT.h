/* $Id: EvenSpacePtsT.h,v 1.5 2006-08-18 18:45:11 tdnguye Exp $ */
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

	/** list of jacobian determinants */
	virtual const dArrayT& Jacobians(const double theta, const C1FunctionT* func);

	/* generate points with the given orientation angle theta */
	virtual const dArray2DT& CirclePoints(double theta);

	/* generate points with the given orientation angle theta */
	virtual const dArrayT& CircleAngles(double theta);

private:

	/* parameters */
	int fNtheta;
			
};

} // namespace Tahoe 
#endif /* _EVENSPACE_PTS_T_H_ */
