/* $Id: GaussPtsT.h,v 1.4 2006-08-18 18:45:11 tdnguye Exp $ */
/* created: paklein (11/02/1997) */
#ifndef _GAUSS_PTS_T_H_
#define _GAUSS_PTS_T_H_

/* base class */
#include "CirclePointsT.h"

namespace Tahoe {

class GaussPtsT: public CirclePointsT
{
public:

	/** constructor */
	GaussPtsT(int n);

	/** generate points with the given orientation angle theta */
	virtual const dArray2DT& CirclePoints(double theta);

	virtual const dArrayT& Jacobians(const double theta, const C1FunctionT* func); 

	/* generate points with the given orientation angle theta */
	virtual const dArrayT& CircleAngles(double theta);
	
private:

	/*
	 * Returns the correct data pointer for the specified number of
	 * integration points
	 */
	void SetCoords(int numint);
	void SetJacobians(int numint);

private:

	/* parameters */
	int	fN;
			
};

} // namespace Tahoe 
#endif /* _GAUSS_PTS_T_H_ */
