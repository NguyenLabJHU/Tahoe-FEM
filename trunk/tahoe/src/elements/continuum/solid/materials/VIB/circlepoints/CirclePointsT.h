/* $Id: CirclePointsT.h,v 1.4 2004-07-15 08:28:03 paklein Exp $ */
/* created: paklein (11/02/1997) */
#ifndef _CIRCLE_PTS_T_H_
#define _CIRCLE_PTS_T_H_

/* direct members */
#include "dArray2DT.h"
#include "dMatrixT.h"
#include "dArrayT.h"

namespace Tahoe {

/** base class for circular integration point generators */
class CirclePointsT
{
public:

	/** constructor */
	CirclePointsT(void);

	/** destructor */
	virtual ~CirclePointsT(void);

	/** enerate points with the given orientation angle theta */
	virtual const dArray2DT& CirclePoints(double theta) = 0;

	/** list of jacobian determinants */
	const dArrayT& Jacobians(void) const;
	
protected:

	/** orient points with given rotations (in degrees) */
	void TransformPoints(double theta);
	
protected:

	/** point table */
	dArray2DT	fPoints;
	
	/** jacobians */
	dArrayT		fJacobians;

private:

	/** tranformation tensor */
	dMatrixT fQ;
			
};

} // namespace Tahoe 
#endif /* _CIRCLE_PTS_T_H_ */
