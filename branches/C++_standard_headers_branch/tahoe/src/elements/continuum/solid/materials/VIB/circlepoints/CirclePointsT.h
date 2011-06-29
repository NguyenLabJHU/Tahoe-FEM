/* $Id: CirclePointsT.h,v 1.5 2006-08-18 18:45:11 tdnguye Exp $ */
/* created: paklein (11/02/1997) */
#ifndef _CIRCLE_PTS_T_H_
#define _CIRCLE_PTS_T_H_

/* direct members */
#include "dArray2DT.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "C1FunctionT.h"

namespace Tahoe {

/** base class for circular integration point generators */
class CirclePointsT
{
public:

	/** constructor */
	CirclePointsT(void);

	/** destructor */
	virtual ~CirclePointsT(void);

	/** generate points with the given orientation angle theta */
	virtual const dArray2DT& CirclePoints(double theta) = 0;

	/** generate points with the given orientation angle theta */
	virtual const dArrayT& CircleAngles(double theta) = 0;

	/** list of jacobian determinants */
	const dArrayT& Jacobians(void) const;

	/** list of jacobian determinants */
	virtual const dArrayT& Jacobians(const double theta, const C1FunctionT* func) = 0;
	
protected:

	/** orient points with given rotations (in degrees) */
	void TransformPoints(double theta);
	
protected:

	/** point table */
	dArray2DT	fPoints;
	dArrayT   fAngles;
	
	/** jacobians */
	dArrayT		fJacobians;

private:

	/** tranformation tensor */
	dMatrixT fQ;
			
};

} // namespace Tahoe 
#endif /* _CIRCLE_PTS_T_H_ */
