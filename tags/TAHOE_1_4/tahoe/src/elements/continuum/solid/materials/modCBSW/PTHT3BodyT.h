/* $Id: PTHT3BodyT.h,v 1.4 2003-12-28 23:37:12 paklein Exp $ */
/* created: paklein (10/11/1997)                                          */

#ifndef _PTHT3_BODY_T_H_
#define _PTHT3_BODY_T_H_

/* base class */
#include "ThreeBodyT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

class PTHT3BodyT: public ThreeBodyT
{
public:

	/** constructor */
	PTHT3BodyT(const dArrayT& lengths, const dArrayT& angles,
		const iArray2DT& bondpairs, const ThermalDilatationT* thermal,
		ifstreamT& in);

	/** destructor */
	virtual ~PTHT3BodyT(void) { };

	/* triggers recomputation */
	virtual void Set(void);

private:

	/* return the 3-body potential */	
	double ComputePhi(double r1, double r2, double c12, double a) const;

	/* compute Gradient of the 3-body potential */	
	void ComputeGradient(dArrayT& grad, double r1, double r2, double c12,
		double a);

	/* compute Hessian of the 3-body potential */
	void ComputeHessian(dMatrixT& hessian, double r1, double r2, double c12,
		double a);

private:

	/* potential parameters */
	double fB;
	double fZ;
	
};

} // namespace Tahoe 
#endif /* _PTHT3_BODY_T_H_ */
