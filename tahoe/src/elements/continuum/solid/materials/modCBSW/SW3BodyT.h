/* $Id: SW3BodyT.h,v 1.1.1.1 2001-01-29 08:20:26 paklein Exp $ */
/* created: paklein (05/22/1997)                                          */

#ifndef _SW3_BODY_T_H_
#define _SW3_BODY_T_H_

/* base class */
#include "ThreeBodyT.h"

/* forward declaration */
class SWDataT;

class SW3BodyT: public ThreeBodyT
{
public:

	/* constructor */
	SW3BodyT(const dArrayT& lengths, const dArrayT& angles,
		const iArray2DT& bondpairs, const ThermalDilatationT* thermal,
		const SWDataT& SW);

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

	/* SW potential parameters */
	const SWDataT& fSW;

};

#endif /* _SW3_BODY_T_H_ */
