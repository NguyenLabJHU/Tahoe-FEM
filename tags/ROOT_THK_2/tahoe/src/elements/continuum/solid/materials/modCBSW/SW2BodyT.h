/* $Id: SW2BodyT.h,v 1.3 2002-07-05 22:28:22 paklein Exp $ */
/* created: paklein (05/20/1997)                                          */

#ifndef _SW2_BODY_T_H_
#define _SW2_BODY_T_H_

/* base class */
#include "TwoBodyT.h"

namespace Tahoe {

/* forward declarations */
class SWDataT;

class SW2BodyT: public TwoBodyT
{
public:

	/* constructor */
	SW2BodyT(const dArrayT& lengths, const ThermalDilatationT* thermal,
		const SWDataT& SW);

	/* set free dof - triggers recomputation */
	virtual void Set(void);

private:

	/* 2 body potential and derivatives */
	double U2body(double r, double a) const;
	double DU2body(double r, double a) const;
	double DDU2body(double r, double a) const;

private:

	/* Stillinger-Weber parameters */
	const SWDataT&	fSW;

};

} // namespace Tahoe 
#endif /* _SW2_BODY_T_H_ */
