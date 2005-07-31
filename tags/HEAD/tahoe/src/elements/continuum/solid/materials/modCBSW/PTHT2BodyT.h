/* $Id: PTHT2BodyT.h,v 1.1.1.1 2001-01-29 08:20:26 paklein Exp $ */
/* created: paklein (10/11/1997)                                          */

#ifndef _PTHT2_BODY_T_H_
#define _PTHT2_BODY_T_H_

/* base class */
#include "TwoBodyT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;

class PTHT2BodyT: public TwoBodyT
{
public:

	/* constructor */
	PTHT2BodyT(const dArrayT& lengths, const ThermalDilatationT* thermal,
		ifstreamT& in);

	/* set free dof - triggers recomputation */
	virtual void Set(void);

private:

	/* 2 body potential and derivatives */
	double U2body(double r, double a) const;
	double DU2body(double r, double a) const;
	double DDU2body(double r, double a) const;

private:

	/* potential parameters */
	double fA;
	double fA1;
	double fA2;
};

#endif /* _PTHT2_BODY_T_H_ */
