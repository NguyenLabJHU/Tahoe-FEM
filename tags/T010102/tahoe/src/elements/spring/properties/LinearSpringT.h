/* $Id: LinearSpringT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (11/20/1996)                                          */

#ifndef _LINEARSPRINGT_H_
#define _LINEARSPRINGT_H_

/* base class */
#include "RodMaterialT.h"

class LinearSpringT: public RodMaterialT
{
public:

	/* constructor */
	LinearSpringT(ifstreamT& in);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* potential function and derivatives */
	virtual double Potential(double rmag, double Rmag) const;
	virtual double DPotential(double rmag, double Rmag) const;
	virtual double DDPotential(double rmag, double Rmag) const;

private:

	double	fSpringConstant;
};

#endif /* _LINEARSPRINGT_H_ */
