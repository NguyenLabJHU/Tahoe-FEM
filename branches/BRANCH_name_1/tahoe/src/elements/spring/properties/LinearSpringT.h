/* $Id: LinearSpringT.h,v 1.1.1.1.10.1 2002-06-27 18:03:54 cjkimme Exp $ */
/* created: paklein (11/20/1996)                                          */

#ifndef _LINEARSPRINGT_H_
#define _LINEARSPRINGT_H_

/* base class */
#include "RodMaterialT.h"


namespace Tahoe {

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

} // namespace Tahoe 
#endif /* _LINEARSPRINGT_H_ */
