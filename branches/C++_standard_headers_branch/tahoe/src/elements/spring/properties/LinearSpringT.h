/* $Id: LinearSpringT.h,v 1.4 2005-11-06 00:37:58 paklein Exp $ */
/* created: paklein (11/20/1996) */
#ifndef _LINEARSPRINGT_H_
#define _LINEARSPRINGT_H_

/* base class */
#include "RodMaterialT.h"

namespace Tahoe {

class LinearSpringT: public RodMaterialT
{
public:

	/* constructor */
	LinearSpringT(double mass, double k);
	
	/* potential function and derivatives */
	virtual double Potential(double rmag, double Rmag) const;
	virtual double DPotential(double rmag, double Rmag) const;
	virtual double DDPotential(double rmag, double Rmag) const;

private:

	double	fSpringConstant;
};

} // namespace Tahoe 
#endif /* _LINEARSPRINGT_H_ */
