/* $Id: HHTalpha.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/14/1996)                                          */

#ifndef _HHT_ALPHA_H_
#define _HHT_ALPHA_H_

#include "Environment.h"

/* base class */
#include "ControllerT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class dArrayT;

/* automatic setting of gamma and beta */
const int	kHHTalphaAuto_O2 = 1;

class HHTalpha: public ControllerT
{
public:

	/* constructor */
	HHTalpha(ifstreamT& in, ostream& out, int auto2ndorder = 0);

protected:

	/* set time integration to single parameter 2nd order */
	void Set2ndOrder(double alpha);
	
protected:

	/* autoset parameters */
	int		fAuto2ndOrder;

	/* time integration parameters */
	double	fgamma;
	double	fbeta;
	double	falpha;
		
};

#endif /* _HHT_ALPHA_H_ */
