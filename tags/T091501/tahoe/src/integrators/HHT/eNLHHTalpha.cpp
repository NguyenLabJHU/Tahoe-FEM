/* $Id: eNLHHTalpha.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/17/1996)                                          */

#include "eNLHHTalpha.h"
#include "ExceptionCodes.h"

/* constructor */
eNLHHTalpha::eNLHHTalpha(ifstreamT& in, ostream& out, int auto2ndorder):
	HHTalpha(in, out, auto2ndorder),
	eLinearHHTalpha(in, out, auto2ndorder)
{

}

/* components of the internal force vector */
int eNLHHTalpha::FormMa(double& constMa) const
{
	constMa = fconstMa;
	return 1;
}

int eNLHHTalpha::FormCv(double& constCv) const
{
	constCv = fconstCv;
	return 1;
}

int eNLHHTalpha::FormKd(double& constKd) const
{
	constKd = fconstKd;
	return 1;
}

/*************************************************************************
* Protected
*************************************************************************/

/* recalculate constants */
void eNLHHTalpha::eComputeParameters(void)
{
	/* inherited */
	eLinearHHTalpha::eComputeParameters();

	/* element residual force coefficients */
	fconstMa = 1.0;
	fconstCv = 1.0 + falpha;
	fconstKd = 1.0 + falpha;
}
