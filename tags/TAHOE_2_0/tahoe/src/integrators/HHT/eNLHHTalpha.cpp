/* $Id: eNLHHTalpha.cpp,v 1.5 2004-07-15 08:30:28 paklein Exp $ */
/* created: paklein (10/17/1996) */
#include "eNLHHTalpha.h"
#include "ExceptionT.h"

using namespace Tahoe;

/* constructor */
eNLHHTalpha::eNLHHTalpha(double alpha):
	HHTalpha(alpha),
	eLinearHHTalpha(alpha)
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
