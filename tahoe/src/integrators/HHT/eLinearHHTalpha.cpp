/* $Id: eLinearHHTalpha.cpp,v 1.5 2004-07-15 08:30:28 paklein Exp $ */
/* created: paklein (10/14/1996) */
#include "eLinearHHTalpha.h"

using namespace Tahoe;

/* constructor */
eLinearHHTalpha::eLinearHHTalpha(double alpha):
	HHTalpha(alpha)
{

}

/* returns 1 if the algorithm requires M, C, or K and sets const equal
* to the coefficient for the linear combination of components in the
* element effective mass matrix */
int eLinearHHTalpha::FormM(double& constM) const
{
	constM = fconstM;
	return 1;
}

int eLinearHHTalpha::FormC(double& constC) const
{
	constC = fconstC;
	return 1;
}

int eLinearHHTalpha::FormK(double& constK) const
{
	constK = fconstK;
	return 1;
}

/* components of the internal force vector */
int eLinearHHTalpha::FormMa(double& constMa) const
{
#pragma unused(constMa)
	return 0;
}

int eLinearHHTalpha::FormCv(double& constCv) const
{
	constCv = 1.0;
	return 1;
}

int eLinearHHTalpha::FormKd(double& constKd) const
{
	constKd = 1.0;
	return 1;
}

/*************************************************************************
* Protected
*************************************************************************/

/* recalculate constants */
void eLinearHHTalpha::eComputeParameters(void)
{
	/* effective mass coefficients */
	fconstM = 1.0;
	fconstC = (1.0 + falpha)*fgamma*fdt;
	fconstK = (1.0 + falpha)*fbeta*fdt*fdt;
}
