/* $Id: eLinearHHTalpha.cpp,v 1.4 2002-10-20 22:48:08 paklein Exp $ */
/* created: paklein (10/14/1996) */

#include "eLinearHHTalpha.h"

#include "Environment.h"
#include "ExceptionT.h"

/* constructor */

using namespace Tahoe;

eLinearHHTalpha::eLinearHHTalpha(ifstreamT& in, ostream& out, bool auto2ndorder):
	HHTalpha(in, out, auto2ndorder)
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
