/* $Header: */
/* created: a-kopacz (08/08/2006) */

#include "eMixed.h"

#include "Environment.h"
#include "ExceptionT.h"

/* constructor */

using namespace Tahoe;

eMixed::eMixed(void) { }

/* returns 1 if the algorithm requires M, C, or K and sets const equal
* to the coefficient for the linear combination of components in the
* element effective mass matrix */
int eMixed::FormM(double& constM) const
{
#pragma unused(constM)
	return 0;
}

int eMixed::FormC(double& constC) const
{
	constC = fconstC;
	return 1;
}

int eMixed::FormK(double& constK) const
{
	constK = fconstK;
	return 1;
}

/* components of the internal force vector */
int eMixed::FormMa(double& constMa) const
{
#pragma unused(constMa)
	return 0;
}

int eMixed::FormCv(double& constCv) const
{
	constCv = 1.0;
	return 1;
}

int eMixed::FormKd(double& constKd) const
{
	constKd = 1.0;
	return 1;
}

/*************************************************************************
* Protected
*************************************************************************/

/* recalculate constants */
void eMixed::eComputeParameters(void)
{
	/* effective mass coefficients */
	fconstC = 1.0;
	fconstK = 0.5*fdt;
}
