/* $Id: eExplicitCD.cpp,v 1.4 2002-10-20 22:48:10 paklein Exp $ */
/* created: paklein (03/23/1997) */

#include "eExplicitCD.h"

#include "Environment.h"
#include "ExceptionT.h"

/* constructor */

using namespace Tahoe;

eExplicitCD::eExplicitCD(void) { }

/* returns 1 if the algorithm requires M, C, or K and sets const equal
* to the coefficient for the linear combination of components in the
* element effective mass matrix */
int eExplicitCD::FormM(double& constM) const
{
	constM = 1.0;
	return 1;
}

int eExplicitCD::FormC(double& constC) const
{
	constC = fconstC;
	return 1;
}

int eExplicitCD::FormK(double& constK) const
{
	constK = 0.0;
	return 0;
}

/* components of the internal force vector */
int eExplicitCD::FormMa(double& constMa) const
{
	constMa = 0.0;
	return 0;
}

int eExplicitCD::FormCv(double& constCv) const
{
	constCv = 1.0;
	return 1;
}

int eExplicitCD::FormKd(double& constKd) const
{
	constKd = 1.0;
	return 1;
}

/*************************************************************************
* Protected
*************************************************************************/

/* recalculate constants */
void eExplicitCD::eComputeParameters(void)
{
	/* effective mass coefficients */
	fconstC = 0.5*fdt;
}
