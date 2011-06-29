/* $Id: eVerlet.cpp,v 1.4 2004-07-15 08:30:57 paklein Exp $ */
#include "eVerlet.h"

using namespace Tahoe;

/* constructor */
eVerlet::eVerlet(void) { }

/* returns 1 if the algorithm requires M, C, or K and sets const equal
* to the coefficient for the linear combination of components in the
* element effective mass matrix */
int eVerlet::FormM(double& constM) const
{
	constM = 1.0;
	return 1;
}

int eVerlet::FormC(double& constC) const
{
	constC = fconstC;
	return 1;
}

int eVerlet::FormK(double& constK) const
{
	constK = 0.0;
	return 0;
}

/* components of the internal force vector */
int eVerlet::FormMa(double& constMa) const
{
	constMa = 0.0;
	return 0;
}

int eVerlet::FormCv(double& constCv) const
{
	constCv = 1.0;
	return 1;
}

int eVerlet::FormKd(double& constKd) const
{
	constKd = 1.0;
	return 1;
}

/*************************************************************************
* Protected
*************************************************************************/

/* recalculate constants */
void eVerlet::eComputeParameters(void)
{
	/* effective mass coefficients */
	fconstC = 0.5*fdt;
}
