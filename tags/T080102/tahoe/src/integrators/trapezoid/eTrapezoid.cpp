/* $Id: eTrapezoid.cpp,v 1.3 2002-07-02 19:55:10 cjkimme Exp $ */
/* created: paklein (10/03/1999)                                          */

#include "eTrapezoid.h"

#include "Environment.h"
#include "ExceptionCodes.h"

/* constructor */

using namespace Tahoe;

eTrapezoid::eTrapezoid(void) { }

/* returns 1 if the algorithm requires M, C, or K and sets const equal
* to the coefficient for the linear combination of components in the
* element effective mass matrix */
int eTrapezoid::FormM(double& constM) const
{
#pragma unused(constM)
	return 0;
}

int eTrapezoid::FormC(double& constC) const
{
	constC = fconstC;
	return 1;
}

int eTrapezoid::FormK(double& constK) const
{
	constK = fconstK;
	return 1;
}

/* components of the internal force vector */
int eTrapezoid::FormMa(double& constMa) const
{
#pragma unused(constMa)
	return 0;
}

int eTrapezoid::FormCv(double& constCv) const
{
	constCv = 1.0;
	return 1;
}

int eTrapezoid::FormKd(double& constKd) const
{
	constKd = 1.0;
	return 1;
}

/*************************************************************************
* Protected
*************************************************************************/

/* recalculate constants */
void eTrapezoid::eComputeParameters(void)
{
	/* effective mass coefficients */
	fconstC = 1.0;
	fconstK = 0.5*fdt;
}
