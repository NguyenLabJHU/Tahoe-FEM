/* $Id: eExplicitCD.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (03/23/1997)                                          */
/* Element controller for an explicit 2nd order                           */
/* accurate, central difference time-stepping algorithm.                  */

#include "eExplicitCD.h"

#include "Environment.h"
#include "ExceptionCodes.h"

/* constructor */
eExplicitCD::eExplicitCD(void) { }

/* time integration parameters */
eControllerT::StatDynFlagT eExplicitCD::StaticDynamic(void) const { return kDynamic; }
eControllerT::ImpExpFlagT eExplicitCD::ImplicitExplicit(void) const { return kExplicit; }

/* return order time discretization */
int eExplicitCD::Order(void) const { return 2; }

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
