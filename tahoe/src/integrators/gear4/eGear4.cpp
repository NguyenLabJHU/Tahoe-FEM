#include "eGear4.h"

#include "Environment.h"
#include "ExceptionT.h"

using namespace Tahoe;

/* constructor */
eGear4::eGear4(void) { }

/* returns 1 if the algorithm requires M, C, or K and sets const equal
* to the coefficient for the linear combination of components in the
* element effective mass matrix */
int eGear4::FormM(double& constM) const
{
	constM = 1.0;
	return 1;
}

int eGear4::FormC(double& constC) const
{
	constC = fconstC;
	return 1;
}

int eGear4::FormK(double& constK) const
{
	constK = 0.0;
	return 0;
}

/* components of the internal force vector */
int eGear4::FormMa(double& constMa) const
{
	constMa = 0.0;
	return 0;
}

int eGear4::FormCv(double& constCv) const
{
	constCv = 1.0;
	return 1;
}

int eGear4::FormKd(double& constKd) const
{
	constKd = 1.0;
	return 1;
}

/*************************************************************************
* Protected
*************************************************************************/

/* recalculate constants */
void eGear4::eComputeParameters(void)
{
	/* effective mass coefficients */
	fconstC = 0.5*fdt;
}
