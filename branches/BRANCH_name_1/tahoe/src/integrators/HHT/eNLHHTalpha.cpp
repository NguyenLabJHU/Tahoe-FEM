/* $Id: eNLHHTalpha.cpp,v 1.2.4.1 2002-06-27 18:02:27 cjkimme Exp $ */
/* created: paklein (10/17/1996) */

#include "eNLHHTalpha.h"
#include "ExceptionCodes.h"

/* constructor */

using namespace Tahoe;

eNLHHTalpha::eNLHHTalpha(ifstreamT& in, ostream& out, bool auto2ndorder):
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
