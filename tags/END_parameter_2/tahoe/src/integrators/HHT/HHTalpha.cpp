/* $Id: HHTalpha.cpp,v 1.5.16.2 2004-03-09 08:55:19 paklein Exp $ */
/* created: paklein (10/14/1996) */
#include "HHTalpha.h"

#include <iostream.h>

#include "fstreamT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
HHTalpha::HHTalpha(double alpha):
	fAuto2ndOrder(true),
	fgamma(0.5),
	fbeta(0.25),
	falpha(alpha)
{
	/* range checks (incomplete) */
	if (falpha > 0.0 || falpha < -1.0) ExceptionT::BadInputValue("HHTalpha::HHTalpha");
}

/*************************************************************************
* Protected
*************************************************************************/

/* set time integration to single parameter 2nd order */
void HHTalpha::Set2ndOrder(double alpha)
{
	/* unconditionally stable, 2nd order */
	falpha = alpha;
	fgamma = 0.5*(1.0 - 2.0*falpha);
	fbeta  = 0.25*(1 - falpha)*(1 - falpha);
}
