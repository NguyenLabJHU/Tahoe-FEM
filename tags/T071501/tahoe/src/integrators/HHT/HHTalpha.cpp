/* $Id: HHTalpha.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/14/1996)                                          */

#include "HHTalpha.h"

#include <iostream.h>

#include "fstreamT.h"
#include "dArrayT.h"

/* constructor */
HHTalpha::HHTalpha(ifstreamT& in, ostream& out, int auto2ndorder):
	fAuto2ndOrder(auto2ndorder),
	fgamma(0.5),
	fbeta(0.25)
{
	if (fAuto2ndOrder != 0 &&
	    fAuto2ndOrder != kHHTalphaAuto_O2) throw eGeneralFail;

	/* read parameter */
	in >> falpha;
		
	/* reduce to single parameter fanmily */
	if (fAuto2ndOrder == kHHTalphaAuto_O2) Set2ndOrder(falpha);

	/* echo */
	out << " HHT-alpha time integration parameters:\n\n";
	out << " gamma . . . . . . . . . . . . . . . . . . . . . = " << fgamma << '\n';
	out << " beta. . . . . . . . . . . . . . . . . . . . . . = " << fbeta  << '\n';
	out << " alpha . . . . . . . . . . . . . . . . . . . . . = " << falpha << endl;

	/* range checks (incomplete) */
	if (falpha < 0.0 || falpha > 1.0) throw eBadInputValue;
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
