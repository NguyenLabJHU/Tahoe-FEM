/* $Id: HHTalpha.cpp,v 1.4 2002-10-20 22:48:08 paklein Exp $ */
/* created: paklein (10/14/1996) */

#include "HHTalpha.h"

#include <iostream.h>

#include "fstreamT.h"
#include "dArrayT.h"

/* constructor */

using namespace Tahoe;

HHTalpha::HHTalpha(ifstreamT& in, ostream& out, bool auto2ndorder):
	fAuto2ndOrder(auto2ndorder),
	fgamma(0.5),
	fbeta(0.25)
{
	/* read parameter */
	in >> falpha;
		
	/* reduce to single parameter fanmily */
	if (fAuto2ndOrder) Set2ndOrder(falpha);

	/* echo */
	out << " HHT-alpha time integration parameters:\n\n";
	out << " gamma . . . . . . . . . . . . . . . . . . . . . = " << fgamma << '\n';
	out << " beta. . . . . . . . . . . . . . . . . . . . . . = " << fbeta  << '\n';
	out << " alpha . . . . . . . . . . . . . . . . . . . . . = " << falpha << endl;

	/* range checks (incomplete) */
	if (falpha < 0.0 || falpha > 1.0) throw ExceptionT::kBadInputValue;
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
