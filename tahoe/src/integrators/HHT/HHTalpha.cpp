/* $Id: HHTalpha.cpp,v 1.6 2004-06-17 07:13:57 paklein Exp $ */
/* created: paklein (10/14/1996) */
#include "HHTalpha.h"

#include <iostream.h>

#include "ifstreamT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
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
