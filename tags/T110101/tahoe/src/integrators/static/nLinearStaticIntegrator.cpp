/* $Id: nLinearStaticIntegrator.cpp,v 1.1 2001-08-27 17:12:15 paklein Exp $ */
/* created: paklein (10/14/1996) */

#include "nLinearStaticIntegrator.h"
#include "dArray2DT.h"

/* constructor */
nLinearStaticIntegrator::nLinearStaticIntegrator(void) { };

/* predictor. Maps ALL degrees of freedom forward. */
void nLinearStaticIntegrator::Predictor(void)
{
#if __option(extended_errorcheck)
	if (fU[0] == NULL)
	{
		cout << "\n nStaticIntegrator::Predictor: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* clear all displacements */
	*fU[0] = 0.0;
}
