/* $Id: ExplicitCDIntegrator.cpp,v 1.1.1.1.10.1 2002-06-27 18:02:30 cjkimme Exp $ */
/* created: paklein (03/23/1997)                                          */
/* Integrator for an explicit 2nd order accurate, central difference      */
/* time-stepping algorithm.                                               */

#include "ExplicitCDIntegrator.h"

#include <iostream.h>

/* constructor */

using namespace Tahoe;

ExplicitCDIntegrator::ExplicitCDIntegrator(ostream& out)
{
	out << "\n Explicit central-difference parameters: NONE" << endl;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void ExplicitCDIntegrator::ComputeParameters(void)
{
	nComputeParameters();
	eComputeParameters();
}
