/* $Id: ExplicitCDIntegrator.cpp,v 1.2.64.1 2004-07-06 06:54:31 paklein Exp $ */
/* created: paklein (03/23/1997) */
#include "ExplicitCDIntegrator.h"

using namespace Tahoe;

/* constructor */
ExplicitCDIntegrator::ExplicitCDIntegrator(void) { }

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void ExplicitCDIntegrator::ComputeParameters(void)
{
	nComputeParameters();
	eComputeParameters();
}
