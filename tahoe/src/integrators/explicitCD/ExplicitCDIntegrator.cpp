/* $Id: ExplicitCDIntegrator.cpp,v 1.3 2004-07-15 08:30:38 paklein Exp $ */
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
