/* $Id: ExplicitCDIntegrator.cpp,v 1.2.56.1 2004-04-08 07:33:40 paklein Exp $ */
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
