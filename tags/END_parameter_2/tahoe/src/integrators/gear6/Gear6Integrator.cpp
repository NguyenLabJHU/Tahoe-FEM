/* $Id: Gear6Integrator.cpp,v 1.1.40.1 2004-01-28 01:34:06 paklein Exp $ */
#include "Gear6Integrator.h"

using namespace Tahoe;

/* constructor */
Gear6Integrator::Gear6Integrator(void) { }

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void Gear6Integrator::ComputeParameters(void)
{
	nComputeParameters();
	eComputeParameters();
}
