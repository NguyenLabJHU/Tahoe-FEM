/* $Id: Gear6Integrator.cpp,v 1.2 2004-07-15 08:30:43 paklein Exp $ */
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
