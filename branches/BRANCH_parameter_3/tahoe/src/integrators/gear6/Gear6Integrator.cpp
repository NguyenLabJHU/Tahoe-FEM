/* $Id: Gear6Integrator.cpp,v 1.1.56.1 2004-04-08 07:33:41 paklein Exp $ */
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
