/* $Id: VerletIntegrator.cpp,v 1.4 2004-07-15 08:30:57 paklein Exp $ */
#include "VerletIntegrator.h"
#include <iostream.h>

using namespace Tahoe;

/* constructor */
VerletIntegrator::VerletIntegrator(void) { }

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void VerletIntegrator::ComputeParameters(void)
{
	nComputeParameters();
	eComputeParameters();
}
