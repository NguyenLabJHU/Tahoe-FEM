/* $Id: VerletIntegrator.cpp,v 1.3.40.1 2004-01-28 01:34:10 paklein Exp $ */
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
