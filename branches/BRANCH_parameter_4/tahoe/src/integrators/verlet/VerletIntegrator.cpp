/* $Id: VerletIntegrator.cpp,v 1.3.64.1 2004-07-06 06:54:36 paklein Exp $ */
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
