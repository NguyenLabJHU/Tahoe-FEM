/* $Id: VerletIntegrator.cpp,v 1.3.56.1 2004-04-08 07:33:46 paklein Exp $ */
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
