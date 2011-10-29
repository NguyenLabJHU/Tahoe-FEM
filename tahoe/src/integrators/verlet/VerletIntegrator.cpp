/* $Id: VerletIntegrator.cpp,v 1.4.38.1 2011-10-29 06:09:10 bcyansfn Exp $ */
#include "VerletIntegrator.h"
#include <iostream>

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
