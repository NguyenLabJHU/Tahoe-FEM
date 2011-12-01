/* $Id: VerletIntegrator.cpp,v 1.5 2011-12-01 21:11:39 bcyansfn Exp $ */
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
