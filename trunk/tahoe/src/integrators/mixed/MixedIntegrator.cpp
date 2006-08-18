/* $Header: /home/cvs/t/ta/tahoe/tahoe/src/integrators/mixed/MixedIntegrator.cpp,v 1.2 2006-08-18 01:15:50 a-kopacz Exp $ */
/* created: a-kopacz (08/08/2006) */

#include "MixedIntegrator.h"

using namespace Tahoe;

/* constructor */
MixedIntegrator::MixedIntegrator(void) { }

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void MixedIntegrator::ComputeParameters(void)
{
	nComputeParameters();
	eComputeParameters();
}
