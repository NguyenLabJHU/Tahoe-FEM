/* $Header: */
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
