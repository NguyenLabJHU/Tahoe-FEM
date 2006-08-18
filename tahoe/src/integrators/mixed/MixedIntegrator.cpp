/* $Header: /home/cvs/t/ta/tahoe/tahoe/src/integrators/mixed/MixedIntegrator.cpp,v 1.4 2006-08-18 21:52:30 tdnguye Exp $ */
/* created: a-kopacz (08/08/2006) */

#include "MixedIntegrator.h"

using namespace Tahoe;

/* constructor */
MixedIntegrator::MixedIntegrator(void) 
{
	/*Initialize values*/
	fNumDOF = 0; 
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void MixedIntegrator::ComputeParameters(void)
{
	/*calling compute parameters*/
	nComputeParameters();
	eComputeParameters();
}
