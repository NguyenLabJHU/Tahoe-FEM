/* $Id: TrapezoidIntegrator.cpp,v 1.3 2004-07-15 08:30:53 paklein Exp $ */
/* created: paklein (10/03/1999) */
#include "TrapezoidIntegrator.h"

using namespace Tahoe;

/* constructor */
TrapezoidIntegrator::TrapezoidIntegrator(void) { }

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void TrapezoidIntegrator::ComputeParameters(void)
{
	nComputeParameters();
	eComputeParameters();
}
