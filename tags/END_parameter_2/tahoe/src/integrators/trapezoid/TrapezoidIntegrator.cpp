/* $Id: TrapezoidIntegrator.cpp,v 1.2.40.1 2004-01-28 01:34:09 paklein Exp $ */
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
