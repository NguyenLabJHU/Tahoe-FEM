/* $Id: TrapezoidIntegrator.cpp,v 1.2.56.1 2004-04-08 07:33:44 paklein Exp $ */
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
