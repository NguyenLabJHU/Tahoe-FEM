/* $Id: TrapezoidIntegrator.cpp,v 1.1 2002-04-02 23:20:42 paklein Exp $ */
/* created: paklein (10/03/1999) */

#include "TrapezoidIntegrator.h"

#include <iostream.h>

/* constructor */
TrapezoidIntegrator::TrapezoidIntegrator(ostream& out)
{
	out << "\n Trapezoid parameters: NONE" << endl;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void TrapezoidIntegrator::ComputeParameters(void)
{
	nComputeParameters();
	eComputeParameters();
}
