/* $Id: TrapezoidIntegrator.cpp,v 1.1.4.1 2002-06-27 18:02:32 cjkimme Exp $ */
/* created: paklein (10/03/1999) */

#include "TrapezoidIntegrator.h"

#include <iostream.h>

/* constructor */

using namespace Tahoe;

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
