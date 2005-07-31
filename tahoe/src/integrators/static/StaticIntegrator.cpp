/* $Id: StaticIntegrator.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/14/1996)                                          */

#include "StaticIntegrator.h"
#include <iostream.h>

/* constructor */
StaticIntegrator::StaticIntegrator(ostream& out)
{
	out << "\n Static controller parameters: NONE" << endl;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void StaticIntegrator::ComputeParameters(void) { }
