/* $Id: LinearStaticIntegrator.cpp,v 1.1 2001-08-27 17:12:14 paklein Exp $ */
/* created: paklein (10/14/1996) */

#include "LinearStaticIntegrator.h"
#include <iostream.h>

/* constructor */
LinearStaticIntegrator::LinearStaticIntegrator(ostream& out)
{
	out << "\n Linear static controller parameters: NONE" << endl;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void LinearStaticIntegrator::ComputeParameters(void) { }
