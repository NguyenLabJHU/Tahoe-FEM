/* $Id: LinearStaticIntegrator.cpp,v 1.1.6.1 2002-06-27 18:02:31 cjkimme Exp $ */
/* created: paklein (10/14/1996) */

#include "LinearStaticIntegrator.h"
#include <iostream.h>

/* constructor */

using namespace Tahoe;

LinearStaticIntegrator::LinearStaticIntegrator(ostream& out)
{
	out << "\n Linear static controller parameters: NONE" << endl;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void LinearStaticIntegrator::ComputeParameters(void) { }
