/* $Id: LinearStaticIntegrator.cpp,v 1.2 2002-07-02 19:55:09 cjkimme Exp $ */
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
