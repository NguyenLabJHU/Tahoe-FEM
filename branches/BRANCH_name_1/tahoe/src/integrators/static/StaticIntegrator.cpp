/* $Id: StaticIntegrator.cpp,v 1.1.1.1.10.1 2002-06-27 18:02:31 cjkimme Exp $ */
/* created: paklein (10/14/1996)                                          */

#include "StaticIntegrator.h"
#include <iostream.h>

/* constructor */

using namespace Tahoe;

StaticIntegrator::StaticIntegrator(ostream& out)
{
	out << "\n Static controller parameters: NONE" << endl;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void StaticIntegrator::ComputeParameters(void) { }
