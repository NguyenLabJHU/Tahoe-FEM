/* $Id: StaticIntegrator.cpp,v 1.2 2002-07-02 19:55:09 cjkimme Exp $ */
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
