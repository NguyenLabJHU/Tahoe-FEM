#include "Gear4Integrator.h"
#include <iostream.h>

using namespace Tahoe;

/* constructor */
Gear4Integrator::Gear4Integrator(ostream& out)
{
	out << "\n Gear4 parameters: NONE" << endl;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void Gear4Integrator::ComputeParameters(void)
{
	nComputeParameters();
	eComputeParameters();
}
