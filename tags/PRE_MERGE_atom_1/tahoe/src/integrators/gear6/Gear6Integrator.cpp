#include "Gear6Integrator.h"
#include <iostream.h>

using namespace Tahoe;

/* constructor */
Gear6Integrator::Gear6Integrator(ostream& out)
{
	out << "\n Gear6 parameters: NONE" << endl;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void Gear6Integrator::ComputeParameters(void)
{
	nComputeParameters();
	eComputeParameters();
}
