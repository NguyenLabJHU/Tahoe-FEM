#include "VerletIntegrator.h"
#include <iostream.h>

using namespace Tahoe;

/* constructor */
VerletIntegrator::VerletIntegrator(ostream& out)
{
	out << "\n Verlet parameters: NONE" << endl;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void VerletIntegrator::ComputeParameters(void)
{
	nComputeParameters();
	eComputeParameters();
}
