#include "VerletIntegrator.h"

#include <iostream.h>

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
