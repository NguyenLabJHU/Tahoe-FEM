/* $Id: IntegratorT.cpp,v 1.6 2003-05-20 10:26:13 paklein Exp $ */
/* created: paklein (10/14/1996) */
#include "IntegratorT.h"
#include "dArrayT.h"
#include <iostream.h>

using namespace Tahoe;

/* constructor */
IntegratorT::IntegratorT(void): fdt(-1.0) { }

/* destructor */
IntegratorT::~IntegratorT(void) { }

/* set the time step size */
void IntegratorT::SetTimeStep(double timestep)
{
	/* check */
	if (timestep < 0.0) throw ExceptionT::kGeneralFail;
	
	/* re-calculate time-stepping parameters */
	fdt = timestep;
	ComputeParameters();
}

/* take control of the time at which the external force vector
* is formed.  Default action to simply call the NodeManagerT's
* FormRHS function */
void IntegratorT::FormNodalForce(NodeManagerT* nodeboss) const
{
#pragma unused(nodeboss)
#pragma message("IntegratorT::FormNodalForce: need this????")
//	nodeboss->FormRHS();
}
