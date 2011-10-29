/* $Id: IntegratorT.cpp,v 1.6.84.1 2011-10-29 06:09:10 bcyansfn Exp $ */
/* created: paklein (10/14/1996) */
#include "IntegratorT.h"
#include "dArrayT.h"
#include <iostream>

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
