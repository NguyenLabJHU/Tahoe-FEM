/* $Id: IntegratorT.cpp,v 1.7 2011-12-01 21:11:39 bcyansfn Exp $ */
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
