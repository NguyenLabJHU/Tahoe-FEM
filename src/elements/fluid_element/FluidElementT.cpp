/* $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/fluid_element/FluidElementT.cpp,v 1.3 2006-07-06 16:42:45 a-kopacz Exp $ */

#include "FluidElementT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "ifstreamT.h"

using namespace Tahoe;

/* constructor */
FluidElementT::FluidElementT(const ElementSupportT& support):
	ContinuumElementT(support)
{
	SetName("fluid");
}

/* destructor */
FluidElementT::~FluidElementT(void)
{

}


/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct the effective mass matrix */
void FluidElementT::LHSDriver(GlobalT::SystemTypeT sys_type)
{

}

void FluidElementT::RHSDriver(void)
{
  
}

/***********************************************************************
 * Private
 ***********************************************************************/
