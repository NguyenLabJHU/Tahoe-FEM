/* $CVSHeader$ */

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
