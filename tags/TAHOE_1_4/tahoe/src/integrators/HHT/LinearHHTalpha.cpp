/* $Id: LinearHHTalpha.cpp,v 1.5 2002-11-14 17:05:49 paklein Exp $ */
/* created: paklein (10/11/1996) */
#include "LinearHHTalpha.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "TimeManagerT.h"

using namespace Tahoe;

/* constructor */
LinearHHTalpha::LinearHHTalpha(TimeManagerT& TM, ifstreamT& in, ostream& out, 
	bool auto2ndorder):
	HHTalpha(in, out, auto2ndorder),
	nLinearHHTalpha(in, out, auto2ndorder),
	eLinearHHTalpha(in, out, auto2ndorder),
	TimeBoss(TM)
{

}

/* take responsibility for forming the nodal contribution
* to the RHS vector:
*
*                     F(t_n+1+alpha)
*/
void LinearHHTalpha::FormNodalForce(NodeManagerT* nodeboss) const
{
	/* shift time back */
	TimeBoss.ShiftTime(fTimeShift);
	
	/* form nodal contribution to RHS */
//	nodeboss->FormRHS();
#pragma unused(nodeboss)
#pragma message("LinearHHTalpha::FormNodalForce: need this????")
	
	/* reset the time */
	TimeBoss.ResetTime();
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void LinearHHTalpha::ComputeParameters(void)
{
	nComputeParameters();
	eComputeParameters();
	
	fTimeShift = falpha*fdt;
}
