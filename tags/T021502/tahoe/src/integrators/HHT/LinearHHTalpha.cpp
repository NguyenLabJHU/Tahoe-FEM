/* $Id: LinearHHTalpha.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/11/1996)                                          */

#include "LinearHHTalpha.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "NodeManagerT.h"
#include "TimeManagerT.h"

/* constructor */
LinearHHTalpha::LinearHHTalpha(TimeManagerT& TM, ifstreamT& in, ostream& out,
	int auto2ndorder):
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
	nodeboss->FormRHS();
	
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
