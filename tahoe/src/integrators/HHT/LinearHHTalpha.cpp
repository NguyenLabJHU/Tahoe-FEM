/* $Id: LinearHHTalpha.cpp,v 1.6 2004-07-15 08:30:27 paklein Exp $ */
/* created: paklein (10/11/1996) */
#include "LinearHHTalpha.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "TimeManagerT.h"

using namespace Tahoe;

/* constructor */
LinearHHTalpha::LinearHHTalpha(double alpha):
	HHTalpha(alpha),
	nLinearHHTalpha(alpha),
	eLinearHHTalpha(alpha),
	fTimeBoss(NULL)
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
	fTimeBoss->ShiftTime(fTimeShift);
	
	/* form nodal contribution to RHS */
//	nodeboss->FormRHS();
#pragma unused(nodeboss)
#pragma message("LinearHHTalpha::FormNodalForce: need this????")
	
	/* reset the time */
	fTimeBoss->ResetTime();
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
