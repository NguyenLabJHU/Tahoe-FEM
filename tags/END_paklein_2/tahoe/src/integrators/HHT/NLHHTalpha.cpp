/* $Id: NLHHTalpha.cpp,v 1.5 2002-09-12 17:49:50 paklein Exp $ */
/* created: paklein (10/11/1996) */

#include "NLHHTalpha.h"

#include <iostream.h>
#include "toolboxConstants.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "TimeManagerT.h"

/* constructor */

using namespace Tahoe;

NLHHTalpha::NLHHTalpha(TimeManagerT& TM, ifstreamT& in, ostream& out,
	bool auto2ndorder):
	HHTalpha(in, out, auto2ndorder),
	nNLHHTalpha(in, out, auto2ndorder),
	eNLHHTalpha(in, out, auto2ndorder),
	TimeBoss(TM)
{
	/* time-shifting not supported for NL */
	if (falpha*falpha > kSmall)
	{
		/* reset parameters */
		Set2ndOrder(0.0);

		/* echo */
		out << " NOTE: (alpha != 0.0) not supported for nonlinear HHT\n\n";
		out << " gamma . . . . . . . . . . . . . . . . . . . . . = " << fgamma << '\n';
		out << " beta. . . . . . . . . . . . . . . . . . . . . . = " << fbeta  << '\n';
		out << " alpha . . . . . . . . . . . . . . . . . . . . . = " << falpha << endl;
	}
}

/* take responsibility for forming the nodal contribution
* to the RHS vector:
*
*                     F(t_n+1+alpha)
*/
void NLHHTalpha::FormNodalForce(NodeManagerT* nodeboss) const
{
	/* shift time back */
	TimeBoss.ShiftTime(fTimeShift);
	
	/* form nodal contribution to RHS */
//	nodeboss->FormRHS();
#pragma unused(nodeboss)
#pragma message("NLHHTalpha::FormNodalForce: need this???")
	
	/* reset the time */
	TimeBoss.ResetTime();
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void NLHHTalpha::ComputeParameters(void)
{
	nComputeParameters();
	eComputeParameters();
	
	fTimeShift = falpha*fdt;
}
