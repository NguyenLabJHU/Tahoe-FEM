/* $Id: NLHHTalpha.cpp,v 1.6 2004-07-15 08:30:28 paklein Exp $ */
/* created: paklein (10/11/1996) */
#include "NLHHTalpha.h"

#include <iostream.h>
#include "toolboxConstants.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "TimeManagerT.h"

using namespace Tahoe;

/* constructor */
NLHHTalpha::NLHHTalpha(double alpha):
	HHTalpha(alpha),
	nNLHHTalpha(alpha),
	eNLHHTalpha(alpha),
	fTimeBoss(NULL)
{
#if 0
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
#endif
}

/* take responsibility for forming the nodal contribution
* to the RHS vector:
*
*                     F(t_n+1+alpha)
*/
void NLHHTalpha::FormNodalForce(NodeManagerT* nodeboss) const
{
	/* shift time back */
	fTimeBoss->ShiftTime(fTimeShift);
	
	/* form nodal contribution to RHS */
//	nodeboss->FormRHS();
#pragma unused(nodeboss)
#pragma message("NLHHTalpha::FormNodalForce: need this???")
	
	/* reset the time */
	fTimeBoss->ResetTime();
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
