/* $Id: DPPrimitiveT.cpp,v 1.10 2004-07-15 08:28:48 paklein Exp $ */
/* created: myip (06/01/1999) */
#include "DPPrimitiveT.h"

#include "dSymMatrixT.h"
#include <math.h>

using namespace Tahoe;

const double sqrt23 = sqrt(2.0/3.0);
const double sqrt32 = sqrt(3.0/2.0);

/* constructor */
DPPrimitiveT::DPPrimitiveT(void): 
	ParameterInterfaceT("DP_primitive"),
	falpha_bar(-1.0),
	ffriction(-1.0),
	fdilation(-1.0),
	fH_prime(0.0)
{

}

/* destructor */
DPPrimitiveT::~DPPrimitiveT(void) { }

/* describe the parameters needed by the interface */
void DPPrimitiveT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT alpha_bar(falpha_bar, "alpha_bar");
	alpha_bar.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_bar);

	ParameterT friction(ffriction, "friction");
	friction.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(friction);

	list.AddParameter(fdilation, "dilation");
	list.AddParameter(fH_prime, "H_prime");
}

/* accept parameter list */
void DPPrimitiveT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	falpha_bar = list.GetParameter("alpha_bar");
	ffriction = list.GetParameter("friction");
	fdilation = list.GetParameter("dilation");
	fH_prime = list.GetParameter("H_prime");
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/*
 * Returns the value of the yield function given the
 * stress vector and state variables, where alpha
 * represents isotropic hardening.
 */
double DPPrimitiveT::YieldCondition(const dSymMatrixT& devstress, 
			const double meanstress, double alpha) const
{
  double kTemp;
  kTemp  = sqrt32*sqrt(devstress.ScalarProduct());
  kTemp += sqrt(3.0)*(-falpha_bar + ffriction*meanstress);
  kTemp += alpha;
  return   kTemp;
}

