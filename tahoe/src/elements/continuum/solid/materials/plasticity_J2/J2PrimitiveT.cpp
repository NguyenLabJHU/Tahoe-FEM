/* $Id: J2PrimitiveT.cpp,v 1.5 2004-07-15 08:28:54 paklein Exp $ */
/* created: paklein (02/17/1997) */
#include "J2PrimitiveT.h"
#include "dSymMatrixT.h"
#include <math.h>

using namespace Tahoe;

const double sqrt23 = sqrt(2.0/3.0);

/* constructor */
J2PrimitiveT::J2PrimitiveT(void):
	ParameterInterfaceT("J2_primitive"),
	fYield(0.0),
	ftheta(-1.0),
	fH_bar(-1.0)
{

}

/* destructor */
J2PrimitiveT::~J2PrimitiveT(void) { }

/* describe the parameters needed by the interface */
void J2PrimitiveT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT yield(fYield, "yield");
	yield.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(yield);

	ParameterT H(fH_bar, "H_bar");
	H.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(H);

	ParameterT t(ftheta, "theta");
	t.AddLimit(0.0, LimitT::LowerInclusive);
	t.AddLimit(1.0, LimitT::UpperInclusive);
	list.AddParameter(t);
}

/* accept parameter list */
void J2PrimitiveT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	fYield = list.GetParameter("yield");
	fH_bar = list.GetParameter("H_bar");
	ftheta = list.GetParameter("theta");
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* returns the value value of the yield function given the
* relative stress vector and state variables, where  alpha
* represents isotropic hardening.  NOTE: the relative stress
* should already contain the correction for any kinematic
* hardening. */
double J2PrimitiveT::YieldCondition(const dSymMatrixT& relstress,
	double alpha) const
{
	return sqrt(relstress.ScalarProduct()) - sqrt23*K(alpha);
}
