
/* $Id: DPPrimitiveLocT.cpp,v 1.5 2004-09-10 01:07:59 cfoster Exp $ */
/* created: myip (06/01/1999)                                             */

/* Base class for Druker-Prager, nonassociative, small strain,
   pressure dependent plasticity model with linear isotropic hardening
   and localization 
   */
#include "DPPrimitiveLocT.h"

#include "dSymMatrixT.h"
#include <math.h>

using namespace Tahoe;

const double sqrt23 = sqrt(2.0/3.0);
const double sqrt32 = sqrt(3.0/2.0);

/* constructor */
DPPrimitiveLocT::DPPrimitiveLocT(void): 
	ParameterInterfaceT("DP_Loc_primitive"),
	falpha_bar(-1.0),
	ffriction(-1.0),
	fdilation(-1.0),
	fH_prime(0.0),
	fH_delta(1.0),
	fEta(-1.0)
{

}

/* destructor */
DPPrimitiveLocT::~DPPrimitiveLocT(void) { }

/* describe the parameters needed by the interface */
void DPPrimitiveLocT::DefineParameters(ParameterListT& list) const
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
	list.AddParameter(fH_delta, "H_delta");

	ParameterT eta(fEta, "fluidity_eta");
	eta.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(eta);
}

/* accept parameter list */
void DPPrimitiveLocT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	falpha_bar = list.GetParameter("alpha_bar");
	ffriction = list.GetParameter("friction");
	fdilation = list.GetParameter("dilation");
	fH_prime = list.GetParameter("H_prime");
	fH_delta = list.GetParameter("H_delta");
	fEta = list.GetParameter("fluidity_eta");
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/*
 * Returns the value of the yield function given the
 * stress vector and state variables, where alpha
 * represents isotropic hardening.
 */
double DPPrimitiveLocT::YieldCondition(const dSymMatrixT& devstress, 
				const double meanstress, double alpha) const
{
	double kTemp;
	kTemp  = sqrt32*sqrt(devstress.ScalarProduct());
	kTemp += sqrt(3.0)*(-falpha_bar + ffriction*meanstress);
	kTemp += alpha;
	return   kTemp;
}

