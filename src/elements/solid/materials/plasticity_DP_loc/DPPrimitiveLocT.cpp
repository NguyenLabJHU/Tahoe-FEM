
/* $Id: DPPrimitiveLocT.cpp,v 1.6 2005-04-08 19:22:46 raregue Exp $ */
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
	fkappa(-1.0),
	ffriction(-1.0),
	fdilation(-1.0),
	fH(0.0),
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

	ParameterT kappa(fkappa, "kappa");
	kappa.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(kappa);

	ParameterT friction(ffriction, "friction");
	friction.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(friction);

	list.AddParameter(fdilation, "dilation");
	
	list.AddParameter(fH, "H");

	ParameterT eta(fEta, "fluidity_eta");
	eta.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(eta);
}

/* accept parameter list */
void DPPrimitiveLocT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	fkappa = list.GetParameter("kappa");
	ffriction = list.GetParameter("friction");
	fdilation = list.GetParameter("dilation");
	fH = list.GetParameter("H");
	fEta = list.GetParameter("fluidity_eta");
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/*
 * Returns the value of the yield function given the
 * stress vector and state variables, where kappa
 * represents isotropic hardening.
 */
double DPPrimitiveLocT::YieldCondition(const dSymMatrixT& devstress, 
				const double meanstress, double kappa) const
{
	double kTemp;
	kTemp  = sqrt32*sqrt(devstress.ScalarProduct());
	kTemp += sqrt(3.0)*(-fkappa + ffriction*meanstress);
	kTemp += kappa;
	return   kTemp;
}

