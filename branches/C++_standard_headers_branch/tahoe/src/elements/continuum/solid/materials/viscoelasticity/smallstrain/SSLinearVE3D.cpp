/* $Id: SSLinearVE3D.cpp,v 1.6 2008-08-09 15:00:37 tdnguye Exp $ */
/* created: TDN (5/31/2001) */
#include "SSLinearVE3D.h"
#include "SSMatSupportT.h"

#include <math.h>
#include <iostream.h>

#include "ExceptionT.h"

using namespace Tahoe;

const double third = 1.0/3.0;
const int kNumOutputVar = 1;
static const char* Labels[kNumOutputVar] = {"Dvisc"};

SSLinearVE3D::SSLinearVE3D(void):
	ParameterInterfaceT("linear_viscoelastic")
{

}	

/** information about subordinate parameter lists */
void SSLinearVE3D::DefineSubs(SubListT& sub_list) const
{
	/*override*/
	SSSolidMatT::DefineSubs(sub_list);
}

/** information about subordinate parameter lists */
ParameterInterfaceT* SSLinearVE3D::NewSub(const StringT& name) const
{
	/*override*/
	SSSolidMatT::NewSub(name);
}

/* describe the parameters needed by the interface */
void SSLinearVE3D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSSolidMatT::DefineParameters(list);

	/* common limit */
	LimitT positive(0.0, LimitT::Lower);

	/* relaxation times */
	ParameterT tau_shear(ParameterT::Double, "tau_shear");
	ParameterT tau_bulk(ParameterT::Double, "tau_bulk");
	tau_shear.AddLimit(positive);
	tau_bulk.AddLimit(positive);
	list.AddParameter(tau_shear);
	list.AddParameter(tau_bulk);

	/* elastic properties */
	ParameterT mu_EQ(ParameterT::Double, "mu_EQ");
	ParameterT kappa_EQ(ParameterT::Double, "kappa_EQ");
	ParameterT mu_NEQ(ParameterT::Double, "mu_NEQ");
	ParameterT kappa_NEQ(ParameterT::Double, "kappa_NEQ");
	mu_EQ.AddLimit(positive);
	kappa_EQ.AddLimit(positive);
	mu_NEQ.AddLimit(positive);
	kappa_NEQ.AddLimit(positive);
	list.AddParameter(mu_EQ);
	list.AddParameter(kappa_EQ);
	list.AddParameter(mu_NEQ);
	list.AddParameter(kappa_NEQ);
}

/* accept parameter list */
void SSLinearVE3D::TakeParameterList(const ParameterListT& list)
{
	/* override */
	SSSolidMatT::TakeParameterList(list);

	/* dimension params */
	fnumprocess = 1;
	fMu.Dimension(fnumprocess+1);
	fKappa.Dimension(fnumprocess+1);
	ftauS.Dimension(fnumprocess);
	ftauB.Dimension(fnumprocess);
	
	/* dimension work space */
	fStrain3D.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));
	fStress.Dimension(NumSD());

	SS_Visc_Support::SetSateVariables();

	/* relaxation times */
	ftauS[0] = list.GetParameter("tau_shear");
	ftauB[0] = list.GetParameter("tau_bulk");

	/* elastic properties */
	fMu[kEquilibrium] = list.GetParameter("mu_EQ");
	fKappa[kEquilibrium] = list.GetParameter("kappa_EQ");
	fMu[kNonEquilibrium] = list.GetParameter("mu_NEQ");
	fKappa[kNonEquilibrium] = list.GetParameter("kappa_NEQ");
}
