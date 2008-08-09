/* $Id: SSLinearVE2D.cpp,v 1.7 2008-08-09 15:00:36 tdnguye Exp $ */
/* created: TDN (5/31/2001) */
#include "SSLinearVE2D.h"
#include "SSMatSupportT.h"

#include <math.h>
#include <iostream.h>

#include "ExceptionT.h"

using namespace Tahoe;

const double third = 1.0/3.0;
const int kNumOutputVar = 1;
static const char* Labels[kNumOutputVar] = {"Dvisc"};

SSLinearVE2D::SSLinearVE2D(void):
	ParameterInterfaceT("linear_viscoelastic_2D")
{
	/*set default*/
	fConstraint = kPlaneStrain;
}	

/** information about subordinate parameter lists */
void SSLinearVE2D::DefineSubs(SubListT& sub_list) const
{
	/*override*/
	SSSolidMatT::DefineSubs(sub_list);
}

/** a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSLinearVE2D::NewSub(const StringT& name) const
{
	/*override*/
	SSSolidMatT::NewSub(name);
}


/* describe the parameters needed by the interface */
void SSLinearVE2D::DefineParameters(ParameterListT& list) const
{
	/* override */
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
void SSLinearVE2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSSolidMatT::TakeParameterList(list);

	/* dimension work space */
	fStrain3D.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));
	fStress.Dimension(NumSD());

	/* dimension params */
	fnumprocess = 1;

	fMu.Dimension(fnumprocess+1);
	fKappa.Dimension(fnumprocess+1);
	ftauS.Dimension(fnumprocess);
	ftauB.Dimension(fnumprocess);
	
	SS_Visc_Support::SetSateVariables();

	ftauS[0] = list.GetParameter("tau_shear");
	ftauB[0] = list.GetParameter("tau_bulk");

	/* elastic properties */
	fMu[kEquilibrium] = list.GetParameter("mu_EQ");
	fKappa[kEquilibrium] = list.GetParameter("kappa_EQ");
	fMu[kNonEquilibrium] = list.GetParameter("mu_NEQ");
	fKappa[kNonEquilibrium] = list.GetParameter("kappa_NEQ");
}
