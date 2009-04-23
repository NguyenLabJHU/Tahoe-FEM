 /* $Id: SSVE_test.cpp,v 1.1 2009-04-23 14:45:29 tdnguye Exp $ */
#include "SSVE_test.h"
#include "ExceptionT.h"
#include "ParameterContainerT.h"
#include "SSMatSupportT.h"

#include <math.h>
#include <iostream.h>

#include "ExceptionT.h"

using namespace Tahoe;

const double third = 1.0/3.0;

using namespace Tahoe;

SSVE_test::SSVE_test(void):
	ParameterInterfaceT("ssve_test")
{
}	

/* information about subordinate parameter lists */
void SSVE_test::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("ssve_test_eq_params", ParameterListT::Once);
	sub_list.AddSub("ssve_test_neq_params", ParameterListT::OnePlus);

}

/* describe the parameters needed by the interface */
ParameterInterfaceT* SSVE_test::NewSub(const StringT& name) const
{

	/* common limit */
	LimitT positive(0.0, LimitT::Lower);

	if (name == "ssve_test_eq_params")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT moduli("ssve_params_eq");
		moduli.AddParameter(ParameterT::Double, "mu_EQ");
		moduli.AddParameter(ParameterT::Double, "kappa");
		choice->AddSub(moduli, ParameterListT::Once, false);
		
		return choice;
	}
	else if (name == "ssve_test_neq_params")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);

		ParameterContainerT moduli("ssve_params_neq");
		moduli.AddParameter(ParameterT::Double, "mu_NEQ");
		moduli.AddParameter(ParameterT::Double, "eta_S");

		choice->AddSub(moduli, ParameterListT::Once, false);

		return choice;
	}
	else return (SSSolidMatT::NewSub(name));
}

/* accept parameter list */
void SSVE_test::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SSVE_test::TakeParameterList";
	/* inherited */
	SSSolidMatT::TakeParameterList(list);

	fnumprocess = list.NumLists("ssve_test_neq_params");

	if (NumSD() == 2 && Constraint() == SolidMaterialT::kPlaneStress && fnumprocess >1)
		ExceptionT::GeneralFail(caller, "Plane stress formulation only implemented for 1 neq process");
	
	fMu.Dimension(fnumprocess+1);
	fKappa.Dimension(fnumprocess+1);
	ftauS.Dimension(fnumprocess);
	ftauB.Dimension(fnumprocess);
	
	const ParameterListT& eq_params = list.GetListChoice(*this, "ssve_test_eq_params");
	fMu[0] = eq_params.GetParameter("mu_EQ");
	fKappa[0] = eq_params.GetParameter("kappa");
	
	for (int i = 0; i < fnumprocess; i++)
	{
		const ParameterListT& neq_params = list.GetListChoice(*this, "ssve_test_neq_params",i);
		fMu[i+1] = neq_params.GetParameter("mu_NEQ");
		fKappa[i+1] = 0.0;
		double etaS = neq_params.GetParameter("eta_S");
		ftauS[i] = etaS/fMu[i+1];
		ftauB[i] = 100000;
	}
	
	/* dimension work space */
	fStrain3D.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));
	fStress.Dimension(NumSD());

	SS_Visc_Support::SetSateVariables();
}
