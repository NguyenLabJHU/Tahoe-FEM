 /* $Id: SSVisco_Optimize.cpp,v 1.1 2009-04-23 03:03:52 thao Exp $ */
#include "SSVisco_Optimize.h"
#include "ExceptionT.h"
#include "ParameterContainerT.h"
//#include "SSLinearVE3D.h"

using namespace Tahoe;

SSVisco_Optimize::SSVisco_Optimize(void):
	ParameterInterfaceT("linear_prony_optimize_ve")
{
}

/* information about subordinate parameter lists */
void SSVisco_Optimize::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("ss_opto_eq_params", ParameterListT::Once);
	sub_list.AddSub("ss_opto_neq_params", ParameterListT::OnePlus);

}

/* describe the parameters needed by the interface */
ParameterInterfaceT* SSVisco_Optimize::NewSub(const StringT& name) const
{
	/* common limit */
	LimitT positive(0.0, LimitT::Lower);

	if (name == "ss_opto_eq_params")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT moduli("opto_params_eq");
		moduli.AddParameter(ParameterT::Double, "mu_EQ");
		moduli.AddParameter(ParameterT::Double, "kappa_EQ");
		choice->AddSub(moduli, ParameterListT::Once, false);
		
		return choice;
	}
	else if (name == "ss_opto_neq_params")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		choice->SetDescription("Assumes elastic bulk response. kappa_NEQ=0.0");

		ParameterContainerT moduli("opto_params_neq");
		moduli.AddParameter(ParameterT::Double, "mu_NEQ");

		ParameterT kappa_neq(ParameterT::Double, "kappa_NEQ");
		kappa_neq.SetDefault(0.0);
		kappa_neq.AddLimit(0.0, LimitT::Only);
		moduli.AddParameter(kappa_neq);

		moduli.AddParameter(ParameterT::Double, "tau_shear");

		ParameterT tau_b(ParameterT::Double, "tau_bulk");
		tau_b.SetDefault(10000.0);
		tau_b.AddLimit(10000.0, LimitT::LowerInclusive);
		moduli.AddParameter(tau_b);

		choice->AddSub(moduli, ParameterListT::Once, false);

		return choice;
	}
	else return (SSSolidMatT::NewSub(name));
}

void SSVisco_Optimize::DefineParameters(ParameterListT& list) const
{
	/*inherited*/
	SSSolidMatT::DefineParameters(list);
	list.SetDescription("parameter order: [kappa, mu_eq, mu_neq_1, tau_s_1, ... , mu_neq_N, tau_s_N]");
}

/* accept parameter list */
void SSVisco_Optimize::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SSVisco_Optimize::TakeParameterList";
	/* inherited */
	SSSolidMatT::TakeParameterList(list);

	fnumprocess = list.NumLists("ss_opto_neq_params");

	if (NumSD() == 2 && Constraint() == SolidMaterialT::kPlaneStress && fnumprocess >1)
		ExceptionT::GeneralFail(caller, "Plane stress formulation only implemented for 1 neq process");
	
	fMu.Dimension(fnumprocess+1);
	fKappa.Dimension(fnumprocess+1);
	ftauS.Dimension(fnumprocess);
	ftauB.Dimension(fnumprocess);
	
	const ParameterListT& eq_params = list.GetListChoice(*this, "ss_opto_eq_params");
	fMu[0] = eq_params.GetParameter("mu_EQ");
	fKappa[0] = eq_params.GetParameter("kappa_EQ");
	
	for (int i = 0; i < fnumprocess; i++)
	{
		const ParameterListT& neq_params = list.GetListChoice(*this, "ss_opto_neq_params",i);
		fMu[i+1] = neq_params.GetParameter("mu_NEQ");
		fKappa[i+1] = neq_params.GetParameter("kappa_NEQ");
		ftauS[i] = neq_params.GetParameter("tau_shear");
		ftauB[i] = neq_params.GetParameter("tau_bulk");
	}
	
	/* dimension work space */
	fStrain3D.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));
	fStress.Dimension(NumSD());

	SS_Visc_Support::SetSateVariables();

	fnum_params = 2; /*equilibrium params*/	
	fnum_params += 2*fnumprocess; /*neq params*/
		
	fParamGrads.Dimension(dSymMatrixT::NumValues(NumSD()), fnum_params);
	fconstraint_grad.Dimension(fnum_params);
	flabels.Dimension(fnum_params);

	fParamGrads = 0.0;
	fconstraint_grad = 0.0;
	
	flabels[0] = "kappa_eq";
	flabels[1] = "mu_eq";
	int dex = 2;
	for (int i = 0; i < fnumprocess; i++)
	{
		StringT mu, tau;
		mu = "mu_neq_";
		tau = "tau_s_";
		int j = i+1;
		flabels[dex++] = mu.Append(j);
		flabels[dex++] = tau.Append(j);
	}
}

const dArray2DT& SSVisco_Optimize::ds_ij_dlambda_q(void)
{
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	double third =1.0/3.0;
	const dSymMatrixT& strain = e();

	fStrain3D = strain;

	double I1 = strain.Trace();
	
	fStrain3D[0] -= third * I1;
	fStrain3D[1] -= third * I1;
	fStrain3D[2] -= third * I1;

	fParamGrads = 0.0;
	fParamGrads(0,0) = I1;
	fParamGrads(1,0) = I1;
	fParamGrads(2,0) = I1;
	  
	fParamGrads(0,1) = 2.0*fStrain3D[0];
	fParamGrads(1,1) = 2.0*fStrain3D[1];
	fParamGrads(2,1) = 2.0*fStrain3D[2];
	fParamGrads(3,1) = 2.0*fStrain3D[3];
	fParamGrads(4,1) = 2.0*fStrain3D[4];
	fParamGrads(5,1) = 2.0*fStrain3D[5];
	
	/*ep_n+1-ep_n*/
	const dSymMatrixT& strain_last = e_last();
	fStrain3D -= strain_last;

	double I1_last = strain_last.Trace();	
	fStrain3D[0] += third * I1_last;
	fStrain3D[1] += third * I1_last;
	fStrain3D[2] += third * I1_last;
	
	/*derivatives wrt to noneq params*/
	double dex = 2;
	for (int i = 0; i < fnumprocess; i++)
	{
		const dSymMatrixT& qn = fdevQ_n[i];

		/*mu_neq*/
		const double mu_neq = fMu[i+1];
		const double tau = ftauS[i];
		double alphaS = exp(-0.5*fdt/tau);
		double betaS = exp(-fdt/tau);
		

		fParamGrads(0,dex) = 2.0*alphaS*fStrain3D[0] + betaS/mu_neq*qn[0];
		fParamGrads(1,dex) = 2.0*alphaS*fStrain3D[1] + betaS/mu_neq*qn[1];
		fParamGrads(2,dex) = 2.0*alphaS*fStrain3D[2] + betaS/mu_neq*qn[2];
		fParamGrads(3,dex) = 2.0*alphaS*fStrain3D[3] + betaS/mu_neq*qn[3];
		fParamGrads(4,dex) = 2.0*alphaS*fStrain3D[4] + betaS/mu_neq*qn[4];
		fParamGrads(5,dex) = 2.0*alphaS*fStrain3D[5] + betaS/mu_neq*qn[5];
		dex++;
		
		/*tau_s*/
		double coeff =fdt/(tau*tau)*alphaS;
				
		fParamGrads(0,dex) = coeff* (mu_neq*fStrain3D[0]+alphaS*qn[0]);
		fParamGrads(1,dex) = coeff* (mu_neq*fStrain3D[1]+alphaS*qn[1]);
		fParamGrads(2,dex) = coeff* (mu_neq*fStrain3D[2]+alphaS*qn[2]);
		fParamGrads(3,dex) = coeff* (mu_neq*fStrain3D[3]+alphaS*qn[3]);
		fParamGrads(4,dex) = coeff* (mu_neq*fStrain3D[4]+alphaS*qn[4]);
		fParamGrads(5,dex) = coeff* (mu_neq*fStrain3D[5]+alphaS*qn[5]);
		dex++;
		
//		cout << "\nqn: "<<qn;
	}

	return(fParamGrads);
}

