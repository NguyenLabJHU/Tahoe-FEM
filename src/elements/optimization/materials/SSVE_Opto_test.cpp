 /* $Id: SSVE_Opto_test.cpp,v 1.1 2009-04-23 03:03:52 thao Exp $ */
#include "SSVE_Opto_test.h"
#include "ExceptionT.h"

using namespace Tahoe;

SSVE_Opto_test::SSVE_Opto_test(void):
	ParameterInterfaceT("ssve_opto_test")
{
}

void SSVE_Opto_test::DefineParameters(ParameterListT& list) const
{
	/*inherited*/
	SSSolidMatT::DefineParameters(list);
	list.SetDescription("parameter order: [kappa, mu_eq, mu_neq_1, eta_s_1, ... , mu_neq_N, eta_s_N]");
}

/* accept parameter list */
void SSVE_Opto_test::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SSVE_Opto_test::TakeParameterList";
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

	fnum_params = 2; /*equilibrium params*/	
	fnum_params += 2*fnumprocess; /*neq params*/
		
	fParamGrads.Dimension(dSymMatrixT::NumValues(NumSD()), fnum_params);
	fconstraint_grad.Dimension(fnum_params);
	fconstraint_grad = 0.0;
	flabels.Dimension(fnum_params);


	flabels[0] = "kappa";
	flabels[1] = "mu_eq";
	int dex = 2;
	for (int i = 0; i < fnumprocess; i++)
	{
		StringT mu, tau;
		mu = "mu_neq_";
		tau = "eta_s_";
		int j = i+1;
		flabels[dex++] = mu.Append(j);
		flabels[dex++] = tau.Append(j);
	}
}

const double SSVE_Opto_test::constraint(void)
{
	double constraint = 0.0;
	
	for (int i = 0; i < fnumprocess; i++)
	{
		const double mu_neq = fMu[i+1];
		const double tau = ftauS[i];		
		double etaS = tau*mu_neq;

		double temp = fdt/tau;
		constraint += temp*temp;
	}
	constraint = -1.0;
	return(constraint);

}

const dArrayT& SSVE_Opto_test::constraint_grad(void)
{
	fconstraint_grad = 0.0;

	double dex = 2;
	for (int i = 0; i < fnumprocess; i++)
	{
		const double mu_neq = fMu[i+1];
		const double tau = ftauS[i];
		double etaS = tau*mu_neq;

		double temp = fdt/tau;
		fconstraint_grad[dex++] += 2.0*temp*temp/mu_neq;
		fconstraint_grad[dex++] -= 2.0*temp*temp/etaS;
	}
	fconstraint_grad = 0.0;
	return(fconstraint_grad);
}

const dArray2DT& SSVE_Opto_test::ds_ij_dlambda_q(void)
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
	
	double dex = 2;
	for (int i = 0; i < fnumprocess; i++)
	{
		const dSymMatrixT& qn = fdevQ_n[i];
		
		const double mu_neq = fMu[i+1];
		const double tau = ftauS[i];
		double alphaS = exp(-0.5*fdt/tau);
		double betaS = exp(-fdt/tau);
		
		double etaS = tau*mu_neq;
		double g1 = alphaS*(2.0-fdt/tau);
		double g2 = betaS*fdt/etaS;
		double g2b = betaS/mu_neq;

		fParamGrads(0,dex) = g1*fStrain3D[0]+(g2b-g2)*qn[0];
		fParamGrads(1,dex) = g1*fStrain3D[1]+(g2b-g2)*qn[1];
		fParamGrads(2,dex) = g1*fStrain3D[2]+(g2b-g2)*qn[2];
		fParamGrads(3,dex) = g1*fStrain3D[3]+(g2b-g2)*qn[3];
		fParamGrads(4,dex) = g1*fStrain3D[4]+(g2b-g2)*qn[4];
		fParamGrads(5,dex) = g1*fStrain3D[5]+(g2b-g2)*qn[5];
		dex++;
		
		double g3 = fdt/(tau*tau)*alphaS;  
		double g4 = fdt/tau/etaS*betaS;
		
		fParamGrads(0,dex) = g3*fStrain3D[0]+g4*qn[0];
		fParamGrads(1,dex) = g3*fStrain3D[1]+g4*qn[1];
		fParamGrads(2,dex) = g3*fStrain3D[2]+g4*qn[2];
		fParamGrads(3,dex) = g3*fStrain3D[3]+g4*qn[3];
		fParamGrads(4,dex) = g3*fStrain3D[4]+g4*qn[4];
		fParamGrads(5,dex) = g3*fStrain3D[5]+g4*qn[5];
		dex++;
		
/*	cout << "paramgrads: \n"<<fParamGrads<<endl;
	cout << "fStrain3D: \n"<<fStrain3D<<endl;
	cout << "qn: \n"<<qn<<endl;*/
	}
	
	return(fParamGrads);
}

