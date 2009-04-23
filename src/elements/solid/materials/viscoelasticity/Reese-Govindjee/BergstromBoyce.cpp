/* $Id: BergstromBoyce.cpp,v 1.1 2009-04-23 15:00:25 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "BergstromBoyce.h"
#include "ifstreamT.h"
#include "ExceptionT.h"
#include <math.h>
#include <iostream.h>
#include <stdlib.h>
#include "ParameterContainerT.h"

#include "MooneyRivlin.h"
#include "NeoHookean.h"
#include "VWPotentialT.h"
#include "ArrudaBoyce.h"

#include "PowerLawT.h"
#include "LinearExponentialT.h"
#include "ScaledSinh.h"

using namespace Tahoe;

const double loge = log10(exp(1.0));
const double third = 1.0/3.0; 
const int kNumOutputVar = 3; 
static const char* Labels[kNumOutputVar] = {"lv_eff", "smag", "etaS"}; 

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
/* constructors */
BergstromBoyce::BergstromBoyce(void):
  ParameterInterfaceT("BergstromBoyce")
{
	fNumProcess = 1;
//	cout << setprecision(16);
}

/******************PROTECTED************************/
double BergstromBoyce::ShearViscosity(const double smag, const double sy, const double lv)
{
	/*calculate the temperature part*/
	double etaS;
	
	if (fViscType == kBergstromBoyce)
	{
		double strain_v = lv - 1.0 + fepsilon;
		etaS = fetaS0*exp(-smag/sy)*(pow(strain_v, fm));
	}
	else if (fViscType == kExpBergstromBoyce)
	{
		double strain_v = lv - 1.0;
		etaS =  fetaS0*exp(-smag/sy)*exp(fm*(lv -1.0));  
	}
	return(etaS);
}

double BergstromBoyce::DVisc_Dlamv(const double smag, const double sy, const double lv)
{
	/*calculate the temperature part*/
	double etaS, detaS;
	
	if (fViscType == kBergstromBoyce)
	{
		double strain_v = lv - 1.0 + fepsilon;
		etaS = fetaS0*exp(-smag/sy)*(pow(strain_v, fm));
		detaS = etaS*fm/strain_v;
	}
	else if (fViscType == kExpBergstromBoyce)
	{
		double strain_v = lv - 1.0;
		etaS = fetaS0*exp(-smag/sy)*exp(fm*(lv -1.0));  
		detaS = etaS*fm;
	}
	return(etaS);
}

/*************************************************************************
*	PUBLIC
**************************************************************************/
int BergstromBoyce::NumOutputVariables() const {return kNumOutputVar;} 

void BergstromBoyce::OutputLabels(ArrayT<StringT>& labels) const 
{ 
     /*allocates space for labels*/
     labels.Dimension(kNumOutputVar); 
  
     /*copy labels*/
     for (int i = 0; i< kNumOutputVar; i++) 
       labels[i] = Labels[i]; 
} 


void BergstromBoyce::ComputeOutput(dArrayT& output)
{
	/*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
	
	/*effective viscous stretch*/
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fC_v[0], false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();
	double lv0 = fEigs[0];
	double lv1 = fEigs[1];
	double lv2 = fEigs[2];
	
	/*magnitude of nonequilibrium stress*/
	const dMatrixT& F = MechanicalDeformation();
	if (NumSD() == 2)
	{
		fF3D[0] = F[0];
		fF3D[1] = F[1];
		fF3D[2] = 0.0;
	    
		fF3D[3] = F[2];
		fF3D[4] = F[3];
		fF3D[5] = 0.0;
	    
		fF3D[6] = 0.0;
		fF3D[7] = 0.0;
		fF3D[8] = 1.0;
	}
	else fF3D = F;
	double J = fF3D.Det();
	
	/*calc elastic stretch*/
	fInverse.Inverse(fC_v[0]);
	fbe.MultQBQT(fF3D, fInverse);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
	fEigs_e = fSpectralDecompSpat.Eigenvalues(); 

	double Je = sqrt(fEigs_e.Product());
	fEigs_dev = fEigs_e;
	fEigs_dev *= pow(Je,-2.0*third);
		
	fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
	ftau_NEQ += fPot[1]->MeanStress(Je);
	ftau_NEQ /= J;
	
	/*effective visc stretch*/
	double lveff = sqrt(third*(lv0 + lv1 + lv2));
	output[0] = lveff;
	
	/*flow stress*/
	double smag = sqrt(0.5*(ftau_NEQ[0]*ftau_NEQ[0] + ftau_NEQ[1]*ftau_NEQ[1] + ftau_NEQ[2]*ftau_NEQ[2]));
	output[1] = smag;
	
	double etaS = ShearViscosity(smag, fsy0, lveff);
	output[2] = etaS;
}

/* information about subordinate parameter lists */
void BergstromBoyce::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	RGViscoelasticityT::DefineSubs(sub_list);

	sub_list.AddSub("rg_eq_potential", ParameterListT::Once);
	sub_list.AddSub("rg_neq_potential", ParameterListT::Once);

	/* choice of viscosity */
	sub_list.AddSub("viscosity_function", ParameterListT::Once);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* BergstromBoyce::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = RGSplitT2::NewSub(name);
	if (sub) 
	{
		return sub;
	}
	else if (name == "viscosity_function")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT bb("bergstrom_boyce");		
		{
			LimitT zero(0.0, LimitT::Lower);
			LimitT one(1.0, LimitT::Lower);
			LimitT positive(0.0, LimitT::LowerInclusive);

			ParameterT etaSR(ParameterT::Double, "etaS_ref");
			ParameterT sy(ParameterT::Double, "activation_stress");
			ParameterT m(ParameterT::Double, "strain_sensitivity");
			ParameterT ep(ParameterT::Double, "regularization_parameter");

			etaSR.AddLimit(zero);
			sy.AddLimit(zero);
			m.AddLimit(0, LimitT::LowerInclusive);
			m.AddLimit(1.0, LimitT::UpperInclusive);
			ep.AddLimit(0,LimitT::Lower);
			ep.AddLimit(1.0,LimitT::UpperInclusive);
		
			bb.AddParameter(etaSR);
			bb.AddParameter(sy);
			bb.AddParameter(m);
			bb.AddParameter(ep);
			
			bb.SetDescription("eta0 (lv_eff - 1 + epsilon)^m exp(s/sy)");
		}
		choice->AddSub(bb);
		
		ParameterContainerT exponential("exponential_bergstrom_boyce");		
		{
			LimitT zero(0.0, LimitT::Lower);
			LimitT one(1.0, LimitT::Lower);
			LimitT positive(0.0, LimitT::LowerInclusive);

			ParameterT etaSR(ParameterT::Double, "etaS_ref");
			ParameterT sy(ParameterT::Double, "activation_stress");
			ParameterT m(ParameterT::Double, "strain_sensitivity");

			etaSR.AddLimit(zero);
			sy.AddLimit(zero);
			m.AddLimit(0, LimitT::LowerInclusive);
		
			exponential.AddParameter(etaSR);
			exponential.AddParameter(sy);
			exponential.AddParameter(m);

			exponential.SetDescription("eta0 exp(m*(lv_eff - 1)) exp(s/sy)");
		}
		choice->AddSub(exponential);
		
		return(choice);
	}
}

void BergstromBoyce::TakeParameterList(const ParameterListT& list)
{
  const char caller[] = "BergstromBoyce::TakeParameterList";

  int num_neq =  list.NumLists("rg_neq_potential");
  int num_shear_visc = list.NumLists("viscosity_function");
  if (num_neq != num_shear_visc)
	 ExceptionT::GeneralFail("BergstromBoyce::TakeParameterList", 
		"number of  viscosity functions does not match number of nonequilibrium potentials");
  fNumProcess = list.NumLists("viscosity_function");

  /* inherited */
  RGViscoelasticityT::TakeParameterList(list);
  
  fPot.Dimension(fNumProcess+1);
		
  const ParameterListT& eq_pot = list.GetListChoice(*this, "rg_eq_potential");
  if(eq_pot.Name() == "neo-hookean")
	 fPot[0] = new NeoHookean;
  else if(eq_pot.Name() == "mooney-rivlin")
	 fPot[0] = new MooneyRivlin;
  else if(eq_pot.Name() == "veronda-westmann")
	 fPot[0] = new VWPotentialT;
  else if(eq_pot.Name() == "arruda-boyce")
	 fPot[0] = new ArrudaBoyce;
  else 
	 ExceptionT::GeneralFail(caller, "no such potential");

  if (!fPot[0]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", eq_pot.Name().Pointer());			
	 fPot[0]->TakeParameterList(eq_pot);
	
  const ParameterListT& pot_neq = list.GetListChoice(*this, "rg_neq_potential",0);
  if(pot_neq.Name() == "mooney-rivlin")
	 fPot[1] = new MooneyRivlin;
  else if(pot_neq.Name() == "neo-hookean")
	 fPot[1] = new NeoHookean;
  else if(pot_neq.Name() == "veronda-westmann")
	 fPot[1] = new VWPotentialT;
  else if(pot_neq.Name() == "arruda-boyce")
	 fPot[1] = new ArrudaBoyce;
  else 
	 ExceptionT::GeneralFail(caller, "no such potential");
  if (!fPot[1]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", pot_neq.Name().Pointer());			
	 fPot[1]->TakeParameterList(pot_neq);
  
	const ParameterListT& visc = list.GetListChoice(*this, "viscosity_function", 0);
	if (visc.Name() == "bergstrom_boyce")
	{		
		fetaS0 = visc.GetParameter("etaS_ref");
		fsy0 = visc.GetParameter("activation_stress");
		fm = visc.GetParameter("strain_sensitivity");
		fepsilon = visc.GetParameter("regularization_parameter");
		
		fViscType = kBergstromBoyce;
	}
	else if (visc.Name() == "exponential_bergstrom_boyce")
	{
		fetaS0 = visc.GetParameter("etaS_ref");
		fsy0 = visc.GetParameter("activation_stress");
		fm = visc.GetParameter("strain_sensitivity");
		
		fViscType = kExpBergstromBoyce;
	}
	  
	/*set dimension of workspaces*/
	Initialize();

	fRes.Dimension(3);
	fDelta.Dimension(3);
	fGAB.Dimension(3,3);
	fDAB.Dimension(3,3);
	fDABbar.Dimension(3);
	fMat.Dimension(3);
}

/***********************************************************************
 * Private
 ***********************************************************************/
/* set inverse of thermal transformation - return true if active */
 void BergstromBoyce::Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
	const double& dtau_m, dMatrixT& Calg, const int type)
 {
		const dMatrixT& F = F_total();
		double iJ = 1.0/F.Det();

		/*temperature and temperature step*/
		/*time step*/
		double dt = fFSMatSupport->TimeStep();

	    const double& le0 = fEigs_e[0];
	    const double& le1 = fEigs_e[1];
	    const double& le2 = fEigs_e[2];
	    const double& l0 = fEigs[0];
	    const double& l1 = fEigs[1];
	    const double& l2 = fEigs[2];
		double lv0 = l0/le0;
		double lv1 = l1/le1;
		double lv2 = l2/le2;

	    double s0 = iJ*tau_dev[0];
	    double s1 = iJ*tau_dev[1];
	    double s2 = iJ*tau_dev[2];
	    
		double c0 = iJ*dtau_dev(0,0);
		double c1 = iJ*dtau_dev(1,1);
		double c2 = iJ*dtau_dev(2,2);

		double c12 = iJ*dtau_dev(1,2);
		double c02 = iJ*dtau_dev(0,2);
		double c01 = iJ*dtau_dev(0,1);
						
	    /*caculate smag*/
	    double smag = sqrt(0.5*(s0*s0 + s1*s1 + s2*s2));

		/*calculate effective viscous stretch*/
		double lveff = sqrt(third*(lv0+lv1+lv2));
		
		double etaS = ShearViscosity(smag, fsy0, lveff);
		double ietaS = 1.0/etaS;
						
		double gamdot = 0.5*smag*ietaS;

		/*calculate stiffness matrix*/
		/*derivative of retardation time wrt to Tf*/
		double detaS_dsmag = -etaS /fsy0;
		double detaS_dlv = DVisc_Dlamv(smag,fsy0,lveff);
//		cout << "\ndetaS_dsmag: "<<detaS_dsmag;
		
		/*initialize*/
		fiKAB = 0.0;
		
		/*K_epA_epB*/
		double coef0 = 0.5*(s0*c0 + s1*c01 + s2*c02);
		double coef1 = 0.5*(s0*c01 + s1*c1 + s2*c12);
		double coef2 = 0.5*(s0*c02 + s1*c12 + s2*c2);
		if (smag > kSmall)
		{
			coef0 /= smag;
			coef1 /= smag;
			coef2 /= smag;
		}
		fiKAB(0,0) = 1.0 + 0.5*ietaS*dt*c0 - 0.5*dt*ietaS*s0* ietaS*(detaS_dsmag*coef0 + detaS_dlv*(-lv0/lveff));
		fiKAB(1,1) = 1.0 + 0.5*ietaS*dt*c1 - 0.5*dt*ietaS*s1* ietaS*(detaS_dsmag*coef1 + detaS_dlv*(-lv1/lveff));
		fiKAB(2,2) = 1.0 + 0.5*ietaS*dt*c2 - 0.5*dt*ietaS*s2* ietaS*(detaS_dsmag*coef2 + detaS_dlv*(-lv2/lveff));
		
		fiKAB(1,2) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s1* ietaS*(detaS_dsmag*coef2 + detaS_dlv*(-lv2/lveff));
		fiKAB(0,2) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s0* ietaS*(detaS_dsmag*coef2 + detaS_dlv*(-lv2/lveff));
		fiKAB(0,1) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s0* ietaS*(detaS_dsmag*coef1 + detaS_dlv*(-lv1/lveff));

		fiKAB(2,1) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s2* ietaS*(detaS_dsmag*coef1 + detaS_dlv*(-lv1/lveff));
		fiKAB(2,0) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s2* ietaS*(detaS_dsmag*coef0 + detaS_dlv*(-lv0/lveff));
		fiKAB(1,0) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s1* ietaS*(detaS_dsmag*coef0 + detaS_dlv*(-lv0/lveff));
//		cout << "\nfiKAB: "<<fiKAB;

		/*inverts KAB*/
		fiKAB.Inverse();

		/*initialize*/
		fGAB = 0.0;
		
		/*G_epeA_epB*/
		fGAB(0,0) = 1.0 + 0.5*dt*ietaS*s0 - 0.5*dt*ietaS *ietaS*s0*(s0 + detaS_dlv*(lv0/lveff));
		fGAB(0,1) = 0.5*dt*ietaS*s0 - 0.5*dt*ietaS *ietaS*s0*(s0 + detaS_dlv*(lv1/lveff));
		fGAB(0,2) = 0.5*dt*ietaS*s0 - 0.5*dt*ietaS *ietaS*s0*(s0 + detaS_dlv*(lv2/lveff));
		
		fGAB(1,0) = 0.5*dt*ietaS*s1 - 0.5*dt*ietaS *ietaS*s1*(s1 + detaS_dlv*(lv0/lveff));
		fGAB(1,1) = 1.0 + 0.5*dt*ietaS*s1 - 0.5*dt*ietaS *ietaS*s1*(s1 + detaS_dlv*(lv1/lveff));
		fGAB(1,2) = 0.5*dt*ietaS*s1 - 0.5*dt*ietaS *ietaS*s1*(s1 + detaS_dlv*(lv2/lveff));

		fGAB(2,0) = 0.5*dt*ietaS*s2 - 0.5*dt*ietaS *ietaS*s2*(s2 + detaS_dlv*(lv0/lveff));
		fGAB(2,1) = 0.5*dt*ietaS*s2 - 0.5*dt*ietaS *ietaS*s2*(s2 + detaS_dlv*(lv1/lveff));
		fGAB(2,2) = 1.0 + 0.5*dt*ietaS*s2 - 0.5*dt*ietaS *ietaS*s2*(s2 + detaS_dlv*(lv2/lveff));
		
		
		/*Calg = dtau/depe*fiKA*fG	*/	
		/*calculating delta_internval_vars = K^-1.G. delta_epsilon*/
		fDAB.MultAB(fiKAB,fGAB);
		/*copy subset*/
	
		for (int i = 0; i< fDABbar.Rows(); i++)
			for (int j = 0; j< fDABbar.Cols(); j++)
				fDABbar(i,j) = fDAB(i,j);
				
		dtau_dev.ToMatrix(fMat);
		fMat(0,0) += dtau_m;
		fMat(1,1) += dtau_m;
		fMat(2,2) += dtau_m;
		
		Calg.MultAB(fMat, fDABbar);

/*		cout << "\nfDAB: "<<fDAB;
		cout << "\nfDAbar: "<<fDABbar;
		cout <<"\ndtau: "<<dtau_dev;
		cout << "\nCalg0: "<<Calg;
*/
		Calg(0,0) -= 2.0* tau_dev[0];
		Calg(1,1) -= 2.0* tau_dev[1];
		Calg(2,2) -= 2.0* tau_dev[2];
}

void BergstromBoyce::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
			     dArrayT& eigenstress, dSymMatrixT& eigenmodulus,  const int type) 
{		
	const double ctol = 1.00e-12;
		
	/*set references to principle stretches*/
     
	double& le0 = eigenstretch_e[0];
	double& le1 = eigenstretch_e[1];
	double& le2 = eigenstretch_e[2];
	
	double l0 = eigenstretch[0];
	double l1 = eigenstretch[1];
	double l2 = eigenstretch[2];
	
	double tol;

	/*initialize principle elastic and trial elastic log strains */
	const double ep_tr0 = 0.5*log(le0);
	const double ep_tr1 = 0.5*log(le1);
	const double ep_tr2 = 0.5*log(le2);
	
	double ep_e0 = ep_tr0;		
	double ep_e1 = ep_tr1;	
	double ep_e2 = ep_tr2;

//	cout << "\neigs: "<<eigenstretch;
//	cout << "\neigs_e: "<<eigenstretch_e;
	/*jacobian*/
	const dMatrixT& F = F_total();
	double iJ = 1.0/F.Det();

	/*time step*/
	double dt = fFSMatSupport->TimeStep();
	
	int maxiteration = 100;
	/*initializes principle viscous stretch*/
	int iteration = 0;
	do 
	{
		iteration ++;

		double lv0 = l0/le0;
		double lv1 = l1/le1;
		double lv2 = l2/le2;
	
	    double Je=sqrt(le0*le1*le2);
	    fEigs_dev = eigenstretch_e;
	    fEigs_dev *= pow(Je,-2.0*third);

//	cout << "\neigenstretch_edev: "<<fEigs_dev;

	    /*calculate stresses and moduli*/
	    fPot[1]->DevStress(fEigs_dev, eigenstress);
	    
	    double s0 = iJ*eigenstress[0];
	    double s1 = iJ*eigenstress[1];
	    double s2 = iJ*eigenstress[2];
	    
	    fPot[1]->DevMod(fEigs_dev,eigenmodulus);
		/*deviatoric values*/
		double c0 = iJ*eigenmodulus(0,0);
		double c1 = iJ*eigenmodulus(1,1);
		double c2 = iJ*eigenmodulus(2,2);

		double c12 = iJ*eigenmodulus(1,2);
		double c02 = iJ*eigenmodulus(0,2);
		double c01 = iJ*eigenmodulus(0,1);
		
//		cout << "\neigenmod: "<<c0<<"\t"<<c1<<"\t"<<c2;
	    		
	    /*caculate smag*/
	    double smag = sqrt(0.5*(s0*s0 + s1*s1 + s2*s2));
//		cout << "\nsmag: "<<smag;

		/*calculate effective shear strain*/
		double lveff = sqrt(third*(lv0 + lv1 + lv2));
//		cout << "\nlveff: "<<lveff;
		
		/*calculate mobilities*/
		double etaS = ShearViscosity(smag, fsy0, lveff);
		double ietaS = 1.0/etaS;
//		cout << "\neatS: "<<etaS;
					
	    /*calculate the residual*/
	    fRes[0] = ep_e0 + 0.5*dt*ietaS*s0 - ep_tr0;
	    fRes[1] = ep_e1 + 0.5*dt*ietaS*s1 - ep_tr1;
	    fRes[2] = ep_e2 + 0.5*dt*ietaS*s2 - ep_tr2;
//		cout << "\nfRes: "<<fRes<<endl;

		/*calculate stiffness matrix*/
		/*derivative of retardation time wrt to Tf*/
		double detaS_dsmag = -etaS /fsy0;
		double detaS_dlv = DVisc_Dlamv(smag,fsy0,lveff);
//		cout << "\ndetaS_dsmag: "<<detaS_dsmag;
		
		/*initialize*/
		fiKAB = 0.0;
		
		/*K_epA_epB*/
		double coef0 = 0.5*(s0*c0 + s1*c01 + s2*c02);
		double coef1 = 0.5*(s0*c01 + s1*c1 + s2*c12);
		double coef2 = 0.5*(s0*c02 + s1*c12 + s2*c2);
		if (smag > kSmall)
		{
			coef0 /= smag;
			coef1 /= smag;
			coef2 /= smag;
		}
		fiKAB(0,0) = 1.0 + 0.5*ietaS*dt*c0 - 0.5*dt*ietaS*s0* ietaS*(detaS_dsmag*coef0 + detaS_dlv*(-lv0/lveff));
		fiKAB(1,1) = 1.0 + 0.5*ietaS*dt*c1 - 0.5*dt*ietaS*s1* ietaS*(detaS_dsmag*coef1 + detaS_dlv*(-lv1/lveff));
		fiKAB(2,2) = 1.0 + 0.5*ietaS*dt*c2 - 0.5*dt*ietaS*s2* ietaS*(detaS_dsmag*coef2 + detaS_dlv*(-lv2/lveff));
		
		fiKAB(1,2) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s1* ietaS*(detaS_dsmag*coef2 + detaS_dlv*(-lv2/lveff));
		fiKAB(0,2) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s0* ietaS*(detaS_dsmag*coef2 + detaS_dlv*(-lv2/lveff));
		fiKAB(0,1) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s0* ietaS*(detaS_dsmag*coef1 + detaS_dlv*(-lv1/lveff));

		fiKAB(2,1) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s2* ietaS*(detaS_dsmag*coef1 + detaS_dlv*(-lv1/lveff));
		fiKAB(2,0) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s2* ietaS*(detaS_dsmag*coef0 + detaS_dlv*(-lv0/lveff));
		fiKAB(1,0) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s1* ietaS*(detaS_dsmag*coef0 + detaS_dlv*(-lv0/lveff));
//		cout << "\nfiKAB: "<<fiKAB;

		/*inverts KAB*/
		fiKAB.Inverse();
	    		
	    /*solve for the principal strain increments*/
		fiKAB.Multx(fRes, fDelta, -1.0);
		
//		cout << "\ndep: "<<fDelta;				
	    ep_e0 += fDelta[0];
	    ep_e1 += fDelta[1];
	    ep_e2 += fDelta[2];
	    		
	    le0 = exp(2.0*ep_e0);
	    le1 = exp(2.0*ep_e1);
	    le2 = exp(2.0*ep_e2);
	    
	    /*Check that the L2 norm of the residual is less than tolerance*/
	    tol = sqrt(dArrayT::Dot(fRes, fRes));
//		cout << "\ntol: "<<tol;
//		cout << "\neps_e: "<<eigenstretch_e;

	}while (tol>ctol && iteration < maxiteration); 
	if (iteration >= maxiteration) 
		ExceptionT::GeneralFail("BergstromBoyce::ComputeEigs_e", 
			"number of iteration exceeds maximum");
//	cout << "\neigs: "<<eigenstretch;
//	cout << "\neigse: "<<eigenstretch_e;
//	cout << "\nlveff: "<<sqrt(third*(eigenstretch[0]/eigenstretch_e[0]+eigenstretch[1]/eigenstretch_e[1]+eigenstretch[2]/eigenstretch_e[2]));
}
