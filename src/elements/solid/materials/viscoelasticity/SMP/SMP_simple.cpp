/* $Id: SMP_simple.cpp,v 1.3 2006-11-14 22:58:06 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "SMP_simple.h"

#include "ifstreamT.h"
#include "ExceptionT.h"
#include <math.h>
#include <iostream.h>
#include <stdlib.h>
#include "ParameterContainerT.h"

using namespace Tahoe;

const double third = 1.0/3.0; 
const int kNumOutputVar = 4; 
static const char* Labels[kNumOutputVar] = {"AbsTemp","thermal_dialation", "etaB", "etaS"}; 

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
/* constructors */
SMP_simple::SMP_simple(void):
  ParameterInterfaceT("SMP_simple")
{

}

int SMP_simple::NumOutputVariables() const {return kNumOutputVar;} 

void SMP_simple::OutputLabels(ArrayT<StringT>& labels) const 
{ 
     /*allocates space for labels*/
     labels.Dimension(kNumOutputVar); 
  
     /*copy labels*/
     for (int i = 0; i< kNumOutputVar; i++) 
       labels[i] = Labels[i]; 
} 

/*viscosity*/
void SMP_simple::Viscosity(double& ietaS, double& ietaB, const int type)
{
	/*get delta temperature*/
	double T = Compute_Temperature()+fRefTemperature;
	
 	switch(fViscType)
	{
		case (SMP_simple::kSimple):
		{
			double etaS_high = fVisc(type,0);
			double etaS_low = fVisc(type, 1);
			double etaB_high = fVisc(type, 2);
			double etaB_low = fVisc(type, 3);
			double a = fVisc(type, 4);
	
			double etaS = 0.5*(1.0 - tanh(a*(T/fTg-1.0)))*etaS_high + etaS_low;
			double etaB = 0.5*(1.0 - tanh(a*(T/fTg-1.0)))*etaB_high + etaB_low;
		
			ietaS = 1.0/etaS;
			ietaB = 1.0/etaB;
			
			break;
		}
		default:
			ExceptionT::GeneralFail("SMP_simple::Viscosity", "Unknown  potential", fViscType);
	}
}

double SMP_simple::Energy(const dArrayT& lambda_bar, const double J, const int type)
{
	double phi;
	switch(fPotType)
	{
		case (SMP_simple::kMooneyRivlin):
		{
			double c1, c2, gamma, beta;
			if (type > -1)
			{
				c1 = fPot(type+1,0);
				c2 = fPot(type+1,1);
				gamma = fPot(type+1,2);
				beta = fPot(type+1,3);
			}
			else
			{
				double T = Compute_Temperature();
				double r = (1.0+T/fRefTemperature);
				c1 = fPot(0,0)*r;
				c2 = fPot(0,1)*r;
				gamma = fPot(0,2)*r;
			beta = fPot(0,3);
			}

			double I1 = lambda_bar[0]+lambda_bar[1]+lambda_bar[2];
			double I2 = lambda_bar[0]*lambda_bar[1]+lambda_bar[1]*lambda_bar[2]+lambda_bar[0]*lambda_bar[2];

			phi = 0.5*c1*(I1-3.0)-0.5*c2*(I2-3.0)+0.25*gamma*(J*J + 1.0/beta*pow(J,-2.0*beta));
			
			break;
		}
		default:	
			ExceptionT::GeneralFail("SMP_simple::Energy", "Unknown  potential", fPotType);
	}

	return(phi);
}
void SMP_simple::DevStress(const dArrayT& lambda_bar,dArrayT& tau, const int type)
{
	int nsd = tau.Length();

	const double& l0 = lambda_bar[0];
	const double& l1 = lambda_bar[1];
	const double& l2 = lambda_bar[2];

	switch(fPotType)
	{
		case (SMP_simple::kMooneyRivlin):
		{
			/*Mooney Rivlin Potential*/
			/*Psi = 0.5*c1 (I1_bar -3) - 0.5*c2 (I2_bar -3)*/

			double c1, c2;
			if (type > -1)
			{
				/*get delta temperature*/
				c1 = fPot(type+1,0);
				c2 = fPot(type+1,1);
			}
			else
			{
				double T = Compute_Temperature();
				double r = (1.0+T/fRefTemperature);
				c1 = fPot(0,0)*r;
				c2 = fPot(0,1)*r;
			}

			/*c1(I1bar - 3.0) eigen kirchhoff stress*/
			tau[0] = c1*third*(2.0*l0 - l1 - l2);
			tau[1] = c1*third*(2.0*l1 - l2 - l0);

			/*c2(I2bar - 3.0) eigen kirchhoff stress*/
			tau[0] += c2*third*(2.0/l0 - 1.0/l1 - 1.0/l2);
			tau[1] += c2*third*(2.0/l1 - 1.0/l2 - 1.0/l0);

			if (nsd == 3)
			{
				tau[2] = c1*third*(2.0*l2 - l0 - l1);
				tau[2] += c2*third*(2.0/l2 - 1.0/l0 - 1.0/l1);
			}
			break;
		}
		default:
			ExceptionT::GeneralFail("SMP_simple::DevStress", "Unknown  potential", fPotType);
	}
}

double SMP_simple::MeanStress(const double J, const int type) 
{
	/* U(J) = 0.25*gamma*(J^2 + 1/beta J^(-2beta))*/
	double gamma, beta;
	if (type > -1)
	{
		gamma = fPot(type+1,2);
		beta = fPot(type+1,3);
	}
	else
	{
		double T = Compute_Temperature();
		double r = (1.0+T/fRefTemperature);
		gamma = fPot(0,2)*r;
		beta = fPot(0,3);
	}

	return(0.5*gamma*(J*J - pow(J, -2.0*beta)));
}

void SMP_simple::DevMod(const dArrayT& lambda_bar, dSymMatrixT& eigenmodulus, const int type)
{
	/*Mooney Rivlin Potential*/
	/*Psi = 0.5*c1 (I1_bar -3) - 0.5*c2 (I2_bar -3)*/

	int nsd = eigenmodulus.Rows();
	double ninth = third*third;
	
 	switch(fPotType)
	{
		case (SMP_simple::kMooneyRivlin):
		{
			double c1, c2;
			if (type > -1)
			{
				c1 = fPot(type+1,0);
				c2 = fPot(type+1,1);
			}
			else
			{
				double T = Compute_Temperature();
				double r = (1.0+T/fRefTemperature);
				c1 = fPot(0,0)*r;
				c2 = fPot(0,1)*r;
			}
	
			const double& l0 = lambda_bar[0];
			const double& l1 = lambda_bar[1];
			const double& l2 = lambda_bar[2];
	
			eigenmodulus[0] = 2.0*c1*ninth*(4.0*l0 + l1 + l2);
			eigenmodulus[1] = 2.0*c1*ninth*(4.0*l1 + l2 + l0);
	
			eigenmodulus[0] += -2.0*c2*ninth*(4.0/l0 + 1.0/l1 + 1.0/l2);
			eigenmodulus[1] += -2.0*c2*ninth*(4.0/l1 + 1.0/l2 + 1.0/l0);

			if (nsd == 2)
			{
				eigenmodulus[2] = 2.0*c1*ninth*(-2.0*l0 - 2.0*l1 + l2);
				eigenmodulus[2] += -2.0*c2*ninth*(-2.0/l0 - 2.0/l1 + 1.0/l2);
			}
			else 
			{
				eigenmodulus[2] = 2.0*c1*ninth*(4.0*l2 + l0 + l1);
				eigenmodulus[3] = 2.0*c1*ninth*(-2.0*l1 - 2.0*l2 + l0);
				eigenmodulus[4] = 2.0*c1*ninth*(-2.0*l0 - 2.0*l2 + l1);
				eigenmodulus[5] = 2.0*c1*ninth*(-2.0*l0 - 2.0*l1 + l2);

				eigenmodulus[2] += -2.0*c2*ninth*(4.0/l2 + 1.0/l0 + 1.0/l1);
				eigenmodulus[3] += -2.0*c2*ninth*(-2.0/l1 - 2.0/l2 + 1.0/l0);
				eigenmodulus[4] += -2.0*c2*ninth*(-2.0/l0 - 2.0/l2 + 1.0/l1);
				eigenmodulus[5] += -2.0*c2*ninth*(-2.0/l0 - 2.0/l1 + 1.0/l2);
			}
			break;
		}
		default:
			ExceptionT::GeneralFail("SMP_simple::DevMod", "Unknown  potential", fPotType);
	}
}

double SMP_simple::MeanMod(const double J, const int type) 
{
	/* U(J) = 0.25*gamma*(J^2 + 1/beta J^(-2beta))*/
	double gamma, beta;
	if (type > -1)
	{
		gamma = fPot(type+1,2);
		beta = fPot(type+1,3);
	}
	else
	{
		double T = Compute_Temperature();
		double r = (1.0+T/fRefTemperature);
		gamma = fPot(0,2)*r;
		beta = fPot(0,3);
	}
	return(gamma*(J*J+beta*pow(J, -2.0*beta)));
}

void SMP_simple::ComputeOutput(dArrayT& output)
{
	output[0] = fRefTemperature + Compute_Temperature();
	const dMatrixT& F_thermal = F_thermal_inverse();
	double dialation = F_thermal.Det();
	output[1] = 1.0/dialation;
	if (fViscType != SMP_simple::kNone)
	{
	  Viscosity(fietaS, fietaB, 0);
	  output[2] = 1.0/fietaB;
	  output[3] = 1.0/fietaS;
	}
	else
	  {
	    output[2] = 0.0;
	    output[3] = 0.0;
	  }

}

/* describe the parameters needed by the interface */
void SMP_simple::DefineParameters(ParameterListT& list) const
{
  /* inherited */
  RGViscoelasticityT::DefineParameters(list);

  /* common limit */
  LimitT positive(0.0, LimitT::Lower);

  ParameterT reftemp(ParameterT::Double, "ref_temperature");
  list.AddParameter(reftemp);

  ParameterT Tg(ParameterT::Double, "glass_transition_temperature");
  list.AddParameter(Tg);
}

/* information about subordinate parameter lists */
void SMP_simple::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("eq_potential_params", ParameterListT::Once);
	sub_list.AddSub("neq_potential_params", ParameterListT::Any);

	/* choice of viscosity */
	sub_list.AddSub("viscosity", ParameterListT::Any);
}


/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SMP_simple::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = RGViscoelasticityT::NewSub(name);
	if (sub) 
	{
		return sub;
	}
	else if (name == "eq_potential_params" || name == "neq_potential_params")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT matrix("Mooney-Rivlin");
		matrix.SetDescription("Psi=0.5*c1(I1_bar-3) - 0.5*c2(I2_bar -3) + 0.25*gamma*(J^2 + 1/beta * J^(-2*beta) )");
		/* bound */
		LimitT lower(0.0, LimitT::Lower);
		LimitT positive(0.0, LimitT::LowerInclusive);
		
		ParameterT c1(ParameterT::Double, "c1");
		ParameterT c2(ParameterT::Double, "c2");
		ParameterT gamma(ParameterT::Double, "gamma");
		ParameterT beta(ParameterT::Double, "beta");

		c1.AddLimit(lower);
		c2.AddLimit(positive);
		gamma.AddLimit(lower);
		beta.AddLimit(positive);
		
		matrix.AddParameter(c1);
		matrix.AddParameter(c2);
		matrix.AddParameter(gamma);
		matrix.AddParameter(beta);

		choice->AddSub(matrix);
		return choice;
	}
	else if (name == "viscosity")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT tanh("temp_dependence_only");		
		/* set the description */
		LimitT lower(0.0, LimitT::Lower);
		LimitT positive(0.0, LimitT::LowerInclusive);
		tanh.SetDescription("eta(T) = (1-tanh(a*(T/Tref - 1)))*eta_high + eta_low;");	

		ParameterT etaS_high(ParameterT::Double, "etaS_high");
		ParameterT etaS_low(ParameterT::Double, "etaS_low");
		ParameterT etaB_high(ParameterT::Double, "etaB_high");
		ParameterT etaB_low(ParameterT::Double, "etaB_low");
		ParameterT a(ParameterT::Double, "a");
		
		tanh.AddParameter(etaS_high);
		tanh.AddParameter(etaS_low);
		tanh.AddParameter(etaB_high);
		tanh.AddParameter(etaB_low);
		tanh.AddParameter(a);
		
		etaS_high.AddLimit(lower);
		etaS_low.AddLimit(lower);
		etaB_high.AddLimit(lower);
		etaB_low.AddLimit(lower);
		a.AddLimit(positive);
		
		choice->AddSub(tanh);
		return(choice);
	}
}

void SMP_simple::TakeParameterList(const ParameterListT& list)
{
  /* inherited */
  RGViscoelasticityT::TakeParameterList(list);

  /*allows one neq process: */
  int num_pot = list.NumLists("neq_potential_params");
  int num_visc = list.NumLists("viscosity");
  
  if (num_pot != num_visc)
		ExceptionT::GeneralFail("SMP_simple::TakeParameterList", 
			"number of viscosity functions does not match number of nonequilibrium potentials");

  fNumProcess = num_pot;
  fPot.Dimension(fNumProcess+1, 4);
  fVisc.Dimension(fNumProcess+1, 5);

  fRefTemperature = list.GetParameter("ref_temperature");
  fTg = list.GetParameter("glass_transition_temperature");
  
  const ParameterListT& eq = list.GetListChoice(*this, "eq_potential_params");
  if (eq.Name() == "Mooney-Rivlin")
  {
		fPotType = SMP_simple::kMooneyRivlin;
		fPot(0,0) = eq.GetParameter("c1");
		fPot(0,1) = eq.GetParameter("c2");
		fPot(0,2) = eq.GetParameter("gamma");
		fPot(0,3) = eq.GetParameter("beta");
  }
  for (int i = 0; i < fNumProcess; i++)
  {
	  const ParameterListT& neq = list.GetListChoice(*this, "neq_potential_params",i);
	  if (neq.Name() == "Mooney-Rivlin")
	  {
		  fPot(i+1,0) = neq.GetParameter("c1");
		  fPot(i+1,1) = neq.GetParameter("c2");
		  fPot(i+1,0) = neq.GetParameter("c1");
		  fPot(i+1,1) = neq.GetParameter("c2");
		  fPot(i+1,2) = neq.GetParameter("gamma");
		  fPot(i+1,3) = neq.GetParameter("beta");
	 }
	 const ParameterListT& visc  = list.GetListChoice(*this, "viscosity", i);
	 if (visc.Name() == "temp_dependence_only")
	 {
		fViscType = SMP_simple::kSimple;
		fVisc(i,0) = visc.GetParameter("etaS_high");
		fVisc(i,1) = visc.GetParameter("etaS_low");
		fVisc(i,2) = visc.GetParameter("etaB_high");
		fVisc(i,3) = visc.GetParameter("etaB_low");
		fVisc(i,4) = visc.GetParameter("a");
	 }
	 else
	   fViscType = SMP_simple::kNone;
  }


  /*dimension rest of work space*/
  Initialize();
}

void SMP_simple::InitStep(void)
{
	/*inherited*/
	RGSplitT::InitStep();
}

/***********************************************************************
 * Private
 ***********************************************************************/
/* set inverse of thermal transformation - return true if active */
 void SMP_simple::Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
	const double& dtau_m, dMatrixT& Calg, const int type)
 {
		double c0 = dtau_dev(0,0);
		double c1 = dtau_dev(1,1);
		double c2 = dtau_dev(2,2);

		double c12 = dtau_dev(1,2);
		double c02 = dtau_dev(0,2);
		double c01 = dtau_dev(0,1);
	    
		double cm = dtau_m;
		
		/*calculates  KAB = 1+dt*D(dWdE_Idev/nD+isostress/nV)/Dep_e*/
		double dt = fFSMatSupport->TimeStep();
		Viscosity(fietaS, fietaB, type);

		fiKAB(0,0) = 1+0.5*fietaS*dt*c0+third*fietaB*dt*cm;
		fiKAB(1,1) = 1+0.5*fietaS*dt*c1+third*fietaB*dt*cm;
		fiKAB(2,2) = 1+0.5*fietaS*dt*c2+third*fietaB*dt*cm;

		fiKAB(1,2) = 0.5*fietaS*dt*c12+third*fietaB*dt*cm;
		fiKAB(0,2) = 0.5*fietaS*dt*c02+third*fietaB*dt*cm;
		fiKAB(0,1) = 0.5*fietaS*dt*c01+third*fietaB*dt*cm;
       
		fiKAB(2,1) = fiKAB(1,2);
		fiKAB(2,0) = fiKAB(0,2);
		fiKAB(1,0) = fiKAB(0,1);
	
		/*inverts KAB*/
		fiKAB.Inverse();

		dSymMatrixT& DAB = fDtauDe_NEQ;
		DAB += cm; 
	
		Calg(0,0) = (c0+cm)*fiKAB(0,0) + (c01+cm)*fiKAB(1,0) + (c02+cm)*fiKAB(2,0) - 2.0*(tau_dev[0]+tau_m);
		Calg(1,0) = (c01+cm)*fiKAB(0,0) + (c1+cm)*fiKAB(1,0) + (c12+cm)*fiKAB(2,0);
		Calg(2,0) = (c02+cm)*fiKAB(0,0) + (c12+cm)*fiKAB(1,0) + (c2+cm)*fiKAB(2,0);
		Calg(0,1) = (c0+cm)*fiKAB(0,1) + (c01+cm)*fiKAB(1,1) + (c02+cm)*fiKAB(2,1);
		Calg(1,1) = (c01+cm)*fiKAB(0,1) + (c1+cm)*fiKAB(1,1) + (c12+cm)*fiKAB(2,1) - 2.0*(tau_dev[1]+tau_m);
		Calg(2,1) = (c02+cm)*fiKAB(0,1) + (c12+cm)*fiKAB(1,1) + (c2+cm)*fiKAB(2,1);
		Calg(0,2) = (c0+cm)*fiKAB(0,2) + (c01+cm)*fiKAB(1,2) + (c02+cm)*fiKAB(2,2);
		Calg(1,2) = (c01+cm)*fiKAB(0,2) + (c1+cm)*fiKAB(1,2) + (c12+cm)*fiKAB(2,2);
		Calg(2,2) = (c02+cm)*fiKAB(0,2) + (c12+cm)*fiKAB(1,2) + (c2+cm)*fiKAB(2,2) - 2.0*(tau_dev[2]+tau_m);
}

void SMP_simple::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
			     dArrayT& eigenstress, dSymMatrixT& eigenmodulus, const int type) 
{		
	const double ctol = 1.00e-14;
		
	/*set references to principle stretches*/
     
	double& le0 = eigenstretch_e[0];
	double& le1 = eigenstretch_e[1];
	double& le2 = eigenstretch_e[2];
  
	double tol;

	/*initialize principle elastic and trial elastic log strains */
	double ep_tr0 = 0.5*log(le0);
	double ep_tr1 = 0.5*log(le1);
	double ep_tr2 = 0.5*log(le2);

	double ep_e0 = ep_tr0;		
	double ep_e1 = ep_tr1;	
	double ep_e2 = ep_tr2;
	

	/*initializes principle viscous stretch*/
	do 
	{
	    double Je=sqrt(le0*le1*le2);
	    fEigs_dev = eigenstretch_e;
	    fEigs_dev *= pow(Je,-2.0*third);
		
	    /*calculate stresses and moduli*/
	    DevStress(fEigs_dev, eigenstress, type);
	    
	    double& s0 = eigenstress[0];
	    double& s1 = eigenstress[1];
	    double& s2 = eigenstress[2];
	    
	    DevMod(fEigs_dev,eigenmodulus, type);
		/*deviatoric values*/
		double& c0 = eigenmodulus(0,0);
		double& c1 = eigenmodulus(1,1);
		double& c2 = eigenmodulus(2,2);

		double& c12 = eigenmodulus(1,2);
		double& c02 = eigenmodulus(0,2);
		double& c01 = eigenmodulus(0,1);
	    
	    /*caculate means*/
	    double sm = MeanStress(Je, type);
	    double cm = MeanMod(Je, type);
	    
		Viscosity(fietaS, fietaB, type);
		double dt = fFSMatSupport->TimeStep();
		fiKAB(0,0) = 1+0.5*fietaS*dt*c0+third*fietaB*dt*cm;
		fiKAB(1,1) = 1+0.5*fietaS*dt*c1+third*fietaB*dt*cm;
		fiKAB(2,2) = 1+0.5*fietaS*dt*c2+third*fietaB*dt*cm;

		fiKAB(1,2) = 0.5*fietaS*dt*c12+third*fietaB*dt*cm;
		fiKAB(0,2) = 0.5*fietaS*dt*c02+third*fietaB*dt*cm;
		fiKAB(0,1) = 0.5*fietaS*dt*c01+third*fietaB*dt*cm;
       
		fiKAB(2,1) = fiKAB(1,2);
		fiKAB(2,0) = fiKAB(0,2);
		fiKAB(1,0) = fiKAB(0,1);
	
		/*inverts KAB*/
		fiKAB.Inverse();
	    
	    /*calculate the residual*/
	    double res0 = ep_e0 + dt*(0.5*fietaS*s0 +
			  third*fietaB*sm) - ep_tr0;
	    double res1 = ep_e1 + dt*(0.5*fietaS*s1 +
			  third*fietaB*sm) - ep_tr1;
	    double res2 = ep_e2 + dt*(0.5*fietaS*s2 +
			  third*fietaB*sm) - ep_tr2;
		
	    /*solve for the principal strain increments*/
	    double dep_e0=-fiKAB(0,0)*res0-fiKAB(0,1)*res1-fiKAB(0,2)*res2;
	    double dep_e1=-fiKAB(1,0)*res0-fiKAB(1,1)*res1-fiKAB(1,2)*res2;
	    double dep_e2=-fiKAB(2,0)*res0-fiKAB(2,1)*res1-fiKAB(2,2)*res2;
	    
	    /*updates principal elastic stretches*/ 
	    ep_e0 += dep_e0;
	    ep_e1 += dep_e1;
	    ep_e2 += dep_e2;
	    
	    le0 = exp(2.0*ep_e0);
	    le1 = exp(2.0*ep_e1);
	    le2 = exp(2.0*ep_e2);
	    
	    /*Check that the L2 norm of the residual is less than tolerance*/
	    tol = sqrt(res0*res0 + res1*res1+res2*res2);
	}while (tol>ctol); 
}
