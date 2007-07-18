/* $Id: SMP_simple.cpp,v 1.8 2007-07-18 03:08:26 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "SMP_simple.h"

#include "PotentialT.h"
#include "NeoHookean.h"
#include "ArrudaBoyce.h"

#include "ifstreamT.h"
#include "ExceptionT.h"
#include <math.h>
#include <iostream.h>
#include <stdlib.h>
#include "ParameterContainerT.h"

using namespace Tahoe;

const double loge = log10(exp(1.0));
const double third = 1.0/3.0; 
const int kNumOutputVar = 4; 
static const char* Labels[kNumOutputVar] = {"thermal_dialation", "deltaneq", "lv_eff", "sy"}; 

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
/* constructors */
SMP_simple::SMP_simple(void):
  ParameterInterfaceT("SMP_simple")
{
	fNumProcess = 1;
//	cout << setprecision(16);
}

double SMP_simple::FictiveTemperature(const double deltaneq)
{
	double Temp = Compute_Temperature();
	return (Temp + ( deltaneq - (falphar-falphag)*(Temp-fT0) )/(falphar-falphag));
}


double SMP_simple::RetardationTime(const double Temperature, const double deltaneq)
{
	/*calculate fictive temperature*/
	double Tf = FictiveTemperature(deltaneq);
	/*Hodge's Model*/
//	double Teff = Temperature*(1.0 - fT2/Tf);
//	double tauR = ftauR0*exp(fQR/Teff);
	double coeff = -fC1/log10(exp(1))*(fC2*(Temperature - Tf) + Temperature*(Tf-fTg))/(Temperature*(fC2 + Tf - fTg));
	double tauR = ftaug*exp(coeff);
	/*check the limits*/
	if (tauR > ftauRH)
		tauR = ftauRH;
	else if (tauR < ftauRL)
		tauR = ftauRL;
	
	return(tauR);
}


double SMP_simple::ShearViscosity(const double Temperature, const double deltaneq, const double smag, const double sy)
{
	/*calculate the temperature part*/
//	double g = RetardationTime(Temperature, deltaneq)/ftauR0;
	double Tf = FictiveTemperature(deltaneq);
	double coeff = -fC1/log10(exp(1.0))*(fC2*(Temperature - Tf) + Temperature*(Tf-fTg))/(Temperature*(fC2 + Tf - fTg));
	double g = exp(coeff);
	double etaS = fetaS0*g*exp(-fQS/Temperature * smag/sy);
	
	return(etaS);
}

const dMatrixT& SMP_simple::ThermalDeformation_Inverse(void)
{
	/*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());

	double Temp = Compute_Temperature();
	
	/*thermal volume deformation*/
	double ThetaT = 1.0 + falphar*(Temp-fT0) + (-(falphar-falphag)*(Temp-fT0)+*fdelneq);
//	cout << "\nTheta: "<<ThetaT;
	fF_T_inv.Identity(pow(ThetaT,-1.0*third));
	return(fF_T_inv);
}

void SMP_simple::ComputeOutput(dArrayT& output)
{
	/*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());

	double Temp = Compute_Temperature();
	
	/*thermal volume deformation*/
	output[0] = 1.0 + falphar*(Temp-fT0) + (*fdelneq);

	/*neq thermal dilatation*/
	output[1] = (*fdelneq);

	/*effective viscous stretch*/
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fC_v[0], false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();
	output[2] = sqrt(third*(fEigs[0]+fEigs[1]+fEigs[2]));
	
	/*yield strength*/
	output[3] = *fsy;
}

/*************************************************************************
*	PUBLIC
**************************************************************************/
/* describe the parameters needed by the interface */
void SMP_simple::DefineParameters(ParameterListT& list) const
{
  /* inherited */
  RGViscoelasticityT::DefineParameters(list);

  /* common limit */
  LimitT positive(0.0, LimitT::Lower);

  ParameterT reftemp(ParameterT::Double, "ref_temperature");
  list.AddParameter(reftemp);

}

/* information about subordinate parameter lists */
void SMP_simple::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("smp_thermal_expansion", ParameterListT::Once);

	sub_list.AddSub("smp_eq_potential", ParameterListT::Once);
	sub_list.AddSub("smp_neq_potential", ParameterListT::Once);

	/* choice of viscosity */
	sub_list.AddSub("smp_retardation_time", ParameterListT::Once);
	sub_list.AddSub("smp_shear_viscosity", ParameterListT::Once);
}


/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SMP_simple::NewSub(const StringT& name) const
{
		LimitT zero(0.0, LimitT::Lower);
		LimitT one(1.0, LimitT::Lower);
		LimitT positive(0.0, LimitT::LowerInclusive);

	PotentialT* pot = NULL;
	if (name == "neo-hookean")
		pot = new NeoHookean;
	else if (name == "arruda-boyce")
		pot = new ArrudaBoyce;
	if (pot)
		return pot;

	/* inherited */
	ParameterInterfaceT* sub = RGViscoelasticityT::NewSub(name);
	if (sub) 
	{
		return sub;
	}
	else if (name == "smp_thermal_expansion")
	{
		ParameterContainerT* CTE = new ParameterContainerT(name);
		
		ParameterT alphar(ParameterT::Double, "high_temp_CTE");
		ParameterT alphag(ParameterT::Double, "low_temp_CTE");

		alphar.AddLimit(zero);
		alphag.AddLimit(zero);

		CTE->AddParameter(alphar);
		CTE->AddParameter(alphag);
		return CTE;
	}
	else if (name == "smp_eq_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
	
		/* choice of parameters */
		choice->AddSub("arruda-boyce");
		return(choice);
	}
	else if (name == "smp_neq_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
	
		/* choice of parameters */
		choice->AddSub("neo-hookean");
		return(choice);
	}
	else if (name == "smp_retardation_time")
	{
		ParameterContainerT* tauR = new ParameterContainerT(name);

		ParameterT tauRg(ParameterT::Double, "tauR_ref");
		ParameterT Tg(ParameterT::Double, "Tg");
		ParameterT C1(ParameterT::Double, "WLF_C1");
		ParameterT C2(ParameterT::Double, "WLF_C2");
		
		tauRg.AddLimit(zero);
		Tg.AddLimit(zero);

		tauR->AddParameter(tauRg);
		tauR->AddParameter(Tg);
		tauR->AddParameter(C1);
		tauR->AddParameter(C2);

		return(tauR);
	}
	else if (name == "smp_shear_viscosity")
	{
		ParameterContainerT* etaS = new ParameterContainerT(name);

		ParameterT etaSR(ParameterT::Double, "etaS_ref");
		ParameterT A(ParameterT::Double, "activation_energy");
		ParameterT sy_0(ParameterT::Double, "init_yield_strength");
		ParameterT sy_ss(ParameterT::Double, "sat_yield_strength");
		ParameterT h(ParameterT::Double, "hardening_modulus");
		
		etaSR.AddLimit(zero);
		A.AddLimit(zero);
		sy_0.AddLimit(zero);
		sy_ss.AddLimit(zero);
		h.AddLimit(positive);

		etaS->AddParameter(etaSR);
		etaS->AddParameter(A);
		etaS->AddParameter(sy_0);
		etaS->AddParameter(sy_ss);
		etaS->AddParameter(h);

		return(etaS);
	}
}

void SMP_simple::TakeParameterList(const ParameterListT& list)
{
  const char caller[] = "SMP_simple::TakeParameterList";
  /* inherited */
  FSSolidMatT::TakeParameterList(list);


  fT0 = list.GetParameter("ref_temperature");

  const ParameterListT* tm = list.List("smp_thermal_expansion");
  if (tm)
  {
	falphar = tm->GetParameter("high_temp_CTE");
	falphag = tm->GetParameter("low_temp_CTE");
  }
  
	fPot.Dimension(2);
		
	const ParameterListT& eq_pot = list.GetListChoice(*this, "smp_eq_potential");
	if(eq_pot.Name() == "arruda-boyce")
		fPot[0] = new ArrudaBoyce;
	else 
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot[0]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", eq_pot.Name().Pointer());			
	fPot[0]->TakeParameterList(eq_pot);
  
	const ParameterListT& neq_pot = list.GetListChoice(*this, "smp_neq_potential");
	if(neq_pot.Name() == "neo-hookean")
		fPot[1] = new NeoHookean;
	else 
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot[1]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", neq_pot.Name().Pointer());			
	fPot[1]->TakeParameterList(neq_pot);
  
  const ParameterListT* tauR = list.List("smp_retardation_time");
  if (tauR)
  {
	fTg = tauR->GetParameter("Tg");
	ftaug = tauR->GetParameter("tauR_ref");
	fC1 = tauR->GetParameter("WLF_C1");
	fC2 = tauR->GetParameter("WLF_C2");
	
	fT2 = fTg - fC2;
	fQR = fC1*fC2/log10(exp(1.0));
	ftauR0 = ftaug*exp(-fC1/log10(exp(1.0)));
	ftauRL = 1.0e-10*ftaug;
	ftauRH = 1.0e+10*ftaug;
  }
  
  const ParameterListT* etaS = list.List("smp_shear_viscosity");
  if (etaS)
  {
	fetaS0 = etaS->GetParameter("etaS_ref");
	fQS = etaS->GetParameter("activation_energy");
	fsy0 = etaS->GetParameter("init_yield_strength");
	fsinf = etaS->GetParameter("sat_yield_strength");
	fh = etaS->GetParameter("hardening_modulus");
  }
  
	/*set dimension of workspaces*/
	Initialize();
}

/*initializes history variable */
void  SMP_simple::PointInitialize(void)
{
	/* allocate element storage */
	ElementCardT& element = CurrentElement();	
	if (CurrIP() == 0)
	{
		ElementCardT& element = CurrentElement();
		element.Dimension(0, fnstatev*NumIP());
	
		/* initialize internal variables to identity*/
		for (int ip = 0; ip < NumIP(); ip++)
		{
		      /* load state variables */
		      Load(element, ip);
			  double Temp = Compute_Temperature();
		      
			  fC_v[0].Identity();
			  fC_vn[0].Identity();
			  *fdelneq = 0.0;
			  *fdelneq_n = 0.0;
			  *fsy = fsy0;
			  *fsy_n = fsy0;

		      /* write to storage */
		      Store(element, ip);
		}
	}
}
 
void SMP_simple::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		fC_vn[0] = fC_v[0];
		*fsy_n = *fsy;
		*fdelneq_n = *fdelneq;
		
		/* write to storage */
		Store(element, ip);
	}
}

void SMP_simple::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		fC_v[0] = fC_vn[0];
		*fsy = *fsy_n;
		*fdelneq = *fdelneq_n;
		
		/* write to storage */
		Store(element, ip);
	}
}

void SMP_simple::InitStep(void)
{
	/*inherited*/
	RGSplitT2::InitStep();
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


/***********************************************************************
 * Protected
 ***********************************************************************/
void SMP_simple::Initialize(void)
{
 /* dimension work space */
  
	/*Dimension workspace*/
	fC_v.Dimension(1);
	fC_vn.Dimension(1);
	
	int ndof = 3;
	int numstress = dSymMatrixT::NumValues(ndof);

	fnstatev = 0;
	fnstatev += numstress;   /*current C_v*/
	fnstatev += numstress;   /*last C_vn*/
	fnstatev ++;			/*current neq thermal strain delta_neq*/ 
	fnstatev ++;			/*last neq thermal strain*/
	fnstatev ++;			/*current yield strength*/
	fnstatev ++;			/*last yield strength*/
	
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
		
	/* assign pointers to current and last blocks of state variable array */
	fC_v[0].Set(ndof, pstatev);
	pstatev += numstress;
	fC_vn[0].Set(ndof, pstatev);
	pstatev += numstress;
	fdelneq = pstatev; 
	pstatev++;
	fdelneq_n = pstatev;
	pstatev++;
	fsy = pstatev;
	pstatev++;
	fsy_n = pstatev;
	pstatev++;


  fF_M.Dimension(ndof);
  fF_T_inv.Dimension(ndof);

  fF3D.Dimension(ndof);
  fInverse.Dimension(ndof);

  fb.Dimension(ndof);
  fbe.Dimension(ndof);
  fb_tr.Dimension(ndof);

  fEigs_dev.Dimension(ndof);
  fEigs.Dimension(ndof);
  fEigs_e.Dimension(ndof);
  fEigs_tr.Dimension(ndof);

  ftau_EQ.Dimension(ndof);
  ftau_NEQ.Dimension(ndof);

  fStress.Dimension(NumSD());
  fStress3D.Dimension(ndof);

  fDtauDe_EQ.Dimension(ndof);
  fDtauDe_NEQ.Dimension(ndof);

   fRes.Dimension(ndof+2);
   fDelta.Dimension(ndof+2);
   
   fiKAB.Dimension(ndof+2);
   fGAB.Dimension(ndof+2,ndof);
   fDAB.Dimension(ndof+2,ndof);
   fDABbar.Dimension(ndof);
   fMat.Dimension(ndof);
   fCalg.Dimension(ndof);
	
  fModulus3D.Dimension(dSymMatrixT::NumValues(ndof));
  fModMat.Dimension(dSymMatrixT::NumValues(ndof));
  fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));

}


/***********************************************************************
 * Private
 ***********************************************************************/
/* set inverse of thermal transformation - return true if active */
 void SMP_simple::Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
	const double& dtau_m, dMatrixT& Calg, const int type)
 {
		const dMatrixT& F = F_total();
		double iJ = 1.0/F.Det();

		const double& delneq = *fdelneq;
		const double& sy = *fsy;

		/*temperature and temperature step*/
		double Temp = Compute_Temperature();
		double Temp_n = Compute_Temperature_last();
		double dT = Temp - Temp_n;
	
		/*time step*/
		double dt = fFSMatSupport->TimeStep();

		double dalpha = falphar-falphag;


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
	    double smag = sqrt(s0*s0 + s1*s1 + s2*s2);

		/*calculate mobilities*/
		double Tf = FictiveTemperature(delneq);

		double tauR = RetardationTime(Temp, delneq);
		double etaS = ShearViscosity(Temp, delneq, smag, sy);
		double itauR = 1.0/tauR;
		double ietaS = 1.0/etaS;
						
		double gamdot = 0.5*smag*ietaS;

		/*calculate stiffness matrix*/
		/*derivative of retardation time wrt to Tf*/
		double dtauR_dTf = tauR*fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Tf - fTg)*(fC2 + Tf - fTg));
		double detaS_dTf = etaS*fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Tf - fTg)*(fC2 + Tf - fTg));
		double detaS_dsmag = -etaS*fQS / (Temp*sy);
		double detaS_dsy = etaS*fQS*smag / (Temp*sy*sy);
		
		/*initialize*/
		fiKAB = 0.0;
		
		/*K_del_del*/
		double hft = itauR*dtauR_dTf/dalpha;
		fiKAB(0,0) =(1.0 + dt*itauR) - dt*itauR*hft*(delneq-dalpha*(Temp-fT0));
		
		/*K_epA_del*/
		fiKAB(1,0) = -0.5*dt*ietaS*s0* 1.0/dalpha*ietaS*detaS_dTf;
		fiKAB(2,0) = -0.5*dt*ietaS*s1* 1.0/dalpha*ietaS*detaS_dTf;
		fiKAB(3,0) = -0.5*dt*ietaS*s2* 1.0/dalpha*ietaS*detaS_dTf;

		/*K_epA_epB*/
		double coef0 = (s0*c0 + s1*c01 + s2*c02);
		double coef1 = (s0*c01 + s1*c1 + s2*c12);
		double coef2 = (s0*c02 + s1*c12 + s2*c2);
		if (smag > kSmall)
		{
			coef0 /= smag;
			coef1 /= smag;
			coef2 /= smag;
		}
		fiKAB(1,1) = 1.0 + 0.5*ietaS*dt*c0 - 0.5*dt*ietaS*s0* ietaS*detaS_dsmag*coef0;
		fiKAB(2,2) = 1.0 + 0.5*ietaS*dt*c1 - 0.5*dt*ietaS*s1* ietaS*detaS_dsmag*coef1;
		fiKAB(3,3) = 1.0 + 0.5*ietaS*dt*c2 - 0.5*dt*ietaS*s2* ietaS*detaS_dsmag*coef2;
		
		fiKAB(2,3) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s1* ietaS*detaS_dsmag*coef2;
		fiKAB(1,3) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s0* ietaS*detaS_dsmag*coef2;
		fiKAB(1,2) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s0* ietaS*detaS_dsmag*coef1;

		fiKAB(3,2) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s2* ietaS*detaS_dsmag*coef1;
		fiKAB(3,1) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s2* ietaS*detaS_dsmag*coef0;
		fiKAB(2,1) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s1* ietaS*detaS_dsmag*coef0;
       
		/*K_epA_sy*/
		fiKAB(1,4) = -0.5*dt*ietaS*s0 *ietaS*detaS_dsy;
		fiKAB(2,4) = -0.5*dt*ietaS*s1 *ietaS*detaS_dsy;
		fiKAB(3,4) = -0.5*dt*ietaS*s2 *ietaS*detaS_dsy;
	
		/*K_sy_del*/ 
		fiKAB(4,0) = 0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*smag*ietaS *1.0/dalpha*detaS_dTf;

		/*K_sy_epB*/
		fiKAB(4,1) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*coef0*(1.0 - smag*ietaS*detaS_dsmag);
		fiKAB(4,2) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*coef1*(1.0 - smag*ietaS*detaS_dsmag);
		fiKAB(4,3) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*coef2*(1.0 - smag*ietaS*detaS_dsmag);
		
		/*K_sy_sy*/
	    fiKAB(4,4) = 1.0 + 0.5*dt*ietaS*fh* (smag/fsinf + (1.0 - sy/fsinf)*smag*ietaS*detaS_dsy);

		/*inverts KAB*/
		fiKAB.Inverse();

		/*initialize*/
		fGAB = 0.0;
		
		/*G_epeA_epB*/
		fGAB(1,0) = 1.0 + 0.5*dt*ietaS*s0 - 0.5*dt*ietaS *ietaS*s0*s0;
		fGAB(1,1) = 0.5*dt*ietaS*s0 - 0.5*dt*ietaS *ietaS*s0*s0;
		fGAB(1,2) = 0.5*dt*ietaS*s0 - 0.5*dt*ietaS *ietaS*s0*s0;
		
		fGAB(2,0) = 0.5*dt*ietaS*s1 - 0.5*dt*ietaS *ietaS*s1*s1;
		fGAB(2,1) = 1.0 + 0.5*dt*ietaS*s1 - 0.5*dt*ietaS *ietaS*s1*s1;
		fGAB(2,2) = 0.5*dt*ietaS*s1 - 0.5*dt*ietaS *ietaS*s1*s1;

		fGAB(3,0) = 0.5*dt*ietaS*s2 - 0.5*dt*ietaS *ietaS*s2*s2;
		fGAB(3,1) = 0.5*dt*ietaS*s2 - 0.5*dt*ietaS *ietaS*s2*s2;
		fGAB(3,2) = 1.0 + 0.5*dt*ietaS*s2 - 0.5*dt*ietaS *ietaS*s2*s2;
		
		/*G_sy_epB*/
//		double coef = (s0*s0 + s1*s1 + s2*s2)/smag;
		fGAB(4,0) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*smag*(1.0 - smag*ietaS*detaS_dsmag);
		fGAB(4,1) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*smag*(1.0 - smag*ietaS*detaS_dsmag);
		fGAB(4,2) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*smag*(1.0 - smag*ietaS*detaS_dsmag);
		
		/*Calg = dtau/depe*fiKA*fG	*/	
		/*calculating delta_internval_vars = K^-1.G. delta_epsilon*/
		fDAB.MultAB(fiKAB,fGAB);
		/*copy subset*/
	
		for (int i = 0; i< fDABbar.Rows(); i++)
			for (int j = 0; j< fDABbar.Cols(); j++)
				fDABbar(i,j) = fDAB(i+1,j);
				
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

void SMP_simple::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
			     dArrayT& eigenstress, dSymMatrixT& eigenmodulus,  const int type) 
{		
	const double ctol = 1.00e-12;
		
	/*set references to principle stretches*/
     
	double& le0 = eigenstretch_e[0];
	double& le1 = eigenstretch_e[1];
	double& le2 = eigenstretch_e[2];
	
	double& delneq = *fdelneq;
	double& sy = *fsy;
  
	double tol;

	/*initialize principle elastic and trial elastic log strains */
	const double ep_tr0 = 0.5*log(le0);
	const double ep_tr1 = 0.5*log(le1);
	const double ep_tr2 = 0.5*log(le2);
	const double syn = *fsy_n;
	const double delneq_n = *fdelneq_n;
	
	double ep_e0 = ep_tr0;		
	double ep_e1 = ep_tr1;	
	double ep_e2 = ep_tr2;

	/*jacobian*/
	const dMatrixT& F = F_total();
	double iJ = 1.0/F.Det();

	/*time step*/
	double dt = fFSMatSupport->TimeStep();

	/*temperature and temperature step*/
	double Temp = Compute_Temperature();
	double Temp_n = Compute_Temperature_last();
	double dT = Temp - Temp_n;
	
	double dalpha = falphar-falphag;
	int maxiteration = 100;
	/*initializes principle viscous stretch*/
	int iteration = 0;
	do 
	{
		iteration ++;
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
	    		
	    /*caculate smag*/
	    double smag = sqrt(s0*s0 + s1*s1 + s2*s2);
//		cout << "\nsmag: "<<smag;

		/*calculate mobilities*/
		double Tf = FictiveTemperature(delneq);
		double tauR = RetardationTime(Temp, delneq);
		double etaS = ShearViscosity(Temp, delneq, smag, sy);
//		cout << "\ntauR: "<<tauR;
//		cout << "\netaS: "<<etaS;
		
		double itauR = 1.0/tauR;
		double ietaS = 1.0/etaS;
						
		double gamdot = 0.5*smag*ietaS;

	    /*calculate the residual*/
		fRes[0] = delneq + dt*itauR*(delneq-dalpha*(Temp-fT0)) - delneq_n;
		
	    fRes[1] = ep_e0 + 0.5*dt*ietaS*s0 - ep_tr0;
	    fRes[2] = ep_e1 + 0.5*dt*ietaS*s1 - ep_tr1;
	    fRes[3] = ep_e2 + 0.5*dt*ietaS*s2 - ep_tr2;
		
		fRes[4] = sy - dt*fh*(1.0-sy/fsinf)*gamdot - syn;
//		cout << "\nfRes: "<<fRes;
		/*calculate stiffness matrix*/
		/*derivative of retardation time wrt to Tf*/
		double dtauR_dTf = tauR*fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Tf - fTg)*(fC2 + Tf - fTg));
		double detaS_dTf = etaS*fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Tf - fTg)*(fC2 + Tf - fTg));
		double detaS_dsmag = -etaS*fQS / (Temp*sy);
		double detaS_dsy = etaS*fQS*smag / (Temp*sy*sy);
//		cout << "\ndtauR_dTf: "<<dtauR_dTf;
//		cout << "\ndetaS_dTf: "<<detaS_dTf;
//		cout << "\ndetaS_dsmag: "<<detaS_dsmag;
//		cout << "\ndetaS_dsy: "<<detaS_dsy;
		
		/*initialize*/
		fiKAB = 0.0;
		
		/*K_del_del*/
		double hft = itauR*dtauR_dTf/dalpha;
		fiKAB(0,0) =(1.0 + dt*itauR) - dt*itauR*hft*(delneq-dalpha*(Temp-fT0));
		
		/*K_epA_del*/
		fiKAB(1,0) = -0.5*dt*ietaS*s0* 1.0/dalpha*ietaS*detaS_dTf;
		fiKAB(2,0) = -0.5*dt*ietaS*s1* 1.0/dalpha*ietaS*detaS_dTf;
		fiKAB(3,0) = -0.5*dt*ietaS*s2* 1.0/dalpha*ietaS*detaS_dTf;

		/*K_epA_epB*/
		double coef0 = (s0*c0 + s1*c01 + s2*c02);
		double coef1 = (s0*c01 + s1*c1 + s2*c12);
		double coef2 = (s0*c02 + s1*c12 + s2*c2);
		if (smag > kSmall)
		{
			coef0 /= smag;
			coef1 /= smag;
			coef2 /= smag;
		}
		fiKAB(1,1) = 1.0 + 0.5*ietaS*dt*c0 - 0.5*dt*ietaS*s0* ietaS*detaS_dsmag*coef0;
		fiKAB(2,2) = 1.0 + 0.5*ietaS*dt*c1 - 0.5*dt*ietaS*s1* ietaS*detaS_dsmag*coef1;
		fiKAB(3,3) = 1.0 + 0.5*ietaS*dt*c2 - 0.5*dt*ietaS*s2* ietaS*detaS_dsmag*coef2;
		
		fiKAB(2,3) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s1* ietaS*detaS_dsmag*coef2;
		fiKAB(1,3) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s0* ietaS*detaS_dsmag*coef2;
		fiKAB(1,2) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s0* ietaS*detaS_dsmag*coef1;

		fiKAB(3,2) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s2* ietaS*detaS_dsmag*coef1;
		fiKAB(3,1) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s2* ietaS*detaS_dsmag*coef0;
		fiKAB(2,1) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s1* ietaS*detaS_dsmag*coef0;
       
		/*K_epA_sy*/
		fiKAB(1,4) = -0.5*dt*ietaS*s0 *ietaS*detaS_dsy;
		fiKAB(2,4) = -0.5*dt*ietaS*s1 *ietaS*detaS_dsy;
		fiKAB(3,4) = -0.5*dt*ietaS*s2 *ietaS*detaS_dsy;
	
		/*K_sy_del*/ 
		fiKAB(4,0) = 0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*smag*ietaS *1.0/dalpha*detaS_dTf;

		/*K_sy_epB*/
		fiKAB(4,1) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*coef0*(1.0 - smag*ietaS*detaS_dsmag);
		fiKAB(4,2) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*coef1*(1.0 - smag*ietaS*detaS_dsmag);
		fiKAB(4,3) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*coef2*(1.0 - smag*ietaS*detaS_dsmag);
		
		/*K_sy_sy*/
	    fiKAB(4,4) = 1.0 + 0.5*dt*ietaS*fh* (smag/fsinf + (1.0 - sy/fsinf)*smag*ietaS*detaS_dsy);

//		cout << "\nfiKAB: "<<fiKAB;
		/*inverts KAB*/
		fiKAB.Inverse();
	    
		
	    /*solve for the principal strain increments*/
		fiKAB.Multx(fRes, fDelta, -1.0);
		
	    /*updates principal elastic stretches*/ 
		delneq += fDelta[0];
		
	    ep_e0 += fDelta[1];
	    ep_e1 += fDelta[2];
	    ep_e2 += fDelta[3];
	    
		sy += fDelta[4];
		
	    le0 = exp(2.0*ep_e0);
	    le1 = exp(2.0*ep_e1);
	    le2 = exp(2.0*ep_e2);
	    
	    /*Check that the L2 norm of the residual is less than tolerance*/
	    tol = sqrt(dArrayT::Dot(fRes, fRes));
/*		cout << "\ntol: "<<tol;
		cout << "\neps_e: "<<eigenstretch_e;
*/
	}while (tol>ctol && iteration < maxiteration); 
	if (iteration >= 10) 
		ExceptionT::GeneralFail("SMP_simple::ComputeEigs_e", 
			"number of iteration exceeds maximum");
}
