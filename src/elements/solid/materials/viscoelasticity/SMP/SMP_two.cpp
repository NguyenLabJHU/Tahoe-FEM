/* $Id: SMP_two.cpp,v 1.1 2009-04-23 02:51:50 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "SMP_two.h"

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
const double small = 1.0e-12;

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
/* constructors */
SMP_two::SMP_two(void):
  ParameterInterfaceT("SMP_two")
{
	fNumProcess = 2;
}

double SMP_two::ShearViscosityConst(const double Temperature, const double deltaneq)
{
	/*calculate the temperature part*/
	double Tf = FictiveTemperature(deltaneq);
	double coeff = -fC1/log10(exp(1.0))*(fC2*(Temperature - Tf) + Temperature*(Tf-fTg))/(Temperature*(fC2 + Tf - fTg));
	double g = exp(coeff);

	double etaS = fetaS1*g;		
	return(etaS);
}


/*************************************************************************
*	PUBLIC
**************************************************************************/
/* describe the parameters needed by the interface */

/* information about subordinate parameter lists */
void SMP_two::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SMP_simple::DefineSubs(sub_list);

	sub_list.AddSub("smp_neq_relax_potential", ParameterListT::Once);

	/* choice of viscosity */
	sub_list.AddSub("smp_shear_relax_viscosity", ParameterListT::Once);
}


/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SMP_two::NewSub(const StringT& name) const
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
	ParameterInterfaceT* sub = SMP_simple::NewSub(name);
	if (sub) 
	{
		return sub;
	}
	else if (name == "smp_neq_relax_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
	
		/* choice of parameters */
		choice->AddSub("neo-hookean");
		return(choice);
	}
	else if (name == "smp_shear_relax_viscosity")
	{
		ParameterContainerT* etaS = new ParameterContainerT(name);

		ParameterT etaSR(ParameterT::Double, "etaS_rubbery");
			
		etaSR.AddLimit(zero);

		etaS->AddParameter(etaSR);

		return(etaS);
	}
}

void SMP_two::TakeParameterList(const ParameterListT& list)
{
  const char caller[] = "SMP_two::TakeParameterList";
  /* inherited */
  FSSolidMatT::TakeParameterList(list);
  fNumProcess = 2;

  fT0 = list.GetParameter("ref_temperature");
  cout << "\nfT0: "<<fT0;
  const ParameterListT* tm = list.List("smp_thermal_expansion");
  if (tm)
  {
	falphar = tm->GetParameter("high_temp_CTE");
	falphag = tm->GetParameter("low_temp_CTE");
  }
  
	fPot.Dimension(3);
		
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
 
  	const ParameterListT& neq_relax_pot = list.GetListChoice(*this, "smp_neq_relax_potential");
	if(neq_pot.Name() == "neo-hookean")
		fPot[2] = new NeoHookean;
	else 
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot[2]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", neq_pot.Name().Pointer());			
	fPot[2]->TakeParameterList(neq_pot);
 
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

  const ParameterListT* etaR = list.List("smp_shear_relax_viscosity");
  if (etaR)
  {
	fetaS1 = etaR->GetParameter("etaS_rubbery");
  }
  
	/*set dimension of workspaces*/
	Initialize();
}

/*initializes history variable */
void  SMP_two::PointInitialize(void)
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
			  fC_v[1].Identity();
			  fC_vn[1].Identity();
			  *fdelneq = 0.0;
			  *fdelneq_n = 0.0;
			  *fsy = fsy0;
			  *fsy_n = fsy0;
		      /* write to storage */
		      Store(element, ip);
		}
		
	}
}
 
void SMP_two::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		fC_vn[0] = fC_v[0];
		fC_vn[1] = fC_v[1];
		*fsy_n = *fsy;
		*fdelneq_n = *fdelneq;
		
		/* write to storage */
		Store(element, ip);
	}
}

void SMP_two::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		fC_v[0] = fC_vn[0];
		fC_v[1] = fC_vn[1];
		*fsy = *fsy_n;
		*fdelneq = *fdelneq_n;
		
		/* write to storage */
		Store(element, ip);
	}
}



/***********************************************************************
 * Protected
 ***********************************************************************/
void SMP_two::Initialize(void)
{
 /* dimension work space */
  
	/*Dimension workspace*/
	fC_v.Dimension(2);
	fC_vn.Dimension(2);
	
	int ndof = 3;
	int numstress = dSymMatrixT::NumValues(ndof);

	fnstatev = 0;
	fnstatev += numstress;   /*current C_v*/
	fnstatev += numstress;   /*last C_vn*/
	fnstatev += numstress;   /*current C_v*/
	fnstatev += numstress;   /*last C_vn*/
	fnstatev ++;			/*current neq thermal strain delta_neq*/ 
	fnstatev ++;			/*last neq thermal strain*/
	fnstatev ++;			/*current yield strength*/
	fnstatev ++;			/*last yield strength*/
	cout << "\nfnstatev: "<<fnstatev;
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
		
	/* assign pointers to current and last blocks of state variable array */
	fC_v[0].Set(ndof, pstatev);
	pstatev += numstress;
	fC_vn[0].Set(ndof, pstatev);
	pstatev += numstress;
	fC_v[1].Set(ndof, pstatev);
	pstatev += numstress;
	fC_vn[1].Set(ndof, pstatev);
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
 void SMP_two::Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
	const double& dtau_m, dMatrixT& Calg, const int type)
 {
		const dMatrixT& F = F_total();
		double iJ = 1.0/F.Det();

		if (type == 0)
		{
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
			double smag = sqrt(0.5*(s0*s0 + s1*s1 + s2*s2));

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
			double detaS_dTf = (etaS)*fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Tf - fTg)*(fC2 + Tf - fTg));
			double x = fQS*smag/(Temp*sy*sqrt(2.0));
			double detaS_dsmag = 0.0;
			double detaS_dsy = (etaS)/sy;
			if (smag > small)
			{
				double cothx = cosh(x)/sinh(x);
				detaS_dsmag = (etaS)*(1.0-x*cothx)/smag;
				detaS_dsy *= (-1.0 + x*cothx);
			}
		
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
			double coef0 = 0.5*(s0*c0 + s1*c01 + s2*c02);
			double coef1 = 0.5*(s0*c01 + s1*c1 + s2*c12);
			double coef2 = 0.5*(s0*c02 + s1*c12 + s2*c2);
			if (smag > small)
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
		}
		else
		{
			const double& delneq = *fdelneq;

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
						
			/*calculate mobilities*/
			double Tf = FictiveTemperature(delneq);

			double tauR = RetardationTime(Temp, delneq);
			double etaS = ShearViscosityConst(Temp, delneq);
			double itauR = 1.0/tauR;
			double ietaS = 1.0/etaS;

			/*calculate stiffness matrix*/
			/*derivative of retardation time wrt to Tf*/
			double dtauR_dTf = tauR*fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Tf - fTg)*(fC2 + Tf - fTg));
			double detaS_dTf = (etaS)*fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Tf - fTg)*(fC2 + Tf - fTg));
		
			/*initialize*/
			fiKAB = 0.0;
		
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
			fiKAB(1,1) = 1.0 + 0.5*ietaS*dt*c0;
			fiKAB(2,2) = 1.0 + 0.5*ietaS*dt*c1;
			fiKAB(3,3) = 1.0 + 0.5*ietaS*dt*c2;
			
			fiKAB(2,3) = 0.5*ietaS*dt*c12;
			fiKAB(1,3) = 0.5*ietaS*dt*c02;
			fiKAB(1,2) = 0.5*ietaS*dt*c01;

			fiKAB(3,2) = 0.5*ietaS*dt*c12;
			fiKAB(3,1) = 0.5*ietaS*dt*c02;
			fiKAB(2,1) = 0.5*ietaS*dt*c01;
					
			/*K_sy_sy*/
			fiKAB(4,4) = 1.0 + 0.5*dt*ietaS;

			/*inverts KAB*/
			fiKAB.Inverse();

			/*initialize*/
			fGAB = 0.0;
		
			/*G_epeA_epB*/
			fGAB(1,0) = 1.0 + 0.5*dt*ietaS*s0;
			fGAB(1,1) = 0.5*dt*ietaS*s0;
			fGAB(1,2) = 0.5*dt*ietaS*s0;
		
			fGAB(2,0) = 0.5*dt*ietaS*s1;
			fGAB(2,1) = 1.0 + 0.5*dt*ietaS*s1;
			fGAB(2,2) = 0.5*dt*ietaS*s1;

			fGAB(3,0) = 0.5*dt*ietaS*s2;
			fGAB(3,1) = 0.5*dt*ietaS*s2;
			fGAB(3,2) = 1.0 + 0.5*dt*ietaS*s2;
		
		}

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

void SMP_two::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
			     dArrayT& eigenstress, dSymMatrixT& eigenmodulus,  const int type) 
{		
	const double ctol = 1.00e-8;
		
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
	double Je=sqrt(le0*le1*le2);
	fEigs_dev = eigenstretch_e;
	fEigs_dev *= pow(Je,-2.0*third);

	/*calculate stresses and moduli*/
	fPot[type+1]->DevStress(fEigs_dev, eigenstress);
	    
	double s0 = iJ*eigenstress[0];
	double s1 = iJ*eigenstress[1];
	double s2 = iJ*eigenstress[2];
		    		
	/*caculate smag*/
	double smag = sqrt(0.5*(s0*s0 + s1*s1 + s2*s2));
//	cout << "\nsmag: "<<smag;
	
	/*calculate mobilities*/
	double Tf = FictiveTemperature(delneq);
	double tauR = RetardationTime(Temp, delneq);
	double etaS;
	double gamdot;
	if(type==0)
	{
		etaS=ShearViscosity(Temp, delneq, smag, sy);
		gamdot = 0.5*smag/etaS;
	}
	else
		etaS=ShearViscosityConst(Temp,delneq);
			
	double itauR = 1.0/tauR;
	double ietaS = 1.0/etaS;
					

	/*calculate the residual*/
	if(type ==0)
	{
		fRes[0] = delneq + dt*itauR*(delneq-dalpha*(Temp-fT0)) - delneq_n;
	
		fRes[1] = ep_e0 + 0.5*dt*ietaS*s0 - ep_tr0;
		fRes[2] = ep_e1 + 0.5*dt*ietaS*s1 - ep_tr1;
		fRes[3] = ep_e2 + 0.5*dt*ietaS*s2 - ep_tr2;
	
		fRes[4] = sy - dt*fh*(1.0-sy/fsinf)*gamdot - syn;
	}
	else
	{
		fRes[0] = 0;
	
		fRes[1] = ep_e0 + 0.5*dt*ietaS*s0 - ep_tr0;
		fRes[2] = ep_e1 + 0.5*dt*ietaS*s1 - ep_tr1;
		fRes[3] = ep_e2 + 0.5*dt*ietaS*s2 - ep_tr2;
	
		fRes[4] = 0;		
//		cout << "\nfRes: "<<fRes;
	}
	tol = sqrt(dArrayT::Dot(fRes, fRes));
/*	cout << "\ntol: "<<tol;
	cout << "\neps_e: "<<eigenstretch_e;
*/
//	cout << "\ntype: "<<type;
	int iteration = 0;
	while (tol>ctol && iteration < maxiteration)
	{
		iteration ++;
//		cout <<"\ntype: "<< type;
//		cout << "\niteration: "<<iteration;
		if(type==0)
		{
			double dtauR_dTf = tauR*fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Tf - fTg)*(fC2 + Tf - fTg));
			double detaS_dTf = (etaS)*fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Tf - fTg)*(fC2 + Tf - fTg));
			double x = fQS*smag/(Temp*sy*sqrt(2.0));
			double detaS_dsmag = 0.0;
			double detaS_dsy = (etaS)/sy;
			if (smag > small)
			{
				double cothx = cosh(x)/sinh(x);
				detaS_dsmag = (etaS)*(1.0-x*cothx)/smag;
				detaS_dsy *= (-1.0 + x*cothx);
			}

			fPot[type+1]->DevMod(fEigs_dev,eigenmodulus);
			/*deviatoric values*/
			double c0 = iJ*eigenmodulus(0,0);
			double c1 = iJ*eigenmodulus(1,1);
			double c2 = iJ*eigenmodulus(2,2);

			double c12 = iJ*eigenmodulus(1,2);
			double c02 = iJ*eigenmodulus(0,2);
			double c01 = iJ*eigenmodulus(0,1);
		
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
			double coef0 = 0.5*(s0*c0 + s1*c01 + s2*c02);
			double coef1 = 0.5*(s0*c01 + s1*c1 + s2*c12);
			double coef2 = 0.5*(s0*c02 + s1*c12 + s2*c2);
			if (smag > small)
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
			fiKAB(4,0) = sqrt(0.5)*dt*ietaS*fh*(1.0 - sy/fsinf)*smag*ietaS *1.0/dalpha*detaS_dTf;

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
	    
			Je=sqrt(le0*le1*le2);
			fEigs_dev = eigenstretch_e;
			fEigs_dev *= pow(Je,-2.0*third);

		//	cout << "\neigenstretch_edev: "<<fEigs_dev;

			/*calculate stresses and moduli*/
			fPot[type+1]->DevStress(fEigs_dev, eigenstress);
	    
			s0 = iJ*eigenstress[0];
			s1 = iJ*eigenstress[1];
			s2 = iJ*eigenstress[2];
	    
					
			/*caculate smag*/
			smag = sqrt(0.5*(s0*s0 + s1*s1 + s2*s2));
	//		cout << "\nsmag: "<<smag;

			/*calculate mobilities*/
			Tf = FictiveTemperature(delneq);
			tauR = RetardationTime(Temp, delneq);
			etaS = ShearViscosity(Temp, delneq, smag, sy);
	//		cout << "\ntauR: "<<tauR;
	//		cout << "\netaS: "<<etaS;
		
			itauR = 1.0/tauR;
			ietaS = 1.0/etaS;
						
			gamdot =0.5*smag*ietaS;
	//		cout << "\ngamdot: "<<gamdot;

			/*calculate the residual*/
			fRes[0] = delneq + dt*itauR*(delneq-dalpha*(Temp-fT0)) - delneq_n;
		
			fRes[1] = ep_e0 + 0.5*dt*ietaS*s0 - ep_tr0;
			fRes[2] = ep_e1 + 0.5*dt*ietaS*s1 - ep_tr1;
			fRes[3] = ep_e2 + 0.5*dt*ietaS*s2 - ep_tr2;
		
			fRes[4] = sy - dt*fh*(1.0-sy/fsinf)*gamdot - syn;
			/*Check that the L2 norm of the residual is less than tolerance*/
			tol = sqrt(dArrayT::Dot(fRes, fRes));
	/*		cout << "\ntol: "<<tol;
			cout << "\neps_e: "<<eigenstretch_e;
	*/
		}
		else
		{
//		cout << "\netaS: "<<etaS;
	
			double dtauR_dTf = tauR*fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Tf - fTg)*(fC2 + Tf - fTg));
			double detaS_dTf = (etaS)*fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Tf - fTg)*(fC2 + Tf - fTg));

			fPot[type+1]->DevMod(fEigs_dev,eigenmodulus);
			/*deviatoric values*/
			double c0 = iJ*eigenmodulus(0,0);
			double c1 = iJ*eigenmodulus(1,1);
			double c2 = iJ*eigenmodulus(2,2);

			double c12 = iJ*eigenmodulus(1,2);
			double c02 = iJ*eigenmodulus(0,2);
			double c01 = iJ*eigenmodulus(0,1);
		
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
			fiKAB(1,1) = 1.0 + 0.5*ietaS*dt*c0;
			fiKAB(2,2) = 1.0 + 0.5*ietaS*dt*c1;
			fiKAB(3,3) = 1.0 + 0.5*ietaS*dt*c2;
			
			fiKAB(2,3) = 0.5*ietaS*dt*c12;
			fiKAB(1,3) = 0.5*ietaS*dt*c02;
			fiKAB(1,2) = 0.5*ietaS*dt*c01;

			fiKAB(3,2) = 0.5*ietaS*dt*c12;
			fiKAB(3,1) = 0.5*ietaS*dt*c02;
			fiKAB(2,1) = 0.5*ietaS*dt*c01;
					
			/*K_sy_sy*/
			fiKAB(4,4) = 1.0 + 0.5*dt*ietaS;

//			cout << "\nfiKAB: "<<fiKAB;
			/*inverts KAB*/
			fiKAB.Inverse();
	    
			
			/*solve for the principal strain increments*/
			fiKAB.Multx(fRes, fDelta, -1.0);
		
			/*updates principal elastic stretches*/ 
			ep_e0 += fDelta[1];
			ep_e1 += fDelta[2];
			ep_e2 += fDelta[3];

			le0 = exp(2.0*ep_e0);
			le1 = exp(2.0*ep_e1);
			le2 = exp(2.0*ep_e2);
	    
			Je=sqrt(le0*le1*le2);
			fEigs_dev = eigenstretch_e;
			fEigs_dev *= pow(Je,-2.0*third);

		//	cout << "\neigenstretch_edev: "<<fEigs_dev;

			/*calculate stresses and moduli*/
			fPot[type+1]->DevStress(fEigs_dev, eigenstress);
	    
			s0 = iJ*eigenstress[0];
			s1 = iJ*eigenstress[1];
			s2 = iJ*eigenstress[2];
					
			/*calculate mobilities*/
			Tf = FictiveTemperature(delneq);
			tauR = RetardationTime(Temp, delneq);
			etaS = ShearViscosityConst(Temp, delneq);
		
			itauR = 1.0/tauR;
			ietaS = 1.0/etaS;
						
			/*calculate the residual*/
			fRes[0] = 0;
		
			fRes[1] = ep_e0 + 0.5*dt*ietaS*s0 - ep_tr0;
			fRes[2] = ep_e1 + 0.5*dt*ietaS*s1 - ep_tr1;
			fRes[3] = ep_e2 + 0.5*dt*ietaS*s2 - ep_tr2;
		
			fRes[4] = 0;
			/*Check that the L2 norm of the residual is less than tolerance*/
			tol = sqrt(dArrayT::Dot(fRes, fRes));
	/*		cout << "\ntol: "<<tol;
			cout << "\neps_e: "<<eigenstretch_e;
	*/
//			cout << "\nRes: "<<fRes;
		}
	}
	if (iteration >= maxiteration) 
		ExceptionT::GeneralFail("SMP_two::ComputeEigs_e", 
			"number of iteration exceeds maximum");
}
