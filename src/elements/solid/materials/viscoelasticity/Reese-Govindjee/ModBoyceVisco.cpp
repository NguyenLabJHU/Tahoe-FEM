/* $Id: ModBoyceVisco.cpp,v 1.1 2007-07-25 06:43:26 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "ModBoyceVisco.h"

#include "PotentialT.h"
#include "MooneyRivlin.h"
#include "NeoHookean.h"
#include "VWPotentialT.h"
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
const int kNumOutputVar = 3; 
static const char* Labels[kNumOutputVar] = {"lv_eff", "smag","sy"}; 

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
/* constructors */
ModBoyceVisco::ModBoyceVisco(void):
  ParameterInterfaceT("ModBoyceVisco")
{
	fNumProcess = 1;
//	cout << setprecision(16);
}


double ModBoyceVisco::ShearViscosity(const double smag, const double sy)
{
	/*calculate the temperature part*/
	double etaS = fetaS0*exp(fQS/fTemperature *(1.0 - smag/sy));
	
	return(etaS);
}

void ModBoyceVisco::ComputeOutput(dArrayT& output)
{
	/*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
	
	/*effective viscous stretch*/
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fC_v[0], false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();
	output[0] = sqrt(third*(fEigs[0]+fEigs[1]+fEigs[2]));
	
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
	
	output[1] = sqrt(ftau_NEQ[0]*ftau_NEQ[0] + ftau_NEQ[1]*ftau_NEQ[1] + ftau_NEQ[2]*ftau_NEQ[2]);

	/*yield strength*/
	output[2] = *fsy;
}

/*************************************************************************
*	PUBLIC
**************************************************************************/
/* describe the parameters needed by the interface */
void ModBoyceVisco::DefineParameters(ParameterListT& list) const
{
  /* inherited */
  RGViscoelasticityT::DefineParameters(list);
  /* common limit */
  LimitT positive(0.0, LimitT::Lower);

  ParameterT reftemp(ParameterT::Double, "temperature");
  list.AddParameter(reftemp);
}

/* information about subordinate parameter lists */
void ModBoyceVisco::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	RGViscoelasticityT::DefineSubs(sub_list);

	sub_list.AddSub("rg_eq_potential", ParameterListT::Once);
	sub_list.AddSub("rg_neq_potential", ParameterListT::Once);

	/* choice of viscosity */
	sub_list.AddSub("modboyce_shear_viscosity", ParameterListT::Once);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ModBoyceVisco::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = RGSplitT2::NewSub(name);
	if (sub) 
	{
		return sub;
	}
	else if (name == "modboyce_shear_viscosity")
	{
		LimitT zero(0.0, LimitT::Lower);
		LimitT one(1.0, LimitT::Lower);
		LimitT positive(0.0, LimitT::LowerInclusive);

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

void ModBoyceVisco::TakeParameterList(const ParameterListT& list)
{
  const char caller[] = "ModBoyceVisco::TakeParameterList";
  /* inherited */
  FSSolidMatT::TakeParameterList(list);

  /*Get reference temperature*/
  fTemperature = list.GetParameter("temperature");

  int num_neq =  list.NumLists("rg_neq_potential");
  int num_shear_visc = list.NumLists("modboyce_shear_viscosity");
  if (num_neq != num_shear_visc)
	 ExceptionT::GeneralFail("ModBoyceVisco::TakeParameterList", 
		"number of  viscosity functions does not match number of nonequilibrium potentials");
  fNumProcess = list.NumLists("modboyce_shear_viscosity");

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
  
	const ParameterListT* etaS = list.List("modboyce_shear_viscosity");
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
void  ModBoyceVisco::PointInitialize(void)
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
		      
			  fC_v[0].Identity();
			  fC_vn[0].Identity();
			  *fsy = fsy0;
			  *fsy_n = fsy0;

		      /* write to storage */
		      Store(element, ip);
		}
	}
}
 
void ModBoyceVisco::UpdateHistory(void)
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
		
		/* write to storage */
		Store(element, ip);
	}
}

void ModBoyceVisco::ResetHistory(void)
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
		
		/* write to storage */
		Store(element, ip);
	}
}

void ModBoyceVisco::InitStep(void)
{
	/*inherited*/
	RGSplitT2::InitStep();
}

int ModBoyceVisco::NumOutputVariables() const {return kNumOutputVar;} 

void ModBoyceVisco::OutputLabels(ArrayT<StringT>& labels) const 
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
void ModBoyceVisco::Initialize(void)
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
	fnstatev ++;			/*current yield strength*/
	fnstatev ++;			/*last yield strength*/
	
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
		
	/* assign pointers to current and last blocks of state variable array */
	fC_v[0].Set(ndof, pstatev);
	pstatev += numstress;
	fC_vn[0].Set(ndof, pstatev);
	pstatev += numstress;
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

	fRes.Dimension(ndof+1);
	fDelta.Dimension(ndof+1);
   
	fiKAB.Dimension(ndof+1);
	fGAB.Dimension(ndof+1,ndof);
	fDAB.Dimension(ndof+1,ndof);
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
 void ModBoyceVisco::Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
	const double& dtau_m, dMatrixT& Calg, const int type)
 {
		const dMatrixT& F = F_total();
		double iJ = 1.0/F.Det();

		const double& sy = *fsy;

		/*temperature and temperature step*/
		/*time step*/
		double dt = fFSMatSupport->TimeStep();

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

		double etaS = ShearViscosity(smag, sy);
		double ietaS = 1.0/etaS;
						
		double gamdot = 0.5*smag*ietaS;

		/*calculate stiffness matrix*/
		/*derivative of retardation time wrt to Tf*/
		double detaS_dsmag = -etaS*fQS / (fTemperature*sy);
		double detaS_dsy = etaS*fQS*smag / (fTemperature*sy*sy);
		
		/*initialize*/
		fiKAB = 0.0;
		
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
		fiKAB(0,0) = 1.0 + 0.5*ietaS*dt*c0 - 0.5*dt*ietaS*s0* ietaS*detaS_dsmag*coef0;
		fiKAB(1,1) = 1.0 + 0.5*ietaS*dt*c1 - 0.5*dt*ietaS*s1* ietaS*detaS_dsmag*coef1;
		fiKAB(2,2) = 1.0 + 0.5*ietaS*dt*c2 - 0.5*dt*ietaS*s2* ietaS*detaS_dsmag*coef2;
		
		fiKAB(1,2) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s1* ietaS*detaS_dsmag*coef2;
		fiKAB(0,2) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s0* ietaS*detaS_dsmag*coef2;
		fiKAB(0,1) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s0* ietaS*detaS_dsmag*coef1;

		fiKAB(2,1) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s2* ietaS*detaS_dsmag*coef1;
		fiKAB(2,0) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s2* ietaS*detaS_dsmag*coef0;
		fiKAB(1,0) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s1* ietaS*detaS_dsmag*coef0;
       
		/*K_epA_sy*/
		fiKAB(0,3) = -0.5*dt*ietaS*s0 *ietaS*detaS_dsy;
		fiKAB(1,3) = -0.5*dt*ietaS*s1 *ietaS*detaS_dsy;
		fiKAB(2,3) = -0.5*dt*ietaS*s2 *ietaS*detaS_dsy;
	
		/*K_sy_epB*/
		fiKAB(3,0) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*coef0*(1.0 - smag*ietaS*detaS_dsmag);
		fiKAB(3,1) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*coef1*(1.0 - smag*ietaS*detaS_dsmag);
		fiKAB(3,2) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*coef2*(1.0 - smag*ietaS*detaS_dsmag);
		
		/*K_sy_sy*/
	    fiKAB(3,3) = 1.0 + 0.5*dt*ietaS*fh* (smag/fsinf + (1.0 - sy/fsinf)*smag*ietaS*detaS_dsy);

//		cout << "\nfiKAB: "<<fiKAB;
		/*inverts KAB*/
		fiKAB.Inverse();

		/*initialize*/
		fGAB = 0.0;
		
		/*G_epeA_epB*/
		fGAB(0,0) = 1.0 + 0.5*dt*ietaS*s0 - 0.5*dt*ietaS *ietaS*s0*s0;
		fGAB(0,1) = 0.5*dt*ietaS*s0 - 0.5*dt*ietaS *ietaS*s0*s0;
		fGAB(0,2) = 0.5*dt*ietaS*s0 - 0.5*dt*ietaS *ietaS*s0*s0;
		
		fGAB(1,0) = 0.5*dt*ietaS*s1 - 0.5*dt*ietaS *ietaS*s1*s1;
		fGAB(1,1) = 1.0 + 0.5*dt*ietaS*s1 - 0.5*dt*ietaS *ietaS*s1*s1;
		fGAB(1,2) = 0.5*dt*ietaS*s1 - 0.5*dt*ietaS *ietaS*s1*s1;

		fGAB(2,0) = 0.5*dt*ietaS*s2 - 0.5*dt*ietaS *ietaS*s2*s2;
		fGAB(2,1) = 0.5*dt*ietaS*s2 - 0.5*dt*ietaS *ietaS*s2*s2;
		fGAB(2,2) = 1.0 + 0.5*dt*ietaS*s2 - 0.5*dt*ietaS *ietaS*s2*s2;
		
		/*G_sy_epB*/
//		double coef = (s0*s0 + s1*s1 + s2*s2)/smag;
		fGAB(3,0) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*smag*(1.0 - smag*ietaS*detaS_dsmag);
		fGAB(3,1) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*smag*(1.0 - smag*ietaS*detaS_dsmag);
		fGAB(3,2) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*smag*(1.0 - smag*ietaS*detaS_dsmag);
		
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

void ModBoyceVisco::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
			     dArrayT& eigenstress, dSymMatrixT& eigenmodulus,  const int type) 
{		
	const double ctol = 1.00e-12;
		
	/*set references to principle stretches*/
     
	double& le0 = eigenstretch_e[0];
	double& le1 = eigenstretch_e[1];
	double& le2 = eigenstretch_e[2];
	
	double& sy = *fsy;
  
	double tol;

	/*initialize principle elastic and trial elastic log strains */
	const double ep_tr0 = 0.5*log(le0);
	const double ep_tr1 = 0.5*log(le1);
	const double ep_tr2 = 0.5*log(le2);
	const double syn = *fsy_n;
	
	double ep_e0 = ep_tr0;		
	double ep_e1 = ep_tr1;	
	double ep_e2 = ep_tr2;

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
		double etaS = ShearViscosity(smag, sy);
		double ietaS = 1.0/etaS;
						
		double gamdot = 0.5*smag*ietaS;

	    /*calculate the residual*/
	    fRes[0] = ep_e0 + 0.5*dt*ietaS*s0 - ep_tr0;
	    fRes[1] = ep_e1 + 0.5*dt*ietaS*s1 - ep_tr1;
	    fRes[2] = ep_e2 + 0.5*dt*ietaS*s2 - ep_tr2;
		
		fRes[3] = sy - dt*fh*(1.0-sy/fsinf)*gamdot - syn;
//		cout << "\nfRes: "<<fRes;
		/*calculate stiffness matrix*/
		/*derivative of retardation time wrt to Tf*/
		double detaS_dsmag = -etaS*fQS / (fTemperature*sy);
		double detaS_dsy = etaS*fQS*smag / (fTemperature*sy*sy);
//		cout << "\ndetaS_dsmag: "<<detaS_dsmag;
//		cout << "\ndetaS_dsy: "<<detaS_dsy;
		
		/*initialize*/
		fiKAB = 0.0;
		
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
		fiKAB(0,0) = 1.0 + 0.5*ietaS*dt*c0 - 0.5*dt*ietaS*s0* ietaS*detaS_dsmag*coef0;
		fiKAB(1,1) = 1.0 + 0.5*ietaS*dt*c1 - 0.5*dt*ietaS*s1* ietaS*detaS_dsmag*coef1;
		fiKAB(2,2) = 1.0 + 0.5*ietaS*dt*c2 - 0.5*dt*ietaS*s2* ietaS*detaS_dsmag*coef2;
		
		fiKAB(1,2) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s1* ietaS*detaS_dsmag*coef2;
		fiKAB(0,2) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s0* ietaS*detaS_dsmag*coef2;
		fiKAB(0,1) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s0* ietaS*detaS_dsmag*coef1;

		fiKAB(2,1) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s2* ietaS*detaS_dsmag*coef1;
		fiKAB(2,0) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s2* ietaS*detaS_dsmag*coef0;
		fiKAB(1,0) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s1* ietaS*detaS_dsmag*coef0;
       
		/*K_epA_sy*/
		fiKAB(0,3) = -0.5*dt*ietaS*s0 *ietaS*detaS_dsy;
		fiKAB(1,3) = -0.5*dt*ietaS*s1 *ietaS*detaS_dsy;
		fiKAB(2,3) = -0.5*dt*ietaS*s2 *ietaS*detaS_dsy;
	
		/*K_sy_epB*/
		fiKAB(3,0) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*coef0*(1.0 - smag*ietaS*detaS_dsmag);
		fiKAB(3,1) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*coef1*(1.0 - smag*ietaS*detaS_dsmag);
		fiKAB(3,2) = -0.5*dt*ietaS*fh*(1.0 - sy/fsinf)*coef2*(1.0 - smag*ietaS*detaS_dsmag);
		
		/*K_sy_sy*/
	    fiKAB(3,3) = 1.0 + 0.5*dt*ietaS*fh* (smag/fsinf + (1.0 - sy/fsinf)*smag*ietaS*detaS_dsy);

//		cout << "\nfiKAB: "<<fiKAB;
		/*inverts KAB*/
		fiKAB.Inverse();
	    
		
	    /*solve for the principal strain increments*/
		fiKAB.Multx(fRes, fDelta, -1.0);
				
	    ep_e0 += fDelta[0];
	    ep_e1 += fDelta[1];
	    ep_e2 += fDelta[2];
	    
		sy += fDelta[3];
		
	    le0 = exp(2.0*ep_e0);
	    le1 = exp(2.0*ep_e1);
	    le2 = exp(2.0*ep_e2);
	    
	    /*Check that the L2 norm of the residual is less than tolerance*/
	    tol = sqrt(dArrayT::Dot(fRes, fRes));
/*		cout << "\ntol: "<<tol;
		cout << "\neps_e: "<<eigenstretch_e;
*/
	}while (tol>ctol && iteration < maxiteration); 
	if (iteration >= maxiteration) 
		ExceptionT::GeneralFail("ModBoyceVisco::ComputeEigs_e", 
			"number of iteration exceeds maximum");
}
