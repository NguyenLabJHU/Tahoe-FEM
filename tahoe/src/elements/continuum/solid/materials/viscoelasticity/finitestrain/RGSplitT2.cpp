/* $Id: RGSplitT2.cpp,v 1.1 2007-04-09 22:27:06 tdnguye Exp $ */
/* created: TDN (01/22/2001) */

#include "RGSplitT2.h"
#include "ParameterContainerT.h"

#include "ifstreamT.h"
#include "ExceptionT.h"
#include <math.h>
#include <iostream.h>
#include <stdlib.h>

#include "MooneyRivlin.h"
#include "NeoHookean.h"
#include "VWPotentialT.h"

#include "LinearExponentialT.h"
#include "ScaledCsch.h"

using namespace Tahoe;

const double third = 1.0/3.0; 
const int kNumOutputVar =1; 
static const char* Labels[kNumOutputVar] = {"Dvisc"}; 

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
RGSplitT2::RGSplitT2(void):
  ParameterInterfaceT("RG_split_general"),
  fSpectralDecompSpat(3)
{}

int RGSplitT2::NumOutputVariables() const {return kNumOutputVar;} 

void RGSplitT2::OutputLabels(ArrayT<StringT>& labels) const 
{ 
     /*allocates space for labels*/
     labels.Dimension(kNumOutputVar); 
  
     /*copy labels*/
     for (int i = 0; i< kNumOutputVar; i++) 
       labels[i] = Labels[i]; 
} 

double RGSplitT2::StrainEnergyDensity(void)
{
	/*calculates equilibrium part*/
	const dMatrixT& F = F_mechanical();
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

	/*elastic stretch*/
	fb.MultAAT(fF3D);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	
	double J = sqrt(fEigs.Product());
	fEigs_dev = fEigs;
	fEigs_dev *= pow(J, -2.0*third);
     
	double energy = 0.0;
	energy = fPot[0]->Energy(fEigs_dev, J);
  
	/*adds nonequilibrium part */
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
  
	for (int i = 0; i < fNumProcess; i++)
	{
		/*calculate be*/
		fInverse.Inverse(fC_v[i]);
		fbe.MultQBQT(fF3D,fInverse);

		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
		fEigs_e = fSpectralDecompSpat.Eigenvalues();
	
		double Je = sqrt(fEigs_e.Product());
		fEigs_dev = fEigs_e;
		fEigs_dev *= pow(Je,-2.0*third);
  
		energy += fPot[i+1]->Energy(fEigs_dev, Je);
	}
	return(energy);
}

/* modulus */
const dMatrixT& RGSplitT2::c_ijkl(void)
{
    
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    
	const dMatrixT& F = F_mechanical();
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

	/*calcualte total stretch*/
    fb.MultAAT(fF3D);
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();
    const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();

	/*calc EQ component of stress and moduli*/
    double J = sqrt(fEigs.Product());
    fEigs_dev = fEigs;
    fEigs_dev *= pow(J, -2.0*third);
	
    fPot[0]->DevStress(fEigs_dev, ftau_EQ);
    ftau_EQ += fPot[0]->MeanStress(J);    
    fPot[0]->DevMod(fEigs_dev,fDtauDe_EQ);
    fDtauDe_EQ += fPot[0]->MeanMod(J);

    dSymMatrixT& Gamma = fDtauDe_EQ;
    Gamma(0,0) -= 2.0*ftau_EQ[0];
    Gamma(1,1) -= 2.0*ftau_EQ[1];
    Gamma(2,2) -= 2.0*ftau_EQ[2];
   
	fModulus3D = fSpectralDecompSpat.EigsToRank4(Gamma);	
  
	double dl, coeff;

    double& l0 = fEigs[0];
    double& l1 = fEigs[1];
    double& l2 = fEigs[2];
	
	dl = l0 - l1;
    if (fabs(dl) > kSmall)
		coeff = (ftau_EQ[0]*l1 - ftau_EQ[1]*l0)/dl;
    else 
		coeff = 0.5*(Gamma(0,0)-Gamma(0,1))-ftau_EQ[0];
    MixedRank4_3D(eigenvectors[0], eigenvectors[1], fModMat);
    fModulus3D.AddScaled(2.0*coeff, fModMat);
    
    dl = l0 - l2;
    if (fabs(dl) > kSmall)
      coeff = (ftau_EQ[0]*l2 - ftau_EQ[2]*l0)/dl;
    else 
      coeff = 0.5*(Gamma(0,0)-Gamma(0,2))-ftau_EQ[2];	
    MixedRank4_3D(eigenvectors[0], eigenvectors[2], fModMat);
    fModulus3D.AddScaled(2.0*coeff, fModMat);
    
    dl = l1 - l2;
   if (fabs(dl) > kSmall)
		coeff  = (ftau_EQ[1]*l2 - ftau_EQ[2]*l1)/dl;
    else
      coeff = 0.5*(Gamma(1,1)-Gamma(1,2))-ftau_EQ[1];	
    MixedRank4_3D(eigenvectors[1], eigenvectors[2], fModMat);
    fModulus3D.AddScaled(2.0*coeff, fModMat);

	/*calc NEQ component of stress and moduli*/
	/*calcualte principal values of elastic stretch*/
	for (int i = 0; i < fNumProcess; i++)
	{
		fInverse.Inverse(fC_vn[i]);
		fb_tr.MultQBQT(fF3D, fInverse);

		fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_tr, false);	
		fEigs_tr = fSpectralDecompSpat.Eigenvalues(); 		

		fInverse.Inverse(fC_v[i]);
		fbe.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
		fEigs_e = fSpectralDecompSpat.Eigenvalues(); 
		const ArrayT<dArrayT>& eigenvectors_e=fSpectralDecompSpat.Eigenvectors();
		
		double Je = sqrt(fEigs_e.Product());
		fEigs_dev = fEigs_e;
		fEigs_dev *= pow(Je,-2.0*third);
    
		/*stresses*/
		fPot[i+1]->DevStress(fEigs_dev, ftau_NEQ);
		double sm =  fPot[i+1]->MeanStress(Je);    

		fPot[i+1]->DevMod(fEigs_dev, fDtauDe_NEQ);
		double cm = fPot[i+1]->MeanMod(Je);
		
		/*Calculate Calg_AB*/
		Compute_Calg(ftau_NEQ, fDtauDe_NEQ, sm, cm, fCalg, i);

		ftau_NEQ += sm;
		fDtauDe_NEQ += cm;
			   
		fModulus3D += fSpectralDecompSpat.NonSymEigsToRank4(fCalg);
    
		double dl_tr;

		double& l0_tr = fEigs_tr[0];
		double& l1_tr = fEigs_tr[1];
		double& l2_tr = fEigs_tr[2];
	
	
		dl_tr = l0_tr - l1_tr;
		if (fabs(dl_tr) > kSmall)
			coeff = (ftau_NEQ[0]*l1_tr - ftau_NEQ[1]*l0_tr)/dl_tr;
		else 
			coeff = 0.5*(fCalg(0,0)-fCalg(0,1))-ftau_NEQ[0];
		MixedRank4_3D(eigenvectors_e[0], eigenvectors_e[1], fModMat);
		fModulus3D.AddScaled(2.0*coeff, fModMat);
    
		dl_tr = l0_tr - l2_tr;
		if (fabs(dl_tr) > kSmall)
			coeff =(ftau_NEQ[0]*l2_tr - ftau_NEQ[2]*l0_tr)/dl_tr;
		else 
			coeff = 0.5*(fCalg(0,0)-fCalg(0,2))-ftau_NEQ[2];	
		MixedRank4_3D(eigenvectors_e[0], eigenvectors_e[2], fModMat);
		fModulus3D.AddScaled(2.0*coeff, fModMat);
    
		dl_tr = l1_tr - l2_tr;
		if (fabs(dl_tr) > kSmall)
			coeff  = (ftau_NEQ[1]*l2_tr - ftau_NEQ[2]*l1_tr)/dl_tr;
		else
			coeff = 0.5*(fCalg(1,1)-fCalg(1,2))-ftau_NEQ[1];	
		MixedRank4_3D(eigenvectors_e[1], eigenvectors_e[2], fModMat);
		fModulus3D.AddScaled(2.0*coeff, fModMat);
    }
	if (NumSD() == 2)
	{
		fModulus[0] = fModulus3D[0];
		fModulus[1] = fModulus3D[1];
		fModulus[2] = fModulus3D[5];

		fModulus[3] = fModulus3D[6];
		fModulus[4] = fModulus3D[7];
		fModulus[5] = fModulus3D[11];
		fModulus[6] = fModulus3D[30];
		fModulus[7] = fModulus3D[31];
		fModulus[8] = fModulus3D[35];
	}
	else fModulus = fModulus3D;

	const dMatrixT& Ftotal = F_total();	
	fModulus *= 1.0/Ftotal.Det();

    return fModulus;
}

/* stresses */
const dSymMatrixT& RGSplitT2::s_ij(void)
{
	const dMatrixT& F = F_mechanical();
	const dMatrixT& F_thermal = F_thermal_inverse();
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

	/*calculate EQ part of the stress*/
	fb.MultAAT(fF3D);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();

	/*jacobian determinant*/
	double J = sqrt(fEigs.Product());
	fEigs_dev = fEigs;
	fEigs_dev *= pow(J,-2.0*third);
    
	fPot[0]->DevStress(fEigs_dev, ftau_EQ);
	ftau_EQ += fPot[0]->MeanStress(J);
/*		const double mu_eq = fPot[0]->GetMu();
		const double kappa_eq = fPot[0]->GetKappa();
		cout << "\neq mu: "<< mu_eq;
		cout << "\neq kappa: "<< kappa_eq;
*/	
	fStress3D = fSpectralDecompSpat.EigsToRank2(ftau_EQ);

    /*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
    {		
		/*calc NEQ component of stress and moduli*/
		for (int i = 0; i < fNumProcess; i++)
		{
			/*calc trial state*/
			fInverse.Inverse(fC_vn[i]);
			fb_tr.MultQBQT(fF3D, fInverse);
			
			fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_tr, false);	
			fEigs_tr = fSpectralDecompSpat.Eigenvalues(); 

			/*calc elastic stretch*/
			fEigs_e = fEigs_tr; /*initial condition*/
			ComputeEigs_e(fEigs, fEigs_e, ftau_NEQ, fDtauDe_NEQ, i);

			double Je = sqrt(fEigs_e.Product());
			fEigs_dev = fEigs_e;
			fEigs_dev *= pow(Je,-2.0*third);
	
			fPot[i+1]->DevStress(fEigs_dev, ftau_NEQ);
			ftau_NEQ += fPot[i+1]->MeanStress(Je);
			fStress3D += fSpectralDecompSpat.EigsToRank2(ftau_NEQ);

/*		const double mu_neq = fPot[i+1]->GetMu();
		const double kappa_neq = fPot[i+1]->GetKappa();
		cout << "\nneq mu: "<< mu_neq;
		cout << "\nneq kappa: "<< kappa_neq;
*/
	
			/*Calculate Cv*/
			fInverse = fSpectralDecompSpat.EigsToRank2(fEigs_e); /*be which is colinear with btr*/
			fInverse.Inverse();
			fC_v[i].MultQTBQ(fF3D, fInverse); 
		}
		Store(element, CurrIP());
	}	
    else 
    {
		/*calc NEQ component of stress and moduli*/
		for (int i = 0; i < fNumProcess; i++)
		{
			/*calc elastic stretch*/
			fInverse.Inverse(fC_v[i]);
			fbe.MultQBQT(fF3D, fInverse);
			fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
			fEigs_e = fSpectralDecompSpat.Eigenvalues(); 

	//		fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	

			double Je = sqrt(fEigs_e.Product());
			fEigs_dev = fEigs_e;
			fEigs_dev *= pow(Je,-2.0*third);
		
			fPot[i+1]->DevStress(fEigs_dev, ftau_NEQ);
			ftau_NEQ += fPot[i+1]->MeanStress(Je);
			fStress3D += fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
		}
    }

	if (NumSD() == 2)
    {
        fStress[0] = fStress3D[0];
        fStress[1] = fStress3D[1];
        fStress[2] = fStress3D[5];
    }
    else fStress = fStress3D;
	
	const dMatrixT& Ftotal = F_total();	
    fStress *= 1.0/Ftotal.Det();
	return fStress;
}

/* material description */
const dMatrixT& RGSplitT2::C_IJKL(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F_total();
  
    /* transform */
    fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl()));
    return fModulus;	
}

const dSymMatrixT& RGSplitT2::S_IJ(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F_total();
  
    /* transform */
    fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij()));
    return fStress;
}

void RGSplitT2::ComputeOutput(dArrayT& output)
{
	const dMatrixT& F = F_mechanical();
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

    fb.MultAAT(fF3D);
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();

	output[0] = 0.0;
	/*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
   
    for (int i = 0; i < fNumProcess; i++)
	{
		/*calc elastic stretch*/
		fInverse.Inverse(fC_v[i]);
		fbe.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
		fEigs_e = fSpectralDecompSpat.Eigenvalues(); 

		/*calc jacobian*/
		double Je = sqrt(fEigs_e.Product()) ;
		fEigs_dev = fEigs_e;
		fEigs_dev *= pow(Je,-2.0*third);
    
		fPot[i+1]->DevStress(fEigs_dev, ftau_NEQ);
		fStress3D = fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
		double sm = fPot[i+1]->MeanStress(Je);

		double stress_mag = sqrt(ftau_NEQ[0]*ftau_NEQ[0] + ftau_NEQ[1]*ftau_NEQ[1] + ftau_NEQ[2]*ftau_NEQ[2]);
		fietaS = 1.0/fVisc_s[i]->Function(stress_mag);
		fietaB = 1.0/fVisc_b[i]->Function(sm);
		
		output[0] += 0.5*(0.5*fietaS*fStress3D.ScalarProduct()+fietaB*sm*sm);
    }
}

/***********************************************************************
 * Protected
 ***********************************************************************/
void RGSplitT2::Initialize(void)
{
 /* dimension work space */
  fInverse.Dimension(3);
  fb.Dimension(3);
  fbe.Dimension(3);
  fb_tr.Dimension(3);
  fF3D.Dimension(3);
  fEigs_dev.Dimension(3);
  fEigs.Dimension(3);
  fEigs_e.Dimension(3);
  fEigs_tr.Dimension(3);
  ftau_EQ.Dimension(3);
  ftau_NEQ.Dimension(3);
  fDtauDe_EQ.Dimension(3);
  fDtauDe_NEQ.Dimension(3);
  fCalg.Dimension(3);
  fModulus3D.Dimension(6);
  fModMat.Dimension(6);
  fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));
  fStress.Dimension(NumSD());
  fStress3D.Dimension(3);
  fiKAB.Dimension(3);
}

/***********************************************************************
 * Private
 ***********************************************************************/
 void RGSplitT2::Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
	const double& dtau_m, dMatrixT& Calg, const int type)
 {
		double s0 = tau_dev[0];
		double s1 = tau_dev[1];
		double s2 = tau_dev[2];
		
		double sm = tau_m;
		
		double c0 = dtau_dev(0,0);
		double c1 = dtau_dev(1,1);
		double c2 = dtau_dev(2,2);

		double c12 = dtau_dev(1,2);
		double c02 = dtau_dev(0,2);
		double c01 = dtau_dev(0,1);
	    
		double cm = dtau_m;
		
		/*calculates  KAB = 1+dt*D(dWdE_Idev/nD+isostress/nV)/Dep_e*/
		double dt = fFSMatSupport->TimeStep();
		
		double stress_mag = sqrt(s0*s0 + s1*s1 + s2*s2);
		fietaS = 1.0/fVisc_s[type]->Function(stress_mag);
		fietaB = 1.0/fVisc_b[type]->Function(tau_m);

		fiKAB(0,0) = 1+0.5*fietaS*dt*c0+third*fietaB*dt*cm;
		fiKAB(1,1) = 1+0.5*fietaS*dt*c1+third*fietaB*dt*cm;
		fiKAB(2,2) = 1+0.5*fietaS*dt*c2+third*fietaB*dt*cm;

		fiKAB(1,2) = 0.5*fietaS*dt*c12+third*fietaB*dt*cm;
		fiKAB(0,2) = 0.5*fietaS*dt*c02+third*fietaB*dt*cm;
		fiKAB(0,1) = 0.5*fietaS*dt*c01+third*fietaB*dt*cm;
       
		fiKAB(2,1) = fiKAB(1,2);
		fiKAB(2,0) = fiKAB(0,2);
		fiKAB(1,0) = fiKAB(0,1);
		
		if (stress_mag > kSmall)
		{
			double coeffs = fietaS*fVisc_s[type]->DFunction(stress_mag);
			double coeffb = fietaB*fVisc_b[type]->DFunction(sm);
			
			fiKAB(0,0) -= 0.5*fietaS*dt*coeffs*(s0*c0+s1*c01+s2*c02)/stress_mag*s0 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(1,1) -= 0.5*fietaS*dt*coeffs*(s0*c01+s1*c1+s2*c12)/stress_mag*s1 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(2,2) -= 0.5*fietaS*dt*coeffs*(s0*c02+s1*c12+s2*c2)/stress_mag*s2 
						- third*fietaB*dt*coeffb*(cm*sm);

			fiKAB(1,2) -= 0.5*fietaS*dt*coeffs*(s0*c02+s1*c12+s2*c2)/stress_mag*s1 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(0,2) -= 0.5*fietaS*dt*coeffs*(s0*c02+s1*c12+s2*c2)/stress_mag*s0 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(0,1) -= 0.5*fietaS*dt*coeffs*(s0*c01+s1*c1+s2*c12)/stress_mag*s0 
						- third*fietaB*dt*coeffb*(cm*sm);
						
			fiKAB(2,1) -= 0.5*fietaS*dt*coeffs*(s0*c01+s1*c1+s2*c12)/stress_mag*s2 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(2,0) -= 0.5*fietaS*dt*coeffs*(s0*c0+s1*c01+s2*c02)/stress_mag*s2 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(1,0) -= 0.5*fietaS*dt*coeffs*(s0*c0+s1*c01+s2*c02)/stress_mag*s1 
						- third*fietaB*dt*coeffb*(cm*sm);
		}
	
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

void RGSplitT2::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
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
	int iteration  = 0;	
	do 
	{
		iteration ++;
	    double Je=sqrt(le0*le1*le2);
	    fEigs_dev = eigenstretch_e;
	    fEigs_dev *= pow(Je,-2.0*third);
		
	    /*calculate stresses and moduli*/
	    fPot[type+1]->DevStress(fEigs_dev, eigenstress);
	    
	    double& s0 = eigenstress[0];
	    double& s1 = eigenstress[1];
	    double& s2 = eigenstress[2];
		
		double stress_mag = sqrt(s0*s0 + s1*s1 + s2*s2);
		fietaS = 1.0/fVisc_s[type]->Function(stress_mag);
		
	    fPot[type+1]->DevMod(fEigs_dev,eigenmodulus);
		/*deviatoric values*/
		double& c0 = eigenmodulus(0,0);
		double& c1 = eigenmodulus(1,1);
		double& c2 = eigenmodulus(2,2);

		double& c12 = eigenmodulus(1,2);
		double& c02 = eigenmodulus(0,2);
		double& c01 = eigenmodulus(0,1);
	    
	    /*caculate means*/
	    double sm = fPot[type+1]->MeanStress(Je);
		fietaB = 1.0/fVisc_b[type]->Function(sm);

	    double cm = fPot[type+1]->MeanMod(Je);
	    
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
	
		if (stress_mag > kSmall)
		{
			double coeffs = fietaS*fVisc_s[type]->DFunction(stress_mag);
			double coeffb = fietaB*fVisc_b[type]->DFunction(sm);
			
			fiKAB(0,0) -= 0.5*fietaS*dt*coeffs*(s0*c0+s1*c01+s2*c02)/stress_mag*s0 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(1,1) -= 0.5*fietaS*dt*coeffs*(s0*c01+s1*c1+s2*c12)/stress_mag*s1 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(2,2) -= 0.5*fietaS*dt*coeffs*(s0*c02+s1*c12+s2*c2)/stress_mag*s2 
						- third*fietaB*dt*coeffb*(cm*sm);

			fiKAB(1,2) -= 0.5*fietaS*dt*coeffs*(s0*c02+s1*c12+s2*c2)/stress_mag*s1 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(0,2) -= 0.5*fietaS*dt*coeffs*(s0*c02+s1*c12+s2*c2)/stress_mag*s0 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(0,1) -= 0.5*fietaS*dt*coeffs*(s0*c01+s1*c1+s2*c12)/stress_mag*s0 
						- third*fietaB*dt*coeffb*(cm*sm);
						
			fiKAB(2,1) -= 0.5*fietaS*dt*coeffs*(s0*c01+s1*c1+s2*c12)/stress_mag*s2 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(2,0) -= 0.5*fietaS*dt*coeffs*(s0*c0+s1*c01+s2*c02)/stress_mag*s2 
						- third*fietaB*dt*coeffb*(cm*sm);
			fiKAB(1,0) -= 0.5*fietaS*dt*coeffs*(s0*c0+s1*c01+s2*c02)/stress_mag*s1 
						- third*fietaB*dt*coeffb*(cm*sm);
		}
	
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
	}while (tol>ctol && iteration < 10); 
	if (iteration >= 10) 
		ExceptionT::GeneralFail("RGSplitT2::ComputeEigs_e", 
			"number of iteration exceeds maximum of 10");
}

/* information about subordinate parameter lists */
void RGSplitT2::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("rg_eq_potential", ParameterListT::Once);
	sub_list.AddSub("rg_neq_potential", ParameterListT::Any);

	/* choice of viscosity */
	sub_list.AddSub("rg_shear_viscosity", ParameterListT::Any);
	sub_list.AddSub("rg_bulk_viscosity", ParameterListT::Any);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* RGSplitT2::NewSub(const StringT& name) const
{
	PotentialT* pot = NULL;
	if (name == "neo-hookean")
		pot = new NeoHookean;
	else if (name == "mooney-rivlin")
		pot = new MooneyRivlin;
	else if (name == "veronda-westmann")
		pot = new VWPotentialT;
	if (pot)
		return pot;

	C1FunctionT* func = NULL;
	if (name == "scaled-csch")
		func = new ScaledCsch;
	else if (name == "linear_exponential")
		func = new LinearExponentialT;

	if (func)
		return func;

	if (name == "rg_eq_potential" || name == "rg_neq_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
	
		/* choice of parameters */
		choice->AddSub("neo-hookean");
		choice->AddSub("mooney-rivlin");
		choice->AddSub("veronda-westmann");
		return(choice);
	}
	else if (name == "rg_shear_viscosity")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetDescription("eta_S(|sig_dev|)");	
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		
		choice->AddSub("scaled-csch");
		choice->AddSub("linear_exponential");
		return(choice);
	}
	else if (name == "rg_bulk_viscosity")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetDescription("eta_B(sig_m)");	
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		
		choice->AddSub("scaled-csch");
		choice->AddSub("linear_exponential");
		return(choice);
	}
	else return(RGViscoelasticityT::NewSub(name));
}

/* accept parameter list */
void RGSplitT2::TakeParameterList(const ParameterListT& list)
{
	StringT caller = "RGSplitT2::TakeParameterList";
	int num_neq =  list.NumLists("rg_neq_potential");
	int num_shear_visc = list.NumLists("rg_shear_viscosity");
	int num_bulk_visc = list.NumLists("rg_bulk_viscosity");
	if (num_neq != num_shear_visc || num_neq != num_bulk_visc)
		ExceptionT::GeneralFail("RGSplitT2::TakeParameterList", 
			"number of matrix viscosity functions does not match number of matrix nonequilibrium potentials");
	fNumProcess = list.NumLists("rg_shear_viscosity");

	/*inherited*/
	RGViscoelasticityT::TakeParameterList(list);

	fPot.Dimension(fNumProcess+1);
	fVisc_s.Dimension(fNumProcess);
	fVisc_b.Dimension(fNumProcess);
		
	const ParameterListT& eq_pot = list.GetListChoice(*this, "rg_eq_potential");
	if(eq_pot.Name() == "neo-hookean")
		fPot[0] = new NeoHookean;
	else if(eq_pot.Name() == "mooney-rivlin")
		fPot[0] = new MooneyRivlin;
	else if(eq_pot.Name() == "veronda-westmann")
		fPot[0] = new VWPotentialT;
	else 
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot[0]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", eq_pot.Name().Pointer());			
	fPot[0]->TakeParameterList(eq_pot);
	

	for (int i = 0; i < fNumProcess; i++)
	{
		const ParameterListT& pot_neq = list.GetListChoice(*this, "rg_neq_potential",i);
		if(pot_neq.Name() == "mooney-rivlin")
			fPot[i+1] = new MooneyRivlin;
		else if(pot_neq.Name() == "neo-hookean")
			fPot[i+1] = new NeoHookean;
		else if(pot_neq.Name() == "veronda-westmann")
			fPot[i+1] = new VWPotentialT;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fPot[i+1]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", pot_neq.Name().Pointer());			
		fPot[i+1]->TakeParameterList(pot_neq);

		const ParameterListT& shear_visc = list.GetListChoice(*this, "rg_shear_viscosity", i);
		if (shear_visc.Name() == "linear_exponential")
			fVisc_s[i] = new LinearExponentialT;
		else if (shear_visc.Name() == "scaled-csch")
			fVisc_s[i] = new ScaledCsch;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fVisc_s[i]) throw ExceptionT::kOutOfMemory;
		fVisc_s[i]->TakeParameterList(shear_visc);

		const ParameterListT& bulk_visc = list.GetListChoice(*this, "rg_bulk_viscosity", i);
		if (bulk_visc.Name() == "linear_exponential")
			fVisc_b[i] = new LinearExponentialT;
		else if (bulk_visc.Name() == "scaled-csch")
			fVisc_b[i] = new ScaledCsch;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fVisc_b[i]) throw ExceptionT::kOutOfMemory;
		fVisc_b[i]->TakeParameterList(bulk_visc);
	}
	
	/*set dimension of workspaces*/
	Initialize();
}

