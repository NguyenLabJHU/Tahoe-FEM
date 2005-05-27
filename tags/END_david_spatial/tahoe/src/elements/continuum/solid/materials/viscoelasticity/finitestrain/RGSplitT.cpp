/* $Id: RGSplitT.cpp,v 1.5 2004-12-01 17:50:21 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "RGSplitT.h"
#include "PotentialT.h"
#include "NeoHookean.h"

#include "ifstreamT.h"
#include "ExceptionT.h"
#include <math.h>
#include <iostream.h>
#include <stdlib.h>

using namespace Tahoe;

const double third = 1.0/3.0; 
const int kNumOutputVar =1; 
static const char* Labels[kNumOutputVar] = {"Dvisc"}; 

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
/* constructors */
RGSplitT::RGSplitT(void):
  ParameterInterfaceT("Reese-Govindjee_split"),
  fSpectralDecompSpat(3),
  fPot_EQ(NULL),
  fPot_NEQ(NULL)
{

}

RGSplitT::~RGSplitT(void)
{
    delete fPot_EQ;
    delete fPot_NEQ;
}

int RGSplitT::NumOutputVariables() const {return kNumOutputVar;} 

void RGSplitT::OutputLabels(ArrayT<StringT>& labels) const 
{ 
     /*allocates space for labels*/
     labels.Dimension(kNumOutputVar); 
  
     /*copy labels*/
     for (int i = 0; i< kNumOutputVar; i++) 
       labels[i] = Labels[i]; 
} 

double RGSplitT::StrainEnergyDensity(void)
{
	/*calculates equilibrium part*/
	Compute_b(fStretch);
	if (NumSD() == 2)
	{
		fb3D[0] = fStretch[0];
		fb3D[1] = fStretch[1];
		fb3D[2] = 1.0;
       
		fb3D[3] = 0.0;
		fb3D[4] = 0.0;
		fb3D[5] = fStretch[2];
	}
	else fb3D = fStretch;

	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb3D, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	
	double J = sqrt(fEigs.Product());
	fEigs_dev = fEigs;
	fEigs_dev *= pow(J, -2.0*third);
     
	double energy = 0.0;
	energy = fPot_EQ->Energy(fEigs_dev, J);
  
	/*adds nonequilibrium part */
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
  
	/*calculate be*/
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
	fInverse = fC_v;
	fInverse.Inverse();
	fbe.MultQBQT(fF3D,fInverse);

	fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
	fEigs_e = fSpectralDecompSpat.Eigenvalues();
	
	double Je = sqrt(fEigs_e.Product());
	fEigs_dev = fEigs_e;
	fEigs_dev *= pow(Je,-2.0*third);
  
	energy += fPot_NEQ->Energy(fEigs_dev, Je);
  
	return(energy);
}

/* modulus */
const dMatrixT& RGSplitT::c_ijkl(void)
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
	
    Compute_b(fStretch);
    if (NumSD() == 2)
    {
      fb3D[0] = fStretch[0];
      fb3D[1] = fStretch[1];
      fb3D[2] = 1.0;
    
      fb3D[3] = 0.0;
      fb3D[4] = 0.0;
      fb3D[5] = fStretch[2];
    }
    else fb3D = fStretch;
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb3D, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();
    const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();

    double J = sqrt(fEigs.Product());
    fEigs_dev = fEigs;
    fEigs_dev *= pow(J, -2.0*third);
	
    fPot_EQ->DevStress(fEigs_dev, ftau_EQ);
    ftau_EQ += fPot_EQ->MeanStress(J);    
    fPot_EQ->DevMod(fEigs_dev, fDtauDe_EQ);
    fDtauDe_EQ += fPot_EQ->MeanMod(J);

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

	fInverse = fC_vn;
	fInverse.Inverse();
	fb_tr.MultQBQT(fF3D, fInverse);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_tr, false);	
	fEigs_tr = fSpectralDecompSpat.Eigenvalues(); 
		
	fInverse = fC_v;
	fInverse.Inverse();
	fbe.MultQBQT(fF3D, fInverse);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
	fEigs_e = fSpectralDecompSpat.Eigenvalues(); 
    const ArrayT<dArrayT>& eigenvectors_e=fSpectralDecompSpat.Eigenvectors();

	double Je = sqrt(fEigs_e.Product());
    fEigs_dev = fEigs_e;
    fEigs_dev *= pow(Je,-2.0*third);
    
    fPot_NEQ->DevStress(fEigs_dev, ftau_NEQ);
    ftau_NEQ += fPot_NEQ->MeanStress(Je);    
    fPot_NEQ->DevMod(fEigs_dev, fDtauDe_NEQ);
    double cm = fPot_NEQ->MeanMod(Je);
	    
	ComputeiKAB(fDtauDe_NEQ, cm);
    dSymMatrixT& DAB = fDtauDe_NEQ;
    DAB += cm; 
	
    fCalg(0,0) = DAB(0,0)*fiKAB(0,0) + DAB(0,1)*fiKAB(1,0) + DAB(0,2)*fiKAB(2,0) - 2.0*ftau_NEQ[0];
    fCalg(1,0) = DAB(1,0)*fiKAB(0,0) + DAB(1,1)*fiKAB(1,0) + DAB(1,2)*fiKAB(2,0);
    fCalg(2,0) = DAB(2,0)*fiKAB(0,0) + DAB(2,1)*fiKAB(1,0) + DAB(2,2)*fiKAB(2,0);
    fCalg(0,1) = DAB(0,0)*fiKAB(0,1) + DAB(0,1)*fiKAB(1,1) + DAB(0,2)*fiKAB(2,1);
    fCalg(1,1) = DAB(1,0)*fiKAB(0,1) + DAB(1,1)*fiKAB(1,1) + DAB(1,2)*fiKAB(2,1) - 2.0*ftau_NEQ[1];
    fCalg(2,1) = DAB(2,0)*fiKAB(0,1) + DAB(2,1)*fiKAB(1,1) + DAB(2,2)*fiKAB(2,1);
    fCalg(0,2) = DAB(0,0)*fiKAB(0,2) + DAB(0,1)*fiKAB(1,2) + DAB(0,2)*fiKAB(2,2);
    fCalg(1,2) = DAB(1,0)*fiKAB(0,2) + DAB(1,1)*fiKAB(1,2) + DAB(1,2)*fiKAB(2,2);
    fCalg(2,2) = DAB(2,0)*fiKAB(0,2) + DAB(2,1)*fiKAB(1,2) + DAB(2,2)*fiKAB(2,2) - 2.0*ftau_NEQ[2];
	   
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

    fModulus *= 1.0/J;

    return fModulus;
}

/* stresses */
const dSymMatrixT& RGSplitT::s_ij(void)
{
    /* stretch tensor */
    Compute_b(fStretch);
    if (NumSD() == 2)
    {
      fb3D[0] = fStretch[0];
      fb3D[1] = fStretch[1];
      fb3D[2] = 1.0;
    
      fb3D[3] = 0.0;
      fb3D[4] = 0.0;
      fb3D[5] = fStretch[2];
    }
    else fb3D = fStretch;
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb3D, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();

    /*jacobian determinant*/
    double J = sqrt(fEigs.Product());
    fEigs_dev = fEigs;
    fEigs_dev *= pow(J,-2.0*third);
    
    fPot_EQ->DevStress(fEigs_dev, ftau_EQ);
	ftau_EQ += fPot_EQ->MeanStress(J);
    
	fStress3D = fSpectralDecompSpat.EigsToRank2(ftau_EQ);
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

    /*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
    {		
		/*calc trial state*/
		fInverse = fC_vn;
		fInverse.Inverse();
		fb_tr.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_tr, false);	
		fEigs_tr = fSpectralDecompSpat.Eigenvalues(); 

		/*calc elastic stretch*/
		fEigs_e = fEigs_tr; /*initial condition*/
		ComputeEigs_e(fEigs, fEigs_e, ftau_NEQ, fDtauDe_NEQ);

		double Je = sqrt(fEigs_e.Product());
		fEigs_dev = fEigs_e;
		fEigs_dev *= pow(Je,-2.0*third);
	
		fPot_NEQ->DevStress(fEigs_dev, ftau_NEQ);
		ftau_NEQ += fPot_NEQ->MeanStress(Je);
		fStress3D += fSpectralDecompSpat.EigsToRank2(ftau_NEQ);

	
		/*Calculate Cv*/
		fInverse = fSpectralDecompSpat.EigsToRank2(fEigs_e); /*be which is colinear with btr*/
		fInverse.Inverse();
		fC_v.MultQTBQ(fF3D, fInverse); 
		Store(element, CurrIP());
	}	
    else 
    {
		/*calc elastic stretch*/
		fInverse = fC_v;
		fInverse.Inverse();
		fbe.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
		fEigs_e = fSpectralDecompSpat.Eigenvalues(); 
		
		double Je = sqrt(fEigs_e.Product());
		fEigs_dev = fEigs_e;
		fEigs_dev *= pow(Je,-2.0*third);
		
		fPot_NEQ->DevStress(fEigs_dev, ftau_NEQ);
		ftau_NEQ += fPot_NEQ->MeanStress(Je);
		fStress3D += fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
    }

    
    if (NumSD() == 2)
    {
        fStress[0] = fStress3D[0];
        fStress[1] = fStress3D[1];
        fStress[2] = fStress3D[5];
    }
    else fStress = fStress3D;
    fStress *= 1.0/J;
	return fStress;
}

/* material description */
const dMatrixT& RGSplitT::C_IJKL(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F();
  
    /* transform */
    fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl()));
    return fModulus;	
}

const dSymMatrixT& RGSplitT::S_IJ(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F();
  
    /* transform */
    fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij()));
    return fStress;
}

void RGSplitT::ComputeOutput(dArrayT& output)
{
    /* spectral decomposition */
    Compute_b(fStretch);
    if (NumSD() == 2)
    {
      fb3D[0] = fStretch[0];
      fb3D[1] = fStretch[1];
      fb3D[2] = 1.0;
    
      fb3D[3] = 0.0;
      fb3D[4] = 0.0;
      fb3D[5] = fStretch[2];
    }
    else fb3D = fStretch;
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb3D, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();

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
	
	/*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
   
	/*calc elastic stretch*/
	fInverse = fC_v;
	fInverse.Inverse();
	fbe.MultQBQT(fF3D, fInverse);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
	fEigs_e = fSpectralDecompSpat.Eigenvalues(); 

    /*calc jacobian*/
    double Je = sqrt(fEigs_e.Product()) ;
    fEigs_dev = fEigs_e;
    fEigs_dev *= pow(Je,-2.0*third);
    
    fPot_NEQ->DevStress(fEigs_dev, ftau_NEQ);
    fStress3D = fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
    double sm = fPot_NEQ->MeanStress(Je);
    
    output[0] = 0.5*(0.5*fietaS*fStress3D.ScalarProduct()+fietaB*sm*sm);
    
}

/***********************************************************************
 * Private
 ***********************************************************************/

void RGSplitT::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
			     dArrayT& eigenstress, dSymMatrixT& eigenmodulus) 
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
	    fPot_NEQ->DevStress(fEigs_dev, eigenstress);
	    
	    double& s0 = eigenstress[0];
	    double& s1 = eigenstress[1];
	    double& s2 = eigenstress[2];
	    
	    fPot_NEQ->DevMod(fEigs_dev,eigenmodulus);
	    
	    /*caculate means*/
	    double sm = fPot_NEQ->MeanStress(Je);
	    double cm = fPot_NEQ->MeanMod(Je);
	    
	    ComputeiKAB(eigenmodulus,cm);
	    
	    /*calculate the residual*/
	    double dt = fFSMatSupport->TimeStep();
	    double res0 = ep_e0 + dt*(0.5*fietaS*s0 +
			  third*fietaB*sm) - ep_tr0;
	    double res1 = ep_e1 + dt*(0.5*fietaS*s1 +
			  third*fietaB*sm) - ep_tr1;
	    double res2 = ep_e2 + dt*(0.5*fietaS*s2 +
			  third*fietaB*sm) - ep_tr2;
		
	    //cout << "\n residual1 "<< res0;
	    /*solve for the principal strain increments*/
	    double dep_e0=-fiKAB(0,0)*res0-fiKAB(0,1)*res1-fiKAB(0,2)*res2;
	    double dep_e1=-fiKAB(1,0)*res0-fiKAB(1,1)*res1-fiKAB(1,2)*res2;
	    double dep_e2=-fiKAB(2,0)*res0-fiKAB(2,1)*res1-fiKAB(2,2)*res2;
	    
	    /*updates principal elastic stretches*/ 
	    ep_e0 += dep_e0;
	    ep_e1 += dep_e1;
	    ep_e2 += dep_e2;
	    //	cout << "\n depsilon1 "<< dep_e0;
	    
	    le0 = exp(2.0*ep_e0);
	    le1 = exp(2.0*ep_e1);
	    le2 = exp(2.0*ep_e2);
	    
	    /*Check that the L2 norm of the residual is less than tolerance*/
	    tol = sqrt(res0*res0 + res1*res1+res2*res2);
	}while (tol>ctol); 
}

void RGSplitT::ComputeiKAB(dSymMatrixT& eigenmodulus, double& bulkmodulus)
{	
        /*inverse viscosities*/

	/*deviatoric values*/
	double& c0 = eigenmodulus(0,0);
	double& c1 = eigenmodulus(1,1);
	double& c2 = eigenmodulus(2,2);

	double& c12 = eigenmodulus(1,2);
	double& c02 = eigenmodulus(0,2);
	double& c01 = eigenmodulus(0,1);
	
	/*mean value*/
	double& cm = bulkmodulus;

	dMatrixT& KAB = fiKAB;
		
	/*calculates  KAB = 1+dt*D(dWdE_Idev/nD+isostress/nV)/Dep_e*/

	double dt = fFSMatSupport->TimeStep();
	KAB(0,0) = 1+0.5*fietaS*dt*c0+third*fietaB*dt*cm;
	KAB(1,1) = 1+0.5*fietaS*dt*c1+third*fietaB*dt*cm;
	KAB(2,2) = 1+0.5*fietaS*dt*c2+third*fietaB*dt*cm;

	KAB(1,2) = 0.5*fietaS*dt*c12+third*fietaB*dt*cm;
	KAB(0,2) = 0.5*fietaS*dt*c02+third*fietaB*dt*cm;
	KAB(0,1) = 0.5*fietaS*dt*c01+third*fietaB*dt*cm;
       
	KAB(2,1) = KAB(1,2);
	KAB(2,0) = KAB(0,2);
	KAB(1,0) = KAB(0,1);
	
	/*inverts KAB*/
	fiKAB.Inverse(KAB);
}

/* describe the parameters needed by the interface */
void RGSplitT::DefineParameters(ParameterListT& list) const
{
  /* inherited */
  RGViscoelasticityT::DefineParameters(list);

  /* common limit */
  LimitT positive(0.0, LimitT::Lower);

  /* viscosities */
  ParameterT eta_shear(ParameterT::Double, "eta_shear");
  ParameterT eta_bulk(ParameterT::Double, "eta_bulk");
  eta_shear.AddLimit(positive);
  eta_bulk.AddLimit(positive);
  list.AddParameter(eta_shear);
  list.AddParameter(eta_bulk);

  /* potentials - could make this a choice but just neo-Hookean for now */
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

void RGSplitT::TakeParameterList(const ParameterListT& list)
{
  /* inherited */
  RGViscoelasticityT::TakeParameterList(list);

  /* dimension work space */
  fStretch.Dimension(NumSD());
  fInverse.Dimension(3);
  fb3D.Dimension(3);
  fbe.Dimension(3);
  fb_tr.Dimension(3);
  fF3D.Dimension(3);
  fEigs_dev.Dimension(3);
  fEigs.Dimension(3);
  fEigs_e.Dimension(3);
  fEigs_v.Dimension(3);
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

  /* viscosities */
  double etaS = list.GetParameter("eta_shear");
  double etaB = list.GetParameter("eta_bulk");
  fietaS = 1.0/etaS;
  fietaB = 1.0/etaB;

  /* potentials - could make this a choice but just neo-Hookean for now */
  double mu_eq = list.GetParameter("mu_EQ");
  double kappa_eq = list.GetParameter("kappa_EQ");
  NeoHookean* pot_eq = new NeoHookean;
  pot_eq->SetKappaMu(kappa_eq, mu_eq);
  fPot_EQ = pot_eq;

  double mu_neq = list.GetParameter("mu_NEQ");
  double kappa_neq = list.GetParameter("kappa_NEQ");
  NeoHookean* pot_neq = new NeoHookean;
  pot_neq->SetKappaMu(kappa_neq, mu_neq);
  fPot_NEQ = pot_neq;
}


