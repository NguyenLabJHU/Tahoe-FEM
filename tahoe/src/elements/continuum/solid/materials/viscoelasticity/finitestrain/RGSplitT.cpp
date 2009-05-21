/* $Id: RGSplitT.cpp,v 1.12 2009-05-21 22:30:27 tdnguye Exp $ */
/* created: TDN (01/22/2001) */

#include "RGSplitT.h"

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
RGSplitT::RGSplitT(void):
  ParameterInterfaceT("Reese-Govindjee_split"),
  fSpectralDecompSpat(3)
{
	fNumProcess = 1;
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
	energy = Energy(fEigs_dev, J, -1);
  
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
  
		energy += Energy(fEigs_dev, Je, i);
	}
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

	/*calcualte total stretch*/
    fb.MultAAT(fF3D);
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();
    const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();

	/*calc EQ component of stress and moduli*/
    double J = sqrt(fEigs.Product());
    fEigs_dev = fEigs;
    fEigs_dev *= pow(J, -2.0*third);
	
    DevStress(fEigs_dev, ftau_EQ, -1);
    ftau_EQ += MeanStress(J, -1);    
    DevMod(fEigs_dev,fDtauDe_EQ, -1);
    fDtauDe_EQ += MeanMod(J, -1);

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
		DevStress(fEigs_dev, ftau_NEQ, i);
		double sm =  MeanStress(Je, i);    

		DevMod(fEigs_dev, fDtauDe_NEQ, i);
		double cm = MeanMod(Je, i);
		
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
const dSymMatrixT& RGSplitT::s_ij(void)
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
    
	DevStress(fEigs_dev, ftau_EQ, -1);
	ftau_EQ += MeanStress(J, -1);
    
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
	
			DevStress(fEigs_dev, ftau_NEQ, i);
			ftau_NEQ += MeanStress(Je, i);
			fStress3D += fSpectralDecompSpat.EigsToRank2(ftau_NEQ);

	
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
		
			DevStress(fEigs_dev, ftau_NEQ, i);
			ftau_NEQ += MeanStress(Je, i);
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
const dMatrixT& RGSplitT::C_IJKL(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F_total();
  
    /* transform */
    fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl()));
    return fModulus;	
}

const dSymMatrixT& RGSplitT::S_IJ(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F_total();
  
    /* transform */
    fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij()));
    return fStress;
}

double RGSplitT::Energy(const dArrayT& lambda_bar, const double J, const int type)
{
  double mu, kappa;
  if (type == -1) {
	mu = fmu_eq;
	kappa = fkappa_eq;
  }
  else {
	mu = fmu_neq;
	kappa = fkappa_neq;
  }
  double I1 = lambda_bar[0]+lambda_bar[1]+lambda_bar[2];
  double phi = 0.5*mu*(I1-3)+0.25*kappa*(J*J-1-2*log(J));
  return(phi);
}
void RGSplitT::DevStress(const dArrayT& lambda_bar,dArrayT& tau, const int type)
{
  int nsd = tau.Length();
  double mu;
  (type == -1) ? mu = fmu_eq : mu = fmu_neq;
  
  const double& l0 = lambda_bar[0];
  const double& l1 = lambda_bar[1];
  const double& l2 = lambda_bar[2];
  
  tau[0] = mu*third*(2.0*l0-l1-l2);
  tau[1] = mu*third*(2.0*l1-l0-l2);
  
  if (nsd == 3)
    tau[2] = mu*third*(2.0*l2-l0-l1);
}

double RGSplitT::MeanStress(const double J, const int type) 
{
	double kappa;
	(type == -1) ? kappa = fkappa_eq : kappa = fkappa_neq;
	return(0.5*kappa*(J*J-1));
}

void RGSplitT::DevMod(const dArrayT& lambda_bar, dSymMatrixT& eigenmodulus, const int type)
{
  int nsd = eigenmodulus.Rows();
  double ninth = third*third;

  double mu;
  (type == -1) ? mu = fmu_eq : mu = fmu_neq;
  
  const double& l0 = lambda_bar[0];
  const double& l1 = lambda_bar[1];
  const double& l2 = lambda_bar[2];
  
  eigenmodulus[0] = 2.0*mu*ninth*(4.0*l0+l1+l2);
  eigenmodulus[1] = 2.0*mu*ninth*(4.0*l1+l2+l0);
  if (nsd == 2)
  {
    eigenmodulus[2] = 2.0*mu*ninth*(-2.0*l0-2.0*l1+l2); 
  }
  else 
  {
    eigenmodulus[2] = 2.0*mu*ninth*(4.0*l2+l0+l1);	
    eigenmodulus[3] = 2.0*mu*ninth*(-2.0*l1-2.0*l2+l0);
    eigenmodulus[4] = 2.0*mu*ninth*(-2.0*l0-2.0*l2+l1);
    eigenmodulus[5] = 2.0*mu*ninth*(-2.0*l0-2.0*l1+l2);
  }
}

double RGSplitT::MeanMod(const double J, const int type) 
{
	double kappa;
	(type == -1) ? kappa = fkappa_eq : kappa = fkappa_neq;
	return(kappa*J*J);
}


void RGSplitT::ComputeOutput(dArrayT& output)
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
    
		DevStress(fEigs_dev, ftau_NEQ, i);
		fStress3D = fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
		double sm = MeanStress(Je, i);
    
		output[0] += 0.5*(0.5*fietaS*fStress3D.ScalarProduct()+fietaB*sm*sm);
    }
}

/***********************************************************************
 * Protected
 ***********************************************************************/
void RGSplitT::Initialize(void)
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
 void RGSplitT::Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
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

/*		dSymMatrixT& DAB = fDtauDe_NEQ;
		DAB += cm; 
*/	
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

void RGSplitT::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
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
	}while (tol>ctol && iteration < 10); 
	if (iteration >= 10) 
		ExceptionT::GeneralFail("RGSplitT::ComputeEigs_e", 
			"number of iteration exceeds maximum of 10");
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
	cout << "\n RGSpliT is no longer supported.  Use the more general implementation RGSplitT2";	

  /* inherited */
  /*allows one neq process: */
  fNumProcess = 1;
  RGViscoelasticityT::TakeParameterList(list);

  /* viscosities */
  double etaS = list.GetParameter("eta_shear");
  double etaB = list.GetParameter("eta_bulk");
  fietaS = 1.0/etaS;
  fietaB = 1.0/etaB;

  /* moduli for neo-hookean potentials */
  fmu_eq = list.GetParameter("mu_EQ");
  fkappa_eq = list.GetParameter("kappa_EQ");

  fmu_neq = list.GetParameter("mu_NEQ");
  fkappa_neq = list.GetParameter("kappa_NEQ");
  
  /*Dimension work space*/
  Initialize();
}


