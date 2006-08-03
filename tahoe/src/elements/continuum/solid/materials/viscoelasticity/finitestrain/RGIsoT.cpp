/* $Id: RGIsoT.cpp,v 1.1 2006-08-03 23:17:28 tdnguye Exp $ */
/* created: TDN (01/22/2001) */

#include "RGIsoT.h"

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
RGIsoT::RGIsoT(void):
  ParameterInterfaceT("Reese-Govindjee_iso"),
  fSpectralDecompSpat(3)
{
}

int RGIsoT::NumOutputVariables() const {return kNumOutputVar;} 

void RGIsoT::OutputLabels(ArrayT<StringT>& labels) const 
{ 
     /*allocates space for labels*/
     labels.Dimension(kNumOutputVar); 
  
     /*copy labels*/
     for (int i = 0; i< kNumOutputVar; i++) 
       labels[i] = Labels[i]; 
} 

/* modulus */
const dMatrixT& RGIsoT::c_ijkl(void)
{
	/*calculate eigenvalues of  stretch tensor b*/
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
				    	
	/*total stretch tensor*/
	fb.MultAAT(fF3D);
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();
    const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();

	/*calculate equilibrium contribution*/
	/*calculate principle 2nd PK principal stress given 2nd PK principle strains*/
	ddWddE(fEigs,  fdWdE_eq, fddWddE_eq, -1);

	/*convert to tau and dtau_de_eq-2.0*tau*eig2 1*/
	fdWdE_eq[0] *= fEigs[0];
	fdWdE_eq[1] *= fEigs[1];
	fdWdE_eq[2] *= fEigs[2];

	fddWddE_eq(0,0) *= fddWddE_eq(0,0)*fEigs[0]*fEigs[0];
	fddWddE_eq(1,1) *= fddWddE_eq(1,1)*fEigs[1]*fEigs[1];
	fddWddE_eq(2,2) *= fddWddE_eq(2,2)*fEigs[2]*fEigs[2];
	fddWddE_eq(0,1) *= fddWddE_eq(2,3)*fEigs[0]*fEigs[1];
	fddWddE_eq(0,2) *= fddWddE_eq(0,2)*fEigs[0]*fEigs[2];
	fddWddE_eq(1,2) *= fddWddE_eq(1,2)*fEigs[1]*fEigs[2];	
   
	fModulus3D = fSpectralDecompSpat.EigsToRank4(fddWddE_eq);	
  
	double dl, coeff;

    double& l0 = fEigs[0];
    double& l1 = fEigs[1];
    double& l2 = fEigs[2];
	
	dl = l0 - l1;
    if (fabs(dl) > kSmall)
		coeff = (fdWdE_eq[0]*l1 - fdWdE_eq[1]*l0)/dl;
    else 
		coeff = 0.5*(fddWddE_eq(0,0)-fddWddE_eq(0,1))-fdWdE_eq[0];
    MixedRank4_3D(eigenvectors[0], eigenvectors[1], fModMat);
    fModulus3D.AddScaled(2.0*coeff, fModMat);
    
	dl = l0 - l2;
	if (fabs(dl) > kSmall)
		coeff = (fdWdE_eq[0]*l2 - fdWdE_eq[2]*l0)/dl;
	else 
		coeff = 0.5*(fddWddE_eq(0,0)-fddWddE_eq(0,2))-fdWdE_eq[2];	
	MixedRank4_3D(eigenvectors[0], eigenvectors[2], fModMat);
	fModulus3D.AddScaled(2.0*coeff, fModMat);
    
	dl = l1 - l2;
	if (fabs(dl) > kSmall)
		coeff  = (fdWdE_eq[1]*l2 - fdWdE_eq[2]*l1)/dl;
	else
		coeff = 0.5*(fddWddE_eq(1,1)-fddWddE_eq(1,2))-fdWdE_eq[1];	
	MixedRank4_3D(eigenvectors[1], eigenvectors[2], fModMat);
	fModulus3D.AddScaled(2.0*coeff, fModMat);


	/*calculate nonequilibrium contribution*/
	/*Load state variables (Cv and Cvn)*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());


	for (int i = 0; i < fNumProcess; i++)
	{
		/*trial stretch tensor*/
		fInverse.Inverse(fC_vn[i]);
		fb_tr.MultQBQT(fF3D, fInverse);

		fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_tr, false);	
		fEigs_tr = fSpectralDecompSpat.Eigenvalues(); 		

		fInverse.Inverse(fC_v[i]);
		fbe.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
		fEigs_e = fSpectralDecompSpat.Eigenvalues(); 
		const ArrayT<dArrayT>& eigenvectors_e=fSpectralDecompSpat.Eigenvectors();
	
		/*calculate principle 2nd PK principal stress given 2nd PK principle strains*/
		ddWddE(fEigs_e,  fdWdE_neq, fddWddE_neq, i);
		/*convert to tau and dtau_de_neq*/

		fdWdE_neq[0] *= fEigs_e[0];
		fdWdE_neq[1] *= fEigs_e[1];
		fdWdE_neq[2] *= fEigs_e[2];

		fddWddE_neq(0,0) *= fddWddE_neq(0,0)*fEigs_e[0]*fEigs_e[0];
		fddWddE_neq(1,1) *= fddWddE_neq(1,1)*fEigs_e[1]*fEigs_e[1];
		fddWddE_neq(2,2) *= fddWddE_neq(2,2)*fEigs_e[2]*fEigs_e[2];
		fddWddE_neq(0,0) += 2.0*fdWdE_neq[0];
		fddWddE_neq(1,1) += 2.0*fdWdE_neq[1];
		fddWddE_neq(2,2) += 2.0*fdWdE_neq[2];
	
		fddWddE_neq(0,1) *= fddWddE_neq(2,3)*fEigs_e[0]*fEigs_e[1];	
		fddWddE_neq(0,2) *= fddWddE_neq(0,2)*fEigs_e[0]*fEigs_e[2];
		fddWddE_neq(1,2) *= fddWddE_neq(1,2)*fEigs_e[1]*fEigs_e[2];	
	    
		ComputeiKAB(fddWddE_neq);
		dSymMatrixT& DAB = fddWddE_neq;
	
		fCalg(0,0) = DAB(0,0)*fiKAB(0,0) + DAB(0,1)*fiKAB(1,0) + DAB(0,2)*fiKAB(2,0) - 2.0*fdWdE_neq[0];
		fCalg(1,0) = DAB(1,0)*fiKAB(0,0) + DAB(1,1)*fiKAB(1,0) + DAB(1,2)*fiKAB(2,0);
		fCalg(2,0) = DAB(2,0)*fiKAB(0,0) + DAB(2,1)*fiKAB(1,0) + DAB(2,2)*fiKAB(2,0);
		fCalg(0,1) = DAB(0,0)*fiKAB(0,1) + DAB(0,1)*fiKAB(1,1) + DAB(0,2)*fiKAB(2,1);
		fCalg(1,1) = DAB(1,0)*fiKAB(0,1) + DAB(1,1)*fiKAB(1,1) + DAB(1,2)*fiKAB(2,1) - 2.0*fdWdE_neq[1];
		fCalg(2,1) = DAB(2,0)*fiKAB(0,1) + DAB(2,1)*fiKAB(1,1) + DAB(2,2)*fiKAB(2,1);
		fCalg(0,2) = DAB(0,0)*fiKAB(0,2) + DAB(0,1)*fiKAB(1,2) + DAB(0,2)*fiKAB(2,2);
		fCalg(1,2) = DAB(1,0)*fiKAB(0,2) + DAB(1,1)*fiKAB(1,2) + DAB(1,2)*fiKAB(2,2);
		fCalg(2,2) = DAB(2,0)*fiKAB(0,2) + DAB(2,1)*fiKAB(1,2) + DAB(2,2)*fiKAB(2,2) - 2.0*fdWdE_neq[2];
	   
		fModulus3D += fSpectralDecompSpat.NonSymEigsToRank4(fCalg);
    
		double dl_tr;

		double& l0_tr = fEigs_tr[0];
		double& l1_tr = fEigs_tr[1];
		double& l2_tr = fEigs_tr[2];
	
	
		dl_tr = l0_tr - l1_tr;
		if (fabs(dl_tr) > kSmall)
			coeff = (fdWdE_neq[0]*l1_tr - fdWdE_neq[1]*l0_tr)/dl_tr;
		else 
			coeff = 0.5*(fCalg(0,0)-fCalg(0,1))-fdWdE_neq[0];
		MixedRank4_3D(eigenvectors_e[0], eigenvectors_e[1], fModMat);
		fModulus3D.AddScaled(2.0*coeff, fModMat);
    
		dl_tr = l0_tr - l2_tr;
		if (fabs(dl_tr) > kSmall)
			coeff =(fdWdE_neq[0]*l2_tr - fdWdE_neq[2]*l0_tr)/dl_tr;
		else 
			coeff = 0.5*(fCalg(0,0)-fCalg(0,2))-fdWdE_neq[2];	
		MixedRank4_3D(eigenvectors_e[0], eigenvectors_e[2], fModMat);
		fModulus3D.AddScaled(2.0*coeff, fModMat);
    
		dl_tr = l1_tr - l2_tr;
		if (fabs(dl_tr) > kSmall)
			coeff  = (fdWdE_neq[1]*l2_tr - fdWdE_neq[2]*l1_tr)/dl_tr;
		else
			coeff = 0.5*(fCalg(1,1)-fCalg(1,2))-fdWdE_neq[1];	
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

	double J = sqrt(fEigs[0]*fEigs[1]*fEigs[2]);
    fModulus *= 1.0/J;

    return fModulus;
}

/* stresses */
const dSymMatrixT& RGIsoT::s_ij(void)
{
	/*calculate eigenvalues of  stretch tensor b*/
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
			
	/*calculates EQ contribution*/
	/*total stretches and eigs*/
	fb.MultAAT(fF3D);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();

	/*calculate principle 2nd PK principal stress given 2nd PK principle strains
	 * for equilibrium process    */
	dWdE(fEigs, fdWdE_eq, -1);

	/*convert to tau*/
	fdWdE_eq[0] *= fEigs[0];
	fdWdE_eq[1] *= fEigs[1];
	fdWdE_eq[2] *= fEigs[2];

	fStress3D = fSpectralDecompSpat.EigsToRank2(fdWdE_eq);

	/*compute NEQ contribution*/
    /*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
    {	
		for (int i = 0; i < fNumProcess; i++)
		{
			/*calc trial state*/
			fInverse.Inverse(fC_vn[i]);
			fb_tr.MultQBQT(fF3D, fInverse);

			fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_tr, false);	
			fEigs_tr = fSpectralDecompSpat.Eigenvalues(); 
			
			/*calc elastic stretch*/
			fEigs_e = fEigs_tr; /*predictor*/
			ComputeEigs_e(fEigs, fEigs_e, i); /*corrector*/

			/*calculate principle 2nd PK principal stress given 2nd PK principle strains
			 * for non-equilibrium process i    */
			dWdE(fEigs_e, fdWdE_neq, i);
	
			fdWdE_neq[0] *= fEigs_e[0];
			fdWdE_neq[1] *= fEigs_e[1];
			fdWdE_neq[2] *= fEigs_e[2];
		
			fStress3D += fSpectralDecompSpat.EigsToRank2(fdWdE_neq);

			/*Calculate Cv*/
			fInverse = fSpectralDecompSpat.EigsToRank2(fEigs_e); /*be which is colinear with b*/
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
		
			/*calculate principle 2nd PK principal stress given 2nd PK principle strains
			 * for non-equilibrium process i    */
			dWdE(fEigs_e, fdWdE_neq, i);

			fdWdE_neq[0] *= fEigs_e[0];
			fdWdE_neq[1] *= fEigs_e[1];
			fdWdE_neq[2] *= fEigs_e[2];
		
			fStress3D += fSpectralDecompSpat.EigsToRank2(fdWdE_neq);
		}
    }
    
    if (NumSD() == 2)
    {
        fStress[0] = fStress3D[0];
        fStress[1] = fStress3D[1];
        fStress[2] = fStress3D[5];
    }
    else fStress = fStress3D;

	double J = sqrt(fEigs[0]*fEigs[1]*fEigs[2]);
    fStress *= 1.0/J;
	return fStress;
}

/* material description */
const dMatrixT& RGIsoT::C_IJKL(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F();
  
    /* transform */
    fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl()));
    return fModulus;	
}

const dSymMatrixT& RGIsoT::S_IJ(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F();
  
    /* transform */
    fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij()));
    return fStress;
}

void RGIsoT::ComputeOutput(dArrayT& output)
{
	/*calculate eigenvalues of  stretch tensor b*/
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

	output[0] = 0.0;
	for (int i = 0; i < fNumProcess; i++)
	{
		/*calc elastic stretch*/
		fInverse.Inverse(fC_v[i]);
		fbe.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
		fEigs_e = fSpectralDecompSpat.Eigenvalues(); 
		
		/*calculate principle 2nd PK principal stress given 2nd PK principle strains*/
		dWdE(fEigs_e, fdWdE_neq, i);


		double sm = third*(fStress3D[0]+fStress3D[1]+fStress3D[2]);
		fStress3D[0] -= sm;
		fStress3D[1] -= sm;
		fStress3D[2] -= sm;

		output[0] += 0.5*(0.5*fietaS*fStress3D.ScalarProduct()+fietaB*sm*sm);
    }
}

/***********************************************************************
 * Private
 ***********************************************************************/

void RGIsoT::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, const int process_index) 
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
		/*calculate principle 2nd PK principal stress given 2nd PK principle strains*/
		ddWddE(eigenstretch_e, fdWdE_neq,  fddWddE_neq, process_index);
	    
		/*convert to tau*/
		fdWdE_neq[0] *= eigenstretch_e[0];
		fdWdE_neq[1] *= eigenstretch_e[1];
		fdWdE_neq[2] *= eigenstretch_e[2];
		
	    double s0 = fdWdE_neq[0];
	    double s1 = fdWdE_neq[1];
	    double s2 = fdWdE_neq[2];
	    	    
	    /*caculate means*/
	    double sm = third*(s0 + s1 + s2);
	    
	    ComputeiKAB(fddWddE_neq);
	    
	    /*calculate the residual*/
	    double dt = fFSMatSupport->TimeStep();
	    double res0 = ep_e0 + dt*(0.5*fietaS*(s0-sm) +
			  third*fietaB*sm) - ep_tr0;
	    double res1 = ep_e1 + dt*(0.5*fietaS*(s1-sm) +
			  third*fietaB*sm) - ep_tr1;
	    double res2 = ep_e2 + dt*(0.5*fietaS*(s2-sm) +
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

void RGIsoT::ComputeiKAB(dSymMatrixT& eigenmodulus_neq)
{	
        /*inverse viscosities*/

	/*deviatoric values*/
	double c0 = eigenmodulus_neq(0,0);
	double c1 = eigenmodulus_neq(1,1);
	double c2 = eigenmodulus_neq(2,2);

	double c12 = eigenmodulus_neq(1,2);
	double c02 = eigenmodulus_neq(0,2);
	double c01 = eigenmodulus_neq(0,1);
	
	double cm0 = third*(c0 + c01 + c02);
	double cm1 = third*(c1 + c01 + c12);
	double cm2 = third*(c2 + c02 + c12);

	/*calculates  KAB = 1+dt*D(dWdE_Idev/nD+isostress/nV)/Dep_e*/

	double dt = fFSMatSupport->TimeStep();
	fiKAB(0,0) = 1+0.5*fietaS*dt*(c0-cm0)+third*fietaB*dt*cm0;
	fiKAB(1,1) = 1+0.5*fietaS*dt*(c1-cm1)+third*fietaB*dt*cm1;
	fiKAB(2,2) = 1+0.5*fietaS*dt*(c2-cm2)+third*fietaB*dt*cm2;

	fiKAB(1,2) = 0.5*fietaS*dt*(c12-cm2)+third*fietaB*dt*cm2;
	fiKAB(0,2) = 0.5*fietaS*dt*(c02-cm2)+third*fietaB*dt*cm2;
	fiKAB(0,1) = 0.5*fietaS*dt*(c01-cm1)+third*fietaB*dt*cm1;
       
	fiKAB(2,1) = 0.5*fietaS*dt*(c12-cm1)+third*fietaB*dt*cm1;
	fiKAB(2,0) = 0.5*fietaS*dt*(c02-cm0)+third*fietaB*dt*cm0;
	fiKAB(1,0) = 0.5*fietaS*dt*(c01-cm0)+third*fietaB*dt*cm0;
	
	/*inverts KAB*/
	fiKAB.Inverse();
}

/* describe the parameters needed by the interface */
void RGIsoT::DefineParameters(ParameterListT& list) const
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

void RGIsoT::TakeParameterList(const ParameterListT& list)
{
  /* inherited */
  RGViscoelasticityT::TakeParameterList(list);

  /* dimension work space */
  fb.Dimension(3);
  fbe.Dimension(3);
  fb_tr.Dimension(3);
  fF3D.Dimension(3);
  fEigs.Dimension(3);
  fEigs_e.Dimension(3);
  fEigs_tr.Dimension(3);
  fdWdE_eq.Dimension(3);
  fddWddE_eq.Dimension(3);
  fdWdE_neq.Dimension(3);
  fddWddE_neq.Dimension(3);
  fCalg.Dimension(3);
  fModulus3D.Dimension(6);
  fModMat.Dimension(6);
  fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));
  fStress.Dimension(NumSD());
  fStress3D.Dimension(3);
  fiKAB.Dimension(3);
  fInverse.Dimension(3);

  /* viscosities */
  double etaS = list.GetParameter("eta_shear");
  double etaB = list.GetParameter("eta_bulk");
  fietaS = 1.0/etaS;
  fietaB = 1.0/etaB;
}


