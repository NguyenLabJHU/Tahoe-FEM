/* $Id: RGSplit3D.cpp,v 1.1 2003-03-19 19:00:55 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "RGSplit3D.h"
#include "PotentialT.h"
#include "NeoHookean.h"

#include "fstreamT.h"
#include "ExceptionT.h"
#include <math.h>
#include <iostream.h>
#include <stdlib.h>

using namespace Tahoe;
 
const int kNumOutputVar =1; 
static const char* Labels[kNumOutputVar] = {"dW_visc"}; 

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
RGSplit3D::RGSplit3D(ifstreamT& in, const FSMatSupportT& support):
  RGBaseT(in, support),
  fb(NumSD()),
  fEigs(NumSD()),
  fEigs_e(NumSD()),
  fEigs_bar(NumSD()),
  fEigs_ebar(NumSD()),
  ftau_EQ(NumSD()),
  ftau_NEQ(NumSD()),
  fDtauDe_EQ(NumSD()),
  fDtauDe_NEQ(NumSD()),
  fModMat(dSymMatrixT::NumValues(NumSD())),
  fModulus(dSymMatrixT::NumValues(NumSD())),
  fStress(NumSD()),
  fiKAB(NumSD(),NumSD()),
  fthird(1.0/3.0)
{
  /*read in potential code*/
  int code;
  in >> code;
  switch(code)
  {
     case PotentialT::kNeoHookean: 
     {
       fPot_EQ = new NeoHookean(in);
       fPot_NEQ = new NeoHookean(in);
       break;
     }
     default:
     {
       throw ExceptionT::kBadInputValue;
     }
  }
  /*read in viscosities*/
  double etaS, etaB;
  in >> etaS;
  fietaS = 1.0/etaS;
  in >> etaB;
  fietaB = 1.0/etaB;
}

RGSplit3D::~RGSplit3D(void)
{
  delete fPot_EQ;
  delete fPot_NEQ;
}

/* print parameters */
void RGSplit3D::Print(ostream& out) const
{
  RGBaseT::Print(out);
  out<<"Equilibrium free energy potential\n";
  fPot_EQ->Print(out);
  out<<"Non Equilibrium free energy potential\n";
  fPot_NEQ->Print(out);
  
  out<<"Constant Viscosity \n";
  out<<"     Shear Viscosity: "<<1.0/fietaS<<'\n';
  out<<"     Bulk Viscosity: "<<1.0/fietaB<<'\n';
}

/* print name */
void RGSplit3D::PrintName(ostream& out) const
{
  /* inherited */
  RGBaseT::PrintName(out);
  out<<"Equilibrium free energy potential\n";
  fPot_EQ->PrintName(out);
  out<<"Non Equilibrium free energy potential\n";
  fPot_NEQ->PrintName(out);
}

int RGSplit3D::NumOutputVariables() const {return kNumOutputVar;} 

void RGSplit3D::OutputLabels(ArrayT<StringT>& labels) const 
{ 
  //allocates space for labels 
  labels.Dimension(kNumOutputVar); 
  
  //copy labels 
  for (int i = 0; i< kNumOutputVar; i++) 
    labels[i] = Labels[i]; 
} 

double RGSplit3D::StrainEnergyDensity(void)
{
  /*calculates deviatoric and volumetric part of the total stretch */
  Compute_b(fb);
  fb.PrincipalValues(fEigs);
  double J = sqrt(fEigs.Product());
  dArrayT eigenstretch_bar = fEigs;
  eigenstretch_bar *= pow(J, -2.0*fthird);
  
  double energy =0.0;
  energy = fPot_EQ->Energy(eigenstretch_bar, J);
  
  /*calculates deviatoric and volumetric part of the elastic stretch */
  ElementCardT& element = CurrentElement();
  Load(element, CurrIP());
  
  fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);	
  dArrayT Eigs_v(3);
  fC_v.PrincipalValues(Eigs_v);
  
  fEigs_e = fEigs;
  fEigs_e /= Eigs_v;
  
  J = sqrt(fEigs_e.Product());
  eigenstretch_bar.SetToScaled(pow(J,-2.0*fthird), fEigs_e);
  
  energy += fPot_NEQ->Energy(eigenstretch_bar, J);
  
  return(energy);
}

/* modulus */
const dMatrixT& RGSplit3D::c_ijkl(void)
{
        /*assumes s_ij() has been called*/

        double iJ=1.0/fJ;
  
	fPot_EQ->DevMod(fEigs_bar, fDtauDe_EQ);
	fDtauDe_EQ += fPot_EQ->MeanMod(fJ);
	dSymMatrixT Gamma = fDtauDe_EQ;
	Gamma *= iJ;

        fPot_NEQ->DevMod(fEigs_ebar, fDtauDe_NEQ);
	double cm = fPot_NEQ->MeanMod(fJe);
	ComputeiKAB(fDtauDe_NEQ, cm);
	dSymMatrixT DAB = fDtauDe_NEQ;
	DAB += cm; 
	dMatrixT Calg(NumSD());
	Calg.MultSymAB(DAB,fiKAB);
	Calg *= iJ;

	const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();
	const ArrayT<dArrayT>& eigenvectors2=fSpectralDecompRef.Eigenvectors();
	double dlamb, coeff;

	/*Assemble moduli*/
	/*axial*/
	Gamma(0,0) -= 2.0*iJ*ftau_EQ[0];
	Gamma(1,1) -= 2.0*iJ*ftau_EQ[1];
	Gamma(2,2) -= 2.0*iJ*ftau_EQ[2];

	Calg(0,0) -= 2.0*iJ*ftau_NEQ[0];
	Calg(1,1) -= 2.0*iJ*ftau_NEQ[1];
	Calg(2,2) -= 2.0*iJ*ftau_NEQ[2];
	
	fModulus = fSpectralDecompSpat.EigsToRank4(Gamma);	
	fModulus += fSpectralDecompSpat.NonSymEigsToRank4(Calg);
	
	double sig1 = iJ*(ftau_EQ[0]+ftau_NEQ[0]);
	double sig2 = iJ*(ftau_EQ[1]+ftau_NEQ[1]);
	double sig3 = iJ*(ftau_EQ[2]+ftau_NEQ[2]);
	
	double& lamb1 = fEigs[0];
	double& lamb2 = fEigs[1];
	double& lamb3 = fEigs[2];
	
	/* 1,2 */
	dlamb = lamb1 - lamb2;
	/* modulus coefficient */
	if (fabs(dlamb) > kSmall)
	  coeff = (sig1*lamb2 - sig2*lamb1)/dlamb;
	else
	  coeff = 0.5*(Gamma(0,0) - Gamma(0,1) + 
		       Calg(0,0) - Calg(0,1));
	MixedRank4_3D(eigenvectors[0], eigenvectors[1], fModMat);
	fModulus.AddScaled(2.0*coeff, fModMat);
	
	/* 1,3 */
	dlamb = lamb1 - lamb3;
	/* modulus coefficient */
	if (fabs(dlamb) > kSmall)
	  coeff = (sig1*lamb3 - sig3*lamb1)/dlamb;
	else
	  coeff = 0.5*(Gamma(0,0) - Gamma(0,2) + 
		       Calg(0,0) - Calg(0,2));	
	MixedRank4_3D(eigenvectors[0], eigenvectors[2], fModMat);
	fModulus.AddScaled(2.0*coeff, fModMat);
	
	/* 2,3 */
	dlamb = lamb2 - lamb3;
	/* modulus coefficient */
	if (fabs(dlamb) > kSmall)
	  coeff = (sig2*lamb3 - sig3*lamb2)/dlamb;
	else
	  coeff = 0.5*(Gamma(1,1) - Gamma(1,2) + 
		       Calg(1,1) - Calg(1,2));	
	MixedRank4_3D(eigenvectors[1], eigenvectors[2], fModMat);
	fModulus.AddScaled(2.0*coeff, fModMat);

	return fModulus;
}

/* stresses */
const dSymMatrixT& RGSplit3D::s_ij(void)
{
	/* stretch tensor */
	Compute_b(fb);

	/* spectral decomposition */
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	/*calculates principal stretch*/
	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	fJ = sqrt(fEigs.Product());
	fEigs_bar.SetToScaled(pow(fJ,-2.0*fthird),fEigs);

	fPot_EQ->DevStress(fEigs_bar, ftau_EQ);
	ftau_EQ += fPot_EQ->MeanStress(fJ);

	/*load the viscoelastic principal stretches from state variable arrays*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	if (fFSMatSupport.RunState() == GlobalT::kFormRHS)
	{
		dSymMatrixT iCvn = fC_vn;
		iCvn.Inverse();

		/*calculate trial state;*/
		dSymMatrixT b_tr(NumSD());
		const dMatrixT& F = F_mechanical();
		b_tr.MultQBQT(F,iCvn);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(b_tr, false);	
		fEigs_e = fSpectralDecompSpat.Eigenvalues();
		ComputeEigs_e(fEigs, fEigs_e, ftau_NEQ, fDtauDe_NEQ);

		fJe = sqrt(fEigs_e.Product());
		fEigs_ebar.SetToScaled(pow(fJe,-2.0*fthird),fEigs_e);
		fPot_NEQ->DevStress(fEigs_ebar, ftau_NEQ);
		ftau_NEQ += fPot_NEQ->MeanStress(fJe);

		/*store material stress in state variable array*/
		if (HasDissipVar())
		{
		  double* pstatev = fstatev.Pointer();
		  pstatev += (fC_v.Length()+fC_vn.Length());
		  dSymMatrixT mat_stress(NumSD(),pstatev);
		  fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
		  mat_stress = fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
		  dSymMatrixT b_e = fSpectralDecompSpat.EigsToRank2(fEigs_e);
		  b_e.Inverse();
		  dSymMatrixT B(b_e.Rows());
		  B.MultAB(mat_stress,b_e);
		  const dMatrixT& F = F_mechanical();
		  mat_stress.MultQTBQ(F,B);
		}

		/*update viscuous stretch tensor*/
		Compute_C(fC_v);
		fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v,false);
		dArrayT Eigs_v = fEigs;
		Eigs_v /= fEigs_e;
		fC_v = fSpectralDecompRef.EigsToRank2(Eigs_v);

		Store(element, CurrIP());
	}	
	else 
	{
		fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);
	        dArrayT Eigs_v = fSpectralDecompRef.Eigenvalues();

		fEigs_e = fEigs;
		fEigs_e /= Eigs_v;

		fJe = sqrt(fEigs_e.Product());
		fEigs_ebar.SetToScaled(pow(fJe,-2.0*fthird),fEigs_e);
		fPot_NEQ->DevStress(fEigs_ebar, ftau_NEQ);
		ftau_NEQ += fPot_NEQ->MeanStress(fJe);
	}

	/*evaluate cauchy stress*/
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	/*	cout <<"\nfEigs: "<<fEigs;*/
	fStress = fSpectralDecompSpat.EigsToRank2(ftau_EQ);
	/*	cout <<"\nfEQ: "<<ftau_EQ;*/
	fStress += fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
	/*	cout <<"\nfNEQ: "<<ftau_NEQ;*/
	/*        cout <<"\nfStress: "<<fStress;*/
	double iJ = 1.0/fJ;
	fStress *= iJ;

	return fStress;
}

/* material description */
const dMatrixT& RGSplit3D::C_IJKL(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(fJ, PullBack(Fmat, c_ijkl()));
	return fModulus;	
}

const dSymMatrixT& RGSplit3D::S_IJ(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(fJ, PullBack(Fmat, s_ij()));
	return fStress;
}

void RGSplit3D::ComputeOutput(dArrayT& output)
{
	/* stretch tensor */
	/* stretch tensor */
	Compute_b(fb);

	/* spectral decomposition */
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	/*calculates principal stretch*/
	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	/*calc jacobian*/
	fJ = sqrt(fEigs.Product()) ;
	fEigs_bar.SetToScaled(pow(fJ,-2.0*fthird),fEigs);

	double iJ = 1.0/fJ;
	/*load the viscoelastic principal stretches from state variable arrays*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);	
	dArrayT Eigs_v = fSpectralDecompRef.Eigenvalues();

	fEigs_e = fEigs;
	fEigs_e /= Eigs_v;
	fJe = sqrt(fEigs_e[0]*fEigs_e[1]*fEigs_e[2]);
	fEigs_ebar.SetToScaled(pow(fJe,-2.0*fthird),fEigs_e);

	dArrayT sig_NEQ(NumSD());
	fPot_NEQ->DevStress(fEigs_ebar, sig_NEQ);
	sig_NEQ *=iJ;
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	dSymMatrixT sdev = fSpectralDecompSpat.EigsToRank2(sig_NEQ);

      	double sm = fPot_NEQ->MeanStress(fJe)*iJ;

	double rate_visc_disp = 0.5*(sdev.ScalarProduct()*0.5*fietaS+
				     sm*sm*fietaB);

	/*put in planestress option*/

	output[0] = rate_visc_disp;
	/*	output[1] = fEigs_e[2];*/
}
/***********************************************************************
 * Private
 ***********************************************************************/
void RGSplit3D::ComputeEigs_e(const dArrayT& eigenstretch, 
			    dArrayT& eigenstretch_e, dArrayT& eigenstress, 
			    dSymMatrixT& eigenmodulus) 
{		
	const double ctol = 1.00e-10;
		
	/*set references to principle stretches*/
	double& l0 = eigenstretch[0];
	double& l1 = eigenstretch[1];
	double& l2 = eigenstretch[2];
      
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

	        fJe=sqrt(le0*le1*le2);
	        fEigs_ebar.SetToScaled(pow(fJe,-2.0*fthird),
					    eigenstretch_e);

	        /*calculate stresses and moduli*/
		fPot_NEQ->DevStress(fEigs_ebar, eigenstress);

		double& s0 = eigenstress[0];
		double& s1 = eigenstress[1];
		double& s2 = eigenstress[2];

		fPot_NEQ->DevMod(fEigs_ebar,eigenmodulus);

		/*caculate means*/
		double sm = fPot_NEQ->MeanStress(fJe);
		double cm = fPot_NEQ->MeanMod(fJe);

 	  	ComputeiKAB(eigenmodulus,cm);
		
	   	/*calculate the residual*/
	   	double dt = fFSMatSupport.TimeStep();
	 	double res0 = ep_e0 + dt*(0.5*fietaS*s0 +
					fthird*fietaB*sm) - ep_tr0;
	 	double res1 = ep_e1 + dt*(0.5*fietaS*s1 +
					fthird*fietaB*sm) - ep_tr1;
	 	double res2 = ep_e2 + dt*(0.5*fietaS*s2 +
					fthird*fietaB*sm) - ep_tr2;
		
		//	cout << "\n residual1 "<< res0;
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

void RGSplit3D::ComputeiKAB(dSymMatrixT& eigenmodulus, double& bulkmodulus)
{	
        /*inverse viscosities*/
        //		cout <<"\n Je: "<<Je;
	//		cout <<"\nfietaS: "<<fietaS;

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

	double dt = fFSMatSupport.TimeStep();
	KAB(0,0) = 1+0.5*fietaS*dt*c0+fthird*fietaB*dt*cm;
	KAB(1,1) = 1+0.5*fietaS*dt*c1+fthird*fietaB*dt*cm;
	KAB(2,2) = 1+0.5*fietaS*dt*c2+fthird*fietaB*dt*cm;

	KAB(1,2) = 0.5*fietaS*dt*c12+fthird*fietaB*dt*cm;
	KAB(0,2) = 0.5*fietaS*dt*c02+fthird*fietaB*dt*cm;
	KAB(0,1) = 0.5*fietaS*dt*c01+fthird*fietaB*dt*cm;
       
	KAB(2,1) = KAB(1,2);
	KAB(2,0) = KAB(0,2);
	KAB(1,0) = KAB(0,1);

	/*inverts KAB*/
	fiKAB.Inverse(KAB);
}

