/* $Id: RGSplit3D.cpp,v 1.3 2003-03-26 22:57:44 thao Exp $ */
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
  fSpectralDecompSpat(3),
  fSpectralDecompRef(3),
  fb(NumSD()),
  fb3D(3),
  fEigs(3),
  fEigs_e(3),
  ftau_EQ(3),
  ftau_NEQ(3),
  fDtauDe_EQ(3),
  fDtauDe_NEQ(3),
  fModMat(6),
  fModulus(dSymMatrixT::NumValues(NumSD())),
  fStress(NumSD()),
  fiKAB(3),
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

void RGSplit3D::Initialize(void)
{
    /*inheritance*/
    RGBaseT::Initialize();
    
    /*fndof initialized in RGBaseT*/
    int numstress = dSymMatrixT::NumValues(fndof);
    double* pstatev = fstatev.Pointer();  
    pstatev +=numstress;  //fC_v
    pstatev +=numstress;  //fC_vn
    fMatInStress.Set(fndof,pstatev);
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
  if (NumSD() == 2)
  {
    fb3D[0] = fb[0];
    fb3D[1] = fb[1];
    fb3D[2] = 1.0;
    
    fb3D[3] = 0.0;
    fb3D[4] = 0.0;
    fb3D[5] = fb[2];
  }
  else fb3D = fb;
  fSpectralDecompSpat.SpectralDecomp_Jacobi(fb3D, false);	
  fEigs = fSpectralDecompSpat.Eigenvalues();
  
  double J = sqrt(fEigs.Product());
  dArrayT eigenstretch_bar = fEigs;
  eigenstretch_bar *= pow(J, -2.0*fthird);
  
  double energy = 0.0;
  energy = fPot_EQ->Energy(eigenstretch_bar, J);
  
  /*calculates deviatoric and volumetric part of the elastic stretch */
  ElementCardT& element = CurrentElement();
  Load(element, CurrIP());
  
  fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);	
  dArrayT Eigs_v = fSpectralDecompRef.Eigenvalues();
  
  fEigs_e = fEigs;
  fEigs_e /= Eigs_v;
  
  double Je = sqrt(fEigs_e.Product());
  dArrayT eigenstretche_bar = fEigs_e;
  eigenstretche_bar *= pow(Je,-2.0*fthird);
  
  energy += fPot_NEQ->Energy(eigenstretche_bar, Je);
  
  return(energy);
}

/* modulus */
const dMatrixT& RGSplit3D::c_ijkl(void)
{
	/*spectral decomposition of stretch tensor*/
    Compute_b(fb);
    if (NumSD() == 2)
    {
      fb3D[0] = fb[0];
      fb3D[1] = fb[1];
      fb3D[2] = 1.0;
    
      fb3D[3] = 0.0;
      fb3D[4] = 0.0;
      fb3D[5] = fb[2];
    }
    else fb3D = fb;
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb3D, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();
	const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();

	/*jacobian determinants*/
	double J = sqrt(fEigs.Product());

	/*deviatoric principal stretches*/
	dArrayT eigenstretch_bar = fEigs;
	eigenstretch_bar *= pow(J, -2.0*fthird);

	/*retrieve viscous stretch tensor from history variables*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	/*calculate elastic principal stretches*/
	fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);
	dArrayT Eigs_v = fSpectralDecompRef.Eigenvalues();
	fEigs_e = fEigs;
	fEigs_e /= Eigs_v;

	/*jacobian determinants*/
	double Je = sqrt(fEigs_e.Product());

	/*deviatoric principal stretches*/
	dArrayT eigenstretche_bar = fEigs_e;
	eigenstretche_bar *= pow(Je,-2.0*fthird);
 
	/*principal components of spatial tangent moduli*/
	fPot_EQ->DevStress(eigenstretch_bar, ftau_EQ);
	ftau_EQ += fPot_EQ->MeanStress(J);

	fPot_EQ->DevMod(eigenstretch_bar, fDtauDe_EQ);
	fDtauDe_EQ += fPot_EQ->MeanMod(J);
	dSymMatrixT Gamma = fDtauDe_EQ;
	Gamma(0,0) -= 2.0*ftau_EQ[0];
	Gamma(1,1) -= 2.0*ftau_EQ[1];
	Gamma(2,2) -= 2.0*ftau_EQ[2];

	fPot_NEQ->DevStress(eigenstretche_bar, ftau_NEQ);
	ftau_NEQ += fPot_NEQ->MeanStress(Je);

    fPot_NEQ->DevMod(eigenstretche_bar, fDtauDe_NEQ);
	double cm = fPot_NEQ->MeanMod(Je);
	ComputeiKAB(fDtauDe_NEQ, cm);
	dSymMatrixT DAB = fDtauDe_NEQ;
	DAB += cm; 

	dMatrixT Calg(3);
	Calg(0,0) = DAB(0,0)*fiKAB(0,0) + DAB(0,1)*fiKAB(1,0) + DAB(0,2)*fiKAB(2,0) - 2.0*ftau_NEQ[0];
	Calg(1,0) = DAB(1,0)*fiKAB(0,0) + DAB(1,1)*fiKAB(1,0) + DAB(1,2)*fiKAB(2,0);
	Calg(2,0) = DAB(2,0)*fiKAB(0,0) + DAB(2,1)*fiKAB(1,0) + DAB(2,2)*fiKAB(2,0);
	Calg(0,1) = DAB(0,0)*fiKAB(0,1) + DAB(0,1)*fiKAB(1,1) + DAB(0,2)*fiKAB(2,1);
	Calg(1,1) = DAB(1,0)*fiKAB(0,1) + DAB(1,1)*fiKAB(1,1) + DAB(1,2)*fiKAB(2,1) - 2.0*ftau_NEQ[1];
	Calg(2,1) = DAB(2,0)*fiKAB(0,1) + DAB(2,1)*fiKAB(1,1) + DAB(2,2)*fiKAB(2,1);
	Calg(0,2) = DAB(0,0)*fiKAB(0,2) + DAB(0,1)*fiKAB(1,2) + DAB(0,2)*fiKAB(2,2);
	Calg(1,2) = DAB(1,0)*fiKAB(0,2) + DAB(1,1)*fiKAB(1,2) + DAB(1,2)*fiKAB(2,2);
	Calg(2,2) = DAB(2,0)*fiKAB(0,2) + DAB(2,1)*fiKAB(1,2) + DAB(2,2)*fiKAB(2,2) - 2.0*ftau_NEQ[2];

	double dl, coeff;
	/*Assemble moduli*/
	/*axial*/
	
	dMatrixT Modulus3D(6);
	Modulus3D = fSpectralDecompSpat.EigsToRank4(Gamma);	
	Modulus3D += fSpectralDecompSpat.NonSymEigsToRank4(Calg);
	
	double s0 = ftau_EQ[0]+ftau_NEQ[0];
	double s1 = ftau_EQ[1]+ftau_NEQ[1];
	double s2 = ftau_EQ[2]+ftau_NEQ[2];
	
	double& l0 = fEigs[0];
	double& l1 = fEigs[1];
	double& l2 = fEigs[2];
	
	/* 1,2 */
	dl = l0 - l1;
	/* modulus coefficient */
	if (fabs(dl) > kSmall)
	  coeff = (s0*l1 - s1*l0)/dl;
	else
	  coeff = 0.5*(Gamma(0,0)-Gamma(0,1)+Calg(0,0)-Calg(0,1))-s0;
	MixedRank4_3D(eigenvectors[0], eigenvectors[1], fModMat);
	Modulus3D.AddScaled(2.0*coeff, fModMat);
	
	/* 1,3 */
	dl = l0 - l2;
	/* modulus coefficient */
	if (fabs(dl) > kSmall)
	  coeff = (s0*l2 - s2*l0)/dl;
	else
	  coeff = 0.5*(Gamma(0,0)-Gamma(0,2)+Calg(0,0)-Calg(0,2))-s2;	
	MixedRank4_3D(eigenvectors[0], eigenvectors[2], fModMat);
	Modulus3D.AddScaled(2.0*coeff, fModMat);
	
	/* 2,3 */
	dl = l1 - l2;
	/* modulus coefficient */
	if (fabs(dl) > kSmall)
	  coeff = (s1*l2 - s2*l1)/dl;
	else
	  coeff = 0.5*(Gamma(1,1)-Gamma(1,2)+Calg(1,1)-Calg(1,2))-s1;	
	MixedRank4_3D(eigenvectors[1], eigenvectors[2], fModMat);
	Modulus3D.AddScaled(2.0*coeff, fModMat);

	if (NumSD() == 2)
    {
      fModulus[0] = Modulus3D[0];
      fModulus[1] = Modulus3D[1];
      fModulus[2] = Modulus3D[5];

      fModulus[3] = Modulus3D[6];
      fModulus[4] = Modulus3D[7];
      fModulus[5] = Modulus3D[11];

      fModulus[6] = Modulus3D[30];
      fModulus[7] = Modulus3D[31];
      fModulus[8] = Modulus3D[35];
    }
	else fModulus = Modulus3D;
	/*	cout << "\n stretch: "<<fEigs;
	cout << "\n estretch: "<<fEigs_e;
	cout << "\n stress eq: "<<ftau_EQ;
	cout << "\n stress neq: "<<ftau_NEQ;
	cout << "\n mod eq: "<<Gamma;
	cout << "\n mod neq: "<<Calg; */
	fModulus *= 1.0/J;
	return fModulus;
}

/* stresses */
const dSymMatrixT& RGSplit3D::s_ij(void)
{
	/* stretch tensor */
    Compute_b(fb);
    if (NumSD() == 2)
    {
      fb3D[0] = fb[0];
      fb3D[1] = fb[1];
      fb3D[2] = 1.0;
    
      fb3D[3] = 0.0;
      fb3D[4] = 0.0;
      fb3D[5] = fb[2];
    }
    else fb3D = fb;
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb3D, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();

	/*jacobian determinant*/
	double J = sqrt(fEigs.Product());

	/*deviatoric principal stretch*/
	dArrayT eigenstretch_bar = fEigs;
	eigenstretch_bar *= pow(J,-2.0*fthird);

	fPot_EQ->DevStress(eigenstretch_bar, ftau_EQ);
	ftau_EQ += fPot_EQ->MeanStress(J);

	/*load the viscoelastic principal stretches from state variable arrays*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	if (fFSMatSupport.RunState() == GlobalT::kFormRHS)
	{
		dSymMatrixT iCvn = fC_vn;
		iCvn.Inverse();

		/*calculate trial state;*/
		SpectralDecompT SpectralDecompTrial(3);
		dSymMatrixT b_tr(3);

		const dMatrixT& F = F_mechanical();
		dMatrixT F3D(3);
		if (NumSD() == 2)
		{
		  F3D[0] = F[0];
		  F3D[1] = F[1];
		  F3D[2] = 0.0;
		  
		  F3D[3] = F[2];
		  F3D[4] = F[3];
		  F3D[5] = 0.0;
		  
		  F3D[6] = 0.0;
		  F3D[7] = 0.0;
		  F3D[8] = 1.0;
		}
		else F3D = F;
		
		b_tr.MultQBQT(F3D,iCvn);
		SpectralDecompTrial.SpectralDecomp_Jacobi(b_tr, false);	
		fEigs_e = SpectralDecompTrial.Eigenvalues();

		ComputeEigs_e(fEigs, fEigs_e, ftau_NEQ, fDtauDe_NEQ);

		double Je = sqrt(fEigs_e.Product());
		dArrayT eigenstretche_bar = fEigs_e;
		eigenstretche_bar *= pow(Je,-2.0*fthird);

		fPot_NEQ->DevStress(eigenstretche_bar, ftau_NEQ);
		ftau_NEQ += fPot_NEQ->MeanStress(Je);

		/*store material stress in state variable array*/
		if (HasDissipVar())
		{
		  fMatInStress = fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
		  dSymMatrixT ib_e = fSpectralDecompSpat.EigsToRank2(fEigs_e);
		  ib_e.Inverse();
		  dSymMatrixT temp(3);
		  temp.MultAB(fMatInStress,ib_e);
		  fMatInStress.MultQTBQ(F3D,temp);
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

		double Je = sqrt(fEigs_e.Product());
		dArrayT eigenstretche_bar = fEigs_e;
		eigenstretche_bar *= pow(Je,-2.0*fthird);
		fPot_NEQ->DevStress(eigenstretche_bar, ftau_NEQ);
		ftau_NEQ += fPot_NEQ->MeanStress(Je);
	}

	/*evaluate cauchy stress*/
	dSymMatrixT Stress3D = fSpectralDecompSpat.EigsToRank2(ftau_EQ);
    Stress3D += fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
    
    if (NumSD() == 2)
    {
        fStress[0] = Stress3D[0];
        fStress[1] = Stress3D[1];
        fStress[2] = Stress3D[5];
    }
    else fStress = Stress3D;
	fStress *= 1.0/J;
	return fStress;
}

/* material description */
const dMatrixT& RGSplit3D::C_IJKL(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl()));
	return fModulus;	
}

const dSymMatrixT& RGSplit3D::S_IJ(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij()));
	return fStress;
}

void RGSplit3D::ComputeOutput(dArrayT& output)
{
	/* spectral decomposition */
	Compute_b(fb);
    if (NumSD() == 2)
    {
      fb3D[0] = fb[0];
      fb3D[1] = fb[1];
      fb3D[2] = 1.0;
    
      fb3D[3] = 0.0;
      fb3D[4] = 0.0;
      fb3D[5] = fb[2];
    }
    else fb3D = fb;
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb3D, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();
	double J = sqrt(fEigs.Product());

	/*load the viscoelastic principal stretches from state variable arrays*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);	
	dArrayT Eigs_v = fSpectralDecompRef.Eigenvalues();

	fEigs_e = fEigs;
	fEigs_e /= Eigs_v;

	/*calc jacobian*/
	double Je = sqrt(fEigs_e.Product()) ;
	dArrayT eigenstretche_bar = fEigs_e;
	eigenstretche_bar *= pow(Je,-2.0*fthird);

	dArrayT tau_NEQ(3);
	fPot_NEQ->DevStress(eigenstretche_bar, tau_NEQ);

	dSymMatrixT sdev = fSpectralDecompSpat.EigsToRank2(tau_NEQ);
	sdev /= J;
    double sm = 1.0/J*(fPot_NEQ->MeanStress(Je));

	double rate_visc_disp = 0.5*(sdev.ScalarProduct()*0.5*fietaS+sm*sm*fietaB);

	/*put in planestress option*/

	output[0] = rate_visc_disp;
}

/***********************************************************************
 * Private
 ***********************************************************************/
void RGSplit3D::ComputeEigs_e(const dArrayT& eigenstretch, 
			    dArrayT& eigenstretch_e, dArrayT& eigenstress, 
			    dSymMatrixT& eigenmodulus) 
{		
	const double ctol = 1.00e-14;
		
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
	    double Je=sqrt(le0*le1*le2);
		dArrayT eigenstretche_bar = eigenstretch_e;
		eigenstretche_bar *= pow(Je,-2.0*fthird);

	        /*calculate stresses and moduli*/
		fPot_NEQ->DevStress(eigenstretche_bar, eigenstress);

		double& s0 = eigenstress[0];
		double& s1 = eigenstress[1];
		double& s2 = eigenstress[2];

		fPot_NEQ->DevMod(eigenstretche_bar,eigenmodulus);

		/*caculate means*/
		double sm = fPot_NEQ->MeanStress(Je);
		double cm = fPot_NEQ->MeanMod(Je);

 	  	ComputeiKAB(eigenmodulus,cm);
		
	   	/*calculate the residual*/
	   	double dt = fFSMatSupport.TimeStep();
	 	double res0 = ep_e0 + dt*(0.5*fietaS*s0 +
					fthird*fietaB*sm) - ep_tr0;
	 	double res1 = ep_e1 + dt*(0.5*fietaS*s1 +
					fthird*fietaB*sm) - ep_tr1;
	 	double res2 = ep_e2 + dt*(0.5*fietaS*s2 +
					fthird*fietaB*sm) - ep_tr2;
		
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

