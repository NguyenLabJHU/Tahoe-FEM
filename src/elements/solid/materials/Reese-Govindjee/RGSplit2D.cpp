/* $Id: RGSplit2D.cpp,v 1.2 2003-03-22 00:40:52 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "RGSplit2D.h"
#include "PotentialT.h"
#include "NeoHookean.h"

#include "fstreamT.h"
#include "ExceptionT.h"
#include <math.h>
#include <iostream.h>
#include <stdlib.h>

using namespace Tahoe;
 
const int kNumOutputVar = 1; 
static const char* Labels[kNumOutputVar] = {"dW_visc"}; 

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
RGSplit2D::RGSplit2D(ifstreamT& in, const FSMatSupportT& support):
  RGBaseT(in, support),
  fb(NumSD()),
  fEigs(NumSD()),
  fEigs_e(NumSD()),
  ftau_EQ(NumSD()),
  ftau_NEQ(NumSD()),
  fDtauDe_EQ(NumSD()),
  fDtauDe_NEQ(NumSD()),
  fModMat(dSymMatrixT::NumValues(NumSD())),
  fModulus(dSymMatrixT::NumValues(NumSD())),
  fStress(NumSD()),
  fiKAB(NumSD()),
  fthird(1.0/3.0)
{
  /*read in potential code*/
  int code;
  in >> code;
  switch(code)
  {
     case PotentialT::kNeoHookean: 
     {
       fPhi_EQ = new NeoHookean(in);
       fPhi_NEQ = new NeoHookean(in);
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

RGSplit2D::~RGSplit2D(void)
{
  delete fPhi_EQ;
  delete fPhi_NEQ;
}
/* print parameters */
void RGSplit2D::Print(ostream& out) const
{
  RGBaseT::Print(out);
  out<<"Equilibrium free energy potential\n";
  fPhi_EQ->Print(out);
  out<<"Non Equilibrium free energy potential\n";
  fPhi_NEQ->Print(out);
  
  out<<"Constant Viscosity \n";
  out<<"     Shear Viscosity: "<<1.0/fietaS<<'\n';
  out<<"     Bulk Viscosity: "<<1.0/fietaB<<'\n';
}

/* print name */
void RGSplit2D::PrintName(ostream& out) const
{
  /* inherited */
  RGBaseT::PrintName(out);
  out<<"       2D PlaneStrain\n";
  out<<"Equilibrium free energy potential\n";
  fPhi_EQ->PrintName(out);
  out<<"Non Equilibrium free energy potential\n";
  fPhi_NEQ->PrintName(out);
}

int RGSplit2D::NumOutputVariables() const {return kNumOutputVar;} 

void RGSplit2D::OutputLabels(ArrayT<StringT>& labels) const 
{ 
  //allocates space for labels 
  labels.Dimension(kNumOutputVar); 
  
  //copy labels 
  for (int i = 0; i< kNumOutputVar; i++) 
    labels[i] = Labels[i]; 
} 

double RGSplit2D::StrainEnergyDensity(void)
{
  /*calculates deviatoric and volumetric part of the total stretch */
  Compute_b(fb);
  fb.PrincipalValues(fEigs);

  double J = sqrt(fEigs.Product());

  dArrayT eigenstretch_bar(3);
  double a = pow(J, -2.0*fthird);
  eigenstretch_bar[0] = a*(fEigs[0]);
  eigenstretch_bar[1] = a*(fEigs[1]);
  eigenstretch_bar[2] = a;
  
  double energy =0.0;
  energy = fPhi_EQ->Energy(eigenstretch_bar, J);
  
  /*calculates deviatoric and volumetric part of the elastic stretch */
  ElementCardT& element = CurrentElement();
  Load(element, CurrIP());
  
  fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);	
  dArrayT Eigs_v(fC_v.Rows());
  fC_v.PrincipalValues(Eigs_v);
  
  fEigs_e = fEigs;
  fEigs_e /= Eigs_v;

  double Je = sqrt(fEigs_e.Product());

  dArrayT eigenstretche_bar(3);
  a = pow(Je, -2.0*fthird);
  eigenstretche_bar[0] = a*(fEigs_e[0]);
  eigenstretche_bar[1] = a*(fEigs_e[1]);
  eigenstretche_bar[2] = a;
  
  energy += fPhi_NEQ->Energy(eigenstretche_bar, Je);
  
  return(energy);
}

/* modulus */
const dMatrixT& RGSplit2D::c_ijkl(void)
{
	Compute_b(fb);

	/* spectral decomposition */
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	

	/*calculates principal stretches*/
	fEigs = fSpectralDecompSpat.Eigenvalues();
	const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();

	/*retrieve viscous stretch tensor from history variables*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	/*calculate elastic principal stretches*/
	fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);
	dArrayT Eigs_v = fSpectralDecompRef.Eigenvalues();
	fEigs_e = fEigs;
	fEigs_e /= Eigs_v;

	/*jacobian determinants*/
	double J = sqrt(fEigs.Product());
	double Je = sqrt(fEigs_e.Product());
  
	dArrayT eigenstretch_bar(3);
       	double a = pow(J, -2.0*fthird);
	eigenstretch_bar[0] = a*(fEigs[0]);
	eigenstretch_bar[1] = a*(fEigs[1]);
	eigenstretch_bar[2] = a;

	dArrayT eigenstretche_bar(3);
       	a = pow(Je, -2.0*fthird);
	eigenstretche_bar[0] = a*(fEigs_e[0]);
	eigenstretche_bar[1] = a*(fEigs_e[1]);
	eigenstretche_bar[2] = a;

	/*principal components of spatial tangent moduli*/
	fPhi_EQ->DevStress(eigenstretch_bar, ftau_EQ);
	ftau_EQ += fPhi_EQ->MeanStress(J);

	fPhi_EQ->DevMod(eigenstretch_bar, fDtauDe_EQ);
	fDtauDe_EQ += fPhi_EQ->MeanMod(J);
	dSymMatrixT Gamma = fDtauDe_EQ;
	Gamma(0,0) -= 2.0*ftau_EQ[0];
	Gamma(1,1) -= 2.0*ftau_EQ[1];
	Gamma /= J;

	fPhi_NEQ->DevStress(eigenstretche_bar, ftau_NEQ);
	ftau_NEQ += fPhi_NEQ->MeanStress(Je);

        fPhi_NEQ->DevMod(eigenstretche_bar, fDtauDe_NEQ);
	double cm = fPhi_NEQ->MeanMod(Je);
	ComputeiKAB(fDtauDe_NEQ, cm);
	dSymMatrixT DAB = fDtauDe_NEQ;
	DAB += cm; 
	dMatrixT Calg(NumSD());
	Calg.MultSymAB(DAB,fiKAB);
	Calg(0,0) = DAB(0,0)*fiKAB(0,0)+DAB(0,1)*fiKAB(1,0) - 2.0*ftau_NEQ[0];
	Calg(1,0) = DAB(1,0)*fiKAB(0,0)+DAB(1,1)*fiKAB(1,0);
	Calg(0,1) = DAB(0,0)*fiKAB(0,1)+DAB(0,1)*fiKAB(1,1);
	Calg(1,1) = DAB(1,0)*fiKAB(0,1)+DAB(1,1)*fiKAB(1,1) - 2.0*ftau_NEQ[1];
	Calg /= J;

	double dlamb, coeff;

	/*Assemble moduli and stresses*/
	/*axial*/
	fModulus = fSpectralDecompSpat.EigsToRank4(Gamma);	
	fModulus += fSpectralDecompSpat.NonSymEigsToRank4(Calg);
	
	double sig1 = (ftau_EQ[0]+ftau_NEQ[0])/J;
	double sig2 = (ftau_EQ[1]+ftau_NEQ[1])/J;
	
	double& lamb1 = fEigs[0];
	double& lamb2 = fEigs[1];
	
	/* 1,2 */
	dlamb = lamb1 - lamb2;
	/* modulus coefficient */
	if (fabs(dlamb) > kSmall)
	  coeff = (sig1*lamb2 - sig2*lamb1)/dlamb;
	else
	  coeff = 0.5*(Gamma(0,0)-Gamma(0,1)+Calg(0,0)-Calg(0,1))-sig1;
	MixedRank4_2D(eigenvectors[0], eigenvectors[1], fModMat);
	fModulus.AddScaled(2.0*coeff, fModMat);

	return fModulus;
}

/* stresses */
const dSymMatrixT& RGSplit2D::s_ij(void)
{
	/* stretch tensor */
	Compute_b(fb);

	/* spectral decomposition */
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	/*calculates principal stretch*/
	
	fEigs = fSpectralDecompSpat.Eigenvalues();

	double J = sqrt(fEigs.Product());

	dArrayT eigenstretch_bar(3);
       	double a = pow(J, -2.0*fthird);
	eigenstretch_bar[0] = a*(fEigs[0]);
	eigenstretch_bar[1] = a*(fEigs[1]);
	eigenstretch_bar[2] = a;

	fPhi_EQ->DevStress(eigenstretch_bar, ftau_EQ);
	ftau_EQ += fPhi_EQ->MeanStress(J);
	
	/*load the viscoelastic principal stretches from state variable arrays*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	if (fFSMatSupport.RunState() == GlobalT::kFormRHS)
	{
		dSymMatrixT iCvn = fC_vn;
		iCvn.Inverse();

		/*calculate trial state;*/
		SpectralDecompT SpectralDecompTrial(NumSD());
		dSymMatrixT b_tr(NumSD());

		const dMatrixT& F = F_mechanical();
		b_tr.MultQBQT(F,iCvn);
		SpectralDecompTrial.SpectralDecomp_Jacobi(b_tr, false);

		/*set initial value of elastic principal stretches 
		 * to trial values*/
		fEigs_e = SpectralDecompTrial.Eigenvalues();
		ComputeEigs_e(fEigs, fEigs_e, ftau_NEQ, fDtauDe_NEQ);

		double Je = sqrt(fEigs_e.Product());

		dArrayT eigenstretche_bar(3);
		double a = pow(Je, -2.0*fthird);
		eigenstretche_bar[0] = a*(fEigs_e[0]);
		eigenstretche_bar[1] = a*(fEigs_e[1]);
		eigenstretche_bar[2] = a;

		fPhi_NEQ->DevStress(eigenstretche_bar, ftau_NEQ);
		ftau_NEQ += fPhi_NEQ->MeanStress(Je);

		/*store material stress in state variable array*/
		if (HasDissipVar())
		{
		  double* pstatev = fstatev.Pointer();
		  pstatev += (fC_v.Length()+fC_vn.Length());

		  dSymMatrixT mat_stress(NumSD(),pstatev);
		  mat_stress = fSpectralDecompSpat.EigsToRank2(ftau_NEQ);

		  dSymMatrixT ib_e = fSpectralDecompSpat.EigsToRank2(fEigs_e);
		  ib_e.Inverse();

		  dSymMatrixT B(NumSD());
		  B.MultAB(mat_stress,ib_e);

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

		double Je = sqrt(fEigs_e.Product());

		dArrayT eigenstretche_bar(3);
		double a = pow(Je, -2.0*fthird);

		eigenstretche_bar[0] = a*(fEigs_e[0]);
		eigenstretche_bar[1] = a*(fEigs_e[1]);
		eigenstretche_bar[2] = a;

		fPhi_NEQ->DevStress(eigenstretche_bar, ftau_NEQ);
		ftau_NEQ += fPhi_NEQ->MeanStress(Je);
	}
	fStress = fSpectralDecompSpat.EigsToRank2(ftau_EQ);
	fStress += fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
	fStress /= J;
	
	return fStress;
}

/* material description */
const dMatrixT& RGSplit2D::C_IJKL(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl()));
	return fModulus;	
}

const dSymMatrixT& RGSplit2D::S_IJ(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij()));
	return fStress;
}

void RGSplit2D::ComputeOutput(dArrayT& output)
{
	/* stretch tensor */
	/* stretch tensor */
	Compute_b(fb);

	/* spectral decomposition */
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	

	/*calculates principal stretch*/
	fEigs = fSpectralDecompSpat.Eigenvalues();

	/*calc jacobian*/
	double J = sqrt(fEigs.Product()) ;

	/*load the viscoelastic principal stretches from state variable arrays*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);	
	dArrayT Eigs_v = fSpectralDecompRef.Eigenvalues();

	fEigs_e = fEigs;
	fEigs_e /= Eigs_v;
	double Je = sqrt(fEigs_e.Product());

	dArrayT eigenstretche_bar(3);
	double a = pow(Je, -2.0*fthird);
	eigenstretche_bar[0] = a*(fEigs_e[0]);
	eigenstretche_bar[1] = a*(fEigs_e[1]);
	eigenstretche_bar[2] = a;

	dArrayT sig_NEQ(NumSD());
	fPhi_NEQ->DevStress(eigenstretche_bar, sig_NEQ);
	sig_NEQ /= J;

	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	dSymMatrixT sdev = fSpectralDecompSpat.EigsToRank2(sig_NEQ);

      	double sm = 1.0/J*fPhi_NEQ->MeanStress(Je);

	double rate_visc_disp = 0.5*(sdev.ScalarProduct()*0.5*fietaS+
				     sm*sm*fietaB);

	/*put in planestress option*/

	output[0] = rate_visc_disp;
}
/***********************************************************************
 * Private
 ***********************************************************************/
void RGSplit2D::ComputeEigs_e(const dArrayT& eigenstretch, 
			    dArrayT& eigenstretch_e, dArrayT& eigenstress, 
			    dSymMatrixT& eigenmodulus) 
{		
	const double ctol = 1.00e-10;
		
	/*set references to principle stretches*/
	double& l0 = eigenstretch[0];
	double& l1 = eigenstretch[1];
      
	double& le0 = eigenstretch_e[0];
	double& le1 = eigenstretch_e[1];
  
	double tol;

	/*initialize principle elastic and trial elastic log strains */
	double ep_tr0 = 0.5*log(le0);
	double ep_tr1 = 0.5*log(le1);

	double ep_e0 = ep_tr0;		
	double ep_e1 = ep_tr1;	

	/*initializes principle viscous stretch*/
	do 
	{

	        double Je=sqrt(le0*le1);

		dArrayT eigenstretche_bar(3);
		double a = pow(Je, -2.0*fthird);
		eigenstretche_bar[0] = a*(eigenstretch_e[0]);
		eigenstretche_bar[1] = a*(eigenstretch_e[1]);
		eigenstretche_bar[2] = a;

	        /*calculate stresses and moduli*/
		fPhi_NEQ->DevStress(eigenstretche_bar, eigenstress);

		double& s0 = eigenstress[0];
		double& s1 = eigenstress[1];

		fPhi_NEQ->DevMod(eigenstretche_bar,eigenmodulus);

		/*caculate means*/
		double sm = fPhi_NEQ->MeanStress(Je);
		double cm = fPhi_NEQ->MeanMod(Je);

 	  	ComputeiKAB(eigenmodulus,cm);
		
	   	/*calculate the residual*/
	   	double dt = fFSMatSupport.TimeStep();
	 	double res0 = ep_e0 + dt*(0.5*fietaS*s0 +
			      fthird*fietaB*sm) - ep_tr0;
	 	double res1 = ep_e1 + dt*(0.5*fietaS*s1 +
			      fthird*fietaB*sm) - ep_tr1;
		//	cout << "\n residual1 "<< res0;
		/*solve for the principal strain increments*/
		double dep_e0=-fiKAB(0,0)*res0-fiKAB(0,1)*res1;
		double dep_e1=-fiKAB(1,0)*res0-fiKAB(1,1)*res1;
	 	
		/*updates principal elastic stretches*/ 
		ep_e0 += dep_e0;
		ep_e1 += dep_e1;
		//	cout << "\n depsilon1 "<< dep_e0;

		le0 = exp(2.0*ep_e0);
		le1 = exp(2.0*ep_e1);
		 
		/*Check that the L2 norm of the residual is less than tolerance*/
		tol = sqrt(res0*res0 + res1*res1);
	}while (tol>ctol); 
}

void RGSplit2D::ComputeiKAB(dSymMatrixT& eigenmodulus, double& bulkmodulus)
{	
        /*inverse viscosities*/
        //		cout <<"\n Je: "<<Je;
	//		cout <<"\nfietaS: "<<fietaS;

	/*deviatoric values*/
	double& c0 = eigenmodulus[0];
	double& c1 = eigenmodulus[1];
	double& c01 = eigenmodulus[2];
	
	/*mean value*/
	double& cm = bulkmodulus;

	dMatrixT& KAB = fiKAB;
		
	/*calculates  KAB = 1+dt*D(dWdE_Idev/nD+isostress/nV)/Dep_e*/

	double dt = fFSMatSupport.TimeStep();
	KAB(0,0) = 1+0.5*fietaS*dt*c0+fthird*fietaB*dt*cm;
	KAB(1,1) = 1+0.5*fietaS*dt*c1+fthird*fietaB*dt*cm;
	KAB(0,1) = 0.5*fietaS*dt*c01+fthird*fietaB*dt*cm;
	KAB(1,0) = KAB(0,1);

	/*inverts KAB*/
	fiKAB.Inverse(KAB);
}

