/* $Id: RGSplitT.cpp,v 1.4 2004-07-15 08:29:30 paklein Exp $ */
/* created: TDN (01/22/2001) */
#include "RGSplitT.h"
#include "PotentialT.h"
#include "NeoHookean.h"

#include "ExceptionT.h"
#include <math.h>
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
	fSpectralDecompSpat(3),
	fSpectralDecompRef(3),
	fSpectralDecompTrial(3),
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
     fEigs_bar = fEigs;
     fEigs_bar *= pow(J, -2.0*third);
     
     double energy = 0.0;
     energy = fPot_EQ->Energy(fEigs_bar, J);
  
     /*calculates deviatoric and volumetric part of the elastic stretch */
     ElementCardT& element = CurrentElement();
     Load(element, CurrIP());
  
     fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);	
     fEigs_v = fSpectralDecompRef.Eigenvalues();
  
     fEigs_e = fEigs;
     fEigs_e /= fEigs_v;
  
     double Je = sqrt(fEigs_e.Product());
     fEigs_ebar = fEigs_e;
     fEigs_ebar *= pow(Je,-2.0*third);
  
     energy += fPot_NEQ->Energy(fEigs_ebar, Je);
  
     return(energy);
}

/* modulus */
const dMatrixT& RGSplitT::c_ijkl(void)
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
    fEigs_bar = fEigs;
    fEigs_bar *= pow(J, -2.0*third);
    
    /*retrieve viscous stretch tensor from history variables*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    
    /*calculate elastic principal stretches*/
    fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);
    fEigs_v = fSpectralDecompRef.Eigenvalues();
    fEigs_e = fEigs;
    fEigs_e /= fEigs_v;
    
    /*jacobian determinants*/
    double Je = sqrt(fEigs_e.Product());
    
    /*deviatoric principal stretches*/
    fEigs_ebar = fEigs_e;
    fEigs_ebar *= pow(Je,-2.0*third);
    
    /*principal components of spatial tangent moduli*/
    fPot_EQ->DevStress(fEigs_bar, ftau_EQ);
    ftau_EQ += fPot_EQ->MeanStress(J);
    
    fPot_EQ->DevMod(fEigs_bar, fDtauDe_EQ);
    fDtauDe_EQ += fPot_EQ->MeanMod(J);
    dSymMatrixT& Gamma = fDtauDe_EQ;
    Gamma(0,0) -= 2.0*ftau_EQ[0];
    Gamma(1,1) -= 2.0*ftau_EQ[1];
    Gamma(2,2) -= 2.0*ftau_EQ[2];
    
    fPot_NEQ->DevStress(fEigs_ebar, ftau_NEQ);
    ftau_NEQ += fPot_NEQ->MeanStress(Je);
    
    fPot_NEQ->DevMod(fEigs_ebar, fDtauDe_NEQ);
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
    
    double dl, coeff;
    /*Assemble moduli*/
    /*axial*/
    
    fModulus3D = fSpectralDecompSpat.EigsToRank4(Gamma);	
    fModulus3D += fSpectralDecompSpat.NonSymEigsToRank4(fCalg);
    
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
      coeff = 0.5*(Gamma(0,0)-Gamma(0,1)+fCalg(0,0)-fCalg(0,1))-s0;
    MixedRank4_3D(eigenvectors[0], eigenvectors[1], fModMat);
    fModulus3D.AddScaled(2.0*coeff, fModMat);
    
    /* 1,3 */
    dl = l0 - l2;
    /* modulus coefficient */
    if (fabs(dl) > kSmall)
      coeff = (s0*l2 - s2*l0)/dl;
    else
      coeff = 0.5*(Gamma(0,0)-Gamma(0,2)+fCalg(0,0)-fCalg(0,2))-s2;	
    MixedRank4_3D(eigenvectors[0], eigenvectors[2], fModMat);
    fModulus3D.AddScaled(2.0*coeff, fModMat);
    
    /* 2,3 */
    dl = l1 - l2;
    /* modulus coefficient */
    if (fabs(dl) > kSmall)
      coeff = (s1*l2 - s2*l1)/dl;
    else
      coeff = 0.5*(Gamma(1,1)-Gamma(1,2)+fCalg(1,1)-fCalg(1,2))-s1;	
    MixedRank4_3D(eigenvectors[1], eigenvectors[2], fModMat);
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
    fEigs_bar = fEigs;
    fEigs_bar *= pow(J,-2.0*third);
    
    fPot_EQ->DevStress(fEigs_bar, ftau_EQ);
    ftau_EQ += fPot_EQ->MeanStress(J);
    
    /*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
    {
        dSymMatrixT& iCvn = fC_vn;
	iCvn.Inverse();
	
	/*calculate trial state;*/
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
	
	fb_tr.MultQBQT(fF3D,iCvn);
	fSpectralDecompTrial.SpectralDecomp_Jacobi(fb_tr, false);	
	fEigs_e = fSpectralDecompTrial.Eigenvalues();
	
	ComputeEigs_e(fEigs, fEigs_e, ftau_NEQ, fDtauDe_NEQ);
	
	double Je = sqrt(fEigs_e.Product());
	fEigs_ebar = fEigs_e;
	fEigs_ebar *= pow(Je,-2.0*third);
	
	fPot_NEQ->DevStress(fEigs_ebar, ftau_NEQ);
	ftau_NEQ += fPot_NEQ->MeanStress(Je);
	
	/*update viscuous stretch tensor*/
	if (NumSD() == 2)
	{
		Compute_C(fC_v_2D);
		fC_v.ExpandFrom2D(fC_v_2D);
		fC_v(2,2) = 1.0; /* no out-of-plane stretch */
	}
	else
		Compute_C(fC_v);
	fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v,false);
	fEigs_v = fEigs;
	fEigs_v /= fEigs_e;
	fC_v = fSpectralDecompRef.EigsToRank2(fEigs_v);
	iCvn.Inverse();

	Store(element, CurrIP());
    }	
    else 
    {
        fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);
	fEigs_v = fSpectralDecompRef.Eigenvalues();
	
	fEigs_e = fEigs;
	fEigs_e /= fEigs_v;
	
	double Je = sqrt(fEigs_e.Product());
	fEigs_ebar = fEigs_e;
	fEigs_ebar *= pow(Je,-2.0*third);
	fPot_NEQ->DevStress(fEigs_ebar, ftau_NEQ);
	ftau_NEQ += fPot_NEQ->MeanStress(Je);
    }

    /*evaluate cauchy stress*/
    fStress3D = fSpectralDecompSpat.EigsToRank2(ftau_EQ);
    fStress3D += fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
    
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
    fEigs_v = fSpectralDecompRef.Eigenvalues();
    
    fEigs_e = fEigs;
    fEigs_e /= fEigs_v;
    
    /*calc jacobian*/
    double Je = sqrt(fEigs_e.Product()) ;
    fEigs_ebar = fEigs_e;
    fEigs_ebar *= pow(Je,-2.0*third);
    
    fPot_NEQ->DevStress(fEigs_ebar, ftau_NEQ);
    
    fStress3D = fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
    double sm = fPot_NEQ->MeanStress(Je);
    
    double rate_visc_disp = 0.5*(0.5*fietaS*fStress3D.ScalarProduct()+fietaB*sm*sm);
    
    /*put in planestress option*/
    
    output[0] = rate_visc_disp;
}

/***********************************************************************
 * Private
 ***********************************************************************/
void RGSplitT::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
			     dArrayT& eigenstress, dSymMatrixT& eigenmodulus) 
{		
	const double ctol = 1.00e-14;
		
	/*set references to principle stretches*/
	const double& l0 = eigenstretch[0];
	const double& l1 = eigenstretch[1];
	const double& l2 = eigenstretch[2];
      
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
	    fEigs_ebar = eigenstretch_e;
	    fEigs_ebar *= pow(Je,-2.0*third);
		
	    /*calculate stresses and moduli*/
	    fPot_NEQ->DevStress(fEigs_ebar, eigenstress);
	    
	    double& s0 = eigenstress[0];
	    double& s1 = eigenstress[1];
	    double& s2 = eigenstress[2];
	    
	    fPot_NEQ->DevMod(fEigs_ebar,eigenmodulus);
	    
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

/* accept parameter list */
void RGSplitT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	RGViscoelasticityT::TakeParameterList(list);

	/* dimension work space */
	fb.Dimension(NumSD());
	fb3D.Dimension(3);
	fbe.Dimension(3);
	fb_tr.Dimension(3);
	fC_v_2D.Dimension(2);
	fF3D.Dimension(3);
	fEigs.Dimension(3);
	fEigs_e.Dimension(3);
	fEigs_bar.Dimension(3);
	fEigs_ebar.Dimension(3);
	fEigs_v.Dimension(3);
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
