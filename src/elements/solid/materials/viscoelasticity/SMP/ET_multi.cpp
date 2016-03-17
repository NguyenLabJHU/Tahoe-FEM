/* $Id: ET_multi.cpp,v 1.1 2016-03-17 13:43:05 tahoe.vickynguyen Exp $ */
/* created: TDN (01/22/2001) */

#include "ET_multi.h"

#include "PotentialT.h"
#include "NeoHookean.h"
#include "ArrudaBoyce.h"

#include "ifstreamT.h"
#include "ExceptionT.h"
#include <cmath>
#include "ofstreamT.h"
#include <cstdlib>
#include "ParameterContainerT.h"

using namespace Tahoe;

const double loge = log10(exp(1.0));
const double pi = 2.0*acos(0);
const double third = 1.0/3.0;
const double small = 1.0e-10;
const double big = 1.0e-10;

const int kNumOutputVar =6;
static const char* Labels[kNumOutputVar] = {"J", "incr_structural_heat",  "incr_latent_heat", "inc_plastic_work","total_incr_heat","heat_deficit"};

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
ET_multi::ET_multi(void):
ParameterInterfaceT("ET_multi")
{
}

double ET_multi::StructuralRelaxationFunc(const double Temperature, const double Sc)
{
    /*calculate entropy at Tg*/
	double ScTg=fBB/fTg/(fdeltac*log(fTg/fT2));
    
	/*Hodge's Model*/
    //	double coeff = (fBB/Temperature/Sc-ScTg)/loge;
	double coeff = (fBB/Temperature/Sc-ScTg);
	double itauR = exp(-coeff);          /*eq.(17)*/
    
	/*check the limits
     if (tauR > big)
     tauR = big;
     else if (tauR < small)
     tauR = small;
     */
	return(itauR);
}

double ET_multi::StressRelaxationFunc(const double Temperature, const double smag)

{
	double itauS = 1.0;
    
	if (smag > kSmall)
	{
		itauS *= sinh(fQS/Temperature * smag);
		itauS *= Temperature/(fQS*smag);
	}
    
	return(itauS);
}

int ET_multi::NumOutputVariables() const {return kNumOutputVar;}

void ET_multi::OutputLabels(ArrayT<StringT>& labels) const
{
    /*allocates space for labels*/
    labels.Dimension(kNumOutputVar);
    
    /*copy labels*/
    for (int i = 0; i< kNumOutputVar; i++)
        labels[i] = Labels[i];
}


/* incremental heat generation */
double ET_multi::IncrementalHeat(void)
{
	/* trust the "current" element is already loaded */
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
//    cout <<"\n:"<< *fHeat;
    return(*fHeat);
}

const dMatrixT& ET_multi::ThermalDeformation_Inverse(void) /*eq.(13)?,eq.(5)TDN*/
{
	fF_T_inv.Identity(1.0);
	return(fF_T_inv);
}

double ET_multi::StrainEnergyDensity(void)
{
	/*calculates equilibrium part*/
	const dMatrixT& F = F_total();
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
	energy = fPot[0]->Energy(fEigs_dev, J) /*eq.(14)TDN*/;
    
    double precoeff=0.0;
    for(int k=0; k<fNumR; k++)
	{
		double phi=fdalpha[k];
        precoeff += (fTfk[k]*phi);
    }
    precoeff /= fT0;
    
    energy *= precoeff;
	/*adds nonequilibrium part */
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
    
	for (int i = 0; i < fNumS; i++)
	{
		/*calculate be*/
		fInverse.Inverse(fC_v[i]);
		fbe.MultQBQT(fF3D,fInverse);
        
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);
		fEigs_e = fSpectralDecompSpat.Eigenvalues();
        
		double Je = sqrt(fEigs_e.Product());
		fEigs_dev = fEigs_e;
		fEigs_dev *= pow(Je,-2.0*third);
		
		energy += fdmu[i]*(fPot[1]->Energy(fEigs_dev, Je));/*eq.(16)TDN*/
	}
	return(energy);
}

/* stresses */
const dSymMatrixT& ET_multi::s_ij(void)
{
    /*load the viscoelastic principal stretches from state variable arraysand calculate NEQ part*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	if (fFSMatSupport->RunState() == GlobalT::kFormRHS)       /*?*/
	{
		/*compute evolution equations for updated effective temperature and viscous deformation tensor*/
		Compute_le(fC_vn, fC_v, fTfk_n, fTfk, *fHeat);
		Store(element, CurrIP());
	}

   	const dMatrixT& F = F_total();
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
	
    double precoeff=0.0;
    for(int k=0; k<fNumR; k++)
	{
		double phi=fdalpha[k];
        precoeff += (fTfk[k]*phi);
    }
    precoeff /= fT0;
    ftau_EQ *= precoeff;
    ftau_EQ += fPot[0]->MeanStress(J);            /*eq.(17)TDN,EQ+p*/
	fStress3D = fSpectralDecompSpat.EigsToRank2(ftau_EQ);
    
/*    if(CurrElementNumber()==0&&CurrIP()==1)
        cout<<setprecision(12)<< "\nEQ stress: "<<fStress3D;
*/
	/*calc NEQ component of stress and moduli*/
	for (int i = 0; i < fNumS; i++)
	{
		/*calc elastic stretch*/
		fInverse.Inverse(fC_v[i]);
		fbe.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);
		fEigs_e = fSpectralDecompSpat.Eigenvalues();
        /*		if(CurrElementNumber()==0&&CurrIP()==0)
         cout << "\nfEigs_e2: "<<fEigs_e;
         */
		double Je = sqrt(fEigs_e.Product());
		fEigs_dev = fEigs_e;
		fEigs_dev *= pow(Je,-2.0*third);
		
		fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
		ftau_NEQ *=fdmu[i];
        /*		if(CurrElementNumber()==0&&CurrIP()==0)
         cout << "\nftau_NEQ: "<<ftau_NEQ;
         */
		fStress3D += fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
	}
/*    if(CurrElementNumber()==0&&CurrIP()==1)
        cout<< "\nTot stress: "<<fStress3D;
*/
	if (NumSD() == 2)
    {
        fStress[0] = fStress3D[0];              /*numbering?*/
        fStress[1] = fStress3D[1];
        fStress[2] = fStress3D[5];
    }
    else fStress = fStress3D;
    /*	if(CurrElementNumber()==0&&CurrIP()==0)
     cout << "\nfStress: "<<fStress;
     */
	const dMatrixT& Ftotal = F_total();
    fStress *= 1.0/Ftotal.Det();
	return fStress;
}


/* modulus */
const dMatrixT& ET_multi::c_ijkl(void)
{
    /*calculate contribution from Tel and epsilone*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    /*Computes modulus using loaded values of Cvj and Tei*/
	Compute_Kneq(fDAB, fDAB2); /*careful since this alters ftau_EQ and DtauEQDe*/
    

	const dMatrixT& F = F_total();
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
//    cout << "\nF: "<<fF3D;
    
	/*calcualte total stretch*/
    fb.MultAAT(fF3D);
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);
    fEigs = fSpectralDecompSpat.Eigenvalues();
    const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();
    
	/*calc EQ component of stress and moduli*/
    double J = sqrt(fEigs.Product());
    fEigs_dev = fEigs;
    fEigs_dev *= pow(J, -2.0*third);
	
    double precoeff=0.0;
    for(int k=0; k<fNumR; k++)
	{
		double phi=fdalpha[k];
        precoeff += (fTfk[k]*phi);
    }
    precoeff /= fT0;
    
    fPot[0]->DevStress(fEigs_dev, ftau_EQ);
    double tau0 = ftau_EQ[0];
    double tau1 = ftau_EQ[1];
    double tau2 = ftau_EQ[2];
    
    ftau_EQ *= precoeff;
	ftau_EQ += fPot[0]->MeanStress(J);
    
	fPot[0]->DevMod(fEigs_dev,fDtauDe_EQ);
    fDtauDe_EQ *= precoeff;
    fDtauDe_EQ += fPot[0]->MeanMod(J);
    
	dSymMatrixT& Gamma = fDtauDe_EQ;
    Gamma(0,0) -= 2.0*ftau_EQ[0];
    Gamma(1,1) -= 2.0*ftau_EQ[1];
    Gamma(2,2) -= 2.0*ftau_EQ[2];
    
    for (int i = 0; i < fNumR; i++)
    {
        double phi=fdalpha[i];
        double kneq10 = tau0/fT0*phi;
        double kneq11 = tau1/fT0*phi;
        double kneq12 = tau2/fT0*phi;
        
        Gamma(0,0) += kneq10*fDAB(i,0);
        Gamma(0,1) += kneq10*fDAB(i,1);
        Gamma(0,2) += kneq10*fDAB(i,2);
        Gamma(1,0) += kneq11*fDAB(i,0);
        Gamma(1,1) += kneq11*fDAB(i,1);
        Gamma(1,2) += kneq11*fDAB(i,2);
        Gamma(2,0) += kneq12*fDAB(i,0);
        Gamma(2,1) += kneq12*fDAB(i,1);
        Gamma(2,2) += kneq12*fDAB(i,2);
    }

    fModulus3D = fSpectralDecompSpat.EigsToRank4(Gamma);	/*nonlinear book p260, eq.6.191, line1,提出J^(-1)*/
	double dl, coeff;
    
    double l0 = fEigs[0];
    double l1 = fEigs[1];
    double l2 = fEigs[2];
	
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
    MixedRank4_3D(eigenvectors[1], eigenvectors[2], fModMat);   /*nonlinear book p260, eq.6.191, line2, eq.6.192*/
    fModulus3D.AddScaled(2.0*coeff, fModMat);
    
/*    if(CurrElementNumber()==0 && CurrIP()==1)
    {
        cout<< "\nfModulus3D: "<<fModulus3D;
    }
 */
    
    /*calculate NEQ*/
//    double* ple = fle.Pointer();
	for (int i = 0; i < fNumS; i++)
	{
        
        fInverse.Inverse(fC_v[i]);
		fbe.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);
		fEigs_e = fSpectralDecompSpat.Eigenvalues();
		double Je = sqrt(fEigs_e.Product());
		fEigs_dev = fEigs_e;
		fEigs_dev *= pow(Je,-2.0*third);

		fInverse.Inverse(fC_vn[i]);
		fb_tr.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_tr, false);
		fEigs_tr = fSpectralDecompSpat.Eigenvalues();
		const ArrayT<dArrayT>& eigenvectors_tr=fSpectralDecompSpat.Eigenvectors();

        /*stresses and moduli*/
		double mu = fdmu[i];
		fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
		ftau_NEQ *= mu;
		
		fPot[1]->DevMod(fEigs_dev, fDtauDe_NEQ);
		fDtauDe_NEQ *=mu;
		
		double a = 3*i;
		fMat(0,0) = fDAB2(a,0);
		fMat(0,1) = fDAB2(a,1);
		fMat(0,2) = fDAB2(a,2);
		fMat(1,0) = fDAB2(a+1,0);
		fMat(1,1) = fDAB2(a+1,1);
		fMat(1,2) = fDAB2(a+1,2);
		fMat(2,0) = fDAB2(a+2,0);
		fMat(2,1) = fDAB2(a+2,1);
		fMat(2,2) = fDAB2(a+2,2);
		
		fCalg(0,0) = fDtauDe_NEQ(0,0)*fMat(0,0) + fDtauDe_NEQ(0,1)*fMat(1,0) + fDtauDe_NEQ(0,2)*fMat(2,0);
		fCalg(0,1) = fDtauDe_NEQ(0,0)*fMat(0,1) + fDtauDe_NEQ(0,1)*fMat(1,1) + fDtauDe_NEQ(0,2)*fMat(2,1);
		fCalg(0,2) = fDtauDe_NEQ(0,0)*fMat(0,2) + fDtauDe_NEQ(0,1)*fMat(1,2) + fDtauDe_NEQ(0,2)*fMat(2,2);
        
		fCalg(1,0) = fDtauDe_NEQ(1,0)*fMat(0,0) + fDtauDe_NEQ(1,1)*fMat(1,0) + fDtauDe_NEQ(1,2)*fMat(2,0);
		fCalg(1,1) = fDtauDe_NEQ(1,0)*fMat(0,1) + fDtauDe_NEQ(1,1)*fMat(1,1) + fDtauDe_NEQ(1,2)*fMat(2,1);
		fCalg(1,2) = fDtauDe_NEQ(1,0)*fMat(0,2) + fDtauDe_NEQ(1,1)*fMat(1,2) + fDtauDe_NEQ(1,2)*fMat(2,2);
		
		fCalg(2,0) = fDtauDe_NEQ(2,0)*fMat(0,0) + fDtauDe_NEQ(2,1)*fMat(1,0) + fDtauDe_NEQ(2,2)*fMat(2,0);
		fCalg(2,1) = fDtauDe_NEQ(2,0)*fMat(0,1) + fDtauDe_NEQ(2,1)*fMat(1,1) + fDtauDe_NEQ(2,2)*fMat(2,1);
		fCalg(2,2) = fDtauDe_NEQ(2,0)*fMat(0,2) + fDtauDe_NEQ(2,1)*fMat(1,2) + fDtauDe_NEQ(2,2)*fMat(2,2);
        
		fCalg(0,0) -= 2.0*ftau_NEQ[0];
		fCalg(1,1) -= 2.0*ftau_NEQ[1];
		fCalg(2,2) -= 2.0*ftau_NEQ[2];
	
		fModulus3D += fSpectralDecompSpat.NonSymEigsToRank4(fCalg);
        
        
		double dl_tr;
        
		double l0_tr = fEigs_tr[0];
		double l1_tr = fEigs_tr[1];
		double l2_tr = fEigs_tr[2];
        
        
		dl_tr = l0_tr - l1_tr;
		if (fabs(dl_tr) > kSmall)
			coeff = (ftau_NEQ[0]*l1_tr -ftau_NEQ[1]*l0_tr)/dl_tr;
		else
			coeff = 0.5*(fCalg(0,0)-fCalg(0,1))-ftau_NEQ[0];
		MixedRank4_3D(eigenvectors_tr[0], eigenvectors_tr[1], fModMat);
		fModulus3D.AddScaled(2.0*coeff, fModMat);
        
		dl_tr = l0_tr - l2_tr;
		if (fabs(dl_tr) > kSmall)
			coeff =(ftau_NEQ[0]*l2_tr -ftau_NEQ[2]*l0_tr)/dl_tr;
		else
			coeff = 0.5*(fCalg(0,0)-fCalg(0,2))-ftau_NEQ[2];
		MixedRank4_3D(eigenvectors_tr[0], eigenvectors_tr[2], fModMat);
		fModulus3D.AddScaled(2.0*coeff, fModMat);
        
		dl_tr = l1_tr - l2_tr;
		if (fabs(dl_tr) > kSmall)
			coeff  = (ftau_NEQ[1]*l2_tr - ftau_NEQ[2]*l1_tr)/dl_tr;
		else
			coeff = 0.5*(fCalg(1,1)-fCalg(1,2))-ftau_NEQ[1];
		MixedRank4_3D(eigenvectors_tr[1], eigenvectors_tr[2], fModMat);
		fModulus3D.AddScaled(2.0*coeff, fModMat);

/*        if(CurrElementNumber()==0 && CurrIP()==1)
        {
            cout<< "\nfModulustot: "<<fModulus3D;
        }
 */
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

void ET_multi::ComputeOutput(dArrayT& output)
{
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
    
    fb.MultAAT(fF3D);
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);
    fEigs = fSpectralDecompSpat.Eigenvalues();
    double l0 = fEigs[0];
    double  l1 = fEigs[1];
    double l2 = fEigs[2];

	double J = sqrt(fEigs.Product());
    fEigs_dev = fEigs;
	fEigs_dev *= pow(J,-2.0*third);
    
    fPot[0]->DevStress(fEigs_dev, ftau_EQ);

    const dMatrixT& F_last = F_total_last();
    if (NumSD() == 2)
    {
        fF3D_last[0] = F_last[0];
        fF3D_last[1] = F_last[1];
        fF3D_last[2] = 0.0;
	    
        fF3D_last[3] = F_last[2];
        fF3D_last[4] = F_last[3];
        fF3D_last[5] = 0.0;
	    
        fF3D_last[6] = 0.0;
        fF3D_last[7] = 0.0;
        fF3D_last[8] = 1.0;
    }
    else fF3D_last = F_last;
    
    fb_last.MultAAT(fF3D_last);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_last, false);
	fEigs_last = fSpectralDecompSpat.Eigenvalues();
    double l0_last = fEigs_last[0];
    double l1_last = fEigs_last[1];
    double l2_last = fEigs_last[2];
    
/*    double depsilon0=(sqrt(l0)-sqrt(l0_last))/sqrt(l0);
    double depsilon1=(sqrt(l1)-sqrt(l1_last))/sqrt(l1);
    double depsilon2=(sqrt(l2)-sqrt(l2_last))/sqrt(l2);
 */
    double depsilon0=(0.5*log(l0)-0.5*log(l0_last));
    double depsilon1=(0.5*log(l1)-0.5*log(l1_last));
    double depsilon2=(0.5*log(l2)-0.5*log(l2_last));

    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
        
    /*calculate incremental heat*/
    const double rho = Density();
    double structural=0.0;
    double coeff_latentheat=0.0;
    for (int k = 0; k < fNumS; k++)
    {
        structural -= rho*fdeltac*fdalpha[k]*(fTfk[k]-fTfk_n[k]);
        coeff_latentheat += (fTfk[k]/fT0)*fdalpha[k];
    }
    /*add contributions from structural relaxation and latent heat*/
    double latentheat = coeff_latentheat*(ftau_EQ[0]*depsilon0+ftau_EQ[1]*depsilon1+ftau_EQ[2]*depsilon2);
    
    double heat = structural + latentheat;
    
    double plasticw =0.0;
    double stress1 = ftau_EQ[0];
    double stress2 = ftau_EQ[1];
    double stress3 = ftau_EQ[2];
    for(int k=0; k<fNumS; k++)
	{
        fInverse.Inverse(fC_v[k]);
		fbe.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);
		fEigs_e = fSpectralDecompSpat.Eigenvalues();
        double epse0 = 0.5*log(fEigs_e[0]);
        double epse1 = 0.5*log(fEigs_e[1]);
        double epse2 = 0.5*log(fEigs_e[2]);
        
 		double Je = sqrt(fEigs_e.Product());
		fEigs_dev = fEigs_e;
		fEigs_dev *= pow(Je,-2.0*third);
		fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
		ftau_NEQ *=fdmu[k];
        
        /*calc trial elastic stretch*/
        fInverse.Inverse(fC_vn[k]);
        fb_tr.MultQBQT(fF3D, fInverse);
        fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_tr, false);
        fEigs_tr = fSpectralDecompSpat.Eigenvalues();
        double epstr0 = 0.5*log(fEigs_tr[0]);
        double epstr1 = 0.5*log(fEigs_tr[1]);
        double epstr2 = 0.5*log(fEigs_tr[2]);
        
        double plasticwk = -(ftau_NEQ[0]*(epse0-epstr0)+ ftau_NEQ[1]*(epse1-epstr1) + ftau_NEQ[2]*(epse2-epstr2));
        plasticw += plasticwk;
        
        stress1 += ftau_NEQ[0];
        stress2 += ftau_NEQ[1];
        stress3 += ftau_NEQ[2];
	}
    
    stress1+=fPot[0]->MeanStress(J);
    stress2+=fPot[0]->MeanStress(J);
    stress3+=fPot[0]->MeanStress(J);

    /*add contribution from plastic work*/
    heat += plasticw;
    const double Temp = Compute_Temperature();
    const double Temp_last = Compute_Temperature_last();

    double dt = fFSMatSupport->TimeStep();

    output[0] = J;
    output[1] = structural/dt;
    output[2] = latentheat/dt;
    output[3] = plasticw/dt;
    output[4] = heat/dt;
    output[5] = rho*fcg*(Temp-Temp_last)/(stress1*depsilon0 + stress2*depsilon1+stress3*depsilon2);
}
/*************************************************************************
 *	PUBLIC
 **************************************************************************/
/* describe the parameters needed by the interface */
void ET_multi::DefineParameters(ParameterListT& list) const
{
    /* inherited */
    RGViscoelasticityT::DefineParameters(list);
    
    /* common limit */
    LimitT positive(0.0, LimitT::LowerInclusive);
    LimitT zero(0.0, LimitT::Lower);
    
    ParameterT inittemp(ParameterT::Double, "initial_temperature_T0");
    ParameterT glasstemp(ParameterT::Double, "glass_trans_temp_Tg");
    ParameterT C1(ParameterT::Double, "WLF_C1");
    ParameterT C2(ParameterT::Double, "WLF_C2");
    
    ParameterT afrac(ParameterT::Double, "fraction_a");
    ParameterT bdecay(ParameterT::Double, "decay_param_b");
    ParameterT Tess(ParameterT::Double, "steady_state_Tess");
    ParameterT BB(ParameterT::Double, "B");
    ParameterT deltac(ParameterT::Double, "delta_c");
    ParameterT cg(ParameterT::Double, "c_g");     /*Glassy heat capacity for computing the heating deficit*/
    ParameterT T2(ParameterT::Double, "T2");
    
    inittemp.AddLimit(positive);
    glasstemp.AddLimit(positive);
    C1.AddLimit(positive);
    C2.AddLimit(positive);
    
    afrac.AddLimit(positive);
    bdecay.AddLimit(positive);
    Tess.AddLimit(positive);
    BB.AddLimit(positive);
    deltac.AddLimit(positive);
    cg.AddLimit(positive);
    T2.AddLimit(positive);
    
    list.AddParameter(inittemp);
    list.AddParameter(glasstemp);
    list.AddParameter(C1);
    list.AddParameter(C2);
    
    list.AddParameter(afrac);
    list.AddParameter(bdecay);
    list.AddParameter(Tess);
    list.AddParameter(BB);
    list.AddParameter(deltac);
    list.AddParameter(cg);
    list.AddParameter(T2);
}

/* information about subordinate parameter lists */
void ET_multi::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);
    
	sub_list.AddSub("ET_multi_eq_potential", ParameterListT::Once);
	sub_list.AddSub("ET_multi_neq_potential", ParameterListT::Once);
    
	/* choice of viscosity */
	sub_list.AddSub("ET_multi_structural_spectrum", ParameterListT::Once);
	sub_list.AddSub("ET_multi_stress_spectrum", ParameterListT::Once);
}


/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ET_multi::NewSub(const StringT& name) const
{
    LimitT zero(0.0, LimitT::Lower);
	LimitT one(1.0, LimitT::Upper);
	LimitT two(2, LimitT::LowerInclusive);
	LimitT positive(0.0, LimitT::LowerInclusive);
    
	PotentialT* pot = NULL;
	if (name == "neo-hookean")
		pot = new NeoHookean;
	else if (name == "arruda-boyce")
		pot = new ArrudaBoyce;
	if (pot)
		return pot;
    
	/* inherited */
	ParameterInterfaceT* sub = RGViscoelasticityT::NewSub(name);
	if (sub)
	{
		return sub;
	}
	else if (name == "ET_multi_eq_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		choice->SetDescription("temperature normalized network stiffness");
		
		/* choice of parameters */
		choice->AddSub("arruda-boyce");
		choice->AddSub("neo-hookean");
		return(choice);
	}
	else if (name == "ET_multi_neq_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
        
		/* choice of parameters */
		choice->AddSub("neo-hookean");
		return(choice);
	}
	else if (name == "ET_multi_structural_spectrum")
	{
		ParameterContainerT* tauR = new ParameterContainerT(name);
		
		ParameterT file(ParameterT::String, "user_input_file");
        
		tauR->AddParameter(file);
        /*
         ParameterT tauRg(ParameterT::Double, "characterisic_relaxation_time_tauRref");
         ParameterT betaR(ParameterT::Double, "stretch_exponential_exponent_beta");
         ParameterT numR(ParameterT::Integer, "number_discrete_relaxation_times");
         ParameterT tauRmin(ParameterT::Double, "min_relaxation_time_spectrum");
         ParameterT tauRmax(ParameterT::Double, "max_relaxation_time_spectrum");
         
         tauRg.AddLimit(zero);
         tauRmin.AddLimit(zero);
         tauRmax.AddLimit(zero);
         betaR.AddLimit(zero);
         betaR.AddLimit(one);
         numR.AddLimit(two);
         
         tauR->AddParameter(tauRg);
         tauR->AddParameter(tauRmin);
         tauR->AddParameter(tauRmax);
         tauR->AddParameter(betaR);
         tauR->AddParameter(numR);
         */
        
		return(tauR);
	}
	else if (name == "ET_multi_stress_spectrum")
	{
		ParameterContainerT* tauS = new ParameterContainerT(name);
        
		ParameterT file(ParameterT::String, "user_input_file");
		tauS->AddParameter(file);
		
		ParameterT A(ParameterT::Double, "activation_volume");
		ParameterT tauY(ParameterT::Double, "characteristic_relaxation_time_tauYref");
        
		A.AddLimit(zero);
        
		tauS->AddParameter(A);
        
        /*
         ParameterT tauSg(ParameterT::Double, "characterisic_relaxation_time_tauSref");
         ParameterT betaS(ParameterT::Double, "stretch_exponential_exponent_beta");
         ParameterT numS(ParameterT::Integer, "number_discrete_relaxation_times");
         ParameterT tauSmin(ParameterT::Double, "min_relaxation_time_spectrum");
         ParameterT tauSmax(ParameterT::Double, "max_relaxation_time_spectrum");
         
         tauSg.AddLimit(zero);
         tauSmin.AddLimit(zero);
         tauSmax.AddLimit(zero);
         betaS.AddLimit(zero);
         betaS.AddLimit(one);
         numS.AddLimit(two);
         
         tauS->AddParameter(tauSg);
         tauS->AddParameter(tauSmin);
         tauS->AddParameter(tauSmax);
         tauS->AddParameter(betaS);
         tauS->AddParameter(numS);
         */
        
		return(tauS);
	}
}

void ET_multi::TakeParameterList(const ParameterListT& list)
{
    const char caller[] = "ET_multi::TakeParameterList";
    /* inherited */
    RGViscoelasticityT::TakeParameterList(list);
    
    
    fT0 = list.GetParameter("initial_temperature_T0");
    fTg = list.GetParameter("glass_trans_temp_Tg");
    fC1 = list.GetParameter("WLF_C1");
    fC2 = list.GetParameter("WLF_C2");
    fafrac = list.GetParameter("fraction_a");
    fb_decay = list.GetParameter("decay_param_b");
    fTess = list.GetParameter("steady_state_Tess");
    fBB = list.GetParameter("B");
    fdeltac = list.GetParameter("delta_c");
    fcg = list.GetParameter("c_g");
    fT2 = list.GetParameter("T2");
    
	fPot.Dimension(2);
    
	const ParameterListT& eq_pot = list.GetListChoice(*this, "ET_multi_eq_potential");
	if(eq_pot.Name() == "arruda-boyce")
		fPot[0] = new ArrudaBoyce;
	else if(eq_pot.Name() == "neo-hookean")
		fPot[0] = new NeoHookean;
	else
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot[0]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", eq_pot.Name().Pointer());
	fPot[0]->TakeParameterList(eq_pot);
	/*set the rubbery modulus*/
	fmur = fPot[0]->GetMu();
    
	const ParameterListT& neq_pot = list.GetListChoice(*this, "ET_multi_neq_potential");
	if(neq_pot.Name() == "neo-hookean")
		fPot[1] = new NeoHookean;
	else
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot[1]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", neq_pot.Name().Pointer());
	fPot[1]->TakeParameterList(neq_pot);
  	/*set the glassy modulus*/
	fmug = fmur + fPot[1]->GetMu();
	/*normalize the neq potential function by setting the shear modulus to one*/
	fPot[1]->SetMu(1.0);
    const ParameterListT* tauR = list.List("ET_multi_structural_spectrum");
    if (tauR)
    {
        fInputR = tauR->GetParameter("user_input_file");
        
        /*	ftauR0 = tauR->GetParameter("characterisic_relaxation_time_tauRref");
         fbetaR = tauR->GetParameter("stretch_exponential_exponent_beta");
         fNumR = tauR->GetParameter("number_discrete_relaxation_times");
         ftauR0min = tauR->GetParameter("min_relaxation_time_spectrum");
         ftauR0max = tauR->GetParameter("max_relaxation_time_spectrum");
         */
    }
    
    const ParameterListT* tauS = list.List("ET_multi_stress_spectrum");
    if (tauS)
    {
        /*	ftauS0 = tauS->GetParameter("characterisic_relaxation_time_tauSref");
         fbetaS = tauS->GetParameter("stretch_exponential_exponent_beta");
         fNumS = tauS->GetParameter("number_discrete_relaxation_times");
         ftauS0min = tauS->GetParameter("min_relaxation_time_spectrum");
         ftauS0max = tauS->GetParameter("max_relaxation_time_spectrum");
         */
        fInputS = tauS->GetParameter("user_input_file");
        
        fQS = tauS->GetParameter("activation_volume");
    }
    
	/*read in structural relaxation spectrum*/
  	ifstreamT inR;
	inR.open(fInputR);
	if(!inR.good())
		ExceptionT::DatabaseFail(caller,
                                 "could not open file \"%s\"",fInputR.Pointer());
	inR >> fNumR;
	ftimesR.Dimension(fNumR);
	fdalpha.Dimension(fNumR);
	for (int i=0; i<fNumR && inR.good(); i++)
	{
		inR >> ftimesR[i];
		inR >> fdalpha[i];
	}
	inR.close();
    
	/*read in stress relaxation spectrum*/
  	ifstreamT inS;
	inS.open(fInputS);
	if(!inS.good())
		ExceptionT::DatabaseFail(caller,
                                 "could not open file \"%s\"",fInputS.Pointer());
	inS >> fNumS;
	ftimesS.Dimension(fNumS);
	fdmu.Dimension(fNumS);
	for (int i=0; i<fNumS && inS.good(); i++)
	{
		inS >> ftimesS[i];
		inS >> fdmu[i];
	}
	inS.close();
 	
	/*Dimension workspace*/
	/*Dimension state variable accessors*/
	fC_v.Dimension(fNumS);
	fC_vn.Dimension(fNumS);
	fl_tr.Dimension(fNumS*3);
	fle.Dimension(fNumS*3);
	fstressk.Dimension(fNumS*3);
	
	int nsd = NumSD();
	int ndof = 3;
	int numstress = dSymMatrixT::NumValues(ndof);
    
    /*Set pointer to internal state variables*/
	fnstatev = 0;
	fnstatev += numstress*fNumS;   /*current C_v*/
	fnstatev += numstress*fNumS;   /*last C_vn*/
	fnstatev += fNumR;			/*current  Tfk*/
	fnstatev += fNumR;			/*last Tfk*/
    fnstatev ++;                /*fHeat incremental heat calculated from dotCv and dotTfk*/
    fnstatev ++;                /*fHeat_n last incremental heat calculated from dotCv and dotTfk*/
	
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
	
	/* assign pointers to current and last blocks of state variable array */
	for (int i = 0; i < fNumS; i++)
	{
		fC_v[i].Set(ndof, pstatev);
		pstatev += numstress;
		fC_vn[i].Set(ndof, pstatev);
		pstatev += numstress;
	}
    
	fTfk.Alias(fNumR,pstatev);
	pstatev += fNumR;
	fTfk_n.Alias(fNumR,pstatev);
	pstatev += fNumR;
    
	fHeat = pstatev;
	pstatev++;
    fHeat_n = pstatev;
    pstatev++;

	fF_M.Dimension(nsd);
	fF_T_inv.Dimension(nsd);
    
	fF3D.Dimension(ndof);
	fInverse.Dimension(ndof);
    fF3D_last.Dimension(ndof);
    //fF3D_last=1.0;
    
	fb.Dimension(ndof);
	fbe.Dimension(ndof);
	fb_tr.Dimension(ndof);
    fb_last.Dimension(ndof);
    //fb_last=1.0;
	
	fEigs_dev.Dimension(ndof);
	fEigs.Dimension(ndof);
	fEigs_e.Dimension(ndof);
	fEigs_tr.Dimension(ndof);
	fEigs_last.Dimension(ndof);
    //fEigs_last=1.0;
	fEigs_e_last.Dimension(ndof);
    //fEigs_e_last=1.0;
    
	ftau_EQ.Dimension(ndof);
	ftau_NEQ.Dimension(ndof);
	ftau.Dimension(ndof);
    
	fStress.Dimension(NumSD());
	fStress3D.Dimension(ndof);
    
	fDtauDe_EQ.Dimension(ndof);
	fDtauDe_NEQ.Dimension(ndof);
	
	fKdel.Dimension(fNumR);
	fRdel.Dimension(fNumR);
    //fRdel=0.0;
    flatenth.Dimension(fNumR);
    //flatenth=0.0;
    fGB0.Dimension(fNumR);
    //fGB0=0.0;
    fGB1.Dimension(fNumR);
    //fGB1=0.0;
    fGB2.Dimension(fNumR);
    //fGB2=0.0;
    fG10.Dimension(fNumR);
    //fG10=0.0;
    fG11.Dimension(fNumR);
    //fG11=0.0;
    fG12.Dimension(fNumR);
    //fG12=0.0;
    
	fKAB.Dimension(fNumS*3);
	fKAB2.Dimension(fNumS*3);
	fRes.Dimension(fNumS*3);
    //fRes=0.0;
    
    fGA0.Dimension(fNumS*3);
    //fGA0=0.0;
    fGA1.Dimension(fNumS*3);
    //fGA1=0.0;
    fGA2.Dimension(fNumS*3);
    //fGA2=0.0;
    fG20.Dimension(fNumS*3);
    //fG20=0.0;
    fG21.Dimension(fNumS*3);
    //fG21=0.0;
    fG22.Dimension(fNumS*3);
    //fG22=0.0;
	
	fCalg.Dimension(3);
	fMat.Dimension(3);
    //fMat=0.0;
	fDAB.Dimension(fNumS,3);
    //fDAB=0.0;
	fDAB2.Dimension(fNumS*3,3);
    //fDAB2=0.0;
    
    fKrTe.Dimension(3*fNumS,fNumR);
    fKTer.Dimension(fNumR,3*fNumS);
    fInverse1.Dimension(fNumR,3*fNumS);
    fInverse2.Dimension(fNumR,3*fNumS);
    fInverse3.Dimension(3*fNumS,fNumR);
    fKdel2.Dimension(fNumR,fNumR);
    fK1.Dimension(fNumR,fNumR);
    fK2.Dimension(3*fNumS,3*fNumS);
    
    fModulus3D.Dimension(dSymMatrixT::NumValues(ndof));
    fModMat.Dimension(dSymMatrixT::NumValues(ndof));
    fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));
    
}

/*initializes history variable */
void  ET_multi::PointInitialize(void)
{
	/* allocate element storage */
    const double Temp = Compute_Temperature();
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
            for (int k = 0; k < fNumS; k++)
            {
				fC_v[k].Identity();
				fC_vn[k].Identity();
            }
            fTfk = Temp;
            fTfk_n = Temp;
            *fHeat = 0.0;
            *fHeat_n = 0.0;
            
            /* write to storage */
            Store(element, ip);
		}
	}
}

void ET_multi::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables */
		Load(element, ip);
        
		/* assign "current" to "last" */
		for (int k = 0; k < fNumS; k++)
			fC_vn[k] = fC_v[k];
		fTfk_n = fTfk;
        *fHeat_n = *fHeat;
		
		/* write to storage */
		Store(element, ip);
	}
}

void ET_multi::ResetHistory(void)
{
    /* current element */
	ElementCardT& element = CurrentElement();
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables*/
		Load(element, ip);
        
		/* assign "last" to "current" */
		for (int k = 0; k < fNumS; k++)
			fC_v[k] = fC_vn[k];
		fTfk = fTfk_n;
        *fHeat = *fHeat_n;

		/* write to storage */
		Store(element, ip);
	}
}

void ET_multi::Compute_Tei(const dArrayT& eigs_n, const dArrayT& eigs, const dArrayT& eigs_tr, const dArrayT& DScDe, const dArrayT& stretche,const dArrayT& stressk, const dArrayT& Tfk_n, dArrayT& Tfk)
{
    double ctol = 1.0e-10;
	/*time step*/
	double dt = fFSMatSupport->TimeStep();
	/*current temperature*/
	double Tn = Compute_Temperature();
    double energy=0.0;
    double plasticw = 0.0;
    double frac1 = 1.0/fdeltac;
   
    double l0 = eigs[0];
    double l1 = eigs[1];
    double l2 = eigs[2];
    double J = sqrt(l0*l1*l2);
    double iJ = 1.0/J;
//    iJ = 1.0;
    
    double l0_last = eigs_n[0];
    double l1_last = eigs_n[1];
    double l2_last = eigs_n[2];
    
    /*    double depsilon0=(sqrt(l0)-sqrt(l0_last))/sqrt(l0);
     double depsilon1=(sqrt(l1)-sqrt(l1_last))/sqrt(l1);
     double depsilon2=(sqrt(l2)-sqrt(l2_last))/sqrt(l2);
     */
    double depsilon0=(0.5*log(l0)-0.5*log(l0_last));
    double depsilon1=(0.5*log(l1)-0.5*log(l1_last));
    double depsilon2=(0.5*log(l2)-0.5*log(l2_last));
	
    double prefactor = 1.0;
	double temp = 0.0;
	double ScTn = 0.0;
	for(int k=0; k<fNumR; k++)
	{
		double phi=fdalpha[k];
        ScTn += fdeltac*phi*log(Tfk[k]/fT2);
        prefactor += phi*fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(Tfk[k]/fTess - 1.0))))*(Tn/Tfk[k]-1.0);
        flatenth[k] =  Tfk[k]/fT0*(DScDe[0]*depsilon0+DScDe[1]*depsilon1+DScDe[2]*depsilon2);
    }
    /*calculate smag and plastic power*/
    const double* pstress = stressk.Pointer();
    const double* ple = stretche.Pointer();
    const double* pltr = eigs_tr.Pointer();
    plasticw =0.0;
    double tau0 = 0.0;
    double tau1 = 0.0;
    double tau2 = 0.0;
	for(int k=0; k<fNumS; k++)
	{
        double tauk0 = *pstress++;
        double tauk1 = *pstress++;
        double tauk2 = *pstress++;
        
        tau0 += tauk0;
        tau1 += tauk1;
        tau2 += tauk2;
        
        double epse0 = 0.5*log(*ple++);
        double epse1 = 0.5*log(*ple++);
        double epse2 = 0.5*log(*ple++);
        
        double epstr0 = 0.5*log(*pltr++);
        double epstr1 = 0.5*log(*pltr++);
        double epstr2 = 0.5*log(*pltr++);
        
		/*viscosity and its derivatives*/
        double plasticwk = -(tauk0*(epse0-epstr0)+ tauk1*(epse1-epstr1) + tauk2*(epse2-epstr2));
        plasticw += plasticwk;
	}
    double smag = prefactor*iJ*sqrt(0.5*(tau0*tau0+tau1*tau1+tau2*tau2));
    
	fEigs_dev = eigs;
	fEigs_dev *= pow(J, -2.0*third);
    energy = fPot[0]->Energy(fEigs_dev, J) /*eq.(14)TDN*/;
    energy -= fPot[0]->MeanEnergy(J);
    ScTn -= energy/fT0;
    
    double itaubar = StructuralRelaxationFunc(Tn, ScTn);

    for(int k=0; k<fNumR; k++)
	{
		double itauRk = itaubar/ftimesR[k];
        double frac2 = fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(Tfk[k]/fTess - 1.0))))/fdeltac;
        /*residual*/
		fRdel[k] = Tfk[k] - Tfk_n[k] + dt*itauRk*(Tfk[k]-Tn)-frac1*flatenth[k]-frac2*plasticw;
		temp += fRdel[k]*fRdel[k];
/*        if(CurrElementNumber()==0 && CurrIP()==1)
       {
            cout<< "\ndepsilon: "<<depsilon0<<"\t"<<depsilon1<<"\t"<<depsilon2;
            cout<< "\nstress: "<<ftau_EQ[0]<<"\t"<<ftau_EQ[1]<<"\t"<<ftau_EQ[2];
            cout << "\nlatenth: "<<flatenth[k]<<"\tScTn: "<<ScTn<<"\titauR: "<<itauRk<<"\tTfk-Tn: "<<Tfk[k]-Tn;
           cout <<"\nitauRk: "<<itauRk;
        cout << "\ndt*itauRk*(Tfk[k]-Tn): "<<dt*itauRk*(Tfk[k]-Tn)<<"\tfrac1*flatenth[k]: "<<frac1*flatenth[k]<<"\tfrac2*plasticw: "<<frac2*plasticw<<"\tTfk-Tn: "<<Tfk[k]-Tn;
        }
 */
	}
    
	double tol = sqrt(temp);
/*    if(CurrElementNumber()==0 && CurrIP()==1)
    {
        cout<<setprecision(12)<< "\nTe 0: "<<tol;
  //      cout << "\nfEigs: "<<eigs;
 //       cout << "\nTfk: "<<Tfk;
    }
 */
    int maxiter=20;
    int iter = 0;
	//while (tol > ctol && reltol>ctol && iter < maxiter)
    while (tol > ctol  && iter < maxiter)
	{
        iter++;
		fKdel = 0.0;
		for(int k=0; k<fNumR; k++)
		{
			/*stiffness matrix*/
			double itauRk = itaubar/ftimesR[k];
            double frac2= fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(Tfk[k]/fTess - 1.0))))/fdeltac;
            double DaDTfk = -fb_decay/fTess/(1.0+exp(-fb_decay*(Tfk[k]/fTess - 1.0))); /*1/afrac * Dafrac/DTfk*/
            fKdel(k,k) = 1.0 + dt*itauRk-frac1*flatenth[k]/Tfk[k] - DaDTfk*frac2*plasticw ;
			for (int l=0; l<fNumR; l++)
            {
                double phi=fdalpha[l];
                double DtaubarDTfk = -fBB/Tn/ScTn/ScTn*(phi*fdeltac/Tfk[l]);  /*1/tauRi * DtauRi/DTfk*/
				fKdel(k,l) += -dt*itauRk*(Tfk[k]-Tn)*DtaubarDTfk;            }
		}
		/*Solve for update*/
		fKdel.LinearSolve(fRdel);
		Tfk -= fRdel;
		
		/*calculate residual*/
        ScTn = 0.0;
        prefactor = 1.0;
        for(int k=0; k<fNumR; k++)
        {
            double phi=fdalpha[k];
            ScTn += fdeltac*phi*log(Tfk[k]/fT2);
            prefactor += phi*fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(Tfk[k]/fTess - 1.0))))*(Tn/Tfk[k]-1.0);
            flatenth[k] = Tfk[k]/fT0*(DScDe[0]*depsilon0+DScDe[1]*depsilon1+DScDe[2]*depsilon2);
        }
        
        fEigs_dev = eigs;
        fEigs_dev *= pow(J, -2.0*third);
        double  energy = fPot[0]->Energy(fEigs_dev, J) /*eq.(14)TDN*/;
        energy -= fPot[0]->MeanEnergy(J);
        ScTn -= energy/fT0;
        
        /*update smag*/
        double smag = prefactor*iJ*sqrt(0.5*(tau0*tau0+tau1*tau1+tau2*tau2));
        itaubar = StructuralRelaxationFunc(Tn, ScTn);  /*1/Adam-Gibbs but without tauRref*/
        
        temp = 0.0;
        for(int k=0; k<fNumR; k++)
        {
            double itauRk = itaubar/ftimesR[k];
            double frac2 = fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(Tfk[k]/fTess - 1.0))))/fdeltac;
            /*residual*/
            fRdel[k] = Tfk[k] + dt*itauRk*(Tfk[k]-Tn) - Tfk_n[k]-frac1*flatenth[k]-frac2*plasticw;
            temp += fRdel[k]*fRdel[k];
 
/*            if(CurrElementNumber()==0 && CurrIP()==1)
            {
                cout<< "\ndepsilon: "<<depsilon0<<"\t"<<depsilon1<<"\t"<<depsilon2;
                cout<< "\nstress: "<<ftau_EQ[0]/fT0<<"\t"<<ftau_EQ[1]/fT0<<"\t"<<ftau_EQ[2]/fT0;
                cout << "\nlatenth: "<<flatenth[k]<<"\tScTn: "<<ScTn<<"\titauR: "<<itauRk<<"\tTfk-Tn: "<<Tfk[k]-Tn;
                cout <<"\nitauRk: "<<itauRk;
                cout << "\ndt*itauRk*(Tfk[k]-Tn): "<<dt*itauRk*(Tfk[k]-Tn)<<"\tfrac1*flatenth[k]: "<<frac1*flatenth[k]<<"\tfrac2*plasticw: "<<frac2*plasticw<<"\tTfk-Tn: "<<Tfk[k]-Tn;
            }
 */
       }
        
        tol = sqrt(temp);

/*        if(CurrElementNumber()==0 && CurrIP()==1)
       {
            cout<< "\nTe: "<<iter<<"\t"<<tol;
 //          cout << "\nTfk: "<<Tfk;
       }
 */
	}
    if (iter >= maxiter)
	{
		cout<<"\n Number of iteration exceeds maximum. tol:  "<<tol;
		cout<< "\nelem: "<<CurrElementNumber()<<"\nIP: "<<CurrIP();
		ExceptionT::GeneralFail("ET_multi::Compute_le", "number of iteration exceeds maximum");
	}
    
}

void ET_multi::Compute_Kneq(dMatrixT& Modulus1, dMatrixT& Modulus2)
{
    double frac1 = 1.0/fdeltac;
    double energy=0.0;
    double ScTn = 0.0;
    
	/*time step*/
	const double dt = fFSMatSupport->TimeStep();
    
	/*temperature and temperature step*/
	const double Tn = Compute_Temperature();
    
    const dMatrixT& F_last = F_total_last();
    
	if (NumSD() == 2)
	{
		fF3D_last[0] = F_last[0];
		fF3D_last[1] = F_last[1];
		fF3D_last[2] = 0.0;
	    
		fF3D_last[3] = F_last[2];
		fF3D_last[4] = F_last[3];
		fF3D_last[5] = 0.0;
	    
		fF3D_last[6] = 0.0;
		fF3D_last[7] = 0.0;
		fF3D_last[8] = 1.0;
	}
	else fF3D_last = F_last;
    
    fb_last.MultAAT(fF3D_last);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_last, false);
	fEigs_last = fSpectralDecompSpat.Eigenvalues();
    
    double l0_last = fEigs_last[0];
    double l1_last = fEigs_last[1];
    double l2_last = fEigs_last[2];
    
	const dMatrixT& F = F_total();
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
    double J = sqrt(fEigs.Product());
	double iJ = 1.0/J;
   // iJ =1.0;
    
	fEigs_dev = fEigs;
	fEigs_dev *= pow(J, -2.0*third);
    
    double l0 = fEigs[0];
    double l1 = fEigs[1];
    double l2 = fEigs[2];
    
	fPot[0]->DevStress(fEigs_dev, ftau_EQ);
    double l0f = -ftau_EQ[0]/fT0;
    double l1f = -ftau_EQ[1]/fT0;
    double l2f = -ftau_EQ[2]/fT0;
    
    fPot[0]->DevMod(fEigs_dev,fDtauDe_EQ);

    /*    double depsilon0=(sqrt(l0)-sqrt(l0_last))/sqrt(l0);
     double depsilon1=(sqrt(l1)-sqrt(l1_last))/sqrt(l1);
     double depsilon2=(sqrt(l2)-sqrt(l2_last))/sqrt(l2);
     */
    double depsilon0=(0.5*log(l0)-0.5*log(l0_last));
    double depsilon1=(0.5*log(l1)-0.5*log(l1_last));
    double depsilon2=(0.5*log(l2)-0.5*log(l2_last));
    
    double prefactor = 1.0;
    for(int k=0; k<fNumR; k++)
	{
		double phi=fdalpha[k];
        ScTn += fdeltac*phi*log(fTfk[k]/fT2);
        prefactor += phi*fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(fTfk[k]/fTess - 1.0))))*(Tn/fTfk[k]-1.0);
        flatenth[k] = fTfk[k]/fT0*(ftau_EQ[0]*depsilon0+ftau_EQ[1]*depsilon1+ftau_EQ[2]*depsilon2);
    }
    
    energy = fPot[0]->Energy(fEigs_dev, J) /*eq.(14)TDN*/;
    energy -= fPot[0]->MeanEnergy(J);
    ScTn -= energy/fT0;
    
    /*calculate  smag*/
    /*update  smag*/
    ftau= 0.0;
    double plasticw = 0.0;
    double* ple = fle.Pointer();
    double* pltr = fl_tr.Pointer();
    for (int i = 0; i < fNumS; i++)
    {
        /*calc trial elastic stretch*/
        fInverse.Inverse(fC_vn[i]);
        fb_tr.MultQBQT(fF3D, fInverse);
        fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_tr, false);
        fEigs_tr = fSpectralDecompSpat.Eigenvalues();
        *pltr++ = fEigs_tr[0];
        *pltr++ = fEigs_tr[1];
        *pltr++ = fEigs_tr[2];
        double epstr0 = 0.5*log(fEigs_tr[0]);
        double epstr1 = 0.5*log(fEigs_tr[1]);
        double epstr2 = 0.5*log(fEigs_tr[2]);
        
        /*calc elastic stretch*/
        /*calc trial elastic stretch*/
        fInverse.Inverse(fC_v[i]);
        fbe.MultQBQT(fF3D, fInverse);
        fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);
        fEigs_e = fSpectralDecompSpat.Eigenvalues();
        *ple++ = fEigs_e[0];
        *ple++ = fEigs_e[1];
        *ple++ = fEigs_e[2];
        double epse0 = 0.5*log(fEigs_e[0]);
        double epse1 = 0.5*log(fEigs_e[1]);
        double epse2 = 0.5*log(fEigs_e[2]);
 
 /*       fEigs_e[0] = *ple++;
        fEigs_e[1] = *ple++;
        fEigs_e[2] = *ple++;
        double epse0 = 0.5*log(fEigs_e[0]);
        double epse1 = 0.5*log(fEigs_e[1]);
        double epse2 = 0.5*log(fEigs_e[2]);
   */
   
   /*           if (CurrElementNumber()==0 && CurrIP()==1)
         {
         cout <<"\neigse"<< fEigs_e;
         cout <<"\neigsdev"<< fEigs_dev;
         cout << "\ntauneq: "<<ftau_NEQ;
         }
         */
        double Je = sqrt(fEigs_e.Product());
        fEigs_dev = fEigs_e;
        fEigs_dev *= pow(Je,-2.0*third);
            
        /*calculate total neq stress*/
        fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
        ftau_NEQ *=fdmu[i];
        
         ftau+= ftau_NEQ;
        double plasticwk = -(ftau_NEQ[0]*(epse0-epstr0)+ ftau_NEQ[1]*(epse1-epstr1) + ftau_NEQ[2]*(epse2-epstr2));
        plasticw += plasticwk;
    }
    double smag = prefactor*iJ*sqrt(0.5*(ftau[0]*ftau[0]+ftau[1]*ftau[1]+ftau[2]*ftau[2]));
 /*  if(CurrElementNumber()==0 && CurrIP()==1)
    {
       cout<< "\nEigse: "<<fle;
        cout << "\nsmag: "<<smag;
        cout << "\nlatenth: "<<flatenth;
    }
*/

	/*update viscosity*/
    double ietabar = StressRelaxationFunc(Tn, smag);   /*eq.(24)*/
    double itaubar = StructuralRelaxationFunc(Tn, ScTn);
	double DetabarDs = 0.0;
    if(smag > kSmall)
    {
        double x = (fQS*smag)/(Tn);
        double cothx = cosh(x)/sinh(x);
        DetabarDs = (1.0 - x*cothx)/smag;
    }
    
    double DetaSDScTn = -fBB/Tn/ScTn/ScTn; /* 1/etaS DetaSDScTn*/
    
    /*******************HERE*********************/
    
	fKAB = 0.0;    /*DRes/Depsilon^e_Bj*/
    fKTer = 0.0;   /*DRdel/Depsilon^e_Bj*/
    fKrTe = 0.0;   /*DRes/DTe_l*/
    fKdel=0.0;     /*DRel/DTe_l*/
    
	fGA0 = 0.0;    /*DRes/Depsilon_1*/
	fGA1 = 0.0;    /*DRes/Depsilon_2*/
	fGA2 = 0.0;    /*DRes/Depsilon_3*/
    
   	fGB0 = 0.0;    /*DRdel/Depsilon_B*/
	fGB1 = 0.0;    /*DRdel/Depsilon_B*/
	fGB2 = 0.0;    /*DRdel/Depsilon_B*/
    
    
	Modulus1 = 0.0;
    Modulus2 = 0.0;

	ple = fle.Pointer();
    pltr = fl_tr.Pointer();
	for(int k=0; k<fNumS; k++)
	{
		/*viscosity and its derivatives*/
        double muk=fdmu[k];
        double itauSk = itaubar/ftimesS[k];
        double ietaS0k = itauSk/muk;
        double ietaSk = ietaS0k*ietabar;
        
        double epstr0 = 0.5*log(*pltr++);
        double epstr1 = 0.5*log(*pltr++);
        double epstr2 = 0.5*log(*pltr++);
        
		fEigs_dev[0] = *ple++;
		fEigs_dev[1] = *ple++;
		fEigs_dev[2] = *ple++;
        double epse0 = 0.5*log(fEigs_dev[0]);
        double epse1 = 0.5*log(fEigs_dev[0]);
        double epse2 = 0.5*log(fEigs_dev[0]);

        double Je = sqrt(fEigs_dev[0]*fEigs_dev[1]*fEigs_dev[2]);
		fEigs_dev *= pow(Je,-2.0*third);
		
		fPot[1]->DevMod(fEigs_dev,fDtauDe_NEQ);
        fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
        ftau_NEQ *=muk;
        
		/*moduli*/
        double dk00 = muk*iJ*fDtauDe_NEQ(0,0);
        double dk11 = muk*iJ*fDtauDe_NEQ(1,1);
        double dk22 = muk*iJ*fDtauDe_NEQ(2,2);
        double dk12 = muk*iJ*fDtauDe_NEQ(1,2);
        double dk02 = muk*iJ*fDtauDe_NEQ(0,2);
        double dk01 = muk*iJ*fDtauDe_NEQ(0,1);
        
        double tauk0 = ftau_NEQ[0];
        double tauk1 = ftau_NEQ[1];
        double tauk2 = ftau_NEQ[2];
        
        double sk0 = tauk0*iJ*prefactor;
        double sk1 = tauk1*iJ*prefactor;
        double sk2 = tauk2*iJ*prefactor;
        
/*        if(CurrElementNumber()==0 && CurrIP()==0)
        {
            cout<< "\nietaSk: "<<ietaSk;
            cout << "\nftau_NEQ: "<<ftau_NEQ;
        }
 */
		/*stiffness matrix*/
		int a = k*3;
        fKAB(a,a) = 1.0 + 0.5*dt*ietaSk*dk00*prefactor;
        fKAB(a,a+1) = 0.5*dt*ietaSk*dk01*prefactor;
        fKAB(a,a+2) = 0.5*dt*ietaSk*dk02*prefactor;
        fKAB(a+1,a) = 0.5*dt*ietaSk*dk01*prefactor;
        fKAB(a+1,a+1) = 1.0 + 0.5*dt*ietaSk*dk11*prefactor;
        fKAB(a+1,a+2) = 0.5*dt*ietaSk*dk12*prefactor;
        fKAB(a+2,a) = 0.5*dt*ietaSk*dk02*prefactor;
        fKAB(a+2,a+1) = 0.5*dt*ietaSk*dk12*prefactor;
        fKAB(a+2,a+2) = 1.0 + 0.5*dt*ietaSk*dk22*prefactor;
        
        double *ple2 = fle.Pointer();
        for (int l=0; l<fNumS; l++)       /*evolution eqs. are coupled in this case because of s (or smag here)*/
        {
            /*moduli*/
            fEigs_dev[0] = *ple2++;
            fEigs_dev[1] = *ple2++;
            fEigs_dev[2] = *ple2++;
            double Je2 = sqrt(fEigs_dev[0]*fEigs_dev[1]*fEigs_dev[2]);
            fEigs_dev *= pow(Je2,-2.0*third);
            fPot[1]->DevMod(fEigs_dev,fDtauDe_NEQ);
            double dl00 = fdmu[l]*iJ*fDtauDe_NEQ(0,0);
            double dl11 = fdmu[l]*iJ*fDtauDe_NEQ(1,1);
            double dl22 = fdmu[l]*iJ*fDtauDe_NEQ(2,2);
            double dl12 = fdmu[l]*iJ*fDtauDe_NEQ(1,2);
            double dl02 = fdmu[l]*iJ*fDtauDe_NEQ(0,2);
            double dl01 = fdmu[l]*iJ*fDtauDe_NEQ(0,1);
            
            double cl0 = 0.0;
            double cl1 = 0.0;
            double cl2 = 0.0;
            
            if(smag >kSmall)
            {
                cl0 = 0.5*(iJ*ftau[0]*dl00 + iJ*ftau[1]*dl01 + iJ*ftau[2]*dl02)*prefactor*prefactor/smag;
                cl1 = 0.5*(iJ*ftau[0]*dl01 + iJ*ftau[1]*dl11 + iJ*ftau[2]*dl12)*prefactor*prefactor/smag;
                cl2 = 0.5*(iJ*ftau[0]*dl02 + iJ*ftau[1]*dl12 + iJ*ftau[2]*dl22)*prefactor*prefactor/smag;
            }
            
            int b = l*3;
            fKAB(a,b) += - 0.5*dt*ietaSk*(DetabarDs*cl0)*sk0;
            fKAB(a,b+1) += - 0.5*dt*ietaSk*(DetabarDs*cl1)*sk0;
            fKAB(a,b+2) += - 0.5*dt*ietaSk*(DetabarDs*cl2)*sk0;
            fKAB(a+1,b) += - 0.5*dt*ietaSk*(DetabarDs*cl0)*sk1;
            fKAB(a+1,b+1) += - 0.5*dt*ietaSk*(DetabarDs*cl1)*sk1;
            fKAB(a+1,b+2) += - 0.5*dt*ietaSk*(DetabarDs*cl2)*sk1;
            fKAB(a+2,b) +=- 0.5*dt*ietaSk*(DetabarDs*cl0)*sk2;
            fKAB(a+2,b+1) += - 0.5*dt*ietaSk*(DetabarDs*cl1)*sk2;
            fKAB(a+2,b+2) += - 0.5*dt*ietaSk*(DetabarDs*cl2)*sk2;
        }
        fGA0[a] = 1.0 + 0.5*dt*ietaSk*sk0*(DetaSDScTn*l0f - DetabarDs*smag + 1.0);
        fGA0[a+1] = 0.5*dt*ietaSk*sk1*(DetaSDScTn*l0f - DetabarDs*smag + 1.0);
        fGA0[a+2] = 0.5*dt*ietaSk*sk2*(DetaSDScTn*l0f - DetabarDs*smag + 1.0);
        
        fGA1[a] = 0.5*dt*ietaSk*sk0*(DetaSDScTn*l1f - DetabarDs*smag + 1.0);
        fGA1[a+1] = 1.0 + 0.5*dt*ietaSk*sk1*(DetaSDScTn*l1f -DetabarDs*smag + 1.0);
        fGA1[a+2] = 0.5*dt*ietaSk*sk2*(DetaSDScTn*l1f - DetabarDs*smag + 1.0);
        
        fGA2[a] = 0.5*dt*ietaSk*sk0*(DetaSDScTn*l2f - DetabarDs*smag + 1.0);
        fGA2[a+1] = 0.5*dt*ietaSk*sk1*(DetaSDScTn*l2f - DetabarDs*smag + 1.0);
        fGA2[a+2] = 1.0 + 0.5*dt*ietaSk*sk2*(DetaSDScTn*l2f - DetabarDs*smag + 1.0);
 
        
        for(int l=0; l<fNumR; l++)
        {
            double frac2 = fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(fTfk[k]/fTess - 1.0))))/fdeltac;
            /*DRdel/Depsilon^e_Bk*/
            fKTer(l,a) = -frac2*(dk00/iJ*(epse0-epstr0) + tauk0 + dk01/iJ*(epse1-epstr1)+dk02/iJ*(epse2-epstr2));
            fKTer(l,a+1) = -frac2*(dk01/iJ*(epse0-epstr0) + tauk1 + dk11/iJ*(epse1-epstr1)+dk12/iJ*(epse2-epstr2));
            fKTer(l,a+2) = -frac2*(dk02/iJ*(epse0-epstr0) + tauk2 + dk12/iJ*(epse1-epstr1)+dk22/iJ*(epse2-epstr2));
            
            double phi = fdalpha[l];
            double coeff1 = -fafrac*phi*(1.0 - 1.0/(1.0+exp(-fb_decay*(fTfk[l]/fTess-1.0))));
            double DprefactorDTfk  =coeff1*(fb_decay*(Tn/fTfk[l]-1.0)/fTess/(1.0+exp(-fb_decay*(fTfk[l]/fTess-1.0))) +Tn/fTfk[l]/fTfk[l]);
            
            double DsDTek = 0.0;
            if (smag>kSmall)
                DsDTek = DprefactorDTfk/prefactor*smag;
            
            /*DRes/DTe_l*/
            fKrTe(a,l) = 0.5*dt*ietaSk*(-fBB/Tn/ScTn/ScTn*(phi*fdeltac/fTfk[l])+DetabarDs*DsDTek)*sk0 -0.5*dt*ietaSk*sk0*DprefactorDTfk/prefactor;
            fKrTe(a+1,l) = 0.5*dt*ietaSk*(-fBB/Tn/ScTn/ScTn*(phi*fdeltac/fTfk[l])+DetabarDs*DsDTek)*sk1 -0.5*dt*ietaSk*sk1*DprefactorDTfk/prefactor;
            fKrTe(a+2,l) = 0.5*dt*ietaSk*(-fBB/Tn/ScTn/ScTn*phi*fdeltac/fTfk[l]+DetabarDs*DsDTek)*sk2 -0.5*dt*ietaSk*sk2*DprefactorDTfk/prefactor;

/*            fKrTe(a,l) = 0.5*dt*ietaSk*(DetaSDScTn*(phi*fdeltac/fTfk[l])+DetabarDs*DsDTek)*sk0 -0.5*dt*ietaSk*iJ*tauk0*DprefactorDTfk;
            fKrTe(a+1,l) = 0.5*dt*ietaSk*(DetaSDScTn*(phi*fdeltac/fTfk[l])+DetabarDs*DsDTek)*sk1 -0.5*dt*iJ*ietaSk*tauk1*DprefactorDTfk;
            fKrTe(a+2,l) = 0.5*dt*ietaSk*(DetaSDScTn*(phi*fdeltac/fTfk[l])+DetabarDs)*sk2 -0.5*dt*iJ*ietaSk*tauk2*DprefactorDTfk;
 */
        }
    }
/*if(CurrElementNumber()==0 && CurrIP()==1)
    {
        cout<< "\nfKAB: "<<fKAB;
        cout << "\nfKTer: "<<fKTer;
        cout << "\nfKrTe: "<<fKrTe;
        cout<< "\nfGA0: "<<fGA0;
        cout<< "\nfGA1: "<<fGA0;
        cout<< "\nfGA2: "<<fGA0;
    }
*/
    for(int k=0; k<fNumR; k++)
    {
        /*stiffness matrix*/
        double itauRk = itaubar/ftimesR[k];
        double frac2= fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(fTfk[k]/fTess - 1.0))))/fdeltac;
        double DaDTfk = -fb_decay/fTess/(1.0+exp(-fb_decay*(fTfk[k]/fTess - 1.0))); /*1/afrac * Dafrac/DTfk*/
        fKdel(k,k) = 1.0 + dt*itauRk-frac1*flatenth[k]/fTfk[k] - DaDTfk*frac2*plasticw ;
        for (int l=0; l<fNumR; l++)
        {
            double phi=fdalpha[l];
            double DtaubarDTfk = -fBB/Tn/ScTn/ScTn*(phi*fdeltac/fTfk[l]);  /*1/tauRi * DtauRi/DTfk*/
            fKdel(k,l) += -dt*itauRk*(fTfk[k]-Tn)*DtaubarDTfk;
        }
    
        fGB0[k] = dt*itauRk*(fTfk[k]-Tn)*DetaSDScTn*l0f+frac1*fTfk[k]/fT0*(fDtauDe_EQ(0,0)*depsilon0+fDtauDe_EQ(1,0)*depsilon1+fDtauDe_EQ(2,0)*depsilon2 + ftau_EQ[0]) + frac2*ftau[0];
        
        fGB1[k] = dt*itauRk*(fTfk[k]-Tn)*DetaSDScTn*l1f+frac1*fTfk[k]/fT0*(fDtauDe_EQ(0,1)*depsilon0+fDtauDe_EQ(1,1)*depsilon1+fDtauDe_EQ(2,1)*depsilon2 + ftau_EQ[1]) + frac2*ftau[1];
        
        fGB2[k] = dt*itauRk*(fTfk[k]-Tn)*DetaSDScTn*l2f+frac1*fTfk[k]/fT0*(fDtauDe_EQ(0,2)*depsilon0+fDtauDe_EQ(1,2)*depsilon1+fDtauDe_EQ(2,2)*depsilon2 + ftau_EQ[2]) + frac2*ftau[2];

    }
   if(CurrElementNumber()==0 && CurrIP()==1)
    {
 //       cout<< "\nfKdel: "<<fKdel;
//        cout<< "\nfGB0: "<<fGB0;
 //       cout<< "\nfGB1: "<<fGB0;
 //       cout<< "\nfGB2: "<<fGB0;
    }

	
    fK1=fKdel;
    fK2=fKAB;
    fKdel.Inverse();
    fKAB.Inverse();
    
    fKdel2.MultABC(fKTer, fKAB, fKrTe);
    fKAB2.MultABC(fKrTe, fKdel, fKTer);
    fK1.AddScaled(-1.0, fKdel2);
    fK2.AddScaled(-1.0, fKAB2);
    
    
    fInverse2.MultAB(fKTer,fKAB);
    fInverse2.Multx(fGA0,fG10);
    fInverse2.Multx(fGA1,fG11);
    fInverse2.Multx(fGA2,fG12);
    fG10 *= -1.0;
    fG11 *= -1.0;
    fG12 *= -1.0;
    fG10 += fGB0;
    fG11 += fGB1;
    fG12 += fGB2;
    
    fInverse3.MultAB(fKrTe,fKdel);
    fInverse3.Multx(fGB0,fG20);
    fInverse3.Multx(fGB1,fG21);
    fInverse3.Multx(fGB2,fG22);
    fG20 *= -1.0;
    fG21 *= -1.0;
    fG22 *= -1.0;
    fG20 += fGA0;
    fG21 += fGA1;
    fG22 += fGA2;
    
    /*Linear Solve changes the object*/
    fKdel=fK1;
    fKAB=fK2;
	/*KAB^-1 GBC*/
	fK1.LinearSolve(fG10);
	fK1 = fKdel;
	fK1.LinearSolve(fG11);
	fK1 = fKdel;
	fK1.LinearSolve(fG12);
	fK1 = fKdel;
    
    fK2.LinearSolve(fG20);
	fK2 = fKAB;
    fK2.LinearSolve(fG21);
	fK2 = fKAB;
    fK2.LinearSolve(fG22);
	fK2 = fKAB;
    
    for (int k =0; k<fNumR; k++)
	{
		Modulus1(k,0) = fG10[k];
		Modulus1(k,1) = fG11[k];
		Modulus1(k,2) = fG12[k];
	}
	for (int k =0; k<3*fNumS; k++)
	{
		Modulus2(k,0) = fG20[k];
		Modulus2(k,1) = fG21[k];
		Modulus2(k,2) = fG22[k];
	}
/*    if(CurrElementNumber()==0 && CurrIP()==1)
    {
        cout<< "\nModulus1: "<<Modulus1;
        cout<< "\nModulus2: "<<Modulus2;
    }
 */
}

void ET_multi::Compute_le(const ArrayT<dSymMatrixT>& C_vn, ArrayT<dSymMatrixT>& C_v, const dArrayT& Tfk_n, dArrayT& Tfk, double& heat) /*Solving eq.(23)*/
{
	
    double ctol = 1.00e-10;
	int maxiter = 20;
    double frac1 = 1.0/fdeltac;
    double ScTn = 0.0;
    double prefactor = 1.0;
    
	/*time step*/
	const double dt = fFSMatSupport->TimeStep();
    
	/*temperature and temperature step*/
	const double Tn = Compute_Temperature();
    
    
    const dMatrixT& F_last = F_total_last();
    	if (NumSD() == 2)
	{
		fF3D_last[0] = F_last[0];
		fF3D_last[1] = F_last[1];
		fF3D_last[2] = 0.0;
	    
		fF3D_last[3] = F_last[2];
		fF3D_last[4] = F_last[3];
		fF3D_last[5] = 0.0;
	    
		fF3D_last[6] = 0.0;
		fF3D_last[7] = 0.0;
		fF3D_last[8] = 1.0;
	}
	else fF3D_last = F_last;
    
    fb_last.MultAAT(fF3D_last);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_last, false);
	fEigs_last = fSpectralDecompSpat.Eigenvalues();
    
    double l0_last = fEigs_last[0];
    double l1_last = fEigs_last[1];
    double l2_last = fEigs_last[2];

    const dMatrixT& F = F_total();
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
    double J = sqrt(fEigs.Product());
	double iJ = 1.0/J;
//    iJ = 1.0;

    fEigs_dev = fEigs;
	fEigs_dev *= pow(J, -2.0*third);
	fPot[0]->DevStress(fEigs_dev, ftau_EQ);
    
    double l0 = fEigs[0];
    double l1 = fEigs[1];
    double l2 = fEigs[2];

    /*    double depsilon0=(sqrt(l0)-sqrt(l0_last))/sqrt(l0);
     double depsilon1=(sqrt(l1)-sqrt(l1_last))/sqrt(l1);
     double depsilon2=(sqrt(l2)-sqrt(l2_last))/sqrt(l2);
     */
    double depsilon0=(0.5*log(l0)-0.5*log(l0_last));
    double depsilon1=(0.5*log(l1)-0.5*log(l1_last));
    double depsilon2=(0.5*log(l2)-0.5*log(l2_last));

	/*calc trial solution l_tr, smag, and platicw*/
    double* pltr = fl_tr.Pointer();
	double* ple = fle.Pointer();
    double* pstress = fstressk.Pointer();
    ftau = 0.0;
	/*calculate trial state*/
	for (int i = 0; i < fNumS; i++)
	{
		/*calc trial elastic stretch*/
		fInverse.Inverse(C_vn[i]);
		fbe.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);
		fEigs_e = fSpectralDecompSpat.Eigenvalues();
		*pltr++ = fEigs_e[0];
		*pltr++ = fEigs_e[1];
		*pltr++ = fEigs_e[2];
        
		/*initial condition be = btr*/
		*ple++ = fEigs_e[0];
		*ple++ = fEigs_e[1];
		*ple++ = fEigs_e[2];
 
		double Je = sqrt(fEigs_e.Product());
		fEigs_dev = fEigs_e;
		fEigs_dev *= pow(Je,-2.0*third);
		
		/*calculate total neq stress*/
		fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
		ftau_NEQ *=fdmu[i];
        
        *pstress++ =ftau_NEQ[0];
        *pstress++ =ftau_NEQ[1];
        *pstress++ =ftau_NEQ[2];

        ftau+= ftau_NEQ;
	}
    double plasticw = 0.0; //initially zero because epe-eptr = 0.0;

    /*compute effective temperatures*/
    Compute_Tei(fEigs_last, fEigs,fl_tr, ftau_EQ, fle, fstressk, Tfk_n, Tfk);

    for(int k=0; k<fNumR; k++)
	{
		double phi=fdalpha[k];
        ScTn += fdeltac*phi*log(Tfk[k]/fT2);
        prefactor += phi*fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(Tfk[k]/fTess - 1.0))))*(Tn/Tfk[k]-1.0);
        flatenth[k] = Tfk[k]/fT0*(ftau_EQ[0]*depsilon0+ftau_EQ[1]*depsilon1+ftau_EQ[2]*depsilon2);
        
    }
	fEigs_dev = fEigs;
	fEigs_dev *= pow(J, -2.0*third);
    double energy = fPot[0]->Energy(fEigs_dev, J) /*eq.(14)TDN*/;
    energy -= fPot[0]->MeanEnergy(J);
    ScTn -= energy/fT0;

    double smag = prefactor*iJ*sqrt(0.5*(ftau[0]*ftau[0]+ftau[1]*ftau[1]+ftau[2]*ftau[2]));

	/*calc viscosity functions*/
	double ietabar = StressRelaxationFunc(Tn, smag);   /*eq.(24)*/
    double itaubar = StructuralRelaxationFunc(Tn, ScTn);
	double DetabarDs = 0.0;
    if(smag > kSmall)
    {
        double x = (fQS*smag)/(Tn);
        double cothx = cosh(x)/sinh(x);
        DetabarDs = (1.0 - x*cothx)/smag;
    }
    
	/*calc residual*/
	/*re-assign pointers*/
	pltr = fl_tr.Pointer();
	ple = fle.Pointer();
    pstress = fstressk.Pointer();
	double* pr = fRes.Pointer();
	fRes = 0.0;
	double r0, r1, r2, tol;
	double temp = 0.0;
	for (int k = 0; k < fNumS; k++)
	{
		double epstr0 = 0.5*log(*pltr++);
		double epstr1 = 0.5*log(*pltr++);
		double epstr2 = 0.5*log(*pltr++);
		
		/*calculate neq stress tauk*/
		fEigs_dev[0] = *ple++;
		fEigs_dev[1] = *ple++;
		fEigs_dev[2] = *ple++;
		double epse0 = 0.5*log(fEigs_dev[0]);
		double epse1 = 0.5*log(fEigs_dev[1]);
		double epse2 = 0.5*log(fEigs_dev[2]);
		
		double Je = sqrt(fEigs_dev.Product());
		fEigs_dev *= pow(Je,-2.0*third);

        ftau_NEQ[0] =  *pstress++;
        ftau_NEQ[1] =  *pstress++;
        ftau_NEQ[2] =  *pstress++;
        double s0 = prefactor*iJ*ftau_NEQ[0];
        double s1 = prefactor*iJ*ftau_NEQ[1];
        double s2 = prefactor*iJ*ftau_NEQ[2];
        
		/*viscosity*/
        double muk=fdmu[k];
        double itauSk = itaubar/ftimesS[k];
        double ietaS0k = itauSk/muk;
        double ietaSk = ietaS0k*ietabar;
		/*calc residual*/
        r0 = epse0 - epstr0 + 0.5*dt*ietaSk*s0;
        r1 = epse1 - epstr1 + 0.5*dt*ietaSk*s1;
        r2 = epse2 - epstr2 + 0.5*dt*ietaSk*s2;
/*       if(CurrElementNumber()==0 && CurrIP()==1)
        {
            cout << "\nitaubar:"<<itaubar;
            cout<< "\nviscosity: "<<ietaSk;
            cout << "\nScTn: "<<ScTn;
            cout<< "\nsmag: "<<smag;
            cout<< "\nres : "<<r0<<"\t"<<r1<<"\t"<<r2;
            cout<< "\nietaS : "<<0.5*dt*ietaSk*s0<<"\t"<<0.5*dt*ietaSk*s1<<"\t"<<0.5*dt*ietaSk*s2;
            cout<< "\nepse0 : "<<epse0<<"\t"<<epse1<<"\t"<<epse0;
            cout<< "\nepstr0 : "<<epstr0<<"\t"<<epstr1<<"\t"<<epstr0;
            cout<< "\nepse0-epstr0 : "<<epse0-epstr0<<"\t"<<epse1-epstr1<<"\t"<<epse0-epstr0;
            cout<<"\nplasticwk: "<<plasticw;
        }
 */
        
		temp += r0*r0 + r1*r1+ r2*r2;
		*pr++ = r0;
		*pr++ = r1;
		*pr++ = r2;
 }
	tol = sqrt(temp);
	double tol0 = tol;
	double reltol =tol0;
	int iter = 0;
    
 /*   if(CurrElementNumber()==0 && CurrIP()==1)
    {
        cout<<setprecision(12)<< "\nle 0: "<<tol;
 //       cout << "\nfle: "<<fle;
    }
 */
 
    while (tol > ctol  && iter < maxiter)
	{
		iter++;
        fKAB2 = 0.0;
		fKAB = 0.0; /*DfRes/Depsilon^e_Bk*/
        fKrTe = 0.0; /*DRes/DTe_l*/
        fKTer = 0.0;/*DRdel/Depsilon^e_Bk*/
		fKdel = 0.0;  /*DRdel/DTel*/
		ple = fle.Pointer();
        pltr = fl_tr.Pointer();
        pstress = fstressk.Pointer();
		double DetabarDs = 0.0;
        if(smag > kSmall)
        {
            double x = (fQS*smag)/(Tn);
            double cothx = cosh(x)/sinh(x);
            DetabarDs = (1.0 - x*cothx)/smag;
        }
       
		for(int k=0; k<fNumS; k++)
		{
            double epstr0 = 0.5*log(*pltr++);
            double epstr1 = 0.5*log(*pltr++);
            double epstr2 = 0.5*log(*pltr++);
            
			fEigs_dev[0] = *ple++;
			fEigs_dev[1] = *ple++;
			fEigs_dev[2] = *ple++;
            
            double epse0 = 0.5*log(fEigs_dev[0]);
            double epse1 = 0.5*log(fEigs_dev[1]);
            double epse2 = 0.5*log(fEigs_dev[2]);

            /*calculate neq stress tauk*/
			double Je = sqrt(fEigs_dev[0]*fEigs_dev[1]*fEigs_dev[2]);
			fEigs_dev *= pow(Je,-2.0*third);
            fPot[1]->DevMod(fEigs_dev,fDtauDe_NEQ);
            
			/*moduli*/
			/*viscosity and its derivatives*/
			double muk=fdmu[k];
            double itauSk = itaubar/ftimesS[k];
			double ietaS0k = itauSk/muk;
			double ietaSk = ietaS0k*ietabar;

			double dk00 = muk*iJ*fDtauDe_NEQ(0,0);
			double dk11 = muk*iJ*fDtauDe_NEQ(1,1);
			double dk22 = muk*iJ*fDtauDe_NEQ(2,2);
			double dk12 = muk*iJ*fDtauDe_NEQ(1,2);
			double dk02 = muk*iJ*fDtauDe_NEQ(0,2);
			double dk01 = muk*iJ*fDtauDe_NEQ(0,1);
			
            double tauk0 =  *pstress++;
            double tauk1 =  *pstress++;
            double tauk2 =  *pstress++;
            double sk0 = (tauk0)*iJ*prefactor;
			double sk1 = (tauk1)*iJ*prefactor;
			double sk2 = (tauk2)*iJ*prefactor;
            
  /*          if(CurrElementNumber()==0 && CurrIP()==1)
            {
                cout<< "\nietaSk: "<<ietaSk;
                cout << "\nfDtauDe_NEQ: "<<fDtauDe_NEQ;
                cout<<"\n"<<0.5*dt*ietaSk*dk22*prefactor;
            }*/
            
			/*stiffness matrix*/
			int a = k*3;
			fKAB(a,a) = 1.0 + 0.5*dt*ietaSk*dk00*prefactor;
			fKAB(a,a+1) = 0.5*dt*ietaSk*dk01*prefactor;
			fKAB(a,a+2) = 0.5*dt*ietaSk*dk02*prefactor;
			fKAB(a+1,a) = 0.5*dt*ietaSk*dk01*prefactor;
			fKAB(a+1,a+1) = 1.0 + 0.5*dt*ietaSk*dk11*prefactor;
			fKAB(a+1,a+2) = 0.5*dt*ietaSk*dk12*prefactor;
			fKAB(a+2,a) = 0.5*dt*ietaSk*dk02*prefactor;
			fKAB(a+2,a+1) = 0.5*dt*ietaSk*dk12*prefactor;
 			fKAB(a+2,a+2) = 1.0 + 0.5*dt*ietaSk*dk22*prefactor;

 /*           if(CurrElementNumber()==0 && CurrIP()==1)
            {
                cout<< "\nietaSk: "<<ietaSk;
                cout << "\nfDtauDe_NEQ: "<<fDtauDe_NEQ;
                gg<<"\n"<<fKAB(a+2,a+2);
            }*/
            
            double *ple2 = fle.Pointer();
			for (int l=0; l<fNumS; l++)       /*evolution eqs. are coupled in this case because of s (or smag here)*/
			{
				/*moduli*/
				fEigs_dev[0] = *ple2++;
				fEigs_dev[1] = *ple2++;
				fEigs_dev[2] = *ple2++;
				double Je2 = sqrt(fEigs_dev[0]*fEigs_dev[1]*fEigs_dev[2]);
				fEigs_dev *= pow(Je2,-2.0*third);
				fPot[1]->DevMod(fEigs_dev,fDtauDe_NEQ);
				double dl00 = fdmu[l]*iJ*fDtauDe_NEQ(0,0);
				double dl11 = fdmu[l]*iJ*fDtauDe_NEQ(1,1);
				double dl22 = fdmu[l]*iJ*fDtauDe_NEQ(2,2);
				double dl12 = fdmu[l]*iJ*fDtauDe_NEQ(1,2);
				double dl02 = fdmu[l]*iJ*fDtauDe_NEQ(0,2);
				double dl01 = fdmu[l]*iJ*fDtauDe_NEQ(0,1);
                
                double cl0 = 0.0;
                double cl1 = 0.0;
                double cl2 = 0.0;
                
                if(smag >kSmall)
                {
                     cl0 = 0.5*(iJ*ftau[0]*dl00 + iJ*ftau[1]*dl01 + iJ*ftau[2]*dl02)*prefactor*prefactor/smag;
                     cl1 = 0.5*(iJ*ftau[0]*dl01 + iJ*ftau[1]*dl11 + iJ*ftau[2]*dl12)*prefactor*prefactor/smag;
                     cl2 = 0.5*(iJ*ftau[0]*dl02 + iJ*ftau[1]*dl12 + iJ*ftau[2]*dl22)*prefactor*prefactor/smag;
                  }
                
				int b = l*3;
				fKAB(a,b) += - 0.5*dt*ietaSk*(DetabarDs*cl0)*sk0;
				fKAB(a,b+1) += - 0.5*dt*ietaSk*(DetabarDs*cl1)*sk0;
				fKAB(a,b+2) += - 0.5*dt*ietaSk*(DetabarDs*cl2)*sk0;
				fKAB(a+1,b) += - 0.5*dt*ietaSk*(DetabarDs*cl0)*sk1;
				fKAB(a+1,b+1) += - 0.5*dt*ietaSk*(DetabarDs*cl1)*sk1;
				fKAB(a+1,b+2)+= - 0.5*dt*ietaSk*(DetabarDs*cl2)*sk1;
				fKAB(a+2,b) += - 0.5*dt*ietaSk*(DetabarDs*cl0)*sk2;
				fKAB(a+2,b+1) += - 0.5*dt*ietaSk*(DetabarDs*cl1)*sk2;
                fKAB(a+2,b+2) += - 0.5*dt*ietaSk*(DetabarDs*cl2)*sk2;
                
/*                if(CurrElementNumber()==0 && CurrIP()==1 && b==a)
                {
                    cout<< "\nsk2: "<<sk2;
                    cout<<"\ncl: "<<cl2;
                    cout<<"\nDeta: "<<ietaS0k*(DetabarDs);
                    cout<<"\nagain"<<- 0.5*dt*ietaS0k*(-DetabarDs*cl2)*sk2<<"\t"<<fKAB(a+2,b+2);
                    
                }
 */
            }
            for(int l=0; l<fNumR; l++)
            {
                double frac2 = fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(Tfk[k]/fTess - 1.0))))/fdeltac;
                /*DRdel/Depsilon^e_Bk*/
                fKTer(l,a) = -frac2*(dk00/iJ*(epse0-epstr0) + tauk0 + dk01/iJ*(epse1-epstr1)+dk02/iJ*(epse2-epstr2));
                fKTer(l,a+1) = -frac2*(dk01/iJ*(epse0-epstr0) + tauk1 + dk11/iJ*(epse1-epstr1)+dk12/iJ*(epse2-epstr2));
                fKTer(l,a+2) = -frac2*(dk02/iJ*(epse0-epstr0) + tauk2 + dk12/iJ*(epse1-epstr1)+dk22/iJ*(epse2-epstr2));
                
/*                double phi = fdalpha[l];
                double coeff1 = fafrac*phi*(1.0 - 1.0/(1.0+exp(-fb_decay*(Tfk[l]/fTess-1.0))));
                double DprefactorDTfk  = coeff1*(fb_decay*(1.0-Tn/Tfk[l])/fTess/(1.0+exp(-fb_decay*(Tfk[l]/fTess-1.0))) -Tn/Tfk[l]/Tfk[l]);
                double DsDTek = 0.0;
                if (smag>kSmall)
                    DsDTek = DprefactorDTfk/prefactor*smag;
*/
                double phi = fdalpha[l];
                double coeff1 = -fafrac*phi*(1.0 - 1.0/(1.0+exp(-fb_decay*(Tfk[l]/fTess-1.0))));
                double DprefactorDTfk  =coeff1*(fb_decay*(Tn/Tfk[l]-1.0)/fTess/(1.0+exp(-fb_decay*(Tfk[l]/fTess-1.0))) +Tn/Tfk[l]/Tfk[l]);
                
                double DsDTek = 0.0;
                if (smag>kSmall)
                    DsDTek = DprefactorDTfk/prefactor*smag;

                /*DRes/DTe_l*/
                fKrTe(a,l) = 0.5*dt*ietaSk*(-fBB/Tn/ScTn/ScTn*(phi*fdeltac/Tfk[l])+DetabarDs*DsDTek)*sk0 -0.5*dt*ietaSk*sk0*DprefactorDTfk/prefactor;
                fKrTe(a+1,l) = 0.5*dt*ietaSk*(-fBB/Tn/ScTn/ScTn*(phi*fdeltac/Tfk[l])+DetabarDs*DsDTek)*sk1 -0.5*dt*ietaSk*sk1*DprefactorDTfk/prefactor;
                fKrTe(a+2,l) = 0.5*dt*ietaSk*(-fBB/Tn/ScTn/ScTn*phi*fdeltac/Tfk[l]+DetabarDs*DsDTek)*sk2 -0.5*dt*ietaSk*sk2*DprefactorDTfk/prefactor;
            }
        }
		for(int k=0; k<fNumR; k++)
		{
			/*stiffness matrix*/
			double itauRk = itaubar/ftimesR[k];
            double frac2= fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(Tfk[k]/fTess - 1.0))))/fdeltac;
            double DaDTfk = -fb_decay/fTess/(1.0+exp(-fb_decay*(Tfk[k]/fTess - 1.0))); /*1/afrac * Dafrac/DTfk*/
            fKdel(k,k) = 1.0 + dt*itauRk-frac1*flatenth[k]/Tfk[k] - DaDTfk*frac2*plasticw ;
			for (int l=0; l<fNumR; l++)
            {
                double phi=fdalpha[l];
                double DtaubarDTfk = -fBB/Tn/ScTn/ScTn*(phi*fdeltac/Tfk[l]);  /*1/tauRi * DtauRi/DTfk*/
				fKdel(k,l) += -dt*itauRk*(Tfk[k]-Tn)*DtaubarDTfk;            }
		}
		fKdel.Inverse();
        fKAB2.MultABC(fKrTe, fKdel, fKTer);
        fKAB.AddScaled(-1.0, fKAB2);
        
 /*       if(CurrElementNumber()==0 && CurrIP()==1)
        {
            cout<<setprecision(12)<< "\nKAB: "<<fKAB;
            //                  cout<<setprecision(12)<< "\nKAB2: "<<fKAB2;
            cout<< "\nKdel: "<<fKdel;
            //                  cout<< "\nKTer: "<<fKrTe;
            //                 cout<< "\nKrTe: "<<fKTer;
            cout << "\nEigse: "<<fle;
        }
  */
       /*Solve  and update*/
		fKAB.LinearSolve(fRes);
        
        /*update Cv and be*/
		ple = fle.Pointer();
		pr = fRes.Pointer();
		for (int k = 0; k < 3*fNumS; k++)
		{
			double dep = -(*pr++);
			*ple++ *=exp(2.0*dep);
		}
        /*update stresses*/
        ple = fle.Pointer();
        pstress = fstressk.Pointer();
        for (int k = 0; k < fNumS; k++)
        {
            fEigs_e[0] = *ple++;
            fEigs_e[1] = *ple++;
            fEigs_e[2] = *ple++;
    
            double Je = sqrt(fEigs_e.Product());
			fEigs_dev[0] = fEigs_e[0]*pow(Je,-2.0*third);
 			fEigs_dev[1] = fEigs_e[1]*pow(Je,-2.0*third);
			fEigs_dev[2] = fEigs_e[2]*pow(Je,-2.0*third);
            fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
            ftau_NEQ *=fdmu[k];

            *pstress++ = ftau_NEQ[0];
            *pstress++ = ftau_NEQ[1];
            *pstress++ = ftau_NEQ[2];
       }
        
		/*update residual*/
        /*update Teff related*/
        Compute_Tei(fEigs_last, fEigs,fl_tr, ftau_EQ, fle,fstressk, Tfk_n, Tfk);
        
        double ScTn = 0.0;
        double prefactor = 1.0;
        for(int k=0; k<fNumR; k++)
        {
            double phi=fdalpha[k];
            ScTn += fdeltac*phi*log(Tfk[k]/fT2);
            prefactor += phi*fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(Tfk[k]/fTess - 1.0))))*(Tn/Tfk[k]-1.0);
//            cout << "\nstrain: "<<-log(l1)<<"\t ratio: "<<Tfk[k]/fTess<<"\tafrac: "<<fafrac*(1.0 - 1.0/(1.0+exp(-fb_decay*(Tfk[k]/fTess - 1.0))));
            flatenth[k] = Tfk[k]/fT0*(ftau_EQ[0]*depsilon0+ftau_EQ[1]*depsilon1+ftau_EQ[2]*depsilon2);
        }
        
        fEigs_dev = fEigs;
        fEigs_dev *= pow(J, -2.0*third);
        double energy = fPot[0]->Energy(fEigs_dev, J) /*eq.(14)TDN*/;
        energy -= fPot[0]->MeanEnergy(J);
        ScTn -= energy/fT0;
        
		/*update  smag*/
		ftau= 0.0;
        plasticw = 0.0;
		ple = fle.Pointer();
        pltr = fl_tr.Pointer();
        pstress = fstressk.Pointer();
		for (int i = 0; i < fNumS; i++)
		{
            double epstr0 = 0.5*log(*pltr++);
            double epstr1 = 0.5*log(*pltr++);
            double epstr2 = 0.5*log(*pltr++);
            
            double epse0 = 0.5*log(*ple++);
            double epse1 = 0.5*log(*ple++);
            double epse2 = 0.5*log(*ple++);
            
            ftau_NEQ[0]= *pstress++;
            ftau_NEQ[1]= *pstress++;
            ftau_NEQ[2]= *pstress++;

            /*           if (CurrElementNumber()==0 && CurrIP()==1)
            {
                cout <<"\neigse"<< fEigs_e;
                cout <<"\neigsdev"<< fEigs_dev;
                cout << "\ntauneq: "<<ftau_NEQ;
            }
  */
            ftau+= ftau_NEQ;
            double plasticwk = -(ftau_NEQ[0]*(epse0-epstr0)+ ftau_NEQ[1]*(epse1-epstr1) + ftau_NEQ[2]*(epse2-epstr2));
            plasticw += plasticwk;
		}
        smag = prefactor*iJ*sqrt(0.5*(ftau[0]*ftau[0]+ftau[1]*ftau[1]+ftau[2]*ftau[2]));
/*        if(CurrElementNumber()==0 && CurrIP()==1)
            cout << "\nsmag: "<<smag;
  */
		/*update viscosity*/
        ietabar = StressRelaxationFunc(Tn, smag);   /*eq.(24)*/
        itaubar = StructuralRelaxationFunc(Tn, ScTn);
        DetabarDs = 0.0;
        if(smag > kSmall)
        {
            double x = (fQS*smag)/(Tn);
            double cothx = cosh(x)/sinh(x);
            DetabarDs = (1.0 - x*cothx)/smag;
        }
        
		pltr = fl_tr.Pointer();
		ple = fle.Pointer();
		pr = fRes.Pointer();
        pstress = fstressk.Pointer();
		fRes = 0.0;
		temp = 0.0;
        /*for calculating KTer DRdel/Depsilon^e_Bk*/
		for (int k = 0; k < fNumS; k++)
		{
			double epstr0 = 0.5*log(*pltr++);
			double epstr1 = 0.5*log(*pltr++);
			double epstr2 = 0.5*log(*pltr++);
			
			/*calculate neq stress tauk*/
			fEigs_dev[0] = *ple++;
			fEigs_dev[1] = *ple++;
			fEigs_dev[2] = *ple++;
			double epse0 = 0.5*log(fEigs_dev[0]);
			double epse1 = 0.5*log(fEigs_dev[1]);
			double epse2 = 0.5*log(fEigs_dev[2]);

            ftau_NEQ[0]= *pstress++;
            ftau_NEQ[1]= *pstress++;
            ftau_NEQ[2]= *pstress++;
            
            double s0 = prefactor*iJ*ftau_NEQ[0];
            double s1 = prefactor*iJ*ftau_NEQ[1];
            double s2 = prefactor*iJ*ftau_NEQ[2];

			/*viscosity*/
            double muk=fdmu[k];
            double itauSk = itaubar/ftimesS[k];
            double ietaS0k = itauSk/muk;
            double ietaSk = ietaS0k*ietabar;
			
			/*calc residual*/
			r0 = epse0 - epstr0 + 0.5*dt*ietaSk*s0;
			r1 = epse1 - epstr1 + 0.5*dt*ietaSk*s1;
			r2 = epse2 - epstr2 + 0.5*dt*ietaSk*s2;
            
/*            if(CurrElementNumber()==0 && CurrIP()==1)
            {
                cout << "\nitaubar:"<<itaubar;
                cout<< "\nviscosity: "<<ietaSk;
                cout << "\nScTn: "<<ScTn;
                cout<< "\nsmag: "<<smag;
                cout<< "\nres : "<<r0<<"\t"<<r1<<"\t"<<r2;
                cout<< "\nietaS : "<<0.5*dt*ietaSk*s0<<"\t"<<0.5*dt*ietaSk*s1<<"\t"<<0.5*dt*ietaSk*s2;
                cout<< "\nepse0 : "<<epse0<<"\t"<<epse1<<"\t"<<epse0;
                cout<< "\nepstr0 : "<<epstr0<<"\t"<<epstr1<<"\t"<<epstr0;
                cout<< "\nepse0-epstr0 : "<<epse0-epstr0<<"\t"<<epse1-epstr1<<"\t"<<epse0-epstr0;
                cout<<"\nplasticwk: "<<plasticw;
           }
 */

			temp += r0*r0 + r1*r1+ r2*r2;
			*pr++ = r0;
			*pr++ = r1;
			*pr++ = r2;

            /*for calculating Dplastic part of KTer DRdel/Depsilon^e_Bk*/
            
            fPot[1]->DevMod(fEigs_dev,fDtauDe_NEQ);
            double dk00 = muk*iJ*fDtauDe_NEQ(0,0);
            double dk11 = muk*iJ*fDtauDe_NEQ(1,1);
            double dk22 = muk*iJ*fDtauDe_NEQ(2,2);
            double dk12 = muk*iJ*fDtauDe_NEQ(1,2);
            double dk02 = muk*iJ*fDtauDe_NEQ(0,2);
            double dk01 = muk*iJ*fDtauDe_NEQ(0,1);
            
            double ck0 = 0.0;
            double ck1 = 0.0;
            double ck2 = 0.0;
            
            if(smag >kSmall)
            {
                ck0 = 0.5*(iJ*ftau[0]*dk00 + iJ*ftau[1]*dk01 + iJ*ftau[2]*dk02)*prefactor*prefactor/smag;
                ck1 = 0.5*(iJ*ftau[0]*dk01 + iJ*ftau[1]*dk11 + iJ*ftau[2]*dk12)*prefactor*prefactor/smag;
                ck2 = 0.5*(iJ*ftau[0]*dk02 + iJ*ftau[1]*dk12 + iJ*ftau[2]*dk22)*prefactor*prefactor/smag;
            }
        }
        tol = sqrt(temp);
        reltol = tol/tol0;
        
/*       if(CurrElementNumber()==0 && CurrIP()==1)
        {
            cout<< "\nle: "<<iter<<"\t"<<tol;
 //           cout << "\nfle: "<<fle;

        }
 */

	}
    if (iter >= maxiter)
	{
		cout<<"\n Number of iteration exceeds maximum. tol:  "<<tol;
		cout<< "\nelem: "<<CurrElementNumber()<<"\nIP: "<<CurrIP();
		ExceptionT::GeneralFail("ET_multi::Compute_le", "number of iteration exceeds maximum");
	}
    
	/*update Cv with converged solution*/
	ple = fle.Pointer();
	for (int k = 0; k < fNumS; k++)
	{
		fInverse.Inverse(C_vn[k]);
		fb_tr.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fb_tr, false);
        
		fEigs_e[0] = *ple++;
		fEigs_e[1] = *ple++;
		fEigs_e[2] = *ple++;
		fbe = fSpectralDecompSpat.EigsToRank2(fEigs_e); /*be which is colinear with btr*/
		fbe.Inverse();
		C_v[k].MultQTBQ(fF3D, fbe);
	}
    
    /*calculate incremental heat*/
    const double rho = Density();
    double structural=0.0;
    double coeff_latentheat=0.0;
    for (int k = 0; k < fNumS; k++)
    {
        structural -= rho*fdeltac*fdalpha[k]*(Tfk[k]-Tfk_n[k]);
        coeff_latentheat += (Tfk[k]/fT0)*fdalpha[k];
    }
    /*add contributions from structural relaxation and latent heat*/
    heat = structural + coeff_latentheat*(ftau_EQ[0]*depsilon0+ftau_EQ[1]*depsilon1+ftau_EQ[2]*depsilon2);
    
    plasticw = 0.0;
    ple = fle.Pointer();
    pltr = fl_tr.Pointer();
    pstress = fstressk.Pointer();
    for (int i = 0; i < fNumS; i++)
    {
        double epstr0 = 0.5*log(*pltr++);
        double epstr1 = 0.5*log(*pltr++);
        double epstr2 = 0.5*log(*pltr++);
        
        double epse0 = 0.5*log(*ple++);
        double epse1 = 0.5*log(*ple++);
        double epse2 = 0.5*log(*ple++);
        
        ftau_NEQ[0]= *pstress++;
        ftau_NEQ[1]= *pstress++;
        ftau_NEQ[2]= *pstress++;
        
        /*           if (CurrElementNumber()==0 && CurrIP()==1)
         {
         cout <<"\neigse"<< fEigs_e;
         cout <<"\neigsdev"<< fEigs_dev;
         cout << "\ntauneq: "<<ftau_NEQ;
         }
         */
        double plasticwk = -(ftau_NEQ[0]*(epse0-epstr0)+ ftau_NEQ[1]*(epse1-epstr1) + ftau_NEQ[2]*(epse2-epstr2));
        plasticw += plasticwk;
    }

    
    /*add contribution from plastic work*/
    heat += plasticw;
}

/*calculates the series approximation of the KWW distribution function - rho(tau/tauKWW)*/
/*Lindsey et al. J. Chem. Phys. 73:3348, 1980*/
double ET_multi::StretchedExponentialSpectrum(const double x, const double tauKWW, const double betakWW)
{
	/*x = tau/taukWW*/
	int nk = 201;
	double rho = 0;
	
	for (int k = 0; k < nk; k++)
	{
		double kfac= fGamma.Function(k+1);
		double kpow = pow(-1.0,k);
		double kb = betakWW*k;
		rho += kpow/kfac * sin(pi*kb)*fGamma.Function(kb+1.0)*pow(x,kb+1.0);
	}
	rho *= -1.0/(pi*tauKWW*x*x);
	return(rho);
	
}

void ET_multi::Compute_discrete_spectrum(const double beta, const double tauWW, const double taumin, const double taumax, const int num, dArrayT& times, dArrayT& spectrum, StringT filename)
{
	cout << "\nCalculating structural relaxation spectrum";
	dArrayT Hcum(num);
	double ratio = taumax/taumin;
	for (int n = 0; n < num; n++)
	{
		double exponent;
		exponent = n/(num-1.0);
		times[n] = taumin*pow(ratio,exponent);
        
		/*set integration step (ds) and initialize  integration intervals Hcum and tau*/
		int nint = 1000;
		double ds, s, Hs;
		if (n==0)
		{
			ds = times[n]/nint;
			Hs = 0.0;
			s = 0.0;
		}
		else
		{
			ds = (times[n]-times[n-1])/nint;
			s= times[n-1];
			Hs = Hcum[n-1];
		}
		for (int i = 0; i < nint; i++)
		{
			s += ds;
			double x = s/tauWW;
			double rho = StretchedExponentialSpectrum(x, tauWW, beta);
			Hs += rho*ds;
		}
		Hcum[n] = Hs;
	}
	/*calculate discrete spectrum by approximating Hcum  Haupt et al. 2000*/
	double cum=0;
	double temp = 0.5*(Hcum[1]+Hcum[0]);
	cum = temp;
	spectrum[0] = temp;
	
	for (int n = 1; n < num-1; n++)
	{
		temp = 0.5*(Hcum[n+1]-Hcum[n-1]);
		cum += temp;
		spectrum[n] = temp;
	}
	spectrum[num-1] = (1.0 - cum);
    
	/*shift spectrums to Tg*/
	double exponent = fC1*(fT0-fTg)/(fC2+fT0-fTg);
	times *= pow(10.0,exponent);
    
	/*write distribution to file*/
	ofstreamT out;
	out.open(filename);
	for (int i=0; i<spectrum.Length(); i++)
		out << times[i]<<"\t"<<Hcum[i]<<"\t"<<spectrum[i]<<"\n";
	out.close();
	cout<< "\nResults written to "<<filename;
    
    
}

