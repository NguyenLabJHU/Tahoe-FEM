/* $Id: SMP_coupled.cpp,v 1.1 2015-09-20 03:48:45 tahoe.vickynguyen Exp $ */
/* created: RX (11/20/2013) */

#include "SMP_coupled.h"
#include "FSThermoMechMatT.h"
#include "PotentialT.h"
#include "NeoHookean.h"
#include "ArrudaBoyce.h"

#include "ifstreamT.h"
#include "ExceptionT.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "ParameterContainerT.h"

using namespace Tahoe;

const double loge = log10(exp(1.0));
const double third = 1.0/3.0; 
const double small = 1.0e-16;
const int kNumOutputVar = 10;
static const char* Labels[kNumOutputVar] = {"thermal_dialation","fictive_temperature", "lm1","lm2","lm3", "lmv1","lmv2","lmv3","T","gamdot"};


/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
/* constructors */
SMP_coupled::SMP_coupled(void):ParameterInterfaceT("SMP_coupled"){ }


/*rewrite this funciton RX */
double SMP_coupled::RetardationTime(const double Temperature, const double Te)
{

	/*Hodge's Model*/
	double coeff = -fC1/log10(exp(1))*(fC2*(Temperature - Te) + Temperature*(Te-fTg))/(Temperature*(fC2 + Te - fTg));
	double tauR = ftaug*exp(coeff);

	/*check the limits*/
	if (tauR > ftauRH)
		tauR = ftauRH;
	else if (tauR < ftauRL)
		tauR = ftauRL;
	
	return(tauR);
}

/*rewrite this funciton RX */
double SMP_coupled::ShearViscosity(const double Temperature, const double Te, const double smag)
{
	double coeff = -fC1/log10(exp(1.0))*(fC2*(Temperature - Te) + Temperature*(Te-fTg))/(Temperature*(fC2 + Te - fTg));
	double g = exp(coeff);
	/*Eyringen model*/
	double etaS;	

	double taumag = smag;
	if (smag/fsy0 > small)
	{
		etaS = fetaS0*g*taumag/sinh(fQS/Temperature * taumag/fsy0);
		etaS *= fQS/(fsy0*Temperature);
	}
	else
	{
		etaS = fetaS0*g;
	}
			
	/*check the limits*/
	if (etaS > fetaSH)
		etaS = fetaSH;
	else if (etaS < fetaSL)
		etaS = fetaSL;

	return(etaS);
}

/*dMatrixT& SMP_coupled::Conductivity(void)
{
    fkij.Identity(fk);
    return fkij;
}

double const SMP_coupled::Capacity(void)
{
    cout<<"\n SMP_coupled the capacity is:"<<fcg*frho<<endl;
    return fcg*frho;
}

double SMP_coupled::Density(void)
{
    return frho;
}
*/
 const dArrayT& SMP_coupled::q_i(void)
{
 const dMatrixT& kij=Conductivity();
    //check this function
    dArrayT fT_i=TemperatureGradient();
    int dim=fT_i.Length();
 //  cout<<"\n the length of TemperatureGradient() "<<fT_i.Length();
    if (dim==2)
    {
    fHeatFlux[0]=kij(0,0)*fT_i[0]+kij(0,1)*fT_i[1];
    fHeatFlux[1]=kij(1,0)*fT_i[0]+kij(1,1)*fT_i[1];
         fHeatFlux[2]=0;
    }
    else
    {
      //  cout<<"\n I came to the 3D place "<<endl;
        fHeatFlux[0]=kij(0,0)*fT_i[0]+kij(0,1)*fT_i[1]+kij(0,2)*fT_i[2];
        fHeatFlux[1]=kij(1,0)*fT_i[0]+kij(1,1)*fT_i[1]+kij(1,2)*fT_i[2];
        fHeatFlux[2]=kij(2,0)*fT_i[0]+kij(2,1)*fT_i[1]+kij(2,2)*fT_i[2];
    }    
    return fHeatFlux;
}

double SMP_coupled::StrainEnergyDensity(void)
{
	/*calculates equilibrium part*/
	double T = Compute_Temperature();
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

	/*elastic stretch*/
	fb.MultAAT(fF3D);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	
	double J = sqrt(fEigs.Product());
	fEigs_dev = fEigs;
	fEigs_dev *= pow(J, -2.0*third);
     
	double energy = 0.0;
	energy = fPot[0]->Energy(fEigs_dev, J, T);
  /* if(CurrElementNumber()==0&&CurrIP()==0)
    {
        cout<<"\n fEigs_dev is :"<<fEigs_dev;
        cout<<"\n J is "<<J;
        cout<<"\n equilibrium strain energy is :"<<fPot[0]->Energy(fEigs_dev, J, T);
    } */

if (fNumProcess > 0)
{
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
      /*  if(CurrElementNumber()==0&&CurrIP()==0)
        {
            cout<<"\n fEigse_dev is :"<<fEigs_dev;
            cout<<"\n Je is "<<Je;
            cout<<"\n nonequilibrium strain energy in SMP_coupled is :"<<fPot[1]->Energy(fEigs_dev, Je);
        }*/
	}
     }
 // if(CurrElementNumber()==0&&CurrIP()==0)
  //  cout<<"\n ip_strain energy in SMP_coupled strain_energy density is :"<<energy;
	return(energy);
}

/*rewrite this funciton RX */
const dMatrixT& SMP_coupled::ThermalDeformation_Inverse(void)
{
	/*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
	double Temp = Compute_Temperature();
	double dalpha=falphar-falphag;

	/*thermal volume deformation*/
    double ThetaT=exp(falphag*(Temp-fT0)+dalpha*(fTe[0]-fT0));
	fF_T_inv.Identity(pow(ThetaT,-1.0*third));
 
	return(fF_T_inv);
}

/* stresses */
const dSymMatrixT& SMP_coupled::s_ijeq(void)
{
	const dMatrixT& F = MechanicalDeformation();
	double T = Compute_Temperature();

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
/*	if(CurrElementNumber()==0&&CurrIP()==0)
		cout <<setprecision(12)<< "\nfEigs: "<<fEigs;
*/
	/*jacobian determinant*/
	double J = sqrt(fEigs.Product());
	
	fEigs_dev = fEigs;
	fEigs_dev *= pow(J,-2.0*third);
	
	fPot[0]->DevStress(fEigs_dev, ftau_EQ, T);	
	ftau_EQ += fPot[0]->MeanStress(J);
	
	fStress3D_EQ = fSpectralDecompSpat.EigsToRank2(ftau_EQ);
    
    if (NumSD() == 2)
    {
        fStress_EQ[0] = fStress3D_EQ[0];
        fStress_EQ[1] = fStress3D_EQ[1];
        fStress_EQ[2] = fStress3D_EQ[5];
    }
    else fStress_EQ = fStress3D_EQ;
    
	const dMatrixT& Ftotal = F_total();	
    fStress_EQ *= 1.0/Ftotal.Det();
  //  cout << "\nfStresseq: "<<fStress_EQ;
	return fStress_EQ;
}

const dSymMatrixT& SMP_coupled::s_ijneq(void)
{
 	const dMatrixT& F = MechanicalDeformation();
	double T = Compute_Temperature();
	
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
    /*load the viscoelastic principal stretches from state variable arrays*/
    if (fNumProcess > 0 )
    {
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
             //  			if(CurrElementNumber()==0&&CurrIP()==0)
               //  cout << "\nfEigs_e: "<<fEigs_e;
                
             //   if(CurrElementNumber()==0&&CurrIP()==0)
               //     cout <<setprecision(12)<< "\nfEigs: "<<fEigs;
                double Je = sqrt(fEigs_e.Product());
                fEigs_dev = fEigs_e;
                fEigs_dev *= pow(Je,-2.0*third);
                
                fPot[i+1]->DevStress(fEigs_dev, ftau_NEQ);
                ftau_NEQ += fPot[i+1]->MeanStress(Je);
               /*			if(CurrElementNumber()==0&&CurrIP()==0)
                 cout << "\nftau_NEQ: "<<ftau_NEQ; */
                fStress3D_NEQ = fSpectralDecompSpat.EigsToRank2(ftau_NEQ); 
                
                /*Calculate Cv*/
                
                fInverse = fSpectralDecompSpat.EigsToRank2(fEigs_e); 
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
                fStress3D_NEQ= fSpectralDecompSpat.EigsToRank2(ftau_NEQ);
            }
        }
    }
    
	if (NumSD() == 2)
    {
        fStress_NEQ[0] = fStress3D_NEQ[0];
        fStress_NEQ[1] = fStress3D_NEQ[1];
        fStress_NEQ[2] = fStress3D_NEQ[5];
    }
    else fStress_NEQ = fStress3D_NEQ;
  	const dMatrixT& Ftotal = F_total();
    fStress_NEQ *= 1.0/Ftotal.Det();
    
 //   cout << "\nfStressneq: "<<fStress_NEQ;
	return fStress_NEQ;
}

const dSymMatrixT& SMP_coupled::s_ij(void)
{
    fStress=s_ijeq();
    fStress +=s_ijneq();
    return fStress;
}
/* modulus */
const dMatrixT& SMP_coupled::c_ijkleq(void)
{
    
	double T = Compute_Temperature();
	const dMatrixT& F = MechanicalDeformation();
   /* if(CurrElementNumber()==0&&CurrIP()==0)
        cout<<"\n deformation gradient at the first interatraiton point is "<<F;
    else if(CurrElementNumber()==0&&CurrIP()==1)
        cout<<"\n deformation gradient at the second interatraiton point is "<<F;*/

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
	
    fPot[0]->DevStress(fEigs_dev, ftau_EQ, T);
	ftau_EQ += fPot[0]->MeanStress(J);    
    
	fPot[0]->DevMod(fEigs_dev,fDtauDe_EQ, T);
    fDtauDe_EQ += fPot[0]->MeanMod(J);

    dSymMatrixT& Gamma = fDtauDe_EQ;
    Gamma(0,0) -= 2.0*ftau_EQ[0];
    Gamma(1,1) -= 2.0*ftau_EQ[1];
    Gamma(2,2) -= 2.0*ftau_EQ[2];
   
	fModulus3D_EQ = fSpectralDecompSpat.EigsToRank4(Gamma);
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
    fModulus3D_EQ.AddScaled(2.0*coeff, fModMat);
    
    dl = l0 - l2;
    if (fabs(dl) > kSmall)
      coeff = (ftau_EQ[0]*l2 - ftau_EQ[2]*l0)/dl;
    else 
      coeff = 0.5*(Gamma(0,0)-Gamma(0,2))-ftau_EQ[2];	
    MixedRank4_3D(eigenvectors[0], eigenvectors[2], fModMat);
    fModulus3D_EQ.AddScaled(2.0*coeff, fModMat);
    
    dl = l1 - l2;
   if (fabs(dl) > kSmall)
		coeff  = (ftau_EQ[1]*l2 - ftau_EQ[2]*l1)/dl;
    else
      coeff = 0.5*(Gamma(1,1)-Gamma(1,2))-ftau_EQ[1];	
    MixedRank4_3D(eigenvectors[1], eigenvectors[2], fModMat);
    fModulus3D_EQ.AddScaled(2.0*coeff, fModMat);
	
	if (NumSD() == 2)
	{
		fModulus_EQ[0] = fModulus3D_EQ[0];
		fModulus_EQ[1] = fModulus3D_EQ[1];
		fModulus_EQ[2] = fModulus3D_EQ[5];

		fModulus_EQ[3] = fModulus3D_EQ[6];
		fModulus_EQ[4] = fModulus3D_EQ[7];
		fModulus_EQ[5] = fModulus3D_EQ[11];
		fModulus_EQ[6] = fModulus3D_EQ[30];
		fModulus_EQ[7] = fModulus3D_EQ[31];
		fModulus_EQ[8] = fModulus3D_EQ[35];
	}
	else fModulus_EQ = fModulus3D_EQ;

	const dMatrixT& Ftotal = F_total();	
	fModulus_EQ *= 1.0/Ftotal.Det();
    
 /* cout << "\nfModuluseq: "<<fModulus_EQ; */
    return fModulus_EQ;
}

const dMatrixT& SMP_coupled::c_ijklneq(void)
{
    
	double T = Compute_Temperature();
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
	/*calc NEQ component of stress and moduli*/
	/*calcualte principal values of elastic stretch*/
    
    if (fNumProcess > 0)
    {
        ElementCardT& element = CurrentElement();
        Load(element, CurrIP());
        
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
            Compute_Calg(ftau_NEQ, fDtauDe_NEQ, sm, cm, fCalg,fdalg,i);
            /* cout << "\nfCalg: "<<fCalg; */
            
            ftau_NEQ += sm;
            fDtauDe_NEQ += cm;
            
            fModulus3D_NEQ = fSpectralDecompSpat.NonSymEigsToRank4(fCalg);
            
            double dl_tr,coeff;
            
            double& l0_tr = fEigs_tr[0];
            double& l1_tr = fEigs_tr[1];
            double& l2_tr = fEigs_tr[2];
            
            
            dl_tr = l0_tr - l1_tr;
            if (fabs(dl_tr) > kSmall)
                coeff = (ftau_NEQ[0]*l1_tr - ftau_NEQ[1]*l0_tr)/dl_tr;
            else
                coeff = 0.5*(fCalg(0,0)-fCalg(0,1))-ftau_NEQ[0];
            MixedRank4_3D(eigenvectors_e[0], eigenvectors_e[1], fModMat);
            fModulus3D_NEQ.AddScaled(2.0*coeff, fModMat);
            
            dl_tr = l0_tr - l2_tr;
            if (fabs(dl_tr) > kSmall)
                coeff =(ftau_NEQ[0]*l2_tr - ftau_NEQ[2]*l0_tr)/dl_tr;
            else
                coeff = 0.5*(fCalg(0,0)-fCalg(0,2))-ftau_NEQ[2];
            MixedRank4_3D(eigenvectors_e[0], eigenvectors_e[2], fModMat);
            fModulus3D_NEQ.AddScaled(2.0*coeff, fModMat);
            
            dl_tr = l1_tr - l2_tr;
            if (fabs(dl_tr) > kSmall)
                coeff  = (ftau_NEQ[1]*l2_tr - ftau_NEQ[2]*l1_tr)/dl_tr;
            else
                coeff = 0.5*(fCalg(1,1)-fCalg(1,2))-ftau_NEQ[1];	
            MixedRank4_3D(eigenvectors_e[1], eigenvectors_e[2], fModMat);
            fModulus3D_NEQ.AddScaled(2.0*coeff, fModMat);
        }
    }
	if (NumSD() == 2)
	{
		fModulus_NEQ[0] = fModulus3D_NEQ[0];
		fModulus_NEQ[1] = fModulus3D_NEQ[1];
		fModulus_NEQ[2] = fModulus3D_NEQ[5];
        
		fModulus_NEQ[3] = fModulus3D_NEQ[6];
		fModulus_NEQ[4] = fModulus3D_NEQ[7];
		fModulus_NEQ[5] = fModulus3D_NEQ[11];
		fModulus_NEQ[6] = fModulus3D_NEQ[30];
		fModulus_NEQ[7] = fModulus3D_NEQ[31];
		fModulus_NEQ[8] = fModulus3D_NEQ[35];
	}
	else fModulus_NEQ = fModulus3D_NEQ;
    
	const dMatrixT& Ftotal = F_total();	
	fModulus_NEQ *= 1.0/Ftotal.Det();
/*    cout << "\nfModulusneq: "<<fModulus_NEQ; */
    return fModulus_NEQ;
}

const dMatrixT& SMP_coupled::c_ijkl(void)
{
    fModulus=c_ijkleq();
    fModulus +=c_ijklneq();
/*    cout << "\n fModulus: "<<fModulus; */
    return fModulus;
}

const dSymMatrixT& SMP_coupled::d_ij(void)
{
    double T = Compute_Temperature();
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
	/*calc NEQ component of stress and moduli*/
	/*calcualte principal values of elastic stretch*/
    
    if (fNumProcess > 0)
    {
        ElementCardT& element = CurrentElement();
        Load(element, CurrIP());
        
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
            Compute_Calg(ftau_NEQ, fDtauDe_NEQ, sm, cm, fCalg,fdalg,i);
            
            fCoupledModulus = fSpectralDecompSpat.EigsToRank2(fdalg);

        }
    }
    
	const dMatrixT& Ftotal = F_total();
	fCoupledModulus *= 1.0/Ftotal.Det();
  //  if(CurrElementNumber()==0&&CurrIP()==0)
   // cout<<"\n fCoupledModulus is"<<fCoupledModulus;
    return fCoupledModulus;
}


double SMP_coupled::a1(const double Temperature, const double Te, const double smag)
{
   	double Temp=Temperature;
    const dMatrixT& F = F_total();
    double iJ = 1.0/F.Det();
    double tauR = RetardationTime(Temp, Te);
	double etaS = ShearViscosity(Temp,Te,smag);
    double dtauR_dTe = fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Te - fTg)*(fC2 + Te - fTg));
    double dtauR_dT=fC1*fC2*Te*log(10)/(Temp*Temp)/(fC2+Te-fTg);
    double detaS_dT=fC1*fC2*Te*log(10)/(Temp*Temp)/(fC2+Te-fTg);
    double detaS_dTe=fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Te - fTg)*(fC2 + Te - fTg));
    double dt = fFSMatSupport->TimeStep();
    if(dt<small)
    {
        dt=small;
    }
	double itauR = 1.0/tauR;
	double ietaS = 1.0/etaS;
    double temp1=itauR-(Te-Temp)*itauR*dtauR_dT+fA*smag*smag*ietaS/iJ*detaS_dT/fdeltac;
    double temp2=1.0/dt+itauR-(Te-Temp)*itauR*dtauR_dTe+fA*smag*smag*ietaS/iJ*detaS_dTe/fdeltac;
    double temp3=temp1/temp2;
    return temp3;
}

double SMP_coupled::a2(const double Temperature, const double Te, const double smag)

{
    double Temp=Temperature;
    const dMatrixT& F = F_total();
   
    double iJ = 1.0/F.Det();
    double tauR = RetardationTime(Temp, Te);
	double etaS = ShearViscosity(Temp,Te,smag);
    double dtauR_dTe = fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Te - fTg)*(fC2 + Te - fTg));
    double detaS_dTe=fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Te - fTg)*(fC2 + Te - fTg));
    double x = fQS*smag/(Temp*fsy0);
    double detaS_dsmag = 0.0;
    if (smag/fsy0 > small)
    {
        double cothx = cosh(x)/sinh(x);
        detaS_dsmag = (1.0-x*cothx)/smag;
    }
    
    double dt = fFSMatSupport->TimeStep();
    if(dt<small)
    {
        dt=small;
    }
	double itauR = 1.0/tauR;
	double ietaS = 1.0/etaS;
    double temp1=fA*(-1.0*smag*smag*detaS_dsmag+2.0*smag)*ietaS/iJ/fdeltac;
    double temp2=1.0/dt+itauR-(Te-Temp)*itauR*dtauR_dTe+fA*smag*smag*ietaS/iJ*detaS_dTe/fdeltac;
    double temp3=temp1/temp2;
    return temp3;
}

double SMP_coupled::b1(void)
{
    double Temp=Compute_Temperature();
    double& Te = fTe[0];
    const dMatrixT& F = F_total();
    double iJ = 1.0/F.Det();
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
    double temp1=0.0;
    /*get smag */
    if (fNumProcess > 0)
    {
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
     for (int i = 0; i < fNumProcess; i++)
  {
    fInverse.Inverse(fC_v[i]);
    fbe.MultQBQT(fF3D, fInverse);
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);
    fEigs_e = fSpectralDecompSpat.Eigenvalues();

        /*calc NEQ component of stress and moduli*/
                    /*calc trial state*/
    double Je = sqrt(fEigs_e.Product());
    fEigs_dev = fEigs_e;
    fEigs_dev *= pow(Je,-2.0*third);
    fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
     double s0 = iJ*ftau_NEQ[0];
    double s1 = iJ*ftau_NEQ[1];
    double s2 = iJ*ftau_NEQ[2];
    
    
    /*caculate smag*/
    double smag = sqrt(0.5*(s0*s0 + s1*s1 + s2*s2));
     /* end */
    double tauR = RetardationTime(Temp, Te);
	double etaS = ShearViscosity(Temp,Te, smag);
    double detaS_dTe = fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Te - fTg)*(fC2 + Te - fTg));
    double detaS_dT=fC1*fC2*Te*log(10)/(Temp*Temp)/(fC2+Te-fTg);
        
    double dt = fFSMatSupport->TimeStep();
      if(dt<small)
      {
          dt=small;
      }
	double itauR = 1.0/tauR;
	double ietaS = 1.0/etaS;
    double fa1=a1(Temp,Te,smag);
    temp1=-smag*smag*ietaS*detaS_dT/iJ+(frho*fdeltac/dt+smag*smag*ietaS*detaS_dTe/iJ)*fa1;
   }
    }
    return temp1;
}
double SMP_coupled::b2(const double Temperature, const double Te, const double smag)
{
    double Temp=Temperature;
    const dMatrixT& F = F_total();
    double iJ = 1.0/F.Det();
    double tauR = RetardationTime(Temp, Te);
	double etaS = ShearViscosity(Temp,Te, smag);
    double detaS_dTe = fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Te - fTg)*(fC2 + Te - fTg));
    double detaS_dsmag = 0.0;
    double x = fQS*smag/(Temp*fsy0);
    if (smag/fsy0 > small)
    {
        double cothx = cosh(x)/sinh(x);
        detaS_dsmag = (1.0-x*cothx)/smag;
    }

    double dt = fFSMatSupport->TimeStep();
    if(dt<small)
    {
        dt=small;
    }
	double itauR = 1.0/tauR;
	double ietaS = 1.0/etaS;
    double fa2=a2(Temp,Te,smag);
    double temp1=(frho*fdeltac/dt+smag*smag*ietaS*detaS_dTe/iJ)*fa2-ietaS*(2.0*smag-smag*smag*detaS_dsmag)/iJ;
    return temp1;
}

const dSymMatrixT& SMP_coupled::h_ij(void)
{
    
    const dSymMatrixT& sijneq=s_ijneq();
        sijneq.ToMatrix(fneqstress);
    const dMatrixT& cijklneq=c_ijklneq();
    Hij.MultAAT(fneqstress);
    Hij *=2.0;
   /*the follwing is addingh sijneq:cijklneq */
    Hij[0] +=sijneq[0]*cijklneq(0,0)+sijneq[1]*cijklneq(1,0)+sijneq[2]*cijklneq(2,0)+2.0*sijneq[3]*cijklneq(3,0)+2.0*sijneq[4]*cijklneq(4,0)+2.0*sijneq[5]*cijklneq(5,0);
    Hij[1] +=sijneq[0]*cijklneq(0,1)+sijneq[1]*cijklneq(1,1)+sijneq[2]*cijklneq(2,1)+2.0*sijneq[3]*cijklneq(3,1)+2.0*sijneq[4]*cijklneq(4,1)+2.0*sijneq[5]*cijklneq(5,1);
    Hij[2] +=sijneq[0]*cijklneq(0,2)+sijneq[1]*cijklneq(1,2)+sijneq[2]*cijklneq(2,2)+2.0*sijneq[3]*cijklneq(3,2)+2.0*sijneq[4]*cijklneq(4,2)+2.0*sijneq[5]*cijklneq(5,2);
    Hij[3] +=sijneq[0]*cijklneq(0,3)+sijneq[1]*cijklneq(1,3)+sijneq[2]*cijklneq(2,3)+2.0*sijneq[3]*cijklneq(3,3)+2.0*sijneq[4]*cijklneq(4,3)+2.0*sijneq[5]*cijklneq(5,3);
    Hij[4] +=sijneq[0]*cijklneq(0,4)+sijneq[1]*cijklneq(1,4)+sijneq[2]*cijklneq(2,4)+2.0*sijneq[3]*cijklneq(3,4)+2.0*sijneq[4]*cijklneq(4,4)+2.0*sijneq[5]*cijklneq(5,4);
    Hij[5] +=sijneq[0]*cijklneq(0,5)+sijneq[1]*cijklneq(1,5)+sijneq[2]*cijklneq(2,5)+2.0*sijneq[3]*cijklneq(3,5)+2.0*sijneq[4]*cijklneq(4,5)+2.0*sijneq[5]*cijklneq(5,5);
    double Temp=Compute_Temperature();
    double& Te = fTe[0];
    /*get smag */
   double smag=0.0;
    const dMatrixT& Ftotal = F_total();
    double iJ = 1.0/Ftotal.Det();
    const dMatrixT& F =  MechanicalDeformation();
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
    if (fNumProcess > 0)
    {
        ElementCardT& element = CurrentElement();
        Load(element, CurrIP());
        for (int i = 0; i < fNumProcess; i++)
        {
            fInverse.Inverse(fC_v[i]);
            fbe.MultQBQT(fF3D, fInverse);
            fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);
            fEigs_e = fSpectralDecompSpat.Eigenvalues();
            
            /*calc NEQ component of stress and moduli*/
            /*calc trial state*/
           double Je = sqrt(fEigs_e.Product());
            fEigs_dev = fEigs_e;
            fEigs_dev *= pow(Je,-2.0*third);
            fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
            double s0 = iJ*ftau_NEQ[0];
            double s1 = iJ*ftau_NEQ[1];
            double s2 = iJ*ftau_NEQ[2];
            
            
            /*caculate smag*/
       smag = sqrt(0.5*(s0*s0 + s1*s1 + s2*s2));
        }
    }
    double fb2=b2(Temp,Te, smag);
    Hij *=1.0/smag*fb2;
    Hij=0.0;
    return Hij;
}


double SMP_coupled::heatres(void)
{
    double Temp=Compute_Temperature();
    double& Te = fTe[0];
    double& Ten=fTe_n[0];

    const dMatrixT& F = F_total();
    double iJ = 1.0/F.Det();
    
    /*get smag */
    double smag=0.0;
    if (fNumProcess > 0)
    {
        ElementCardT& element = CurrentElement();
        Load(element, CurrIP());
        for (int i = 0; i < fNumProcess; i++)
        {
            fInverse.Inverse(fC_v[i]);
    fbe.MultQBQT(fF3D, fInverse);
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);
    fEigs_e = fSpectralDecompSpat.Eigenvalues();
    
    /*calc NEQ component of stress and moduli*/
    /*calc trial state*/
    double Je = sqrt(fEigs_e.Product());
    fEigs_dev = fEigs_e;
    fEigs_dev *= pow(Je,-2.0*third);
    fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
    double s0 = iJ*ftau_NEQ[0];
    double s1 = iJ*ftau_NEQ[1];
    double s2 = iJ*ftau_NEQ[2];
    
    
    /*caculate smag*/
   smag = sqrt(0.5*(s0*s0 + s1*s1 + s2*s2));
        }
    }
    /* end */
	double etaS = ShearViscosity(Temp,Te, smag);
    double dt = fFSMatSupport->TimeStep();
    if(dt<small)
    {
        dt=small;
    }
	double ietaS = 1.0/etaS;
    double temp1=frho*fdeltac*(Te-Ten)/dt-smag*smag*ietaS/iJ;

    return temp1;
}

void SMP_coupled::ComputeOutput(dArrayT& output)
{
	/*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    
	double Temp = Compute_Temperature();
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
	
	/*thermal volume deformation*/
    double dalpha=falphar-falphag;
	output[0] = exp(falphag*(Temp-fT0)+dalpha*(fTe[0]-fT0));
    
	/*neq thermal dilatation*/
	output[1] = (fTe[0]);
    
	/*maximum eigenvalue of mechanical stretch*/
	fInverse.MultATA(fF3D);
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fInverse, false);
    fEigs = fSpectralDecompSpat.Eigenvalues();
	
	output[2] = sqrt(fEigs[0]);
	output[3] = sqrt(fEigs[1]);
	output[4] = sqrt(fEigs[2]);
	
	/*maximum eigenvalue of viscous stretch*/
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fC_v[0], false);
    fEigs = fSpectralDecompSpat.Eigenvalues();
	output[5] = sqrt(fEigs[0]);
	output[6] = sqrt(fEigs[1]);
	output[7] = sqrt(fEigs[2]);
	
	
	/*yield strength*/
	output[8] = Temp;
    
	/*calc elastic stretch*/
	fInverse.Inverse(fC_v[0]);
	fbe.MultQBQT(fF3D, fInverse);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);
	fEigs_e = fSpectralDecompSpat.Eigenvalues();
    
    
	double Je = sqrt(fEigs_e.Product());
	fEigs_dev = fEigs_e;
	fEigs_dev *= pow(Je,-2.0*third);
	fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
	
	const dMatrixT& Ftot = F_total();
	double J = Ftot.Det();
	
	double s0 = ftau_NEQ[0]/J;
	double s1 = ftau_NEQ[1]/J;
	double s2 = ftau_NEQ[2]/J;
	double smag = sqrt(0.5*(s0*s0 + s1*s1 + s2*s2));
    
	/*calculate mobilities*/
	
    double etaS = ShearViscosity(Temp,fTe[0], smag);
	double ietaS = 1.0/etaS;
    
	double gamdot = 0.5*smag*ietaS;
	output[9] = gamdot;
}

/*************************************************************************
*	PUBLIC
**************************************************************************/
/* describe the parameters needed by the interface */
void SMP_coupled::DefineParameters(ParameterListT& list) const
{
  /* inherited */
  RGViscoelasticityT::DefineParameters(list);

  /* common limit */
  LimitT positive(0.0, LimitT::Lower);

  ParameterT reftemp(ParameterT::Double, "ref_temperature");
  list.AddParameter(reftemp);

}

/* information about subordinate parameter lists */
void SMP_coupled::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("smp_thermal_expansion", ParameterListT::Once);
    sub_list.AddSub("smp_heat_conduction", ParameterListT::Once);

	sub_list.AddSub("smp_eq_potential", ParameterListT::Once);
	sub_list.AddSub("smp_neq_potential", ParameterListT::Once);

	/* choice of viscosity */
	sub_list.AddSub("smp_retardation_time", ParameterListT::Once);
	sub_list.AddSub("smp_coupled_shear_viscosity", ParameterListT::Once);
}


/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SMP_coupled::NewSub(const StringT& name) const
{
		LimitT zero(0.0, LimitT::Lower);
		LimitT one(1.0, LimitT::Lower);
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
	else if (name == "smp_thermal_expansion")
	{
		ParameterContainerT* CTE = new ParameterContainerT(name);
		
		ParameterT alphar(ParameterT::Double, "high_temp_CTE");
		ParameterT alphag(ParameterT::Double, "low_temp_CTE");

		alphar.AddLimit(zero);
		alphag.AddLimit(zero);

		CTE->AddParameter(alphar);
		CTE->AddParameter(alphag);
		return CTE;
	}
    else if (name == "smp_heat_conduction")
	{
		ParameterContainerT* heat = new ParameterContainerT(name);
		
		ParameterT rho(ParameterT::Double, "polymer_density");
		ParameterT cg(ParameterT::Double, "instant_heat_capacity");
        ParameterT deltac(ParameterT::Double, "time_dependent_heat_capacity");
        ParameterT kconductivity(ParameterT::Double, "thermal_conductivity");
        ParameterT fraction(ParameterT::Double, "energy_fraction");
        
        
		heat->AddParameter(rho);
		heat->AddParameter(cg);
        heat->AddParameter(deltac);
		heat->AddParameter(kconductivity);
        heat->AddParameter(fraction);
		return heat;
	}
	else if (name == "smp_eq_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		choice->SetDescription("temperature normalized network stiffness");
		
		/* choice of parameters */
		choice->AddSub("arruda-boyce");
		return(choice);
	}
	else if (name == "smp_neq_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
	
		/* choice of parameters */
		choice->AddSub("neo-hookean");
		return(choice);
	}
	else if (name == "smp_retardation_time")
	{
		ParameterContainerT* tauR = new ParameterContainerT(name);

		ParameterT tauRg(ParameterT::Double, "tauR_ref");
		ParameterT Tg(ParameterT::Double, "Tg");
		ParameterT C1(ParameterT::Double, "WLF_C1");
		ParameterT C2(ParameterT::Double, "WLF_C2");
		
		tauRg.AddLimit(zero);
		Tg.AddLimit(zero);

		tauR->AddParameter(tauRg);
		tauR->AddParameter(Tg);
		tauR->AddParameter(C1);
		tauR->AddParameter(C2);

		return(tauR);
	}
	else if (name == "smp_coupled_shear_viscosity")
	{
		ParameterContainerT* etaS = new ParameterContainerT(name);

		ParameterT etaSR(ParameterT::Double, "etaS_ref");
		ParameterT A(ParameterT::Double, "activation_energy");
		ParameterT sy_0(ParameterT::Double, "init_yield_strength");
			
		etaSR.AddLimit(zero);
		A.AddLimit(zero);
		sy_0.AddLimit(zero);

		etaS->AddParameter(etaSR);
		etaS->AddParameter(A);
		etaS->AddParameter(sy_0);

		return(etaS);
	}
}
/*rewrite this funciton RX */
void SMP_coupled::TakeParameterList(const ParameterListT& list)
{
  const char caller[] = "SMP_coupled::TakeParameterList";
  /* inherited */
   // cout << "\n first time SMP_coupled:::TakeParmeterList()"<<endl;

  FSSolidMatT::TakeParameterList(list);
  // cout << "\n second time SMP_coupled:::TakeParmeterList()"<<endl;

  fT0 = list.GetParameter("ref_temperature");

  const ParameterListT* tm = list.List("smp_thermal_expansion");
  if (tm)
  {
	falphar = tm->GetParameter("high_temp_CTE");
	falphag = tm->GetParameter("low_temp_CTE");
  }
    
   const ParameterListT* heat = list.List("smp_heat_conduction");
    if (heat)
    {
        frho=heat->GetParameter("polymer_density");
        fcg=heat->GetParameter("instant_heat_capacity");
        fdeltac=heat->GetParameter("time_dependent_heat_capacity");
        fk=heat->GetParameter("thermal_conductivity");
        fA=heat->GetParameter("energy_fraction");
     //   fSpecificHeat=fcg;
      //  fkij.Identity(fk);
    }
    
	fPot.Dimension(2);
		
	const ParameterListT& eq_pot = list.GetListChoice(*this, "smp_eq_potential");
	if(eq_pot.Name() == "arruda-boyce")
		fPot[0] = new ArrudaBoyce;
	else 
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot[0]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", eq_pot.Name().Pointer());			
	fPot[0]->TakeParameterList(eq_pot);
  
	const ParameterListT& neq_pot = list.GetListChoice(*this, "smp_neq_potential");
	if(neq_pot.Name() == "neo-hookean")
		fPot[1] = new NeoHookean;
	else 
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot[1]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", neq_pot.Name().Pointer());			
	fPot[1]->TakeParameterList(neq_pot);
  
  const ParameterListT* tauR = list.List("smp_retardation_time");
  if (tauR)
  {
	fTg = tauR->GetParameter("Tg");
	ftaug = tauR->GetParameter("tauR_ref");
	fC1 = tauR->GetParameter("WLF_C1");
	fC2 = tauR->GetParameter("WLF_C2");
	
	ftauRL = 1.0e-10*ftaug;
	ftauRH = 1.0e+10*ftaug;
  }
  
  const ParameterListT* etaS = list.List("smp_coupled_shear_viscosity");
  if (etaS)
  {
	fetaS0 = etaS->GetParameter("etaS_ref");
	fetaSL = 1.0e-10*fetaS0;
	fetaSH = 1.0e+10*fetaS0;
	fQS = etaS->GetParameter("activation_energy");
	fsy0 = etaS->GetParameter("init_yield_strength");
  }
  
	/*set dimension of workspaces*/
	Initialize();
  //  cout << "\n done with SMP_coupled:::TakeParmeterList()"<<endl;

}

/*retrieves and returns the effective temperature at the current integration point*/
const dArrayT&  SMP_coupled::Effective_Temperature(void)
{
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    
    return(fTe);
}


const dArrayT&  SMP_coupled::Effective_Temperature_last(void)
{
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    
    return(fTe_n);
}

/*initializes history variable */
void  SMP_coupled::PointInitialize(void)
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
			  double Temp = Compute_Temperature();
		      
			  fC_v[0].Identity();
			  fC_vn[0].Identity();

             fTe[0]=fT0;
            fTe_n[0]=fT0;

            /* write to storage */
		      Store(element, ip);
		}
	}
}
 
void SMP_coupled::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		fC_vn[0] = fC_v[0];

		 fTe_n[0]=fTe[0];
		/* write to storage */
		Store(element, ip);
	}
}

void SMP_coupled::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		fC_v[0] = fC_vn[0];
		fTe[0]=fTe_n[0];

		/* write to storage */
		Store(element, ip);
	}
}


int SMP_coupled::NumOutputVariables() const {return kNumOutputVar;} 

void SMP_coupled::OutputLabels(ArrayT<StringT>& labels) const
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
void SMP_coupled::Initialize(void)
{
 /* dimension work space */
  
	/*Dimension workspace*/
	fC_v.Dimension(1);
	fC_vn.Dimension(1);
	
	int nsd = NumSD();
	int ndof = 3;
	int numstress = dSymMatrixT::NumValues(ndof);

	fnstatev = 0;
	fnstatev += numstress;   /*current C_v*/
	fnstatev += numstress;   /*last C_vn*/
	fnstatev ++;			/*current fictive temperature*/ 
	fnstatev ++;			/*last fictive temperature*/
	//fnstatev ++;			/*current yield strength*/
	//fnstatev ++;			/*last yield strength*/
	
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
		
	/* assign pointers to current and last blocks of state variable array */
	fC_v[0].Set(ndof, pstatev);
	pstatev += numstress;
	fC_vn[0].Set(ndof, pstatev);
	pstatev += numstress;
    fTe.Set(1, pstatev);
    pstatev++;
    fTe_n.Set(1, pstatev);
    pstatev++;


  fF_M.Dimension(nsd);
  fF_T_inv.Dimension(nsd);

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
    
  fStress_EQ.Dimension(NumSD());
  fStress3D_EQ.Dimension(ndof);
    
  fStress_NEQ.Dimension(NumSD());
  fneqstress.Dimension(NumSD());
  fStress3D_NEQ.Dimension(ndof);

  fStress.Dimension(NumSD());
  fStress3D.Dimension(ndof);

  fDtauDe_EQ.Dimension(ndof);
  fDtauDe_NEQ.Dimension(ndof);

   fRes.Dimension(ndof+1);
   fDelta.Dimension(ndof+1);
    fdalg.Dimension(ndof);
   
   fiKAB.Dimension(ndof+1);
   fGAB.Dimension(ndof+1,ndof+1);
   fDAB.Dimension(ndof+1,ndof+1);
   fDABbar.Dimension(ndof);
   fDABbar2.Dimension(ndof);
   fMat.Dimension(ndof);
    fMat2.Dimension(ndof);
   fCalg.Dimension(ndof);
    fCoupledModulus.Dimension(nsd);
    Hij.Dimension(nsd);
    
   fModulus3D_EQ.Dimension(dSymMatrixT::NumValues(ndof));
   fModulus_EQ.Dimension(dSymMatrixT::NumValues(NumSD()));
    
    fModulus3D_NEQ.Dimension(dSymMatrixT::NumValues(ndof));
    fModulus_NEQ.Dimension(dSymMatrixT::NumValues(NumSD()));
    
  fModulus3D.Dimension(dSymMatrixT::NumValues(ndof));
  fModMat.Dimension(dSymMatrixT::NumValues(ndof));
  fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));
    fHeatFlux.Dimension(nsd);
    fkij.Dimension(nsd);

}


/***********************************************************************
 * Private
 ***********************************************************************/
/* set inverse of thermal transformation - return true if active */
 void SMP_coupled::Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
	const double& dtau_m, dMatrixT& Calg, dArrayT& dalg,const int type)
 {
		const dMatrixT& F = F_total();
		double iJ = 1.0/F.Det();

        double& Te = fTe[0];
        double& Ten=fTe_n[0];

		/*temperature and temperature step*/
		double Temp = Compute_Temperature();
		double Temp_n = Compute_Temperature_last();
		double dT = Temp - Temp_n;
	
		/*time step*/
		double dt = fFSMatSupport->TimeStep();
        if(dt<small)
        {
         dt=small;
        }
		double dalpha = falphar-falphag;


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
	    double smag = sqrt(0.5*(s0*s0 + s1*s1 + s2*s2));

		/*calculate mobilities*/

		double tauR = RetardationTime(Temp, Te);
		double etaS = ShearViscosity(Temp, Te, smag);
		double itauR = 1.0/tauR;
		double ietaS = 1.0/etaS;
						
		double gamdot = 0.5*smag*ietaS;

		/*calculate stiffness matrix*/
		/*derivative of retardation time wrt to Te*/
        double dtauR_dTe = fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Te - fTg)*(fC2 + Te - fTg));
        double dtauR_dT=fC1*fC2*Te*log(10)/(Temp*Temp)/(fC2+Te-fTg);
        double detaS_dTe = fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Te - fTg)*(fC2 + Te - fTg));
        double detaS_dT=fC1*fC2*Te*log(10)/(Temp*Temp)/(fC2+Te-fTg);
		double x = fQS*smag/(Temp*fsy0);
		double detaS_dsmag = 0.0;
	/*	double detaS_dsy = (etaS)/fsy0; */
		if (smag/fsy0 > small)
		{
			double cothx = cosh(x)/sinh(x);
			detaS_dsmag = (1.0-x*cothx)/smag;
		}
		
     /*initialize*/
     fiKAB = 0.0;
     
     /*K_del_del*/
     fiKAB(0,0)=(1.0+dt*itauR)-dt*itauR*dtauR_dTe*(Te-Temp)+fA*smag*smag*dt*ietaS/iJ*detaS_dTe/fdeltac+fA/iJ*ietaS*(sqrt(2.0)*smag-smag*smag*detaS_dsmag)*dalpha/fdeltac;
     /*K_epA_del*/
     fiKAB(1,0) = -0.5*dt*ietaS*s0*detaS_dTe+0.5*dt*ietaS*s0*detaS_dsmag*dalpha*smag-0.5*dt*ietaS*s0*dalpha;
     fiKAB(2,0) = -0.5*dt*ietaS*s1*detaS_dTe+0.5*dt*ietaS*s1*detaS_dsmag*dalpha*smag-0.5*dt*ietaS*s1*dalpha;
     fiKAB(3,0) = -0.5*dt*ietaS*s2*detaS_dTe+0.5*dt*ietaS*s2*detaS_dsmag*dalpha*smag-0.5*dt*ietaS*s2*dalpha;
     
     /*K_epA_epB*/
     double coef0 = 0.5*(s0*c0 + s1*c01 + s2*c02);
     double coef1 = 0.5*(s0*c01 + s1*c1 + s2*c12);
     double coef2 = 0.5*(s0*c02 + s1*c12 + s2*c2);
     if (smag/fsy0 > small)
     {
         coef0 /= smag;
         coef1 /= smag;
         coef2 /= smag;
     }
     fiKAB(1,1) = 1.0 + 0.5*ietaS*dt*c0 - 0.5*dt*ietaS*s0*detaS_dsmag*coef0;
     fiKAB(2,2) = 1.0 + 0.5*ietaS*dt*c1 - 0.5*dt*ietaS*s1*detaS_dsmag*coef1;
     fiKAB(3,3) = 1.0 + 0.5*ietaS*dt*c2 - 0.5*dt*ietaS*s2*detaS_dsmag*coef2;
     
     fiKAB(2,3) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s1*detaS_dsmag*coef2;
     fiKAB(1,3) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s0*detaS_dsmag*coef2;
     fiKAB(1,2) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s0*detaS_dsmag*coef1;
     
     fiKAB(3,2) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s2*detaS_dsmag*coef1;
     fiKAB(3,1) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s2*detaS_dsmag*coef0;
     fiKAB(2,1) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s1*detaS_dsmag*coef0;
     /*check the below derivitions*/
     fiKAB(0,1)=-2.0*fA/iJ*dt*ietaS*coef0*smag/fdeltac+fA/iJ*dt*ietaS*smag*smag*detaS_dsmag*coef0/fdeltac;
     fiKAB(0,2)=-2.0*fA/iJ*dt*ietaS*coef1*smag/fdeltac+fA/iJ*dt*ietaS*smag*smag*detaS_dsmag*coef1/fdeltac;
     fiKAB(0,3)=-2.0*fA/iJ*dt*ietaS*coef2*smag/fdeltac+fA/iJ*dt*ietaS*smag*smag*detaS_dsmag*coef2/fdeltac;
     
 /*     cout << "\nfiKAB: "<<fiKAB; */

	 /*inverts KAB*/
		fiKAB.Inverse();

		/*initialize*/
		fGAB = 0.0;
     fGAB(0,0)=dt*itauR-dt*(Te-Temp)*itauR*dtauR_dT+fA*smag*smag*dt*ietaS/iJ*detaS_dT/fdeltac-fA/iJ*ietaS*(sqrt(2.0)*smag-smag*smag*detaS_dsmag)*falphag/fdeltac;;
     fGAB(0,1)=(-iJ*fA*ietaS*smag*smag*dt+iJ*fA*ietaS*smag*smag*smag*detaS_dsmag*dt)/fdeltac;
     fGAB(0,2)=(-iJ*fA*ietaS*smag*smag*dt+iJ*fA*ietaS*smag*smag*smag*detaS_dsmag*dt)/fdeltac;
     fGAB(0,3)=(-iJ*fA*ietaS*smag*smag*dt+iJ*fA*ietaS*smag*smag*smag*detaS_dsmag*dt)/fdeltac;
     
		/*G_epeA_epB*/
        fGAB(1,0)=-0.5*dt*ietaS*s0*detaS_dT-0.5*dt*ietaS*s0*detaS_dsmag*falphag*smag+0.5*dt*ietaS*s0*falphag;
		fGAB(1,1) = 1.0 + 0.5*dt*ietaS*s0 - 0.5*dt*ietaS*s0*detaS_dsmag*smag;
		fGAB(1,2) = 0.5*dt*ietaS*s0 - 0.5*dt*ietaS*s0*detaS_dsmag*smag;
		fGAB(1,3) = 0.5*dt*ietaS*s0 - 0.5*dt*ietaS*s0*detaS_dsmag*smag;
		
        fGAB(2,0)=-0.5*dt*ietaS*s1*detaS_dT-0.5*dt*ietaS*s1*detaS_dsmag*falphag*smag+0.5*dt*ietaS*s1*falphag;
		fGAB(2,1) = 0.5*dt*ietaS*s1 - 0.5*dt*ietaS*s1*detaS_dsmag*smag;
		fGAB(2,2) = 1.0 + 0.5*dt*ietaS*s1 - 0.5*dt*ietaS*s1*detaS_dsmag*smag;
		fGAB(2,3) = 0.5*dt*ietaS*s1 - 0.5*dt*ietaS*s1*detaS_dsmag*smag;

		 fGAB(3,0)=-0.5*dt*ietaS*s2*detaS_dT-0.5*dt*ietaS*s2*detaS_dsmag*falphag*smag+0.5*dt*ietaS*s2*falphag;
         fGAB(3,1) = 0.5*dt*ietaS*s2 - 0.5*dt*ietaS*s2*detaS_dsmag*smag;
		 fGAB(3,2) = 0.5*dt*ietaS*s2 - 0.5*dt*ietaS*s2*detaS_dsmag*smag;
		 fGAB(3,3) = 1.0 + 0.5*dt*ietaS*s2 - 0.5*dt*ietaS*s2*detaS_dsmag*smag;

	 /*Calg = dtau/depe*fiKA*fG	*/
		/*calculating delta_internval_vars = K^-1.G. delta_epsilon*/
		fDAB.MultAB(fiKAB,fGAB);

     /*copy subset*/
/*	 cout << "\nfDAB: "<<fDAB; */
		for (int i = 0; i< fDABbar.Rows(); i++)
			for (int j = 0; j< fDABbar.Cols(); j++)
				fDABbar(i,j) = fDAB(i+1,j+1);
     for (int j = 0; j< fDABbar.Cols(); j++)
         fDABbar2[j]=fDAB(j+1,0);
     
	/*	cout << "\nfDABbar: "<<fDABbar; */
     
		dtau_dev.ToMatrix(fMat);
        dtau_dev.ToMatrix(fMat2);
		fMat(0,0) += dtau_m;
		fMat(1,1) += dtau_m;
		fMat(2,2) += dtau_m;

        Calg.MultAB(fMat, fDABbar);
/*     cout << "\nCalg: "<<Calg; */
     /*get dalg 
       dalg.MultAB(fMat2,fDABbar2); */
     dalg[0]=fMat2(0,0)*fDABbar2[0]+fMat2(0,1)*fDABbar2[1]+fMat2(0,2)*fDABbar2[2];
     dalg[1]=fMat2(1,0)*fDABbar2[0]+fMat2(1,1)*fDABbar2[1]+fMat2(1,2)*fDABbar2[2];
     dalg[2]=fMat2(2,0)*fDABbar2[0]+fMat2(2,1)*fDABbar2[1]+fMat2(2,2)*fDABbar2[2];
		Calg(0,0) -= 2.0* tau_dev[0];
		Calg(1,1) -= 2.0* tau_dev[1];
		Calg(2,2) -= 2.0* tau_dev[2];
 /*     cout << "\nCalg: "<<Calg; */
}

void SMP_coupled::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
			     dArrayT& eigenstress, dSymMatrixT& eigenmodulus,  const int type) 
{		
	const double ctol = 1.00e-10;
		
	/*set references to principle stretches*/
     
	double& le0 = eigenstretch_e[0];
	double& le1 = eigenstretch_e[1];
	double& le2 = eigenstretch_e[2];
	
	double& Te = fTe[0];
    double& Ten = fTe_n[0];
    
	double tol;

	/*initialize principle elastic and trial elastic log strains */
	const double ep_tr0 = 0.5*log(le0);
	const double ep_tr1 = 0.5*log(le1);
	const double ep_tr2 = 0.5*log(le2);
	double ep_e0 = ep_tr0;		
	double ep_e1 = ep_tr1;	
	double ep_e2 = ep_tr2;

	/*jacobian*/
	const dMatrixT& F = F_total();
	double iJ = 1.0/F.Det();

	/*time step*/
	double dt = fFSMatSupport->TimeStep();
    if(dt<small)
    {
        dt=small;
    }
	/*temperature and temperature step*/
	double Temp = Compute_Temperature();
	double Temp_n = Compute_Temperature_last();
	double dT = Temp - Temp_n;

	double dalpha = falphar-falphag;
	int maxiteration = 100;

	/*initializes principle viscous stretch*/
	double Je=sqrt(le0*le1*le2);
	fEigs_dev = eigenstretch_e;
	fEigs_dev *= pow(Je,-2.0*third);

	/*calculate stresses and moduli*/
	fPot[1]->DevStress(fEigs_dev, eigenstress);
	    
	double s0 = iJ*eigenstress[0];
	double s1 = iJ*eigenstress[1];
	double s2 = iJ*eigenstress[2];
		    		
	/*caculate smag*/
	double smag = sqrt(0.5*(s0*s0 + s1*s1 + s2*s2));
	
	/*calculate mobilities*/
	double tauR = RetardationTime(Temp, Te);
	double etaS = ShearViscosity(Temp,Te, smag);
		
	double itauR = 1.0/tauR;
	double ietaS = 1.0/etaS;
					
	double gamdot = 0.5*smag*ietaS;

	/*calculate the residual*/
    fRes[0]=Te-Ten+dt*itauR*(Te-Temp)-fA*smag*smag*dt*ietaS/iJ/fdeltac;
	
	fRes[1] = ep_e0 + 0.5*dt*ietaS*s0 - ep_tr0;
	fRes[2] = ep_e1 + 0.5*dt*ietaS*s1 - ep_tr1;
	fRes[3] = ep_e2 + 0.5*dt*ietaS*s2 - ep_tr2;

	tol = sqrt(dArrayT::Dot(fRes, fRes));
	int iteration = 0;
	while (tol>ctol && iteration < maxiteration)
	{
		iteration ++;
		/*calculate stiffness matrix*/
		/*derivative of retardation time wrt to Te*/
		double dtauR_dTe = fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Te - fTg)*(fC2 + Te - fTg));
		double detaS_dTe = fC1*fC2*(fC2 - fTg)*log(10)/(Temp*(fC2 + Te - fTg)*(fC2 + Te - fTg));
		double x = fQS*smag/(Temp*fsy0);
		double detaS_dsmag = 0.0;
		if (smag/fsy0 > small)
		{
			double cothx = cosh(x)/sinh(x);
			detaS_dsmag = (1.0-x*cothx)/smag;
		}

		fPot[1]->DevMod(fEigs_dev,eigenmodulus);
		/*deviatoric values*/
		double c0 = iJ*eigenmodulus(0,0);
		double c1 = iJ*eigenmodulus(1,1);
		double c2 = iJ*eigenmodulus(2,2);

		double c12 = iJ*eigenmodulus(1,2);
		double c02 = iJ*eigenmodulus(0,2);
		double c01 = iJ*eigenmodulus(0,1);
		
		/*initialize*/
		fiKAB = 0.0;
		
		/*K_del_del*/
		//fiKAB(0,0)=(1.0+dt*itauR)-dt*itauR*dtauR_dTe*(Te-Temp)+fA*smag*smag*dt*ietaS/iJ*detaS_dTe/fdeltac+fA*smag*smag*dt*ietaS/iJ/fdeltac*dalpha;
        fiKAB(0,0)=(1.0+dt*itauR)-dt*itauR*dtauR_dTe*(Te-Temp)+fA*smag*smag*dt*ietaS/iJ*detaS_dTe/fdeltac;

		/*K_epA_del*/
		fiKAB(1,0) = -0.5*dt*ietaS*s0*detaS_dTe;
		fiKAB(2,0) = -0.5*dt*ietaS*s1*detaS_dTe;
		fiKAB(3,0) = -0.5*dt*ietaS*s2*detaS_dTe;

		/*K_epA_epB*/
		double coef0 = 0.5*(s0*c0 + s1*c01 + s2*c02);
		double coef1 = 0.5*(s0*c01 + s1*c1 + s2*c12);
		double coef2 = 0.5*(s0*c02 + s1*c12 + s2*c2);
		if (smag/fsy0 > small)
		{
			coef0 /= smag;
			coef1 /= smag;
			coef2 /= smag;
		}
        fiKAB(1,1) = 1.0 + 0.5*ietaS*dt*c0 - 0.5*dt*ietaS*s0*detaS_dsmag*coef0;
        fiKAB(2,2) = 1.0 + 0.5*ietaS*dt*c1 - 0.5*dt*ietaS*s1*detaS_dsmag*coef1;
        fiKAB(3,3) = 1.0 + 0.5*ietaS*dt*c2 - 0.5*dt*ietaS*s2*detaS_dsmag*coef2;
        
        fiKAB(2,3) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s1*detaS_dsmag*coef2;
        fiKAB(1,3) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s0*detaS_dsmag*coef2;
        fiKAB(1,2) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s0*detaS_dsmag*coef1;
        
        fiKAB(3,2) = 0.5*ietaS*dt*c12 - 0.5*dt*ietaS*s2*detaS_dsmag*coef1;
        fiKAB(3,1) = 0.5*ietaS*dt*c02 - 0.5*dt*ietaS*s2*detaS_dsmag*coef0;
        fiKAB(2,1) = 0.5*ietaS*dt*c01 - 0.5*dt*ietaS*s1*detaS_dsmag*coef0;
        /*check the below derivitions*/
        fiKAB(0,1)=-2.0*fA/iJ*dt*ietaS*coef0*smag/fdeltac+fA/iJ*dt*ietaS*smag*smag*detaS_dsmag*coef0/fdeltac;
        fiKAB(0,2)=-2.0*fA/iJ*dt*ietaS*coef1*smag/fdeltac+fA/iJ*dt*ietaS*smag*smag*detaS_dsmag*coef1/fdeltac;
        fiKAB(0,3)=-2.0*fA/iJ*dt*ietaS*coef2*smag/fdeltac+fA/iJ*dt*ietaS*smag*smag*detaS_dsmag*coef2/fdeltac;
		/*inverts KAB*/
		fiKAB.Inverse();
	    
		
	    /*solve for the principal strain increments*/
		fiKAB.Multx(fRes, fDelta, -1.0);
		
	    /*updates principal elastic stretches*/ 
		Te += fDelta[0];
		
	    ep_e0 += fDelta[1];
	    ep_e1 += fDelta[2];
	    ep_e2 += fDelta[3];
	    
		
	    le0 = exp(2.0*ep_e0);
	    le1 = exp(2.0*ep_e1);
	    le2 = exp(2.0*ep_e2);
	    
		Je=sqrt(le0*le1*le2);
	    fEigs_dev = eigenstretch_e;
	    fEigs_dev *= pow(Je,-2.0*third);

	    /*calculate stresses and moduli*/
	    fPot[1]->DevStress(fEigs_dev, eigenstress);
	    
	    s0 = iJ*eigenstress[0];
	    s1 = iJ*eigenstress[1];
	    s2 = iJ*eigenstress[2];
	    
	    		
	    /*caculate smag*/
	    smag = sqrt(0.5*(s0*s0 + s1*s1 + s2*s2));

		/*calculate mobilities*/
		tauR = RetardationTime(Temp, Te);
		etaS = ShearViscosity(Temp, Te, smag);

		
		itauR = 1.0/tauR;
		ietaS = 1.0/etaS;
						
		gamdot =0.5*smag*ietaS;

	    /*calculate the residual*/
		fRes[0] = Te-Ten+dt*itauR*(Te-Temp)-fA*smag*smag*dt*ietaS/iJ/fdeltac;
		
	    fRes[1] = ep_e0 + 0.5*dt*ietaS*s0 - ep_tr0;
	    fRes[2] = ep_e1 + 0.5*dt*ietaS*s1 - ep_tr1;
	    fRes[3] = ep_e2 + 0.5*dt*ietaS*s2 - ep_tr2;
		
		
	    /*Check that the L2 norm of the residual is less than tolerance*/
	    tol = sqrt(dArrayT::Dot(fRes, fRes));
	}
 	if (iteration >= maxiteration)
	{
		ExceptionT::GeneralFail("SMP_coupled::ComputeEigs_e", 
			"number of iteration exceeds maximum");
	}
}
