/* $Id: SMP_multi.cpp,v 1.1 2010-12-06 20:17:42 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "SMP_multi.h"

#include "PotentialT.h"
#include "NeoHookean.h"
#include "ArrudaBoyce.h"

#include "ifstreamT.h"
#include "ExceptionT.h"
#include <math.h>
#include "ofstreamT.h"
#include <stdlib.h>
#include "ParameterContainerT.h"

using namespace Tahoe;

const double loge = log10(exp(1.0));
const double pi = 2.0*acos(0);
const double third = 1.0/3.0; 
const double small = 1.0e-10;
const double big = 1.0e-10;

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
SMP_multi::SMP_multi(void):
  ParameterInterfaceT("SMP_multi")
{
}

double SMP_multi::FictiveTemperature(const double deltaneq)
{
	return (deltaneq/(falphar-falphag) + fT0);
}

double SMP_multi::StructuralRelaxationFunc(const double Temperature, const double deltaneq)
{
	/*calculate fictive temperature*/
	double Tf = FictiveTemperature(deltaneq);

	/*Hodge's Model*/
	double coeff = fC1/log10(exp(1.0))*(fC2*(Temperature - Tf) + Temperature*(Tf-fTg))/(Temperature*(fC2 + Tf - fTg));
	double itauR = exp(coeff);

	/*check the limits
	if (tauR > big)
		tauR = big;
	else if (tauR < small)
		tauR = small;
	*/	
	return(itauR);
}

double SMP_multi::StressRelaxationFunc(const double Temperature, const double deltaneq, const double smag, const double sy)
{
	/*calculate the temperature part*/
	double Tf = FictiveTemperature(deltaneq);
	double coeff = fC1/log10(exp(1.0))*(fC2*(Temperature - Tf) + Temperature*(Tf-fTg))/(Temperature*(fC2 + Tf - fTg));
	double itauS = exp(coeff);

	if (smag > kSmall)
	{
		itauS *= sinh(fQS/Temperature * smag/sy);
		itauS *= (sy*Temperature)/(fQS*smag);
	}
	else itauS = 1;
	/*check the limits
	if (tauS > big)
		tauS = big;
	else if (tauS < small)
		tauS = small;
*/
	return(itauS);
}

double SMP_multi::YieldRelaxationFunc(const double Temperature, const double deltaneq, const double smag, const double sy)
{
	/*calculate the temperature part*/
	double Tf = FictiveTemperature(deltaneq);
	double coeff = fC1/log10(exp(1.0))*(fC2*(Temperature - Tf) + Temperature*(Tf-fTg))/(Temperature*(fC2 + Tf - fTg));
	double itauY = exp(coeff);
	
	itauY /= fQS/Temperature;
	itauY *= sinh(fQS/Temperature * smag/sy);
			
	return(itauY);
}


const dMatrixT& SMP_multi::ThermalDeformation_Inverse(void)
{
	/*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
	double Temp = Compute_Temperature();
	
	/*thermal volume deformation*/
	double delneq_total = fdelneq.Sum();
	double ThetaT = 1.0 + falphar*(Temp-fT0) + (-(falphar-falphag)*(Temp-fT0)+delneq_total);
	fF_T_inv.Identity(pow(ThetaT,-1.0*third));
	return(fF_T_inv);
}

double SMP_multi::StrainEnergyDensity(void)
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
		
		energy += fdmu[i]*(fPot[1]->Energy(fEigs_dev, Je));
	}
	return(energy);
}

/* stresses */
const dSymMatrixT& SMP_multi::s_ij(void)
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
	/*jacobian determinant*/
	double J = sqrt(fEigs.Product());

/*	if(CurrElementNumber()==0&&CurrIP()==0)
		cout <<setprecision(12)<< "\nfEigs: "<<fEigs;
*/
	fEigs_dev = fEigs;
	fEigs_dev *= pow(J,-2.0*third);
	
	fPot[0]->DevStress(fEigs_dev, ftau_EQ, T);	
	ftau_EQ += fPot[0]->MeanStress(J);
	fStress3D = fSpectralDecompSpat.EigsToRank2(ftau_EQ);

    /*load the viscoelastic principal stretches from state variable arraysand calculate NEQ part*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
	{	
		/*update thermal deformation*/
		Compute_delneq(fdelneq_n, fdelneq);		
		Compute_le(fC_vn, fC_v, *fsy_n, *fsy);		
		Store(element, CurrIP());
	}	
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
	if (NumSD() == 2)
    {
        fStress[0] = fStress3D[0];
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
const dMatrixT& SMP_multi::c_ijkl(void)
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

	/*calculate NEQ*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    
	Compute_Kneq(fDAB);
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
		fMat(0,0) = fDAB(a,0);
		fMat(0,1) = fDAB(a,1);
		fMat(0,2) = fDAB(a,2);
		fMat(1,0) = fDAB(a+1,0);
		fMat(1,1) = fDAB(a+1,1);
		fMat(1,2) = fDAB(a+1,2);
		fMat(2,0) = fDAB(a+2,0);
		fMat(2,1) = fDAB(a+2,1);
		fMat(2,2) = fDAB(a+2,2);
		
		fCalg(0,0) = fDtauDe_NEQ(0,0)*fMat(0,0) + fDtauDe_NEQ(0,1)*fMat(1,0) + fDtauDe_NEQ(0,2)*fMat(2,0);
		fCalg(0,1) = fDtauDe_NEQ(0,0)*fMat(0,1) + fDtauDe_NEQ(0,1)*fMat(1,1) + fDtauDe_NEQ(0,2)*fMat(2,1);
		fCalg(0,2) = fDtauDe_NEQ(0,0)*fMat(0,2) + fDtauDe_NEQ(0,1)*fMat(1,2) + fDtauDe_NEQ(0,2)*fMat(2,2);

		fCalg(1,0) = fDtauDe_NEQ(1,0)*fMat(0,0) + fDtauDe_NEQ(1,1)*fMat(1,0) + fDtauDe_NEQ(1,2)*fMat(2,0);
		fCalg(1,1) = fDtauDe_NEQ(1,0)*fMat(0,1) + fDtauDe_NEQ(1,1)*fMat(1,1) + fDtauDe_NEQ(1,2)*fMat(2,1);
		fCalg(1,2) = fDtauDe_NEQ(1,0)*fMat(0,2) + fDtauDe_NEQ(1,1)*fMat(1,2) + fDtauDe_NEQ(1,2)*fMat(2,2);
		
		fCalg(2,0) = fDtauDe_NEQ(2,0)*fMat(0,0) + fDtauDe_NEQ(2,1)*fMat(1,0) + fDtauDe_NEQ(2,2)*fMat(2,0);
		fCalg(2,1) = fDtauDe_NEQ(2,0)*fMat(0,1) + fDtauDe_NEQ(2,1)*fMat(1,1) + fDtauDe_NEQ(2,2)*fMat(2,1);
		fCalg(2,2) = fDtauDe_NEQ(2,0)*fMat(0,2) + fDtauDe_NEQ(2,1)*fMat(1,2) + fDtauDe_NEQ(2,2)*fMat(2,2);

		fCalg2 += fCalg;
		
		fCalg(0,0) -= 2.0*ftau_NEQ[0];
		fCalg(1,1) -= 2.0*ftau_NEQ[1];
		fCalg(2,2) -= 2.0*ftau_NEQ[2];
		
		fModulus3D += fSpectralDecompSpat.NonSymEigsToRank4(fCalg);
    
		double dl_tr;

		double& l0_tr = fEigs_tr[0];
		double& l1_tr = fEigs_tr[1];
		double& l2_tr = fEigs_tr[2];
	
	
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


/*************************************************************************
*	PUBLIC
**************************************************************************/
/* describe the parameters needed by the interface */
void SMP_multi::DefineParameters(ParameterListT& list) const
{
  /* inherited */
  RGViscoelasticityT::DefineParameters(list);

  /* common limit */
  LimitT positive(0.0, LimitT::LowerInclusive);
  LimitT zero(0.0, LimitT::Lower);

  ParameterT alphar(ParameterT::Double, "rubbery_CTE");
  ParameterT alphag(ParameterT::Double, "glassy_CTE");

  ParameterT inittemp(ParameterT::Double, "initial_temperature_T0");
  ParameterT glasstemp(ParameterT::Double, "glass_trans_temp_Tg");
  ParameterT C1(ParameterT::Double, "WLF_C1");
  ParameterT C2(ParameterT::Double, "WLF_C2");
  
  alphar.AddLimit(zero);
  alphag.AddLimit(zero);
  inittemp.AddLimit(positive);
  glasstemp.AddLimit(positive);
  C1.AddLimit(positive);
  C2.AddLimit(positive);

  list.AddParameter(alphar);
  list.AddParameter(alphag);
  list.AddParameter(inittemp);
  list.AddParameter(glasstemp);
  list.AddParameter(C1);
  list.AddParameter(C2);

}

/* information about subordinate parameter lists */
void SMP_multi::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	sub_list.AddSub("smp_multi_eq_potential", ParameterListT::Once);
	sub_list.AddSub("smp_multi_neq_potential", ParameterListT::Once);

	/* choice of viscosity */
	sub_list.AddSub("smp_multi_structural_spectrum", ParameterListT::Once);
	sub_list.AddSub("smp_multi_stress_spectrum", ParameterListT::Once);
}


/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SMP_multi::NewSub(const StringT& name) const
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
	else if (name == "smp_multi_eq_potential")
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
	else if (name == "smp_multi_neq_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
	
		/* choice of parameters */
		choice->AddSub("neo-hookean");
		return(choice);
	}
	else if (name == "smp_multi_structural_spectrum")
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
	else if (name == "smp_multi_stress_spectrum")
	{
		ParameterContainerT* tauS = new ParameterContainerT(name);

		ParameterT file(ParameterT::String, "user_input_file");
		tauS->AddParameter(file);
		
		ParameterT A(ParameterT::Double, "activation_energy");
		ParameterT sy_0(ParameterT::Double, "init_yield_strength");
		ParameterT sy_ss(ParameterT::Double, "sat_yield_strength");
		ParameterT tauY(ParameterT::Double, "characteristic_relaxation_time_tauYref");
			
		A.AddLimit(zero);
		sy_0.AddLimit(zero);
		sy_ss.AddLimit(zero);
		tauY.AddLimit(positive);

		tauS->AddParameter(A);
		tauS->AddParameter(sy_0);
		tauS->AddParameter(sy_ss);
		tauS->AddParameter(tauY);

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

void SMP_multi::TakeParameterList(const ParameterListT& list)
{
  const char caller[] = "SMP_multi::TakeParameterList";
  /* inherited */
  RGViscoelasticityT::TakeParameterList(list);


  fT0 = list.GetParameter("initial_temperature_T0");
  fTg = list.GetParameter("glass_trans_temp_Tg");
  fC1 = list.GetParameter("WLF_C1");
  fC2 = list.GetParameter("WLF_C2");
  falphar = list.GetParameter("rubbery_CTE");
  falphag = list.GetParameter("glassy_CTE");

	fPot.Dimension(2);
		
	const ParameterListT& eq_pot = list.GetListChoice(*this, "smp_multi_eq_potential");
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
  
	const ParameterListT& neq_pot = list.GetListChoice(*this, "smp_multi_neq_potential");
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
  const ParameterListT* tauR = list.List("smp_multi_structural_spectrum");
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
  
  const ParameterListT* tauS = list.List("smp_multi_stress_spectrum");
  if (tauS)
  {
/*	ftauS0 = tauS->GetParameter("characterisic_relaxation_time_tauSref");
	fbetaS = tauS->GetParameter("stretch_exponential_exponent_beta");
	fNumS = tauS->GetParameter("number_discrete_relaxation_times");
	ftauS0min = tauS->GetParameter("min_relaxation_time_spectrum");
	ftauS0max = tauS->GetParameter("max_relaxation_time_spectrum");
*/
	fInputS = tauS->GetParameter("user_input_file");

	fQS = tauS->GetParameter("activation_energy");
	fsy0 = tauS->GetParameter("init_yield_strength");
	fsinf = tauS->GetParameter("sat_yield_strength");
	ftauY0 = tauS->GetParameter("characteristic_relaxation_time_tauYref");
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
	
	int nsd = NumSD();
	int ndof = 3;
	int numstress = dSymMatrixT::NumValues(ndof);

	fnstatev = 0;
	fnstatev += numstress*fNumS;   /*current C_v*/
	fnstatev += numstress*fNumS;   /*last C_vn*/
	fnstatev += fNumR;			/*current  delta_neq*/ 
	fnstatev += fNumR;			/*last delta_neq_n*/
	fnstatev ++;			/*current yield strength*/
	fnstatev ++;			/*last yield strength*/
	
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

	fdelneq.Alias(fNumR,pstatev); 
	pstatev += fNumR;
	fdelneq_n.Alias(fNumR,pstatev); 
	pstatev += fNumR;
	
	fsy = pstatev;
	pstatev++;
	fsy_n = pstatev;
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
	fsig.Dimension(ndof);
	
	fStress.Dimension(NumSD());
	fStress3D.Dimension(ndof);

	fDtauDe_EQ.Dimension(ndof);
	fDtauDe_NEQ.Dimension(ndof);
	
	fKdel.Dimension(fNumR);
	fRdel.Dimension(fNumR);
   
	fKAB.Dimension(fNumS*3+1);
	fKAB2.Dimension(fNumS*3+1);
	fRes.Dimension(fNumS*3+1);

   fGA0.Dimension(fNumS*3+1);
   fGA1.Dimension(fNumS*3+1);
   fGA2.Dimension(fNumS*3+1);
	
	fCalg.Dimension(3);
	fCalg2.Dimension(3);
	fCalg2=0.0;
	fMat.Dimension(3);
	fDAB.Dimension(fNumS*3,3);
	
  fModulus3D.Dimension(dSymMatrixT::NumValues(ndof));
  fModMat.Dimension(dSymMatrixT::NumValues(ndof));
  fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));

}

/*initializes history variable */
void  SMP_multi::PointInitialize(void)
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
		      
			  for (int k = 0; k < fNumS; k++)
			  {
				fC_v[k].Identity();
				fC_vn[k].Identity();
			  }
			  fdelneq = 0.0;
			  fdelneq_n = 0.0;
			  *fsy = fsy0;
			  *fsy_n = fsy0;

		      /* write to storage */
		      Store(element, ip);
		}
	}
}
 
void SMP_multi::UpdateHistory(void)
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
		*fsy_n = *fsy;
		fdelneq_n = fdelneq;
		
		/* write to storage */
		Store(element, ip);
	}
}

void SMP_multi::ResetHistory(void)
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
		*fsy = *fsy_n;
		fdelneq = fdelneq_n;
		
		/* write to storage */
		Store(element, ip);
	}
	
}

void SMP_multi::Compute_delneq(const dArrayT& delneq_n, dArrayT& delneq)
{
	double ctol = 1.0e-10;
	/*time step*/
	double dt = fFSMatSupport->TimeStep();
	/*current temperature*/
	double Tn = Compute_Temperature();
	
	double delT = delneq.Sum();
	double Tf = FictiveTemperature(delT);
	double itaubar = StructuralRelaxationFunc(Tn, delT);
	double temp = 0.0;
	for(int k=0; k<fNumR; k++)
	{
		double itauRk = itaubar/ftimesR[k];
		/*residual*/
		fRdel[k] = delneq[k] + dt*itauRk*(delneq[k]-fdalpha[k]*(Tn-fT0)) - delneq_n[k];
		temp += fRdel[k]*fRdel[k];
	}
	double tol = sqrt(temp);
	
	while (tol > ctol)
	{
		fKdel = 0.0;
		for(int k=0; k<fNumR; k++)
		{	
			/*stiffness matrix*/
			double itauRk = itaubar/ftimesR[k];
			double DtaubarDTf = (fC1*fC2*(fC2-fTg)*log(10.0))/(Tn*(fC2+Tf-fTg)*(fC2+Tf-fTg));
			fKdel(k,k) = 1.0 + dt*itauRk;
			for (int l=0; l<fNumR; l++)
				fKdel(k,l) -= dt*itauRk*(delneq[k] - fdalpha[k]*(Tn-fT0))*DtaubarDTf/(falphar-falphag);
		}
		/*Solve for update*/
		fKdel.LinearSolve(fRdel);
		delneq -= fRdel;
		
		/*calculate residual*/
		delT = delneq.Sum();
		itaubar = StructuralRelaxationFunc(Tn, delT);
		Tf = FictiveTemperature(delT);
		double temp = 0.0;
		for(int k=0; k<fNumR; k++)
		{
			double itauRk = itaubar/ftimesR[k];
			/*residual*/
			fRdel[k] = delneq[k] + dt*itauRk*(delneq[k]-fdalpha[k]*(Tn-fT0)) - delneq_n[k];
			temp += fRdel[k]*fRdel[k];
		}
		tol = sqrt(temp);
	}
	
}

void SMP_multi::Compute_Kneq(dMatrixT& Modulus)
{
	/*calculate  smag*/
	const dMatrixT& Ftotal = F_total();	
	double iJ = 1.0/Ftotal.Det();

	/*time step*/
	const double dt = fFSMatSupport->TimeStep();

	/*temperature and temperature step*/
	const double Tn = Compute_Temperature();

	/*neq thermal def*/
	double delT = fdelneq.Sum();
	
	/*yield*/
	double sy = *fsy;
	fsig = 0.0;
	double* ple = fle.Pointer();
	for (int i = 0; i < fNumS; i++)
	{
		fInverse.Inverse(fC_v[i]);
		fbe.MultQBQT(fF3D, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
		fEigs_e = fSpectralDecompSpat.Eigenvalues(); 
		*ple++ = fEigs_e[0];
		*ple++ = fEigs_e[1];
		*ple++ = fEigs_e[2];
		
		fEigs_dev = fEigs_e;
		double Je = sqrt(fEigs_dev.Product());
		fEigs_dev *= pow(Je,-2.0*third);
		
		/*calculate total neq stress*/
		fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
		ftau_NEQ *=fdmu[i];
		fsig += ftau_NEQ;
	}
	fsig *= iJ;
	double smag = sqrt(0.5*(fsig[0]*fsig[0]+fsig[1]*fsig[1]+fsig[2]*fsig[2]));

	/*update viscosity*/
	double ietabar = StressRelaxationFunc(Tn,	delT, smag, sy);
	double itauy = YieldRelaxationFunc(Tn,delT, smag, sy);

	fKAB = 0.0;
	fGA0 = 0.0;
	fGA1 = 0.0;
	fGA2 = 0.0;
	Modulus = 0.0;
	ple = fle.Pointer();
	double DetabarDs,DetabarDsy,DtauybarDs,DtauybarDsy;

	double itauy0 = 1.0/ftauY0;
	itauy *= itauy0;
	for(int k=0; k<fNumS; k++)
	{	
		/*viscosity and its derivatives*/
		double mu=fdmu[k];
		double ietaS0k = 1.0/(mu*ftimesS[k]);
		double ietaSk = ietaS0k*ietabar;
		
		if(smag > kSmall)
		{
			double x = (fQS*smag)/(sy*Tn);
			double cothx = cosh(x)/sinh(x);
			DetabarDs = (1.0 - x*cothx)/smag;
			DetabarDsy = -(1.0 - x*cothx)/sy;
			DtauybarDs = -x*cothx/smag;
			DtauybarDsy = x*cothx/sy;
		}
		else
		{
			DetabarDs = 0.0;
			DetabarDsy = 0.0;
			DtauybarDs = 0.0;
			DtauybarDsy = 0.0;
		}
		/*moduli*/
		fEigs_dev[0] = *ple++;
		fEigs_dev[1] = *ple++;
		fEigs_dev[2] = *ple++;
		double Je = sqrt(fEigs_dev[0]*fEigs_dev[1]*fEigs_dev[2]);
		fEigs_dev *= pow(Je,-2.0*third);
		
		fPot[1]->DevMod(fEigs_dev,fDtauDe_NEQ);
		
		fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
		double dk00 = mu*iJ*fDtauDe_NEQ(0,0);
		double dk11 = mu*iJ*fDtauDe_NEQ(1,1);
		double dk22 = mu*iJ*fDtauDe_NEQ(2,2);
		double dk12 = mu*iJ*fDtauDe_NEQ(1,2);
		double dk02 = mu*iJ*fDtauDe_NEQ(0,2);
		double dk01 = mu*iJ*fDtauDe_NEQ(0,1);
		
		double sk0 = mu*iJ*ftau_NEQ[0];
		double sk1 = mu*iJ*ftau_NEQ[1];
		double sk2 = mu*iJ*ftau_NEQ[2];
			
		/*stiffness matrix*/
		int a = k*3;
				
		fKAB(a,a) = 1.0 + 0.5*dt*ietaSk*dk00;
		fKAB(a,a+1) = 0.5*dt*ietaSk*dk01;
		fKAB(a,a+2) = 0.5*dt*ietaSk*dk02;
		fKAB(a+1,a) = 0.5*dt*ietaSk*dk01;
		fKAB(a+1,a+1) = 1.0 + 0.5*dt*ietaSk*dk11;
		fKAB(a+1,a+2) = 0.5*dt*ietaSk*dk12;
		fKAB(a+2,a) = 0.5*dt*ietaSk*dk02;
		fKAB(a+2,a+1) = 0.5*dt*ietaSk*dk12;
		fKAB(a+2,a+2) = 1.0 + 0.5*dt*ietaSk*dk22;

		double *ple2 = fle.Pointer();
		for (int l=0; l<fNumS; l++)
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
		
			double cl0 = (fsig[0]*dl00 + fsig[1]*dl01 + fsig[2]*dl02);
			double cl1 = (fsig[0]*dl01 + fsig[1]*dl11 + fsig[2]*dl12);
			double cl2 = (fsig[0]*dl02 + fsig[1]*dl12 + fsig[2]*dl22);
			
			if(smag >kSmall)
			{
				cl0 /= smag;
				cl1 /= smag;
				cl2 /= smag;
			}
			
			int b = l*3;
			fKAB(a,b) = fKAB(a,b) - 0.5*dt*ietaSk*0.5*DetabarDs*sk0*cl0;
			fKAB(a,b+1) = fKAB(a,b+1) - 0.5*dt*ietaSk*0.5*DetabarDs*sk0*cl1;
			fKAB(a,b+2) = fKAB(a,b+2) - 0.5*dt*ietaSk*0.5*DetabarDs*sk0*cl2;
			fKAB(a+1,b) = fKAB(a+1,b) - 0.5*dt*ietaSk*0.5*DetabarDs*sk1*cl0;
			fKAB(a+1,b+1) = fKAB(a+1,b+1) - 0.5*dt*ietaSk*0.5*DetabarDs*sk1*cl1;
			fKAB(a+1,b+2) = fKAB(a+1,b+2) - 0.5*dt*ietaSk*0.5*DetabarDs*sk1*cl2;
			fKAB(a+2,b) = fKAB(a+2,b) - 0.5*dt*ietaSk*0.5*DetabarDs*sk2*cl0;
			fKAB(a+2,b+1) = fKAB(a+2,b+1) - 0.5*dt*ietaSk*0.5*DetabarDs*sk2*cl1;
			fKAB(a+2,b+2) = fKAB(a+2,b+2) - 0.5*dt*ietaSk*0.5*DetabarDs*sk2*cl2;
			
			if(k==fNumS-1)
			{
				fKAB(3*fNumS,b) =  1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*DtauybarDs*sy*0.5*cl0;
				fKAB(3*fNumS,b+1) =  1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*DtauybarDs*sy*0.5*cl1;
				fKAB(3*fNumS,b+2) =  1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*DtauybarDs*sy*0.5*cl2;
			}
		}
		fKAB(a,3*fNumS) = -0.5*dt*ietaSk*DetabarDsy*sk0;
		fKAB(a+1,3*fNumS) = -0.5*dt*ietaSk*DetabarDsy*sk1;
		fKAB(a+2,3*fNumS) = -0.5*dt*ietaSk*DetabarDsy*sk2;

		fGA0[a] = 1.0 + 0.5*dt*ietaSk*(1.0-DetabarDs*smag)*sk0;
		fGA0[a+1] = 0.5*dt*ietaSk*(1.0-DetabarDs*smag)*sk1;
		fGA0[a+2] = 0.5*dt*ietaSk*(1.0-DetabarDs*smag)*sk2;

		fGA1[a] = 0.5*dt*ietaSk*(1.0-DetabarDs*smag)*sk0;
		fGA1[a+1] = 1.0 + 0.5*dt*ietaSk*(1.0-DetabarDs*smag)*sk1;
		fGA1[a+2] = 0.5*dt*ietaSk*(1.0-DetabarDs*smag)*sk2;

		fGA2[a] = 0.5*dt*ietaSk*(1.0-DetabarDs*smag)*sk0;
		fGA2[a+1] = 0.5*dt*ietaSk*(1.0-DetabarDs*smag)*sk1;
		fGA2[a+2] = 1.0 + 0.5*dt*ietaSk*(1.0-DetabarDs*smag)*sk2;

	}
	fKAB(3*fNumS,3*fNumS) = 1.0+1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*DtauybarDsy*sy - 1.0/sqrt(2.0)*dt*itauy*(1.0-2.0*sy/fsinf);

	fGA0[3*fNumS] = 1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*DtauybarDs*sy*smag;
	fGA1[3*fNumS] = 1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*DtauybarDs*sy*smag;
	fGA2[3*fNumS] = 1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*DtauybarDs*sy*smag;
	
	/*set copy because fKAB is overwritten by LinearSolve*/
	fKAB2=fKAB;
	
	/*KAB^-1 GBC*/
	fKAB.LinearSolve(fGA0);
	fKAB = fKAB2;
	fKAB.LinearSolve(fGA1);
	fKAB = fKAB2;
	fKAB.LinearSolve(fGA2);
	/*condense out Dsy term*/
	for (int k =0; k<3*fNumS; k++)
	{
		Modulus(k,0) = fGA0[k];
		Modulus(k,1) = fGA1[k];
		Modulus(k,2) = fGA2[k];
	}
}

void SMP_multi::Compute_le(const ArrayT<dSymMatrixT>& C_vn, ArrayT<dSymMatrixT>& C_v, const double& sy_n, double& sy)
{		
	double ctol = 1.00e-9;
	int maxiter = 20;
		
	/*time step*/
	const double dt = fFSMatSupport->TimeStep();

	/*temperature and temperature step*/
	const double Tn = Compute_Temperature();

	/*total jacobian*/
	const dMatrixT& Ftotal = F_total();	
	double iJ = 1.0/Ftotal.Det();
	
	double delT = fdelneq.Sum();

	/*calc trial solution l_tr and smag*/
	fsig = 0.0;
	double* pltr = fl_tr.Pointer();
	double* ple = fle.Pointer();
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
		fsig += ftau_NEQ;
	}
	fsig*=iJ;
	double smag = sqrt(0.5*(fsig[0]*fsig[0]+fsig[1]*fsig[1]+fsig[2]*fsig[2]));
	
	/*calc viscosity functions*/
	double ietabar = StressRelaxationFunc(Tn,	delT, smag, sy);
	double itauy = YieldRelaxationFunc(Tn,delT, smag, sy);
	double itauy0 = 1.0/ftauY0;
	itauy *= itauy0;
	/*calc residual*/
	/*re-assign pointers*/
	pltr = fl_tr.Pointer();
	ple = fle.Pointer();
	double* pr = fRes.Pointer();
	fRes = 0.0;
	double epse0, epse1, epse2, epstr0, epstr1, epstr2;
	double r0, r1, r2, ry, tol;
	double temp = 0.0;
	for (int k = 0; k < fNumS; k++)
	{
		epstr0 = 0.5*log(*pltr++);
		epstr1 = 0.5*log(*pltr++);
		epstr2 = 0.5*log(*pltr++);
		
		/*calculate neq stress tauk*/
		fEigs_dev[0] = *ple++;
		fEigs_dev[1] = *ple++;
		fEigs_dev[2] = *ple++;
		epse0 = 0.5*log(fEigs_dev[0]);
		epse1 = 0.5*log(fEigs_dev[1]);
		epse2 = 0.5*log(fEigs_dev[2]);
		
		double Je = sqrt(fEigs_dev.Product());
		fEigs_dev *= pow(Je,-2.0*third);
		fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
		double s0 = fdmu[k]*iJ*ftau_NEQ[0];
		double s1 = fdmu[k]*iJ*ftau_NEQ[1];
		double s2 = fdmu[k]*iJ*ftau_NEQ[2];
		
		/*viscosity*/
		double ietaSk = 1.0/(fdmu[k]*ftimesS[k]);
		ietaSk *= ietabar;

		/*calc residual*/
		r0 = epse0 + 0.5*dt*ietaSk*s0 - epstr0;
		r1 = epse1 + 0.5*dt*ietaSk*s1 - epstr1;
		r2 = epse2 + 0.5*dt*ietaSk*s2 - epstr2;
		temp += r0*r0 + r1*r1+ r2*r2;
		*pr++ = r0;
		*pr++ = r1;
		*pr++ = r2;
	}
	ry = sy - 1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*sy - sy_n;
	*pr++ = ry;
	temp += ry*ry;
	tol = sqrt(temp);
	double tol0 = tol;
	double reltol =tol0;
	int iter = 0;
	while (tol > 10*ctol && reltol>ctol && iter < maxiter)
	{
		iter++;
		fKAB = 0.0;
		ple = fle.Pointer();
		double DetabarDs,DetabarDsy,DtauybarDs,DtauybarDsy;
		double x = (fQS*smag)/(sy*Tn);
		double cothx = cosh(x)/sinh(x);
		DetabarDs = (1.0 - x*cothx)/smag;
		DetabarDsy = -(1.0 - x*cothx)/sy;
		DtauybarDs = -x*cothx/smag;
		DtauybarDsy = x*cothx/sy;
		for(int k=0; k<fNumS; k++)
		{	
			/*viscosity and its derivatives*/
			double mu=fdmu[k];
			double ietaS0k = 1.0/(mu*ftimesS[k]);
			double ietaSk = ietaS0k*ietabar;
			
			/*moduli*/
			fEigs_dev[0] = *ple++;
			fEigs_dev[1] = *ple++;
			fEigs_dev[2] = *ple++;
			double Je = sqrt(fEigs_dev[0]*fEigs_dev[1]*fEigs_dev[2]);
			fEigs_dev *= pow(Je,-2.0*third);
			
			fPot[1]->DevMod(fEigs_dev,fDtauDe_NEQ);
			fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
			double dk00 = mu*iJ*fDtauDe_NEQ(0,0);
			double dk11 = mu*iJ*fDtauDe_NEQ(1,1);
			double dk22 = mu*iJ*fDtauDe_NEQ(2,2);
			double dk12 = mu*iJ*fDtauDe_NEQ(1,2);
			double dk02 = mu*iJ*fDtauDe_NEQ(0,2);
			double dk01 = mu*iJ*fDtauDe_NEQ(0,1);
			
			double sk0 = mu*iJ*ftau_NEQ[0];
			double sk1 = mu*iJ*ftau_NEQ[1];
			double sk2 = mu*iJ*ftau_NEQ[2];
			
			/*stiffness matrix*/
			int a = k*3;
			fKAB(a,a) = 1.0 + 0.5*dt*ietaSk*dk00;
			fKAB(a,a+1) = 0.5*dt*ietaSk*dk01;
			fKAB(a,a+2) = 0.5*dt*ietaSk*dk02;
			fKAB(a+1,a) = 0.5*dt*ietaSk*dk01;
			fKAB(a+1,a+1) = 1.0 + 0.5*dt*ietaSk*dk11;
			fKAB(a+1,a+2) = 0.5*dt*ietaSk*dk12;
			fKAB(a+2,a) = 0.5*dt*ietaSk*dk02;
			fKAB(a+2,a+1) = 0.5*dt*ietaSk*dk12;
			fKAB(a+2,a+2) = 1.0 + 0.5*dt*ietaSk*dk22;
			double *ple2 = fle.Pointer();
			for (int l=0; l<fNumS; l++)
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
			
				double cl0 = (fsig[0]*dl00 + fsig[1]*dl01 + fsig[2]*dl02)/smag;
				double cl1 = (fsig[0]*dl01 + fsig[1]*dl11 + fsig[2]*dl12)/smag;
				double cl2 = (fsig[0]*dl02 + fsig[1]*dl12 + fsig[2]*dl22)/smag;

				int b = l*3;
				fKAB(a,b) = fKAB(a,b) - 0.5*dt*ietaSk*0.5*DetabarDs*sk0*cl0;
				fKAB(a,b+1) = fKAB(a,b+1) - 0.5*dt*ietaSk*0.5*DetabarDs*sk0*cl1;
				fKAB(a,b+2) = fKAB(a,b+2) - 0.5*dt*ietaSk*0.5*DetabarDs*sk0*cl2;
				fKAB(a+1,b) = fKAB(a+1,b) - 0.5*dt*ietaSk*0.5*DetabarDs*sk1*cl0;
				fKAB(a+1,b+1) = fKAB(a+1,b+1) - 0.5*dt*ietaSk*0.5*DetabarDs*sk1*cl1;
				fKAB(a+1,b+2) = fKAB(a+1,b+2) - 0.5*dt*ietaSk*0.5*DetabarDs*sk1*cl2;
				fKAB(a+2,b) = fKAB(a+2,b) - 0.5*dt*ietaSk*0.5*DetabarDs*sk2*cl0;
				fKAB(a+2,b+1) = fKAB(a+2,b+1) - 0.5*dt*ietaSk*0.5*DetabarDs*sk2*cl1;
				fKAB(a+2,b+2) = fKAB(a+2,b+2) - 0.5*dt*ietaSk*0.5*DetabarDs*sk2*cl2;

				if(k==fNumS-1)
				{
					fKAB(3*fNumS,b) =  1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*DtauybarDs*sy*0.5*cl0;
					fKAB(3*fNumS,b+1) =  1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*DtauybarDs*sy*0.5*cl1;
					fKAB(3*fNumS,b+2) =  1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*DtauybarDs*sy*0.5*cl2;
				}
			}
			fKAB(a,3*fNumS) = -0.5*dt*ietaSk*DetabarDsy*sk0;
			fKAB(a+1,3*fNumS) = -0.5*dt*ietaSk*DetabarDsy*sk1;
			fKAB(a+2,3*fNumS) = -0.5*dt*ietaSk*DetabarDsy*sk2;
		}
		fKAB(3*fNumS,3*fNumS) = 1.0+1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*DtauybarDsy*sy - 1.0/sqrt(2.0)*dt*itauy*(1.0-2.0*sy/fsinf);
		
		/*Solve  and update*/
		fKAB.LinearSolve(fRes);	


		ple = fle.Pointer();
		pr = fRes.Pointer();
		for (int k = 0; k < 3*fNumS; k++)
		{
			double dep = -(*pr++);
			*ple++ *=exp(2.0*dep);
		}
		sy -= *pr++;
		
		/*update residual*/
		/*update  smag*/
		fsig = 0.0;
		double* ple = fle.Pointer();
		for (int i = 0; i < fNumS; i++)
		{
			fEigs_dev[0] = *ple++;
			fEigs_dev[1] = *ple++;
			fEigs_dev[2] = *ple++;
			double Je = sqrt(fEigs_dev.Product());
			fEigs_dev *= pow(Je,-2.0*third);
		
			/*calculate total neq stress*/
			fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
			ftau_NEQ *=fdmu[i];
			fsig += ftau_NEQ;
		}
		fsig *= iJ;
		smag = sqrt(0.5*(fsig[0]*fsig[0]+fsig[1]*fsig[1]+fsig[2]*fsig[2]));
		/*update viscosity*/
		ietabar = StressRelaxationFunc(Tn,	delT, smag, sy);
		itauy = YieldRelaxationFunc(Tn,delT, smag, sy);
		itauy /= ftauY0;

		pltr = fl_tr.Pointer();
		ple = fle.Pointer();
		pr = fRes.Pointer();
		fRes = 0.0;
		temp = 0.0;				
		for (int k = 0; k < fNumS; k++)
		{
			epstr0 = 0.5*log(*pltr++);
			epstr1 = 0.5*log(*pltr++);
			epstr2 = 0.5*log(*pltr++);
			
			/*calculate neq stress tauk*/
			fEigs_dev[0] = *ple++;
			fEigs_dev[1] = *ple++;
			fEigs_dev[2] = *ple++;
			epse0 = 0.5*log(fEigs_dev[0]);
			epse1 = 0.5*log(fEigs_dev[1]);
			epse2 = 0.5*log(fEigs_dev[2]);
			
			double Je = sqrt(fEigs_dev.Product());
			fEigs_dev *= pow(Je,-2.0*third);
			fPot[1]->DevStress(fEigs_dev, ftau_NEQ);
			double mu = fdmu[k];
			double s0 = mu*iJ*ftau_NEQ[0];
			double s1 = mu*iJ*ftau_NEQ[1];
			double s2 = mu*iJ*ftau_NEQ[2];
			
			/*viscosity*/
			double ietaSk0 = 1.0/(mu*ftimesS[k]);
			double ietaSk = ietaSk0*ietabar;
			
			/*calc residual*/
			r0 = epse0 + 0.5*dt*ietaSk*s0 - epstr0;
			r1 = epse1 + 0.5*dt*ietaSk*s1 - epstr1;
			r2 = epse2 + 0.5*dt*ietaSk*s2 - epstr2;
			temp += r0*r0 + r1*r1+ r2*r2;
			*pr++ = r0;
			*pr++ = r1;
			*pr++ = r2;
		}
		 ry = sy - 1.0/sqrt(2.0)*dt*itauy*(1.0-sy/fsinf)*sy - sy_n;
		*pr++ = ry;
		temp += ry*ry;
		 tol = sqrt(temp);
		 reltol = tol/tol0;
			
	} /*while loop*/
	if (iter >= maxiter) 
	{
		cout<<"\n Number of iteration exceeds maximum. tol0:  "<<tol0<<"\ttol: "<<tol<<"\treltol: "<<reltol;
		cout<< "\nelem: "<<CurrElementNumber(); 
		ExceptionT::GeneralFail("SMP_multi::Compute_le", "number of iteration exceeds maximum");
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
}

/*calculates the series approximation of the KWW distribution function - rho(tau/tauKWW)*/
/*Lindsey et al. J. Chem. Phys. 73:3348, 1980*/
double SMP_multi::StretchedExponentialSpectrum(const double x, const double tauKWW, const double betakWW)
{
	/*x = tau/taukWW*/
	int nk = 201;
	double rho = 0;
	
	for (int k = 0; k < nk; k++)
	{
		double kfac= fGamma.Function(k+1);
		double kpow = pow(-1,k);
		double kb = betakWW*k;
		rho += kpow/kfac * sin(pi*kb)*fGamma.Function(kb+1.0)*pow(x,kb+1.0);
	}
	rho *= -1.0/(pi*tauKWW*x*x);
	return(rho);
	
}

void SMP_multi::Compute_discrete_spectrum(const double beta, const double tauWW, const double taumin, const double taumax, const int num, dArrayT& times, dArrayT& spectrum, StringT filename)
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




