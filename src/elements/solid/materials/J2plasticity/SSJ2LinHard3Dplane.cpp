/* $Id: SSJ2LinHard3Dplane.cpp,v 1.4 2004-07-15 08:28:34 paklein Exp $ */
/* created: paklein (02/12/1997)                                          */
/* Plane Strain linearly                */
/* isotropically elasto plastic material model subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */

#include "SSJ2LinHard3Dplane.h"

#include <iostream.h>
#include <math.h>

#include "SSMatSupportT.h"
#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"

/* class constants */

using namespace Tahoe;

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {"alpha", "VM", "press"};


/* constructor */
SSJ2LinHard3Dplane::SSJ2LinHard3Dplane(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("SSJ2LinHard3Dplane"),
	SSJ2LinHardBaseT(in, support),
        fStress(4),
	fModulus(dSymMatrixT::NumValues(4)){}


/*computes free energy density*/
double SSJ2LinHard3Dplane::StrainEnergyDensity(void)
{
    double energy = 0.0;
    fDevStrain = e();

    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    
    double I1 = fDevStrain.Trace();
    fDevStrain.PlusIdentity(-fthird*I1);
    
    fDevStrain -= fPlasticStrain;
    
    fStress = fDevStrain;
    fStress *= 2.0*fMu;
    fStress.PlusIdentity(fKappa*I1);

    fDevStrain.PlusIdentity(fthird*I1);

    energy = 0.5*fStress.ScalarProduct(fDevStrain);
    energy += 0.5*falpha[0]*dK()*falpha[0];
    energy += 0.5*fBeta.ScalarProduct(fPlasticStrain);
    return energy;		
}

const dMatrixT& SSJ2LinHard3Dplane::c_ijkl(void)
{
        s_ij();
     
	/*elastic moduli*/
	fModulus= 0.0;
	fModulus(0,0) = fModulus(1,1) = fModulus(3,3) = fKappa+4.0*fMu*fthird;
	fModulus(0,1) = fModulus(0,3) = fKappa-2.0*fMu*fthird;
	fModulus(1,0) = fModulus(1,3) = fKappa-2.0*fMu*fthird;
	fModulus(3,0) = fModulus(3,1) = fKappa-2.0*fMu*fthird;
	fModulus(2,2) = fMu;

       	int iteration = fSSMatSupport->IterationNumber();
       	if (iteration > -1) /* elastic iteration */
        fModulus += ModuliCorrection();

	/*     	const dMatrixT modcorr = ModuliCorrection();
	
       	fModulus(0,0) += modcorr(0,0); fModulus(1,1) += modcorr(1,1);
	fModulus(3,3) += modcorr(2,2);
	fModulus(0,1) += modcorr(0,1); fModulus(0,3) += modcorr(0,2);
	fModulus(1,0) += modcorr(1,0); fModulus(1,3) += modcorr(1,2);
	fModulus(3,0) += modcorr(2,0); fModulus(3,1) += modcorr(2,1);
	fModulus(2,2) += modcorr(5,5);*/

	return fModulus;
}

double SSJ2LinHard3Dplane::Pressure(void) const
{
    /*make sure to call s_ij() before calling Pressure;*/
    return (fthird*fStress.Trace());
}

/* stress */
const dSymMatrixT& SSJ2LinHard3Dplane::s_ij(void)
{
    fDevStrain = e();
    double I1 = fDevStrain.Trace();
    fDevStrain.PlusIdentity(-fthird*I1);
        
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());	

    fDevStrain -= fPlasticStrain_n;
    
    /* deviatoric part of trial stress */
    fStress = fDevStrain;
    fStress *= 2.0*fMu;

    /* modify trial stress (return mapping) */
    int iteration = fSSMatSupport->IterationNumber();
    if (iteration > -1) /* elastic iteration */
      fStress += StressCorrection(fStress);
    Store(element,CurrIP());
	
    /*add volumetric component*/
    fStress.PlusIdentity(fKappa*I1);
    
    return fStress;	
}

/*Note to be called only during post processing*/
const dArrayT& SSJ2LinHard3Dplane::InternalStrainVars(void)
{
        /*non-equilibrium components*/
        ElementCardT& element = CurrentElement();
        Load(element, CurrIP());

	/*evaluate viscous strains*/
	double* pvar = fInternalStrainVars.Pointer();

	*pvar++ = falpha[0];

      	*pvar++ = fPlasticStrain[0];
	*pvar++ = fPlasticStrain[1];
	*pvar++ = fPlasticStrain[2];
	*pvar++ = fPlasticStrain[3];
        
      	*pvar++ = fPlasticStrain[0];
	*pvar++ = fPlasticStrain[1];
	*pvar++ = fPlasticStrain[2];
	*pvar++ = fPlasticStrain[3];

        return(fInternalStrainVars);
}

const dArrayT& SSJ2LinHard3Dplane::InternalStressVars(void)
{
        /*non-equilibrium components*/
        ElementCardT& element = CurrentElement();
        Load(element, CurrIP());
        
	fDevStrain = e();
	double I1 = fDevStrain.Trace();
	fDevStrain.PlusIdentity(-fthird*I1);
	fDevStrain -= fPlasticStrain;

	fStress = fDevStrain;
	fStress *= 2.0*fMu;
	fStress.PlusIdentity(fKappa*I1);
	
	double* pvar = fInternalStressVars.Pointer();

	*pvar++ = -dK()*falpha[0];

	*pvar++ = -fBeta[0];
	*pvar++ = -fBeta[1];
	*pvar++ = -fBeta[2];
	*pvar++ = -fBeta[3];

	*pvar++ = fStress[0];
	*pvar++ = fStress[1];
	*pvar++ = fStress[2];
	*pvar++ = fStress[3];

        return(fInternalStressVars);
}

int SSJ2LinHard3Dplane::NumOutputVariables(void) const  { return kNumOutput; }
void SSJ2LinHard3Dplane::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void SSJ2LinHard3Dplane::ComputeOutput(dArrayT& output)
{

	/* stress tensor (loads element data and sets fStress) */
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());	

        fDevStrain = e();
	double I1 = fDevStrain.Trace();
	output[2] = fKappa*I1;

	fDevStrain.PlusIdentity(-fthird*I1);
 	fDevStrain -= fPlasticStrain;

	fStress = fDevStrain;
	fStress *= 2.0*fMu;

	/* deviatoric Von Mises stress */
	double J2 = 0.5*fStress.ScalarProduct();
	output[1] = sqrt(3.0*J2);
        output[0] = falpha[0];
}


