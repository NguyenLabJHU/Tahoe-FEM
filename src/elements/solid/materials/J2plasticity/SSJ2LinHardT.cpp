/* $Id: SSJ2LinHardT.cpp,v 1.5 2004-07-15 08:28:34 paklein Exp $ */
/* created: paklein (02/12/1997)                                          */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */

#include "SSJ2LinHardT.h"

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
SSJ2LinHardT::SSJ2LinHardT(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("SSJ2LinHardT"),
	SSJ2LinHardBaseT(in, support),
	fStress(3),
	fModulus(dSymMatrixT::NumValues(3)){}

/*computes free energy density*/
double SSJ2LinHardT::StrainEnergyDensity(void)
{
    double energy = 0.0;
    const dSymMatrixT& strain = e();

    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());	    
    
    fDevStrain = strain;
    double I1 = fDevStrain.Trace();
    fDevStrain[0] -= fthird*I1;
    fDevStrain[1] -= fthird*I1;
    fDevStrain[2] -= fthird*I1;
    fDevStrain -= fPlasticStrain;

    fStress = fDevStrain;
    fStress *= 2.0*fMu;
    fStress.PlusIdentity(fKappa*I1);
  
    energy = 0.5*fStress.ScalarProduct(strain);
    energy += 0.5*falpha[0]*dK()*falpha[0];
    energy += 0.5*fBeta.ScalarProduct(fPlasticStrain);    
    return energy;		
}

double SSJ2LinHardT::Pressure(void) const
{
  /*call s_ij before calling this function*/
  return (fthird*fStress.Trace());
}
const dMatrixT& SSJ2LinHardT::c_ijkl(void)
{
  s_ij();
  
  /*deviatoric part*/
  fModulus = 0.0;
  fModulus(0,0) = fModulus(1,1) =  fModulus(2,2) = 2.0*fMu*(1.0 - fthird);
  fModulus(3,3) = fModulus(4,4) =  fModulus(5,5) = fMu;
  fModulus(0,1) = fModulus(0,2) =  fModulus(1,2) = -2.0*fMu*fthird;
  fModulus(1,0) = fModulus(2,0) =  fModulus(2,1) = -2.0*fMu*fthird;
  
  /*volumetric part*/
  fModulus(0,0) += fKappa; fModulus(1,1) += fKappa; fModulus(2,2) += fKappa;
  fModulus(0,1) += fKappa; fModulus(0,2) += fKappa; fModulus(1,2) += fKappa;
  fModulus(1,0) += fKappa; fModulus(2,0) += fKappa; fModulus(2,1) += fKappa;
  
  int iteration = fSSMatSupport->IterationNumber();
  if (iteration > -1) /* elastic iteration */
    fModulus += ModuliCorrection();
  return fModulus;
}

/* stress */
const dSymMatrixT& SSJ2LinHardT::s_ij(void)
{
  ElementCardT& element = CurrentElement();
  Load(element, CurrIP());	

  fDevStrain = e();
  double I1 = fDevStrain.Trace();
  fDevStrain.PlusIdentity(-fthird*I1);  
  fDevStrain -= fPlasticStrain_n;
  
  /*evaluate deviatoric part of trial elastic stress*/
  fStress = fDevStrain;
  fStress *= 2.0*fMu;
  
  /* modify deviatoric part of trial elastic stress (return mapping) */
  int iteration = fSSMatSupport->IterationNumber();
  if (iteration > -1) /* elastic iteration */
    fStress += StressCorrection(fStress);
  Store(element,CurrIP());
  
  /*add volumetric part*/
  fStress.PlusIdentity(fKappa*I1);
  return fStress;	
}

/*Note to be called only during post processing*/
const dArrayT& SSJ2LinHardT::InternalStrainVars(void)
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
	*pvar++ = fPlasticStrain[4];
	*pvar++ = fPlasticStrain[5];
        
      	*pvar++ = fPlasticStrain[0];
	*pvar++ = fPlasticStrain[1];
	*pvar++ = fPlasticStrain[2];
	*pvar++ = fPlasticStrain[3];
	*pvar++ = fPlasticStrain[4];
	*pvar = fPlasticStrain[5];
        
        return(fInternalStrainVars);
}

const dArrayT& SSJ2LinHardT::InternalStressVars(void)
{
        /*non-equilibrium components*/
        ElementCardT& element = CurrentElement();
        Load(element, CurrIP());
        
	fDevStrain = e();
	double I1 = fDevStrain.Trace();
	fDevStrain[0] -= fthird*I1;
	fDevStrain[1] -= fthird*I1;
	fDevStrain[2] -= fthird*I1;
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
	*pvar++ = -fBeta[4];
	*pvar++ = -fBeta[5];

	*pvar++ = fStress[0];
	*pvar++ = fStress[1];
	*pvar++ = fStress[2];
	*pvar++ = fStress[3];
	*pvar++ = fStress[4];
	*pvar = fStress[5];
       
        return(fInternalStressVars);
}

int SSJ2LinHardT::NumOutputVariables(void) const  { return kNumOutput; }
void SSJ2LinHardT::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void SSJ2LinHardT::ComputeOutput(dArrayT& output)
{

        ElementCardT& element = CurrentElement();
	Load(element, CurrIP());	    

	fDevStrain = e();
	double I1 = fDevStrain.Trace();
	output[2] = fKappa*I1;

	fDevStrain[0] -= fthird*I1;
	fDevStrain[1] -= fthird*I1;
	fDevStrain[2] -= fthird*I1;
	fDevStrain -= fPlasticStrain;

	fStress = fDevStrain;
	fStress *= 2.0*fMu;
	
	/* deviatoric Von Mises stress */
	double J2 = fStress.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[1] = sqrt(3.0*J2);

 	output[0] = falpha[0];
}
