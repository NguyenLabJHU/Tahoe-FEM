/* $Id: SSJ2LinHard2D.cpp,v 1.10 2004-07-15 08:28:34 paklein Exp $ */
/* created: paklein (02/12/1997)                                          */
/* Plane Strain linearly                */
/* isotropically elasto plastic material model subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */

#include "SSJ2LinHard2D.h"

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
SSJ2LinHard2D::SSJ2LinHard2D(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("SSJ2LinHard2D"),
	SSJ2LinHardBase2D(in, support),
	fStress(2),
	fModulus(dSymMatrixT::NumValues(2)),
	fStress3D(3),
	fModulus3D(dSymMatrixT::NumValues(3)){}


/*computes free energy density*/
double SSJ2LinHard2D::StrainEnergyDensity(void)
{
    double energy = 0.0;
    const dSymMatrixT& strain = e();

    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    
    double I1 = strain.Trace();
    fDevStrain = 0.0;
    fDevStrain[0] = strain[0] -fthird*I1;
    fDevStrain[1] = strain[1] -fthird*I1;
    fDevStrain[2] -= fthird*I1;
    fDevStrain[5] = strain[2];
    
    fDevStrain -= fPlasticStrain;
    
    energy += 0.5*(2.0*fMu*fDevStrain.ScalarProduct()+fKappa*I1*I1);
    energy += 0.5*falpha[0]*dK()*falpha[0];
    energy += 0.5*fBeta.ScalarProduct(fPlasticStrain);

    return energy;		
}

const dMatrixT& SSJ2LinHard2D::c_ijkl(void)
{
        s_ij();
     
    /*deviatoric part*/
	fModulus3D = 0.0;
	fModulus3D(0,0) = fModulus3D(1,1) =  fModulus3D(2,2) = 2.0*fMu*(1.0 - fthird);
	fModulus3D(3,3) = fModulus3D(4,4) =  fModulus3D(5,5) = fMu;
	fModulus3D(0,1) = fModulus3D(0,2) =  fModulus3D(1,2) = -2.0*fMu*fthird;
	fModulus3D(1,0) = fModulus3D(2,0) =  fModulus3D(2,1) = -2.0*fMu*fthird;

    /*volumetric part*/
	fModulus3D(0,0) += fKappa; fModulus3D(1,1) += fKappa; fModulus3D(2,2) += fKappa;
	fModulus3D(0,1) += fKappa; fModulus3D(0,2) += fKappa; fModulus3D(1,2) += fKappa;
	fModulus3D(1,0) += fKappa; fModulus3D(2,0) += fKappa; fModulus3D(2,1) += fKappa;

	int iteration = fSSMatSupport->IterationNumber();
	if (iteration > -1) /* elastic iteration */
	  fModulus3D += ModuliCorrection();
	
	fModulus = 0.0;
	fModulus(0,0) = fModulus3D(0,0); fModulus(1,1) = fModulus3D(1,1);
	fModulus(0,1) = fModulus3D(0,1); fModulus(1,0) = fModulus3D(1,0);	
	fModulus(2,2) = fModulus3D(5,5);
	fModulus(0,2) = fModulus3D(0,5); fModulus(1,2) = fModulus3D(1,5);
	fModulus(2,0) = fModulus3D(5,0); fModulus(2,1) = fModulus3D(5,1);

	return fModulus;
}

/* stress */
const dSymMatrixT& SSJ2LinHard2D::s_ij(void)
{
    const dSymMatrixT& strain = e();
    double I1 = strain.Trace();
    fDevStrain = 0.0;
    fDevStrain[0] = strain[0] - fthird*I1;
    fDevStrain[1] = strain[1] - fthird*I1;
    fDevStrain[2] = -fthird*I1;
    fDevStrain[5] = strain[2];
        
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());	

    fDevStrain -= fPlasticStrain_n;
    
    /* deviatoric part of trial stress */
    fStress3D = fDevStrain;
    fStress3D *= 2.0*fMu;
    /* modify trial stress (return mapping) */

    int iteration = fSSMatSupport->IterationNumber();
    if (iteration > -1)  
      fStress3D += StressCorrection(fStress3D);
    Store(element,CurrIP());
	
    /*add volumetric component*/
    fStress3D.PlusIdentity(fKappa*I1);
    
    fStress[0] = fStress3D[0]; fStress[1] = fStress3D[1]; fStress[2] = fStress3D[5];

    return fStress;	
}

/*Note to be called only during post processing*/
const dArrayT& SSJ2LinHard2D::InternalStrainVars(void)
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

const dArrayT& SSJ2LinHard2D::InternalStressVars(void)
{
        /*non-equilibrium components*/
        ElementCardT& element = CurrentElement();
        Load(element, CurrIP());
        
	const dSymMatrixT& strain = e();
	double I1 = strain.Trace();
	fDevStrain = 0.0;
	fDevStrain[0] = strain[0]-fthird*I1;
	fDevStrain[1] = strain[1]-fthird*I1;
	fDevStrain[2] -= fthird*I1;
	fDevStrain[5] = strain[2];
	fDevStrain -= fPlasticStrain;

	fStress3D = fDevStrain;
	fStress3D *= 2.0*fMu;
	fStress3D.PlusIdentity(fKappa*I1);
	
	double* pvar = fInternalStressVars.Pointer();

	*pvar++ = -dK()*falpha[0];

	*pvar++ = -fBeta[0];
	*pvar++ = -fBeta[1];
	*pvar++ = -fBeta[2];
	*pvar++ = -fBeta[3];
	*pvar++ = -fBeta[4];
	*pvar++ = -fBeta[5];

	*pvar++ = fStress3D[0];
	*pvar++ = fStress3D[1];
	*pvar++ = fStress3D[2];
	*pvar++ = fStress3D[3];
	*pvar++ = fStress3D[4];
	*pvar = fStress3D[5];

        return(fInternalStressVars);
}

int SSJ2LinHard2D::NumOutputVariables(void) const  { return kNumOutput; }
void SSJ2LinHard2D::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void SSJ2LinHard2D::ComputeOutput(dArrayT& output)
{
	/* stress tensor (loads element data and sets fStress) */
        /*non-equilibrium components*/
        ElementCardT& element = CurrentElement();
        Load(element, CurrIP());
        
	const dSymMatrixT& strain = e();
	double I1 = strain.Trace();
	fDevStrain = 0.0;
	fDevStrain[0] = strain[0]-fthird*I1;
	fDevStrain[1] = strain[1]-fthird*I1;
	fDevStrain[2] -= fthird*I1;
	fDevStrain[5] = strain[2];
	fDevStrain -= fPlasticStrain;

	fStress3D = fDevStrain;
	fStress3D *= 2.0*fMu;
	fStress3D.PlusIdentity(fKappa*I1);	

	/* pressure */
	output[2] = fStress3D.Trace()/3.0;
	
	/* deviatoric Von Mises stress */
	fStress3D.Deviatoric();
	double J2 = fStress3D.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[1] = sqrt(3.0*J2);

        output[0] = falpha[0];
}
