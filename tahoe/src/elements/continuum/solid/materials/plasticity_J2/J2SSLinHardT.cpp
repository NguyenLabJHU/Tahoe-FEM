/* $Id: J2SSLinHardT.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (02/12/1997)                                          */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */

#include "J2SSLinHardT.h"

#include <iostream.h>
#include <math.h>

#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"

/* class constants */
const int    kNumInternal = 4; // number of internal variables
const double sqrt23       = sqrt(2.0/3.0);
const double kYieldTol    = 1.0e-10;

const int kNSD = 3;

/* constructor */
J2SSLinHardT::J2SSLinHardT(ifstreamT& in, int num_ip, double mu):
	J2PrimitiveT(in),
	fNumIP(num_ip),
	fmu(mu),
	fElasticStrain(kNSD),
	fStressCorr(kNSD),
	fModuliCorr(dSymMatrixT::NumValues(kNSD)),
	fRelStress(kNSD),
	fDevStrain(kNSD),
	fTensorTemp(dSymMatrixT::NumValues(kNSD))
{

}

/* returns elastic strain */
const dSymMatrixT& J2SSLinHardT::ElasticStrain(const dSymMatrixT& totalstrain,
	const ElementCardT& element, int ip)
{	
	/* remove plastic strain */
	if (element.IsAllocated())
	{
		/* load internal variables */
		LoadData(element, ip);

		/* compute elastic strain */
		fElasticStrain.DiffOf(totalstrain, fPlasticStrain);
	
		return fElasticStrain;
	}	
	/* no plastic strain */
	else	
		return totalstrain;
}

/* return the correction to stress vector computed by the mapping the
* stress back to the yield surface, if needed */
const dSymMatrixT& J2SSLinHardT::StressCorrection(const dSymMatrixT& trialstrain,
	ElementCardT& element, int ip)
{
	/* check consistency and initialize plastic element */
	if (PlasticLoading(trialstrain, element, ip) &&
	    !element.IsAllocated())
	{
		/* new element */
		AllocateElement(element);
					
		/* initialize element data */
		PlasticLoading(trialstrain, element, ip);
	}

	/* initialize */
	fStressCorr = 0.0;
	
	if (element.IsAllocated())
	{		
		/* fetch data */
		double  ftrial = fInternal[kftrial];
		double& dgamma = fInternal[kdgamma];
		
		/* return mapping (single step) */
		if (ftrial > kYieldTol)
		{
			/* plastic increment */
			dgamma = 0.5*ftrial/(fmu + fH_bar/3.0);
	
			/* plastic increment stress correction */
			fStressCorr.SetToScaled(-2.0*fmu*dgamma, fUnitNorm);
		}
		else
			dgamma = 0.0;			
	}
		
	return fStressCorr;
}	

/* return the correction to moduli due to plasticity (if any)
*
* Note: Return mapping occurs during the call to StressCorrection.
*       The element passed in is already assumed to carry current
*       internal variable values */
const dMatrixT& J2SSLinHardT::ModuliCorrection(const ElementCardT& element,
	int ip)
{
	/* initialize */
	fModuliCorr = 0.0;

	if (element.IsAllocated() &&
	   (element.IntegerData())[ip] == kIsPlastic)
	{
		/* load internal variables */
		LoadData(element,ip);
		
		/* compute constants */
		double alpha    = fInternal[kalpha];
		double thetahat = 2.0*fmu*fInternal[kdgamma]/
		                          fInternal[kstressnorm];
		double thetabar = (3.0/(3.0 + (dK(alpha) + dH(alpha))/fmu)) - thetahat;
		
		/* moduli corrections */
		fTensorTemp.ReducedIndexDeviatoric();
		fModuliCorr.AddScaled(-2.0*fmu*thetahat, fTensorTemp);
					
		fTensorTemp.Outer(fUnitNorm,fUnitNorm);
		fModuliCorr.AddScaled(-2.0*fmu*thetabar, fTensorTemp);
	}

	return fModuliCorr;
}	
	 	
/* return a pointer to a new plastic element object constructed with
* the data from element */
void J2SSLinHardT::AllocateElement(ElementCardT& element)
{
	/* determine storage */
	int i_size = 0;
	i_size += fNumIP; //fFlags

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fPlasticStrain
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fUnitNorm
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fBeta
	d_size += kNumInternal*fNumIP;          //fInternal

	/* construct new plastic element */
	element.Allocate(i_size, d_size);
	
	/* initialize values */
	element.IntegerData() = kIsElastic;
	element.DoubleData()  = 0.0;
}

/***********************************************************************
* Protected
***********************************************************************/

void J2SSLinHardT::PrintName(ostream& out) const
{
	/* inherited */
	J2PrimitiveT::PrintName(out);

	out << "    Small Strain\n";
}

/* element level data */
void J2SSLinHardT::Update(ElementCardT& element)
{
	/* get flags */
	iArrayT& flags = element.IntegerData();

	/* check if reset state (is same for all ip) */
	if (flags[0] == kReset)
	{
		flags = kIsElastic; //don't update again
		return;
	}

	/* update plastic variables */
	for (int ip = 0; ip < fNumIP; ip++)
	{
		/* fetch element data */
		LoadData(element, ip);
		
		/* plastic increment */
		double& dgamma = fInternal[kdgamma];
		
		/* internal variables */
		fInternal[kalpha] += sqrt23*dgamma;
	
		/* kinematic hardening */
		fBeta.AddScaled(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
				
		/* dev plastic strain increment	*/
		fPlasticStrain.AddScaled(dgamma, fUnitNorm);			
	}
}

/* resets to the last converged solution */
void J2SSLinHardT::Reset(ElementCardT& element)
{
	/* flag not to update again */
	(element.IntegerData()) = kReset;
}

/***********************************************************************
* Private
***********************************************************************/

/* load element data for the specified integration point */
void J2SSLinHardT::LoadData(const ElementCardT& element, int ip)
{
	/* check */
	if (!element.IsAllocated()) throw eGeneralFail;

	/* fetch arrays */
	dArrayT& d_array = element.DoubleData();
	
	/* decode */
	int stressdim = dSymMatrixT::NumValues(kNSD);
	int offset    = stressdim*fNumIP;
	int dex       = ip*stressdim;
	
	fPlasticStrain.Set(        kNSD, &d_array[           dex]);
	     fUnitNorm.Set(        kNSD, &d_array[  offset + dex]);
	         fBeta.Set(        kNSD, &d_array[2*offset + dex]);
	     fInternal.Set(kNumInternal, &d_array[3*offset + ip*kNumInternal]);     	
}

/* returns 1 if the trial elastic strain state lies outside of the
* yield surface */
int J2SSLinHardT::PlasticLoading(const dSymMatrixT& trialstrain,
	const ElementCardT& element, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated())
		return( YieldCondition(RelativeStress(trialstrain, element), 0.0)
					> kYieldTol );
	/* already plastic */
	else
	{
		/* get flags */
		iArrayT& Flags = element.IntegerData();
			
		/* load internal variables */
		LoadData(element, ip);
		
		fInternal[kftrial] = YieldCondition(RelativeStress(trialstrain,element),
						    fInternal[kalpha]);

		/* plastic */
		if (fInternal[kftrial] > kYieldTol)
		{		
			/* compute unit normal */
			double& norm = fInternal[kstressnorm];

			norm = sqrt(fRelStress.ScalarProduct());
			fUnitNorm.SetToScaled(1.0/norm, fRelStress);
		
			/* set flag */
			Flags[ip] = kIsPlastic;
	
			return 1;
		}
		/* elastic */
		else
		{
			/* set flag */
			Flags[ip] = kIsElastic;
			
			return 0;
		}
	}
}	

/* computes the relative stress corresponding for the given element
* and elastic strain.  The functions returns a reference to the
* relative stress in fRelStress */
dSymMatrixT& J2SSLinHardT::RelativeStress(const dSymMatrixT& trialstrain,
	const ElementCardT& element)
{
	/* deviatoric strain */
	fDevStrain.Deviatoric(trialstrain);

	/* element is plastic - element data loaded */
	if (element.IsAllocated())
	{
		/* deviatoric part of elastic stress */
		fRelStress.SetToScaled(2.0*fmu, fDevStrain);

		/* kinematic hardening */
		fRelStress -= fBeta;
	}
	else
		fRelStress.SetToScaled(2.0*fmu, fDevStrain);

	return fRelStress;
}
