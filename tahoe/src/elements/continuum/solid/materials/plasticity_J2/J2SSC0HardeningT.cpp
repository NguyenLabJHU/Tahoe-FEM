/* $Id: J2SSC0HardeningT.cpp,v 1.7 2005-04-06 23:32:45 thao Exp $ */
#include "J2SSC0HardeningT.h"
#include "SSMatSupportT.h"
#include "ElementCardT.h"

#include "iArrayT.h"
#include "ElementCardT.h"

#include <math.h>

using namespace Tahoe;

/* class constants */
const int    kNumInternal = 4; // number of internal variables
const double sqrt23       = sqrt(2.0/3.0);
const double kYieldTol    = 1.0e-10;

const int kNSD = 3;

/* constructor */
J2SSC0HardeningT::J2SSC0HardeningT(void):
	ParameterInterfaceT("J2_small_strain_hardening"),
	fElasticStrain(kNSD),
	fStressCorr(kNSD),
	fModuliCorr(dSymMatrixT::NumValues(kNSD)),
	fRelStress(kNSD),
	fDevStrain(kNSD),
	fTensorTemp(dSymMatrixT::NumValues(kNSD))
{
  int dofcount = 0;
  fInternalDOF.Dimension(3);
  fInternalDOF[0] = 1; //alpha, -K*alpha
  dofcount++;
  fInternalDOF[1] = dSymMatrixT::NumValues(kNSD); //plasticstrain, -beta
  dofcount += dSymMatrixT::NumValues(kNSD);
  fInternalDOF[2] = dSymMatrixT::NumValues(kNSD); //plasticstrain, stress
  dofcount += dSymMatrixT::NumValues(kNSD);
  
  fInternalStressVars.Dimension(dofcount);
  fInternalStressVars = 0.0;
  
  fInternalStrainVars.Dimension(dofcount);
  fInternalStrainVars = 0.0;
}

/* returns elastic strain */
const dSymMatrixT& J2SSC0HardeningT::ElasticStrain(const dSymMatrixT& totalstrain,
	const ElementCardT& element, int nip, int ip)
{	
	/* remove plastic strain */
	if (element.IsAllocated())
	{
		/* load internal variables */
		LoadData(element, nip, ip);

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
const dSymMatrixT& J2SSC0HardeningT::StressCorrection(const dSymMatrixT& trialstrain,
	ElementCardT& element, double mu, int nip, int ip)
{
	const char caller[] = "J2SSC0HardeningT::StressCorrection";

	/* check consistency and initialize plastic element */
	if (PlasticLoading(trialstrain, element, mu, nip, ip) &&
	    !element.IsAllocated())
	{
		/* new element */
		AllocateElement(element, nip);
					
		/* initialize element data */
		PlasticLoading(trialstrain, element, mu, nip, ip);
	}

	/* initialize */
	fStressCorr = 0.0;
	
	if (element.IsAllocated())
	{		
		/* fetch data */
		double  ftrial = fInternal[kftrial];
		double& dgamma = fInternal[kdgamma];
		double alpha   = fInternal[kalpha];
		
		/* return mapping (single step) */
		if (ftrial > kYieldTol)
		{
			if (fIsLinear)
			{
				/* plastic increment */
				dgamma = 0.5*ftrial/(mu + dK(alpha)/3.0);
			}
			else /* local Newton iteration */
			{
				double s_tr = ftrial + sqrt23*K(alpha);
				double f =-ftrial;
			
				dgamma = 0.0;
				int max_iteration = 15;
				int count = 0;
				while (fabs(f) > kYieldTol && ++count <= max_iteration)
				{
					/* stiffness */
					double df = 2.0*(dK(alpha + sqrt23*dgamma)/3.0 + mu);
					if (df < kSmall)
						ExceptionT::GeneralFail(caller, "yield function is nonconvex");
				
					/* increment update */
					dgamma -= f/df;
					
					/* update condition */
					f = sqrt23*K(alpha + sqrt23*dgamma) - s_tr + 2.0*mu*dgamma;
				}
				
				/* check for failure */
				if (count == max_iteration)
					ExceptionT::GeneralFail(caller, "local iteration failed after %d iterations", max_iteration);
			}
	
			/* plastic increment stress correction */
			fStressCorr.SetToScaled(-2.0*mu*dgamma, fUnitNorm);
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
const dMatrixT& J2SSC0HardeningT::ModuliCorrection(const ElementCardT& element,
	double mu, int nip, int ip)
{
	/* initialize */
	fModuliCorr = 0.0;

	if (element.IsAllocated() &&
	   (element.IntegerData())[ip] == kIsPlastic)
	{
		/* load internal variables */
		LoadData(element, nip, ip);
		
		/* compute constants */
		double alpha    = fInternal[kalpha];
		double thetahat = 2.0*mu*fInternal[kdgamma]/
		                          fInternal[kstressnorm];
		double thetabar = (3.0/(3.0 + (dK(alpha) + dH(alpha))/mu)) - thetahat;
		
		/* moduli corrections */
		fTensorTemp.ReducedIndexDeviatoric();
		fModuliCorr.AddScaled(-2.0*mu*thetahat, fTensorTemp);
					
		fTensorTemp.Outer(fUnitNorm,fUnitNorm);
		fModuliCorr.AddScaled(-2.0*mu*thetabar, fTensorTemp);
	}

	return fModuliCorr;
}	
	 	
/* return a pointer to a new plastic element object constructed with
* the data from element */
void J2SSC0HardeningT::AllocateElement(ElementCardT& element, int nip)
{
	/* determine storage */
	int i_size = 0;
	i_size += nip; //fFlags

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*nip; //fPlasticStrain
	d_size += dSymMatrixT::NumValues(kNSD)*nip; //fUnitNorm
	d_size += dSymMatrixT::NumValues(kNSD)*nip; //fBeta
	d_size += kNumInternal*nip;          //fInternal

	/* construct new plastic element */
	element.Dimension(i_size, d_size);
	
	/* initialize values */
	element.IntegerData() = kIsElastic;
	element.DoubleData()  = 0.0;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* element level data */
void J2SSC0HardeningT::Update(ElementCardT& element, int nip)
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
	for (int ip = 0; ip < nip; ip++)
		if (flags[ip] ==  kIsPlastic) /* plastic update */
		{
			/* do not repeat if called again */
			flags[ip] = kReset;
			/* NOTE: ComputeOutput writes the updated internal variables
			 *       for output even during iteration output, which is
			 *       called before UpdateHistory */

			/* fetch element data */
			LoadData(element, nip, ip);
		
			/* plastic increment */
			double& dgamma = fInternal[kdgamma];
		
			/* internal variables */
			fInternal[kalpha] += sqrt23*dgamma;
	
			/* kinematic hardening */
			//fBeta.AddScaled(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
				
			/* dev plastic strain increment	*/
			fPlasticStrain.AddScaled(dgamma, fUnitNorm);			
		}
}

/* resets to the last converged solution */
void J2SSC0HardeningT::Reset(ElementCardT& element, int nip)
{
#pragma unused(nip)
	/* flag not to update again */
	(element.IntegerData()) = kReset;
}

/* Access History Variables */
const dSymMatrixT& J2SSC0HardeningT::Get_PlasticStrain(const ElementCardT& element, int nip, int ip)
{
         LoadData(element, nip, ip);
	 return (fPlasticStrain);
}

/* load element data for the specified integration point */
const dSymMatrixT& J2SSC0HardeningT::Get_Beta(const ElementCardT& element, int nip, int ip)
{
         LoadData(element, nip, ip);
	 return (fBeta);
}

/* load element data for the specified integration point */
const dArrayT& J2SSC0HardeningT::Get_Internal(const ElementCardT& element, int nip, int ip)
{
         LoadData(element, nip, ip);
	 return (fInternal);
}

/***********************************************************************
* Private
***********************************************************************/

/* load element data for the specified integration point */
void J2SSC0HardeningT::LoadData(const ElementCardT& element, int nip, int ip)
{
	/* check */
  if (!element.IsAllocated()) { 
    cout << "\nElement is not allocated: ";
    throw ExceptionT::kGeneralFail;
  }

	/* fetch arrays */
	const dArrayT& d_array = element.DoubleData();
	
	/* decode */
	dSymMatrixT::DimensionT dim = dSymMatrixT::int2DimensionT(kNSD);
	int stressdim = dSymMatrixT::NumValues(kNSD);
	int offset    = stressdim*nip;
	int dex       = ip*stressdim;
	
	fPlasticStrain.Alias(         dim, &d_array[           dex]);
	     fUnitNorm.Alias(         dim, &d_array[  offset + dex]);
	         fBeta.Alias(         dim, &d_array[2*offset + dex]);
	     fInternal.Alias(kNumInternal, &d_array[3*offset + ip*kNumInternal]);     	
}

/* returns 1 if the trial elastic strain state lies outside of the
* yield surface */
int J2SSC0HardeningT::PlasticLoading(const dSymMatrixT& trialstrain,
	ElementCardT& element, double mu, int nip, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated())
		return( YieldCondition(RelativeStress(trialstrain, element, mu), 0.0)
					> kYieldTol );
	/* already plastic */
	else
	{
		/* get flags */
		iArrayT& Flags = element.IntegerData();
			
		/* load internal variables */
		LoadData(element, nip, ip);
		
		fInternal[kftrial] = YieldCondition(RelativeStress(trialstrain, element, mu),
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
dSymMatrixT& J2SSC0HardeningT::RelativeStress(const dSymMatrixT& trialstrain,
	const ElementCardT& element, double mu)
{
	/* deviatoric strain */
	fDevStrain.Deviatoric(trialstrain);

	/* element is plastic - element data loaded */
	if (element.IsAllocated())
	{
		/* deviatoric part of elastic stress */
		fRelStress.SetToScaled(2.0*mu, fDevStrain);

		/* kinematic hardening */
		fRelStress -= fBeta;
	}
	else
		fRelStress.SetToScaled(2.0*mu, fDevStrain);

	return fRelStress;
}
