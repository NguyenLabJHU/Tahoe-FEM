/* $Id: J2SSC0Hardening1DT.cpp,v 1.1 2004-07-20 23:38:26 rdorgan Exp $ */
#include "J2SSC0Hardening1DT.h"

#include "iArrayT.h"
#include "ElementCardT.h"

#include <math.h>

using namespace Tahoe;

/* class constants */
const int    kNumInternal = 4; // number of internal variables
const double kYieldTol    = 1.0e-10;

const int kNSD = 1;

/* constructor */
J2SSC0Hardening1DT::J2SSC0Hardening1DT(void):
	ParameterInterfaceT("J2_small_strain_hardening_1D"),
	fElasticStrain(kNSD),
	fStressCorr(kNSD),
	fModuliCorr(dSymMatrixT::NumValues(kNSD)),
	fRelStress(kNSD)
{

}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* returns elastic strain */
const dSymMatrixT& J2SSC0Hardening1DT::ElasticStrain(const dSymMatrixT& totalstrain,
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
const dSymMatrixT& J2SSC0Hardening1DT::StressCorrection(const dSymMatrixT& trialstrain,
	ElementCardT& element, double young, int nip, int ip)
{
	const char caller[] = "J2SSC0Hardening1DT::StressCorrection";

	/* check consistency and initialize plastic element */
	if (PlasticLoading(trialstrain, element, young, nip, ip) &&
	    !element.IsAllocated())
	{
		/* new element */
		AllocateElement(element, nip);
					
		/* initialize element data */
		PlasticLoading(trialstrain, element, young, nip, ip);
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
				dgamma = ftrial/(young + dK(alpha));
			}
			else /* local Newton iteration */
			{
				double s_tr = ftrial + K(alpha);
				double f =-ftrial;
			
				dgamma = 0.0;
				int max_iteration = 15;
				int count = 0;
				while (fabs(f) > kYieldTol && ++count <= max_iteration)
				{
					/* stiffness */
					double df = dK(alpha + dgamma) + young;
					if (df < kSmall)
						ExceptionT::GeneralFail(caller, "yield function is nonconvex");
				
					/* increment update */
					dgamma -= f/df;
					
					/* update condition */
					f = K(alpha + dgamma) - s_tr + young*dgamma;
				}
				
				/* check for failure */
				if (count == max_iteration)
					ExceptionT::GeneralFail(caller, "local iteration failed after %d iterations", max_iteration);
			}
	
			/* plastic increment stress correction */
			fStressCorr.SetToScaled(-young*dgamma, fUnitNorm);
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
const dMatrixT& J2SSC0Hardening1DT::ModuliCorrection(const ElementCardT& element,
	double young, int nip, int ip)
{
#pragma unused(element, young, nip, ip)
	/* initialize */
	fModuliCorr = 0.0;
	return fModuliCorr;
}	
	 	
/* return a pointer to a new plastic element object constructed with
* the data from element */
void J2SSC0Hardening1DT::AllocateElement(ElementCardT& element, int nip)
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

/* element level data */
void J2SSC0Hardening1DT::Update(ElementCardT& element, int nip)
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
			fInternal[kalpha] += dgamma;
	
			/* plastic strain increment	*/
			fPlasticStrain.AddScaled(dgamma, fUnitNorm);			
		}
}

/* resets to the last converged solution */
void J2SSC0Hardening1DT::Reset(ElementCardT& element, int nip)
{
#pragma unused(nip)
	/* flag not to update again */
	(element.IntegerData()) = kReset;
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* load element data for the specified integration point */
void J2SSC0Hardening1DT::LoadData(const ElementCardT& element, int nip, int ip)
{
	/* check */
	if (!element.IsAllocated()) throw ExceptionT::kGeneralFail;

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
int J2SSC0Hardening1DT::PlasticLoading(const dSymMatrixT& trialstrain,
	ElementCardT& element, double young, int nip, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated())
		return( YieldCondition(RelativeStress(trialstrain, element, young), 0.0)
					> kYieldTol );
	/* already plastic */
	else
	{
		/* get flags */
		iArrayT& Flags = element.IntegerData();
			
		/* load internal variables */
		LoadData(element, nip, ip);
		
		fInternal[kftrial] = YieldCondition(RelativeStress(trialstrain, element, young),
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
dSymMatrixT& J2SSC0Hardening1DT::RelativeStress(const dSymMatrixT& trialstrain,
	const ElementCardT& element, double young)
{
	/* element is plastic - element data loaded */
	if (element.IsAllocated())
	{
		/* elastic stress */
		fRelStress.SetToScaled(young, trialstrain);

		/* kinematic hardening */
		fRelStress -= fBeta;
	}
	else
		fRelStress.SetToScaled(young, trialstrain);

	return fRelStress;
}
