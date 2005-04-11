/* $Id: J2SSC0HardeningT.cpp,v 1.5.30.3 2005-04-11 19:40:37 thao Exp $ */
#include "J2SSC0HardeningT.h"

#include <iostream.h>
#include <math.h>

#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"
#include "ifstreamT.h"

/* hardening functions */
#include "CubicSplineT.h"
#include "LinearExponentialT.h"

/* class constants */

using namespace Tahoe;

const int    kNumInternal = 4; // number of internal variables
const double sqrt23       = sqrt(2.0/3.0);
const double kYieldTol    = 1.0e-10;

const int kNSD = 3;

/* constructor */
J2SSC0HardeningT::J2SSC0HardeningT(ifstreamT& in, int num_ip, double mu):
	fNumIP(num_ip),
	fmu(mu),
	ftheta(1.0),
	fK(NULL),
	fElasticStrain(kNSD),
	fStressCorr(kNSD),
	fModuliCorr(dSymMatrixT::NumValues(kNSD)),
	fRelStress(kNSD),
	fDevStrain(kNSD),
	fTensorTemp(dSymMatrixT::NumValues(kNSD))
{
	/* construct hardening function from stream */
	ConstructHardeningFunction(in);

	fYield = K(0.0);

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

/* destructor */
J2SSC0HardeningT::~J2SSC0HardeningT(void) { delete fK; };

/* returns elastic strain */
const dSymMatrixT& J2SSC0HardeningT::ElasticStrain(const dSymMatrixT& totalstrain,
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
const dSymMatrixT& J2SSC0HardeningT::StressCorrection(const dSymMatrixT& trialstrain,
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
		double alpha   = fInternal[kalpha];
		
		/* return mapping (single step) */
		if (ftrial > kYieldTol)
		{
			if (fType == kLinear)
			{
				/* plastic increment */
				dgamma = 0.5*ftrial/(fmu + dK(alpha)/3.0);
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
					double df = 2.0*(dK(alpha + sqrt23*dgamma)/3.0 + fmu);
					if (df < kSmall)
					{
						cout << "\n J2SSC0HardeningT::StressCorrection: consistency function is nonconvex" << endl;
						throw ExceptionT::kGeneralFail;
					}
				
					/* increment update */
					dgamma -= f/df;
					
					/* update condition */
					f = sqrt23*K(alpha + sqrt23*dgamma) - s_tr + 2.0*fmu*dgamma;
				}
				
				/* check for failure */
				if (count == max_iteration)
				{
					cout << "\n J2SSC0HardeningT::StressCorrection: local iteration failed after " 
					     << max_iteration << " iterations" << endl;
					throw ExceptionT::kGeneralFail;
				}
			}
	
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
const dMatrixT& J2SSC0HardeningT::ModuliCorrection(const ElementCardT& element,
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
void J2SSC0HardeningT::AllocateElement(ElementCardT& element)
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
	element.Dimension(i_size, d_size);
	
	/* initialize values */
	element.IntegerData() = kIsElastic;
	element.DoubleData()  = 0.0;
}

/***********************************************************************
* Protected
***********************************************************************/

double J2SSC0HardeningT::YieldCondition(const dSymMatrixT& relstress,
	double alpha) const
{
	return sqrt(relstress.ScalarProduct()) - sqrt23*K(alpha);
}

/* write parameters */
void J2SSC0HardeningT::Print(ostream& out) const
{
	/* hardening function parameters */
	out << " Hardening function:\n";
	fK->Print(out);
	
	/* print out spline coefficients */
	if (fType == kCubicSpline)
	{
#ifndef __NO_RTTI__
		const CubicSplineT* spline = dynamic_cast<const CubicSplineT*>(fK);
		if (spline)
		{
			/* spline coefficients */
			const dArray2DT& coefficients = spline->Coefficients();
			out << " Spline coefficients:\n";
			coefficients.WriteNumbered(out);
		}
		else
			out << " Error: could not cast hardening function to cubic spline" << endl;
#else /* __NO_RTTI__ */
		out << " Note: RTTI not available. Cannot write spline coefficients" << endl;
#endif /* __NO_RTTI__ */
	}
}

void J2SSC0HardeningT::PrintName(ostream& out) const
{
	out << "    J2 Isotropic/Kinematic\n";
	out << "    Hardening with Radial Return\n";
	out << "    Small Strain\n";
}

/* element level data */
void J2SSC0HardeningT::Update(ElementCardT& element)
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
		if (flags[ip] ==  kIsPlastic) /* plastic update */
		{
			/* do not repeat if called again */
			flags[ip] = kReset;
			/* NOTE: ComputeOutput writes the updated internal variables
			 *       for output even during iteration output, which is
			 *       called before UpdateHistory */

			/* fetch element data */
			LoadData(element, ip);
		
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
void J2SSC0HardeningT::Reset(ElementCardT& element)
{
	/* flag not to update again */
	(element.IntegerData()) = kReset;
}

/* Access History Variables */
const dSymMatrixT& J2SSC0HardeningT::Get_PlasticStrain(const ElementCardT& element, int ip)
{
         LoadData(element, ip);
         return fPlasticStrain;
}

/* load element data for the specified integration point */
const dSymMatrixT& J2SSC0HardeningT::Get_Beta(const ElementCardT& element, int ip)
{
         LoadData(element, ip);
         return fBeta;
}

/* load element data for the specified integration point */
const dArrayT& J2SSC0HardeningT::Get_Internal(const ElementCardT& element, int ip)
{
         LoadData(element, ip);
         return fInternal;
}

/***********************************************************************
* Private
***********************************************************************/

/* load element data for the specified integration point */
void J2SSC0HardeningT::LoadData(const ElementCardT& element, int ip)
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
	int offset    = stressdim*fNumIP;
	int dex       = ip*stressdim;
	
	fPlasticStrain.Alias(         dim, &d_array[           dex]);
	     fUnitNorm.Alias(         dim, &d_array[  offset + dex]);
	         fBeta.Alias(         dim, &d_array[2*offset + dex]);
	     fInternal.Alias(kNumInternal, &d_array[3*offset + ip*kNumInternal]);     	
}

/* returns 1 if the trial elastic strain state lies outside of the
* yield surface */
int J2SSC0HardeningT::PlasticLoading(const dSymMatrixT& trialstrain,
	ElementCardT& element, int ip)
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
dSymMatrixT& J2SSC0HardeningT::RelativeStress(const dSymMatrixT& trialstrain,
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

/* construct isotropic hardening function */
void J2SSC0HardeningT::ConstructHardeningFunction(ifstreamT& in)
{
	/* construct hardening function */
	int type;
	in >> type;
	switch (type)
	{
		case kLinear:
		{
			fType = kLinear;
			
			/* parameters */
			double yield = -1;
			double dK = -1;
			in >> yield >> dK;
			if (yield < 0) throw ExceptionT::kBadInputValue;
			
			dArray2DT points(2,2);
			points(0,0) = 0.0;
			points(0,1) = yield;
			points(1,0) = 1.0;
			points(1,1) = yield + dK;
			
			/* construct spline */
			fK = new CubicSplineT(points, CubicSplineT::kFreeRun);
			break;
		}
		case kLinearExponential:
		{
			fType = kLinearExponential;
			double a, b, c, d;
			a = b = c = d = 0.0;
			in >> a >> b >> c >> d;
			
			/* construct function */
			fK = new LinearExponentialT(a, b, c, d);
			break;
		}
		case kCubicSpline:
		{
			fType = kCubicSpline;
			
			/* point data */
			int num_points = -1;
			in >> num_points;
			if (num_points < 2) 
			{
				cout << "\n J2SSC0HardeningT::ConstructHardeningFunction: expecting at least 2 spline points:"
				     << num_points << endl;
				throw ExceptionT::kBadInputValue;
			}
			dArray2DT points(num_points, 2);
			in >> points;
			
			/* construct spline */
			fK = new CubicSplineT(points, CubicSplineT::kFreeRun);
			break;
		}
		default:
			cout << "\n J2SSC0HardeningT::ConstructHardeningFunction: unknown hardening function type: " 
			     << type << endl;
			throw ExceptionT::kBadInputValue;
	}
}
