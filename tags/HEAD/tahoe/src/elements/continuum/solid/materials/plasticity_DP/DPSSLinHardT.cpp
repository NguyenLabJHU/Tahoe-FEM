/* $Id: DPSSLinHardT.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: myip (06/01/1999)                                             */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with isotropic hardening                           */

#include "DPSSLinHardT.h"
#include <iostream.h>
#include <math.h>
#include "ElasticT.h"
#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"

/* class constants */
const int    kNumInternal = 5; // number of internal variables
const double sqrt23       = sqrt(2.0/3.0);
const double sqrt32       = sqrt(3.0/2.0);
const double kYieldTol    = 1.0e-10;
const int kNSD = 3;

/* constructor */  //**mien**//
DPSSLinHardT::DPSSLinHardT(ifstreamT& in, int num_ip, double mu, double lambda):
	DPPrimitiveT(in),
	fNumIP(num_ip),
	fmu(mu),
	flambda(lambda),
	fkappa(flambda + (2.0/3.0*fmu)),
	fX_H(3.0*(fmu+ffriction*fdilation*fkappa) + (3.0*fdilation*fdilation*fK_prime+fH_prime)),
	fElasticStrain(kNSD),
	fStressCorr(kNSD),
	fModuliCorr(dSymMatrixT::NumValues(kNSD)),
	fDevStress(kNSD),
	fMeanStress(0.0 ),
	fDevStrain(kNSD),
	fTensorTemp(dSymMatrixT::NumValues(kNSD)),
IdentityTensor2(kNSD), //**mien
	One(kNSD)
{
/* initialize constant tensor */
One.Identity();
}

/* returns elastic strain */
const dSymMatrixT& DPSSLinHardT::ElasticStrain(const dSymMatrixT& totalstrain,
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
const dSymMatrixT& DPSSLinHardT::StressCorrection(const dSymMatrixT& trialstrain,
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
		        dgamma = ftrial/fX_H;

			/* plastic increment stress correction */
			fStressCorr.PlusIdentity(-sqrt(3.0)*fdilation*fkappa*dgamma);
			fStressCorr.AddScaled(-sqrt(6.0)*fmu*dgamma, fUnitNorm);

			//TEMP
			//construct the trial stress
			dSymMatrixT stress(DeviatoricStress(trialstrain, element));
			stress.PlusIdentity(MeanStress(trialstrain,element));
			
			//corrected stress and internal variables
			stress += fStressCorr;
			double a_dev = fInternal[kalpha_dev] - fH_prime*dgamma;
			double a_vol = fInternal[kalpha_vol] - sqrt(3.0)*fdilation*fK_prime*dgamma ;

			// evaluate plastic consistency
			double p = stress.Trace()/3.0;
			//cout << " pressure part of stress   = " << p << endl;
			double f = YieldCondition(stress.Deviatoric(), p, a_dev, a_vol);
			//cout << " yield function = " << f << endl;
			//TEMP
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
const dMatrixT& DPSSLinHardT::ModuliCorrection(const ElementCardT& element,
	int ip)
{
	/* initialize */
	fModuliCorr = 0.0;

	if (element.IsAllocated() &&
	   (element.IntegerData())[ip] == kIsPlastic)
	{
		/* load internal variables */
	  	LoadData(element,ip);
		
		double c1  = -3.0*ffriction*fdilation*fkappa*fkappa/fX_H;
		       c1 += (4.0/3.0)*sqrt32*fmu*fmu*fInternal[kdgamma]/fInternal[kstressnorm];
		double c2  = -sqrt(6.0)*fmu*fmu*fInternal[kdgamma]/fInternal[kstressnorm];
		double c3  = -(3.0/2.0)/fX_H + sqrt32*fInternal[kdgamma]/fInternal[kstressnorm];
		       c3 *= 4.0*fmu*fmu;
		double c4  = -3.0*sqrt(2.0)*fkappa*fmu/fX_H;

		fTensorTemp.Outer(One, One);
		fModuliCorr.AddScaled(c1, fTensorTemp);

	        fTensorTemp.ReducedIndexI();
		fModuliCorr.AddScaled(2.0*c2, fTensorTemp);

		fTensorTemp.Outer(fUnitNorm,fUnitNorm);
		fModuliCorr.AddScaled(c3, fTensorTemp);

		fTensorTemp.Outer(One,fUnitNorm);
		fModuliCorr.AddScaled(fdilation*c4, fTensorTemp);

		fTensorTemp.Outer(fUnitNorm,One);
		fModuliCorr.AddScaled(ffriction*c4, fTensorTemp);
	}
	return fModuliCorr;
}	
	 	
/* return a pointer to a new plastic element object constructed with
* the data from element */
void DPSSLinHardT::AllocateElement(ElementCardT& element)
{
	/* determine storage */
	int i_size = 0;
	i_size += fNumIP; //fFlags

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fPlasticStrain
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fUnitNorm
	//	d_size += NumValues(kNSD)*fNumIP; //fBeta
	d_size += kNumInternal*fNumIP;          //fInternal

	/* construct new plastic element */
	element.Allocate(i_size, d_size);
	
	/* initialize values */
	element.IntegerData() = kIsElastic;
	element.DoubleData()  = 0.0;  // initialize all double types to 0.0
}

/***********************************************************************
* Protected
***********************************************************************/

void DPSSLinHardT::PrintName(ostream& out) const
{
	/* inherited */
	DPPrimitiveT::PrintName(out);

	out << "    Small Strain\n";
}

/* element level data */
void DPSSLinHardT::Update(ElementCardT& element)
{
	/* get flags */
	iArrayT& Flags = element.IntegerData();

	/* check if reset state */
	if (Flags[0] == kReset)
	{
		Flags = kIsElastic; //don't update again
		return;
	}

	/* update plastic variables */
	for (int ip = 0; ip < fNumIP; ip++)
	{
	/* fetch element data */
		LoadData(element, ip);
	
	/* plastic increment */
		double& dgamma = fInternal[kdgamma];
		//cout << fInternal[kdgamma] << endl;
		
	/* internal variables */
		fInternal[kalpha_dev] -= fH_prime*dgamma;
		fInternal[kalpha_vol] -= sqrt(3.0)*fdilation*fK_prime*dgamma;

	/* dev plastic strain increment	*/
		fPlasticStrain.AddScaled( sqrt32*dgamma, fUnitNorm );
		fPlasticStrain.AddScaled( fdilation*dgamma/sqrt(3.0), One );
	}
}

/* resets to the last converged solution */
void DPSSLinHardT::Reset(ElementCardT& element)
{
	/* flag not to update again */
	(element.IntegerData()) = kReset;
}

/***********************************************************************
* Private
***********************************************************************/

/* load element data for the specified integration point */
void DPSSLinHardT::LoadData(const ElementCardT& element, int ip)
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
	     fInternal.Set(kNumInternal, &d_array[2*offset + ip*kNumInternal]);
}

/* returns 1 if the trial elastic strain state lies outside of the
* yield surface */
int DPSSLinHardT::PlasticLoading(const dSymMatrixT& trialstrain,
	const ElementCardT& element, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated())
		return( YieldCondition(DeviatoricStress(trialstrain,element),
				       MeanStress(trialstrain,element), 0.0, 0.0) > kYieldTol );  //**mien
/* already plastic */
	else
	{
		/* get flags */
		iArrayT& Flags = element.IntegerData();
			
		/* load internal variables */
		LoadData(element, ip);
		
		fInternal[kftrial] = YieldCondition(DeviatoricStress(trialstrain,element),
		    MeanStress(trialstrain,element),fInternal[kalpha_dev],fInternal[kalpha_vol]);  //**mien

		/* plastic */
		if (fInternal[kftrial] > kYieldTol)
		{		
			/* compute unit normal */
			double& norm = fInternal[kstressnorm];

			norm = sqrt(fDevStress.ScalarProduct());  //**mien
			fUnitNorm.SetToScaled(1.0/norm, fDevStress);  //**mien
		
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
* relative stress in fDevStress */
dSymMatrixT& DPSSLinHardT::DeviatoricStress(const dSymMatrixT& trialstrain,
	const ElementCardT& element)
{
#pragma unused(element)

	/* deviatoric strain */
	fDevStrain.Deviatoric(trialstrain);

	/* compute deviatoric elastic stress */
	fDevStress.SetToScaled(2.0*fmu,fDevStrain);  //**mien

	return fDevStress;
}

/* computes the hydrostatic (mean) stress. */    //**mien**//
double DPSSLinHardT::MeanStress(const dSymMatrixT& trialstrain,
	const ElementCardT& element)
{
#pragma unused(element)

fMeanStress = fkappa*trialstrain.Trace();
return fMeanStress;
}
