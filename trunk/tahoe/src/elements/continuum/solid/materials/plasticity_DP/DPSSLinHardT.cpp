/* $Id: DPSSLinHardT.cpp,v 1.19 2006-04-25 18:15:53 regueiro Exp $ */
/* created: myip (06/01/1999) */
#include "DPSSLinHardT.h"
#include <iostream.h>
#include <math.h>

#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* class constants */
const int    kNumInternal = 5; // number of internal state variables
const double sqrt23       = sqrt(2.0/3.0);
const double sqrt32       = sqrt(3.0/2.0);
const double kYieldTol    = 1.0e-10;
const int    kNSD         = 3;

/* constructor */
DPSSLinHardT::DPSSLinHardT(int num_ip, double mu, double lambda):
	fNumIP(num_ip),
	fmu(mu),
	flambda(lambda),
	fkappa(flambda + (2.0/3.0*fmu)),
	fMeanStress(0.0)
{
	SetName("DP_SS_linear_hardening");
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

/* return correction to stress vector computed by mapping the
 * stress back to the yield surface, if needed */
const dSymMatrixT& DPSSLinHardT::StressCorrection(
        		const dSymMatrixT& trialstrain, 
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
		double& dgamma2 = fInternal[kdgamma2];

		/* return mapping (single step) */
		if (ftrial > kYieldTol)
		{

		/* plastic increment */
		dgamma = ftrial/fX_H;

		//construct the trial deviatoric stress and 
		//calculate norm of deviatoric trial stress
		dSymMatrixT stress(DeviatoricStress(trialstrain, element));
		double devstressnorm=sqrt(stress.ScalarProduct());
	
		if ( devstressnorm!=0 && 1-sqrt(6.0)*fmu*dgamma/devstressnorm >= 0)
		{
			/* plastic increment stress correction */
			dgamma2=0.0;
			fStressCorr.PlusIdentity(-sqrt(3.0)*fdilation*fkappa*dgamma);
			fStressCorr.AddScaled(-sqrt(6.0)*fmu*dgamma, fUnitNorm);
		}
		else
		{
			dgamma=devstressnorm/(sqrt(6.0)*fmu);
			//cout << "dgamma = " << dgamma << endl;
			double totalcohesion=sqrt(3.0)*falpha_bar-fInternal[kalpha] + fH_prime*dgamma;	    
			dgamma2=(MeanStress(trialstrain,element)-totalcohesion/ffriction
		    		-sqrt(3.0)*fkappa*fdilation*dgamma)/fkappa;
		    		
			fStressCorr.PlusIdentity(-1*MeanStress(trialstrain,element));
			fStressCorr.PlusIdentity(totalcohesion/(sqrt(3.0)*ffriction));
			fStressCorr.AddScaled(-1.0, stress);
		}

		//augment stress to full trial state
		stress.PlusIdentity(MeanStress(trialstrain,element));
		
		//corrected stress and internal state variables
		stress += fStressCorr;

		//cout << " stress = \n";
		//cout << stress << endl;
		/*
		double a = fInternal[kalpha] - fH_prime*dgamma;

		// evaluate plastic consistency
		double p = stress.Trace()/3.0;
		//cout << " pressure  = " << p << endl;
		
		double f = YieldCondition(stress.Deviatoric(), p, a);
		//cout << " yield function = " << f << endl;
		*/
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
		/* load internal state variables */
		LoadData(element,ip);
		
		if (fInternal[kdgamma2]==0.0)
		{
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

			//cout << "fModuliCorr = " << fModuliCorr << endl;
		}
		else
		{
		    //return fModuliCorr;

		    //cout << "In vertex region \n";

		    //cout << " kstressnorm = " << fInternal[kstressnorm] << endl;
		    double c1 = -fkappa;
		    c1 += 2.0/3.0*fmu;
		    
		    double c2 = -fmu;
		    
		    //double c3  = sqrt32*fInternal[kdgamma]/fInternal[kstressnorm];
		    //c3 *= 4.0*fmu*fmu;
		    //c3 += -2.0*fmu;
		    
		    double c4 = sqrt(2.0)*fH_prime/(3*ffriction);

		    fTensorTemp.Outer(One, One);
		    fModuliCorr.AddScaled(c1, fTensorTemp);
		
		    fTensorTemp.ReducedIndexI();
		    fModuliCorr.AddScaled(2.0*c2, fTensorTemp);
		      
		    //cout << "fModuliCorr = " << fModuliCorr << endl;
	
			//fTensorTemp.Outer(fUnitNorm,fUnitNorm);
		    //fModuliCorr.AddScaled(c3, fTensorTemp);
	
			if ( fInternal[kdgamma] != 0.0)
			{
				fTensorTemp.Outer(One, fUnitNorm);
				fModuliCorr.AddScaled(c4, fTensorTemp);
			}
   
			//cout << "fModuliCorr = " << fModuliCorr << endl;
    
		}
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
	d_size += kNumInternal*fNumIP;        //fInternal

	/* construct new plastic element */
	element.Dimension(i_size, d_size);
	
	/* initialize values */
	element.IntegerData() = kIsElastic;
	element.DoubleData()  = 0.0;  // initialize all double types to 0.0
}

/* accept parameter list */
void DPSSLinHardT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	DPPrimitiveT::TakeParameterList(list);

	/* dimension work space */
	fElasticStrain.Dimension(kNSD);
	fStressCorr.Dimension(kNSD);
	fModuliCorr.Dimension(dSymMatrixT::NumValues(kNSD));
	fDevStress.Dimension(kNSD);
	fDevStrain.Dimension(kNSD); 
	fTensorTemp.Dimension(dSymMatrixT::NumValues(kNSD));
	IdentityTensor2.Dimension(kNSD);
	One.Dimension(kNSD);

	/* constant */
	fX_H = 3.0*(fmu+ffriction*fdilation*fkappa) + fH_prime;
	
	/* initialize constant tensor */
	One.Identity();
}

/***********************************************************************
 * Protected
 ***********************************************************************/

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
		if (Flags[ip] == kIsPlastic) /* plastic update */
		{
			/* do not repeat if called again. */
			Flags[ip] = kIsElastic;
			/* NOTE: ComputeOutput writes the updated internal variables
			 *       for output even during iteration output, which is
			 *       called before UpdateHistory */

			/* fetch element data */
			LoadData(element, ip);
	
			/* plastic increment */
			double& dgamma = fInternal[kdgamma];
			//cout << "kdgamma = " << fInternal[kdgamma] << endl;
		
			/* internal state variable */
			fInternal[kalpha] -= fH_prime*dgamma;

			/* dev plastic strain increment	*/
			fPlasticStrain.AddScaled( sqrt32*dgamma, fUnitNorm);
        	
        	/* vol plastic strain increment	*/
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
	if (!element.IsAllocated()) throw ExceptionT::kGeneralFail;

	/* fetch arrays */
	const dArrayT& d_array = element.DoubleData();
	
	/* decode */
	dSymMatrixT::DimensionT dim = dSymMatrixT::int2DimensionT(kNSD);
	int stressdim = dSymMatrixT::NumValues(kNSD);
	int offset    = stressdim*fNumIP;
	int dex       = ip*stressdim;
	
	fPlasticStrain.Alias(         dim, &d_array[           dex]);
	fUnitNorm.Alias(         dim, &d_array[  offset + dex]);     
	fInternal.Alias(kNumInternal, &d_array[2*offset + ip*kNumInternal]);
}

/* returns 1 if the trial elastic strain state lies outside of the 
 * yield surface */
int DPSSLinHardT::PlasticLoading(const dSymMatrixT& trialstrain, 
	ElementCardT& element, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated()) 
		return( YieldCondition(DeviatoricStress(trialstrain,element),
			       MeanStress(trialstrain,element), 0.0) > kYieldTol );
        /* already plastic */
	else 
	{
		/* get flags */
		iArrayT& Flags = element.IntegerData();
		
		/* load internal variables */
		LoadData(element, ip);
		
		fInternal[kftrial] = YieldCondition(DeviatoricStress(trialstrain,element),
							MeanStress(trialstrain,element),fInternal[kalpha]);

		/* plastic */
		if (fInternal[kftrial] > kYieldTol)
		{		
			/* compute unit normal */
			double& norm = fInternal[kstressnorm];

			norm = sqrt(fDevStress.ScalarProduct());
			fUnitNorm.SetToScaled(1.0/norm, fDevStress);
		
			/* set flag */
			Flags[ip] = kIsPlastic;
	
			return 1;
		}
		/* elastic */
		else
		{
			/* set flag */
		    Flags[ip] = kIsElastic; //removed to avoid resetting 7/01
			return 0;
		}
	}
}	

/* Computes the stress corresponding to the given element
 * and elastic strain.  The function returns a reference to the
 * stress in fDevStress */
dSymMatrixT& DPSSLinHardT::DeviatoricStress(const dSymMatrixT& trialstrain,
	const ElementCardT& element)
{
#pragma unused(element)

	/* deviatoric strain */
	fDevStrain.Deviatoric(trialstrain);

	/* compute deviatoric elastic stress */
	fDevStress.SetToScaled(2.0*fmu,fDevStrain);

	return fDevStress;
}

/* computes the hydrostatic (mean) stress */
double DPSSLinHardT::MeanStress(const dSymMatrixT& trialstrain,
	const ElementCardT& element)
{
#pragma unused(element)

  fMeanStress = fkappa*trialstrain.Trace();
  return fMeanStress;
}
