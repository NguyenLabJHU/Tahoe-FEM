/* $Id: J2SimoLinHardT.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/19/1997)                                          */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */

#include "J2SimoLinHardT.h"

#include <iostream.h>
#include <math.h>

#include "ElasticT.h"
#include "iArrayT.h"
#include "ElementCardT.h"

/* flags */
const int kNumFlags = 2;
const int kEP   = 0;
const int kInit = 1;
	//not used
	
/* elastic/plastic flag values */
const int kNotInit   =-1;
const int kIsPlastic = 0;
const int kIsElastic = 1;

/* init flag values */
const int kIsNotInit = 0;
const int kIsInit    = 1;
	
/* class constants */
const int	kNumInternal   = 6;
const int	kalpha         = 0; /* isotropic hardening         */
const int	kstressnorm    = 1; /* norm of the relative stress */
const int	kdgamma        = 2; /* consistency parameter       */
const int	kftrial        = 3; /* yield function value        */
const int	kI_bar         = 4; /* trace of b_e_bar            */
const int	kDetF_tot      = 5; /* determinant of total F      */

const double sqrt23        = sqrt(2.0/3.0);
const double kYieldTol     = 1.0e-10;

const int kNSD = 3;

/* constructor */
J2SimoLinHardT::J2SimoLinHardT(ifstreamT& in, int num_ip, double mu):
	J2PrimitiveT(in),
	fNumIP(num_ip),
	fmu(mu),
	fbelastic(kNSD),
	fStressCorr(kNSD),
	fModuliCorr(dSymMatrixT::NumValues(kNSD)),
	fRelStress(kNSD),
	fDev_b(kNSD),
	fUnitNorm(kNSD),
	f_f_bar(kNSD),
	f_b_trial(kNSD),
	fMatrixTemp1(kNSD),
	fMatrixTemp2(kNSD),
	fRed2Temp(kNSD),
	fRed4Temp1(dSymMatrixT::NumValues(kNSD)),
	fRed4Temp2(dSymMatrixT::NumValues(kNSD))
{

}

/* returns isochoric elastic stretch */
const dSymMatrixT& J2SimoLinHardT::ElasticStretch(const dMatrixT& F_total,
	const dMatrixT& f_relative, ElementCardT& element, int ip)
{
	/* compute left Cauchy-Green */
	if (element.IsAllocated()) /* element is plastic */
	{
		/* load internal variables */
		LoadData(element, ip);
	
		/* check intermediate state */
		iArrayT& Flags = element.IntegerData();
		if (Flags[ip] == kNotInit)
		{
			InitIntermediate(F_total, f_relative);
			Flags[ip] = kIsElastic;
		}
	
		/* isochoric */
		f_f_bar.SetToScaled( pow(f_relative.Det(),-1.0/3.0), f_relative);

		fb_bar.ToMatrix(fMatrixTemp1);		
		fMatrixTemp2.MultAB(f_f_bar, fMatrixTemp1);
		f_b_trial.MultABT(fMatrixTemp2, f_f_bar);
	}
	else /* element is elastic */
	{
		f_b_trial.MultABT(F_total, F_total);
		f_b_trial *= pow( F_total.Det(), -2.0/3.0);
	}

	/* make reduced index */
	fbelastic.FromMatrix(f_b_trial);
	
	return fbelastic;
}

/* return the correction to stress vector computed by the mapping the
* stress back to the yield surface, if needed */
const dSymMatrixT& J2SimoLinHardT::StressCorrection(const dMatrixT& F_total,
	const dMatrixT& f_relative, ElementCardT& element, int ip)
{
	/* check consistency and initialize plastic element */
	if (PlasticLoading(F_total, f_relative, element, ip) &&
	    !element.IsAllocated())
	{
		/* new element */
		AllocateElement(element);

		/* initialize element data */
		PlasticLoading(F_total, f_relative, element, ip);
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
			/* constants */
			double I_bar = fInternal[kI_bar];

			/* update internal variables */
			dgamma = (3.0*ftrial)/(6.0*fmu*I_bar + 2.0*fH_bar);
	
			/* plastic correction */
			fStressCorr.SetToScaled(-2.0*fmu*I_bar*dgamma,fUnitNorm);
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
const dMatrixT& J2SimoLinHardT::ModuliCorrection(ElementCardT& element, int ip)
{
	/* initialize */
	fModuliCorr = 0.0;

	if (element.IsAllocated() &&
	   (element.IntegerData())[ip] == kIsPlastic)
	{
		/* load internal variables */
		LoadData(element, ip);

		/* fetch element data */
		double stressnorm = fInternal[kstressnorm];
		double dgamma     = fInternal[kdgamma];		
		double I_bar      = fInternal[kI_bar];
		double mu_bar     = fmu*I_bar;

		/* the betas */
		double beta0 = 1.0 + fH_bar/(3.0*mu_bar);
		double beta1 = 2.0*mu_bar*dgamma/stressnorm;
//note: beta2 = 0.0 for perfect plastic
		double beta2 =(1.0 - 1.0/beta0)*(2.0*stressnorm*dgamma)/(3.0*mu_bar);
//		double beta2 =(1.0 - 1.0/beta0)*(2.0*dgamma)/(3.0*mu_bar);
		double beta3 = 1.0/beta0 - beta1 + beta2;
		double beta4 =(1.0/beta0 - beta1)*(stressnorm/mu_bar);
		
		/* assemble deviatoric moduli contribution */
		fRed4Temp1.ReducedIndexDeviatoric();
		fModuliCorr.AddScaled(-beta1*2.0*mu_bar,fRed4Temp1);

		fRed2Temp.Identity();
		fRed4Temp1.Outer(fUnitNorm,fRed2Temp);
		fRed4Temp1.Symmetrize();
		fModuliCorr.AddScaled(beta1*(4.0/3.0)*stressnorm,fRed4Temp1);
		
		fRed4Temp1.Outer(fUnitNorm,fUnitNorm);
		fModuliCorr.AddScaled(-2.0*mu_bar*beta3,fRed4Temp1);

		fRed2Temp.MultAB(fUnitNorm,fUnitNorm);
		fRed2Temp.Deviatoric();
		fRed4Temp1.Outer(fUnitNorm,fRed2Temp);
		//fRed4Temp1.Symmetrize();
		fModuliCorr.AddScaled(-2.0*mu_bar*beta4,fRed4Temp1);
		//fModuliCorr.AddScaled(-4.0*mu_bar*beta4,fRed4Temp1);

		//make symmetric
		//fModuliCorr.Symmetrize();

		/* J factor */
		fModuliCorr /= fInternal[kDetF_tot];		
	}

	return fModuliCorr;
}	
	 	
/* allocate element storage */
void J2SimoLinHardT::AllocateElement(ElementCardT& element)
{
	/* determine storage */
	int i_size = 0;
	i_size += fNumIP; //flags

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fb_bar
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fUnitNorm
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fBeta
	d_size += kNumInternal*fNumIP;          //fInternal

	/* construct new plastic element */
	element.Allocate(i_size, d_size);
	
	/* initialize values */
	element.IntegerData() = kNotInit;
	element.DoubleData()  = 0.0;
//	for (int ip = 0; ip < fNumInt; ip++)
//	{
//		/* fetch element data */
//		LoadData(newplastic,ip);
//		
//		/* initialize intermediate config */			
//		fb_bar.Identity();
//	}
}

/***********************************************************************
* Private
***********************************************************************/

/* set element for next time step. Usually involves computing
* all the converged updated from the previous timestep */
void J2SimoLinHardT::Update(ElementCardT& element)
{
	//disable the material - look at J2QL2DLinHard2DT to
	//verify that this is OK, esp. for successive calls to
	//Update without advancing the simulation
	if (!element.IsAllocated()) throw eGeneralFail;

	/* get flags */
	iArrayT& Flags = element.IntegerData();

	/* update plastic variables */
	for (int ip = 0; ip < fNumIP; ip++)
	{
		/* fetch element data */
		LoadData(element,ip);
		
		/* elastic update */
		fb_bar.SetToScaled(fInternal[kstressnorm]/fmu,fUnitNorm);
		fb_bar.PlusIdentity(fInternal[kI_bar]);
	
		/* plastic update */
		if (Flags[ip] ==  kIsPlastic)
		{
			/* internal variables */
			fInternal[kalpha] += sqrt23*fInternal[kdgamma];
	
			/* kinematic hardening */
			//fBeta.AddScaled(2.0*(1.0-ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
			//??What's the finite deform version??
							
			fb_bar.AddScaled(-2.0*fInternal[kI_bar]*fInternal[kdgamma],
			                  fUnitNorm);
			
		}
	}
}

/* resets to the last converged solution */
void J2SimoLinHardT::Reset(ElementCardT& element)
{
	//disable the material - look at J2QL2DLinHard2DT to
	//verify that this is OK, esp. for successive calls to
	//Update without advancing the simulation
	if (!element.IsAllocated()) throw eGeneralFail;

	for (int ip = 0; ip < fNumIP; ip++)
	{
		/* fetch element data */
		LoadData(element, ip);

		/* clear plastic increment */		
		fInternal[kdgamma] = 0.0;
	}
}

/* initialize intermediate state from F_n */
void J2SimoLinHardT::InitIntermediate(const dMatrixT& F_total,
	const dMatrixT& f_relative)
{
	/* compute F_n */
	fMatrixTemp1.Inverse(f_relative);
	fMatrixTemp2.MultAB(fMatrixTemp1, F_total);
	
	/* b */
	f_b_trial.MultABT(fMatrixTemp2, fMatrixTemp2);
	
	/* b_bar */
	f_b_trial *= pow(f_b_trial.Det(), -1.0/3.0);
	fb_bar.FromMatrix(f_b_trial);
}

/* load element data */
void J2SimoLinHardT::LoadData(const ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* decode */
	int stressdim = dSymMatrixT::NumValues(kNSD);
	int offset    = stressdim*fNumIP;
	int dex       = ip*stressdim;

	        fb_bar.Set(        kNSD, &d_array[           dex]);
	     fUnitNorm.Set(        kNSD, &d_array[  offset + dex]);
	         fBeta.Set(        kNSD, &d_array[2*offset + dex]);
	     fInternal.Set(kNumInternal, &d_array[3*offset + ip*kNumInternal]);     	
}

/*
* Returns 1 if the trial elastic strain state lies outside of the
* yield surface.
*
* NOTE: pass (element = NULL) if element is not yet plastic, ie. has
*       no stored internal variables.
*/
int J2SimoLinHardT::PlasticLoading(const dMatrixT& F_total,
	const dMatrixT& f_relative, ElementCardT& element, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated())
		return(YieldCondition(RelativeStress(F_total, f_relative, element),
					0.0) > kYieldTol);
	/* already plastic */
	else
	{
		/* get flags */
		iArrayT& Flags = element.IntegerData();

		/* load internal variables */
		LoadData(element,ip);
	
		/* initialize intermediate state */
		if (Flags[ip] == kNotInit)
		{
			InitIntermediate(F_total, f_relative);
			Flags[ip] = kIsElastic;
		}
	
		/* set internal variables */		
		fInternal[kftrial] = YieldCondition(
								RelativeStress(F_total, f_relative, element),
								fInternal[kalpha] );
		fInternal[kstressnorm] = sqrt(fRelStress.ScalarProduct());
		fInternal[kI_bar]      = f_b_trial.Trace()/3.0;		
		fInternal[kDetF_tot]   = F_total.Det();
		
		/* compute unit normal */
		fUnitNorm.SetToScaled(1.0/fInternal[kstressnorm], fRelStress);
		
		/* plastic */
		if (fInternal[kftrial] > kYieldTol)
		{
			/* compute unit normal */
//			double& norm = fInternal[kstressnorm];

//			norm = sqrt(fRelStress.ScalarProduct());
//			fUnitNorm.SetToScaled(1.0/norm, fRelStress);
		
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

/* computes the relative stress and returns a reference to the
* relative stress in fRelStress */
dSymMatrixT& J2SimoLinHardT::RelativeStress(const dMatrixT& F_total,
	const dMatrixT& f_relative, const ElementCardT& element)
{
	/* compute left Cauchy-Green */
	if (element.IsAllocated()) /* element is plastic */
	{
		/* isochoric */
		f_f_bar.SetToScaled( pow(f_relative.Det(),-1.0/3.0), f_relative);

		fb_bar.ToMatrix(fMatrixTemp1);		
		fMatrixTemp2.MultAB(f_f_bar, fMatrixTemp1);
		f_b_trial.MultABT(fMatrixTemp2, f_f_bar);
	}
	else /* element is elastic */
	{
		f_b_trial.MultABT(F_total, F_total);
		f_b_trial *= pow( F_total.Det(), -2.0/3.0);
	}

	/* deviatoric trial stress */
	fDev_b.FromMatrix(f_b_trial);
	fRelStress.SetToScaled(fmu,fDev_b.Deviatoric());
	
	return fRelStress;
}
