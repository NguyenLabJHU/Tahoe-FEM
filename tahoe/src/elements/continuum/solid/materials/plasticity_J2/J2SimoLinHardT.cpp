/* $Id: J2SimoLinHardT.cpp,v 1.10 2002-10-20 22:49:05 paklein Exp $ */
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
#include "iArrayT.h"
#include "ElementCardT.h"

/* flags */

using namespace Tahoe;

const int kNumFlags = 2;
const int kEP   = 0;
const int kInit = 1;
	//not used

const double sqrt23    = sqrt(2.0/3.0);
const double kYieldTol = 1.0e-10;

const int kNSD = 3;

/* static variables */
const int J2SimoLinHardT::kNumInternal = 7;

/* constructor */
J2SimoLinHardT::J2SimoLinHardT(ifstreamT& in, int num_ip, double mu):
	J2PrimitiveT(in),
	fNumIP(num_ip),
	fmu(mu),
	fStressCorr(kNSD),
	fModuliCorr(dSymMatrixT::NumValues(kNSD)),
	fRelStress(kNSD),
	f_f_bar(kNSD),
	fb_bar_trial(kNSD),
	fbeta_bar_trial(kNSD),
	fMatrixTemp1(kNSD),
	fMatrixTemp2(kNSD),
	fRed2Temp(kNSD),
	fRed4Temp1(dSymMatrixT::NumValues(kNSD)),
	fRed4Temp2(dSymMatrixT::NumValues(kNSD))
{

}

/** compute trial elastic state - return reference to isochoric,
 * trial elastic stretch */
const dSymMatrixT& J2SimoLinHardT::TrialElasticState(const dMatrixT& F_total,
	const dMatrixT& f_relative, ElementCardT& element, int ip)
{
	/* compute left Cauchy-Green */
	if (element.IsAllocated()) /* element has been plastic */
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
	
		/* isochoric relative deformation gradient */
		f_f_bar.SetToScaled(pow(f_relative.Det(),-1.0/3.0), f_relative);

		/* push stretch forward */
		fb_bar_trial.MultQBQT(f_f_bar, fb_bar);

		/* push kinematic hardening forward */
		fbeta_bar_trial.MultQBQT(f_f_bar, fbeta_bar);
		ftrace_beta_trial = fbeta_bar_trial.Trace();
		fbeta_bar_trial.PlusIdentity(-ftrace_beta_trial/3.0); /* deviatoric part */
		
		/* save */
		fInternal[kDetF_tot] = F_total.Det();
		fb_bar_trial_ = fb_bar_trial;
		fbeta_bar_trial_ = fbeta_bar_trial;
	}
	else /* element is elastic */
	{
		/* trial stretch */
		fb_bar_trial.MultAAT(F_total);
		fb_bar_trial *= pow(F_total.Det(), -2.0/3.0);
		
		/* trial kinematic hardening */
		fbeta_bar_trial = 0.0;
		ftrace_beta_trial = 0.0;
	}

	/* make reduced index */
	return fb_bar_trial;
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
		
		/* set trial state and load data */
		TrialElasticState(F_total, f_relative, element, ip);

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
			double alpha  = fInternal[kalpha];
			double mu_bar = fInternal[kmu_bar];
			double mu_bar_bar = fInternal[kmu_bar_bar];
			
			/* update internal variables */
			dgamma = ftrial/(2.0*mu_bar_bar)/
				(1.0 + (dH(alpha)/3.0/fmu) + (dK(alpha)/3.0/mu_bar_bar));
	
			/* plastic correction */
			fStressCorr.SetToScaled(-2.0*mu_bar_bar*dgamma, fUnitNorm);
			
			/* debugging - check the results of the return mapping */
			bool check_map = false;
			if (check_map)
			{
				/* internal variables */
				double alpha = fInternal[kalpha];
				double mu_bar_bar = fInternal[kmu_bar_bar];
				dSymMatrixT beta_bar = fbeta_bar;

				/* update variables */
				double k = 2.0*mu_bar_bar*dgamma/3.0/fmu;
				alpha += sqrt23*dgamma;
				beta_bar.AddScaled(k*dH(alpha), fUnitNorm);
				
				/* compute relative stress */
				dSymMatrixT relative_stress(3);
				relative_stress.Deviatoric(fb_bar_trial);
				relative_stress *= fmu;
				relative_stress -= fbeta_bar_trial;

				/* return mapping */
				relative_stress += fStressCorr;

				/* compute update yield condition */
				double f = YieldCondition(relative_stress, alpha);
				double t = sqrt(relative_stress.ScalarProduct())/sqrt23;
				cout << "\n J2SimoLinHardT::StressCorrection: check\n"
				     <<   "          f = " << f << '\n'
				     <<   " ||dev[t]|| = " << t << endl;
			}			
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
		double alpha      = fInternal[kalpha];		
		double mu_bar     = fInternal[kmu_bar];
		double mu_bar_bar = fInternal[kmu_bar_bar];
			
		/* scaling factors */
		double f0 = 2.0*mu_bar*dgamma/stressnorm;
		double d0 = 1.0 + dH(alpha)/3.0/fmu + dK(alpha)/3.0/mu_bar_bar;

		double f1 = 1.0/d0 - f0;
		double d1 = 2.0*mu_bar_bar*f1 - (4.0/3.0)*dgamma*((1.0 + dH(alpha)/3.0/fmu)/d0 - 1.0);
		
		double d2 = 2.0*stressnorm*f1;

		/* assemble deviatoric moduli contribution */
		fRed4Temp1.ReducedIndexDeviatoric();
		fModuliCorr.AddScaled(-2.0*mu_bar_bar*f0, fRed4Temp1);

		fRed2Temp.Identity();
		fRed4Temp1.Outer(fUnitNorm, fRed2Temp);
		fRed4Temp1.Symmetrize();
		fModuliCorr.AddScaled(f0*(4.0/3.0)*stressnorm, fRed4Temp1);
		
		fRed4Temp1.Outer(fUnitNorm, fUnitNorm);
		fModuliCorr.AddScaled(-d1, fRed4Temp1);

		fRed2Temp.MultAB(fUnitNorm, fUnitNorm);
		fRed2Temp.Deviatoric();
		fRed4Temp1.Outer(fUnitNorm, fRed2Temp);
		fModuliCorr.AddScaled(-d2, fRed4Temp1);

		/* J factor: Simo's Kirchhoff to spatial tangent modulus */
		//fModuliCorr /= fInternal[kDetF_tot];		
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
	d_size += kNumInternal*fNumIP;                 //fInternal
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fb_bar
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fBeta

	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fUnitNorm
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fb_bar_trial_
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fbeta_bar_trial_

	/* construct new plastic element */
	element.Dimension(i_size, d_size);
	
	/* initialize values */
	element.IntegerData() = kNotInit;
	element.DoubleData()  = 0.0;
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
	if (!element.IsAllocated()) throw ExceptionT::kGeneralFail;

	/* get flags */
	iArrayT& Flags = element.IntegerData();

	/* update plastic variables */
	for (int ip = 0; ip < fNumIP; ip++)
	{
		/* fetch element data */
		LoadData(element,ip);
		
		/* elastic update */
		fb_bar = fb_bar_trial_;
		fbeta_bar = fbeta_bar_trial_;
	
		/* plastic update */
		if (Flags[ip] ==  kIsPlastic)
		{
			/* mark internal state as up to date */
			Flags[ip] = kIsElastic;
			/* NOTE: ComputeOutput writes the updated internal variables
			 *       for output even during iteration output, which is
			 *       called before UpdateHistory */

			/* factors */
			double alpha = fInternal[kalpha];
			double dgamma = fInternal[kdgamma];
			double mu_bar_bar = fInternal[kmu_bar_bar];
			double k = 2.0*mu_bar_bar*dgamma/fmu;
		
			/* internal variables */
			fInternal[kalpha] += sqrt23*dgamma;
			fbeta_bar.AddScaled(k*dH(alpha)/3.0, fUnitNorm);
				
			/* update intermediate configuration */
			fb_bar.AddScaled(-k, fUnitNorm);
		}
	}
}

/* resets to the last converged solution */
void J2SimoLinHardT::Reset(ElementCardT& element)
{
	//disable the material - look at J2QL2DLinHard2DT to
	//verify that this is OK, esp. for successive calls to
	//Update without advancing the simulation
	if (!element.IsAllocated()) throw ExceptionT::kGeneralFail;

	/* get flags */
	iArrayT& Flags = element.IntegerData();

	for (int ip = 0; ip < fNumIP; ip++)
	{
		/* fetch element data */
		LoadData(element, ip);
		
		/* reset loading state */
		Flags[ip] = kIsElastic;

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
	fb_bar.MultAAT(fMatrixTemp2);
	
	/* b_bar - isochoric b */
	fb_bar *= pow(fb_bar.Det(), -1.0/3.0);
	
	/* initialize kinematic hardening */
	fbeta_bar = 0.0;
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

	fb_bar.Set(kNSD, &d_array[dex]);
    fUnitNorm.Set(kNSD, &d_array[offset + dex]);
    fbeta_bar.Set(kNSD, &d_array[2*offset + dex]);
    fb_bar_trial_.Set(kNSD, &d_array[3*offset + dex]);
    fbeta_bar_trial_.Set(kNSD, &d_array[4*offset + dex]);
    fInternal.Set(kNumInternal, &d_array[5*offset + ip*kNumInternal]);     	
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
#pragma unused(F_total)
#pragma unused(f_relative)

	/* compute relative stress */
	fRed2Temp.Deviatoric(fb_bar_trial);
	fRelStress.SetToCombination(fmu, fRed2Temp, -1.0, fbeta_bar_trial);

	/* not yet plastic */
	if (!element.IsAllocated()) 
		return YieldCondition(fRelStress, 0.0) > kYieldTol;
	else /* already plastic */
	{
		/* get flags */
		iArrayT& Flags = element.IntegerData();

		/* initialize intermediate state */
		if (Flags[ip] == kNotInit)
		{
			cout << "\n J2SimoLinHardT::PlasticLoading: should not arrive here\n"
			     <<   "     with uninitialized state" << endl;
			throw ExceptionT::kGeneralFail;
		}
#if 0
		{
			InitIntermediate(F_total, f_relative);
			Flags[ip] = kIsElastic;
		}
#endif
	
		/* set internal variables */		
		fInternal[kftrial]     = YieldCondition(fRelStress, fInternal[kalpha]);
		fInternal[kstressnorm] = sqrt(fRelStress.ScalarProduct());
		fInternal[kmu_bar]     = fmu*fb_bar_trial.Trace()/3.0;
		fInternal[kmu_bar_bar] = fInternal[kmu_bar] - ftrace_beta_trial/3.0;
		
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
		else /* elastic */
		{
			/* set flag */
			Flags[ip] = kIsElastic;
			return 0;
		}
	}
}	
