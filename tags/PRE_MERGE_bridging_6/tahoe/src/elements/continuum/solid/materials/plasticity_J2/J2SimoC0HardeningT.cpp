/* $Id: J2SimoC0HardeningT.cpp,v 1.12 2004-01-27 19:11:40 paklein Exp $ */
/* created: paklein (05/01/2001) */
#include "J2SimoC0HardeningT.h"

#include <iostream.h>
#include <math.h>

#include "ifstreamT.h"
#include "iArrayT.h"
#include "ElementCardT.h"

/* hardening functions */
#include "CubicSplineT.h"
#include "LinearExponentialT.h"
#include "PowerLawT.h"

using namespace Tahoe;

const double sqrt23    = sqrt(2.0/3.0);
const double kYieldTol = 1.0e-10;

const int kNSD = 3;

/* static variables */
const int J2SimoC0HardeningT::kNumInternal = 8;

/* constructor */
J2SimoC0HardeningT::J2SimoC0HardeningT(ifstreamT& in, int num_ip, double mu):
	fNumIP(num_ip),
	fmu(mu),
	fK(NULL),
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
	/* construct hardening function from stream */
	ConstructHardeningFunction(in);
}

/** destructor */
J2SimoC0HardeningT::~J2SimoC0HardeningT(void) { delete fK; }

/***********************************************************************
* Protected
***********************************************************************/

/* write parameters */
void J2SimoC0HardeningT::Print(ostream& out) const
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

void J2SimoC0HardeningT::PrintName(ostream& out) const
{
	out << "    J2 Isotropic/Kinematic\n";
	out << "    Hardening with Radial Return\n";
}

/* compute trial elastic state - return reference to isochoric,
 * trial elastic stretch */
const dSymMatrixT& J2SimoC0HardeningT::TrialElasticState(const dMatrixT& F_mechanical,
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
			InitIntermediate(F_mechanical, f_relative);
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
		fInternal[kDetF_tot] = F_mechanical.Det();
		fb_bar_trial_ = fb_bar_trial;
		fbeta_bar_trial_ = fbeta_bar_trial;
	}
	else /* element is elastic */
	{
		/* trial stretch */
		fb_bar_trial.MultAAT(F_mechanical);
		fb_bar_trial *= pow(F_mechanical.Det(), -2.0/3.0);
		
		/* trial kinematic hardening */
		fbeta_bar_trial = 0.0;
		ftrace_beta_trial = 0.0;
	}

	/* make reduced index */
	return fb_bar_trial;
}

/* determine elastic or plastic loading for the current step */
int J2SimoC0HardeningT::PlasticLoading(ElementCardT& element, int ip)
{
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

		/* should not get here uninitialized */
		if (Flags[ip] == kNotInit)
		{
			cout << "\n J2SimoC0HardeningT::PlasticLoading: should not arrive here\n"
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
		fInternal[kHeatIncr]   = 0.0;
		
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

/* return the correction to stress vector computed by the mapping the
 * stress back to the yield surface, if needed */
const dSymMatrixT& J2SimoC0HardeningT::StressCorrection(ElementCardT& element, int ip)
{
#pragma unused(ip)

	/* initialize */
	fStressCorr = 0.0;

	if (element.IsAllocated())
	{				
		/* fetch data */
		double  ftrial = fInternal[kftrial];
		double& dgamma = fInternal[kdgamma];
		double& heat_incr = fInternal[kHeatIncr];
		
		/* return mapping (single step) */
		if (ftrial > kYieldTol)
		{	
			/* constants */
			double alpha  = fInternal[kalpha];
			double& mu_bar = fInternal[kmu_bar];
			double& mu_bar_bar = fInternal[kmu_bar_bar];
			
			/* single step */
			if (fType == kLinear)
			{
				/* update internal variables */
				dgamma = ftrial/(2.0*mu_bar_bar)/
					(1.0 + (dH(alpha)/3.0/fmu) + (dK(alpha)/3.0/mu_bar_bar));
			}
			else /* local Newton iteration */
			{
				double x_tr = ftrial + sqrt23*K(alpha);
				double f_hat =-ftrial;
				double k = (1.0 + (dH(alpha)/3.0/fmu))*2.0*mu_bar_bar; // no kinematic hardening
				
				dgamma = 0.0;
				int max_iteration = 15;
				int count = 0;
				while (fabs(f_hat) > kYieldTol && ++count <= max_iteration)
				{
					/* stiffness */
					double df_hat = 2.0*dK(alpha + sqrt23*dgamma)/3.0 + k;
					if (df_hat < kSmall)
					{
						cout << "\n J2SimoC0HardeningT::StressCorrection: consistency function is nonconvex" << endl;
						throw ExceptionT::kGeneralFail;
					}
				
					/* increment update */
					dgamma -= f_hat/df_hat;
					
					/* update condition */
					f_hat = sqrt23*K(alpha + sqrt23*dgamma) - x_tr + k*dgamma;
				}
				
				/* check for failure */
				if (count == max_iteration)
				{
					cout << "\n J2SimoC0HardeningT::StressCorrection: local iteration failed after " 
					     << max_iteration << " iterations" << endl;
					throw ExceptionT::kGeneralFail;
				}
			}

			/* plastic correction - to Cauchy stress */
			fStressCorr.SetToScaled(-2.0*mu_bar_bar*dgamma/fInternal[kDetF_tot], fUnitNorm);

			/* incremental heat generation - 90% of plastic work */
			heat_incr = 0.9*sqrt23*dgamma*fInternal[kstressnorm];
				
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

				/* return mapping - Kirchhoff stress */
				relative_stress.AddScaled(fInternal[kDetF_tot], fStressCorr);

				/* compute update yield condition */
				double f = YieldCondition(relative_stress, alpha);
				double t = sqrt(relative_stress.ScalarProduct())/sqrt23;
				cout << "\n J2SimoC0HardeningT::StressCorrection: check\n"
				     <<   "          f = " << f << '\n'
				     <<   " ||dev[t]|| = " << t << endl;
			}			
		}
		else {
			dgamma = 0.0;
			heat_incr = 0.0;
		}
	}
		
	return fStressCorr;
}	

/* return the correction to moduli due to plasticity (if any)
*
* Note: Return mapping occurs during the call to StressCorrection.
*       The element passed in is already assumed to carry current
*       internal variable values */
const dMatrixT& J2SimoC0HardeningT::ModuliCorrection(ElementCardT& element, int ip)
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
		fModuliCorr /= fInternal[kDetF_tot];		
	}

	return fModuliCorr;
}	
	 	
/* allocate element storage */
void J2SimoC0HardeningT::AllocateElement(ElementCardT& element)
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
void J2SimoC0HardeningT::Update(ElementCardT& element)
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
void J2SimoC0HardeningT::Reset(ElementCardT& element)
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
void J2SimoC0HardeningT::InitIntermediate(const dMatrixT& F_mechanical,
	const dMatrixT& f_relative)
{
	/* compute F_n */
	fMatrixTemp1.Inverse(f_relative);
	fMatrixTemp2.MultAB(fMatrixTemp1, F_mechanical);
	
	/* b */
	fb_bar.MultAAT(fMatrixTemp2);
	
	/* b_bar - isochoric b */
	fb_bar *= pow(fb_bar.Det(), -1.0/3.0);
	
	/* initialize kinematic hardening */
	fbeta_bar = 0.0;
}

/* load element data */
void J2SimoC0HardeningT::LoadData(const ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	const dArrayT& d_array = element.DoubleData();

	/* decode */
	int stressdim = dSymMatrixT::NumValues(kNSD);
	dSymMatrixT::DimensionT dim = dSymMatrixT::int2DimensionT(kNSD);
	int offset = stressdim*fNumIP;
	int dex = ip*stressdim;
	
	/* already set */
	if (fb_bar.Pointer() == d_array.Pointer(dex)) 
		return;
	else
	{
		/* set pointers */
		fb_bar.Alias(dim, &d_array[dex]);
		fUnitNorm.Alias(dim, &d_array[offset + dex]);
		fbeta_bar.Alias(dim, &d_array[2*offset + dex]);
		fb_bar_trial_.Alias(dim, &d_array[3*offset + dex]);
		fbeta_bar_trial_.Alias(dim, &d_array[4*offset + dex]);
		fInternal.Alias(kNumInternal, &d_array[5*offset + ip*kNumInternal]);
	}
}

/* returns the value value of the yield function given the
 * relative stress vector and state variables, where  alpha
 * represents isotropic hardening.  NOTE: the relative stress
 * should already contain the correction for any kinematic
 * hardening. */
double J2SimoC0HardeningT::YieldCondition(const dSymMatrixT& stress, 
	double alpha) const
{
	return sqrt(stress.ScalarProduct()) - sqrt23*K(alpha);
}

/* construct isotropic hardening function */
void J2SimoC0HardeningT::ConstructHardeningFunction(ifstreamT& in)
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
				cout << "\n J2SimoC0HardeningT::ConstructHardeningFunction: expecting at least 2 spline points:"
				     << num_points << endl;
				throw ExceptionT::kBadInputValue;
			}
			dArray2DT points(num_points, 2);
			in >> points;
			
			/* construct spline */
			fK = new CubicSplineT(points, CubicSplineT::kFreeRun);
			break;
		}
		case kPowerLaw:
		{
			fType = kPowerLaw;
			double a, b, n;
			a = b = n = 0.0;
			in >> a >> b >> n;
			
			/* construct function */
			fK = new PowerLawT(a, b, n);
			break;
		}
		default:
			cout << "\n J2SimoC0HardeningT::ConstructHardeningFunction: unknown hardening function type: " 
			     << type << endl;
			throw ExceptionT::kBadInputValue;
	}
}
