/* $Id: J2QL2DLinHardT.cpp,v 1.13 2003-11-21 22:46:48 paklein Exp $ */
/* created: paklein (06/29/1997) */
#include "J2QL2DLinHardT.h"

#include <iostream.h>
#include <math.h>

#include "toolboxConstants.h"

#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* flags */
const int kNumFlags = 2;
const int kEP   = 0;
const int kInit = 1;
	//not used
	
/* elastic/plastic flag values */
const int kE0        =-2; // elastic but intermediate state not initialized
const int kP0        =-1; // hit plastic but intermediate state not initialized
const int kNotInit   = 0;
const int kIsPlastic = 1;
const int kIsElastic = 2;
const int kReset     = 3; // indicate not to repeat update

/* init flag values */
const int kIsNotInit = 0;
const int kIsInit    = 1;
	
/* class constants */
const int	kNumInternal = 5;
const int	kalpha       = 0; /* isotropic hardening         */
const int	kstressnorm  = 1; /* norm of the relative stress */
const int	kdgamma      = 2; /* consistency parameter       */
const int	kftrial      = 3; /* yield function value        */
const int	kDetF_tot    = 4; /* determinant of total F      */

const double sqrt23      = sqrt(2.0/3.0);
const double kYieldTol   = 1.0e-10;

const int kNSD = 3;

/* element output data */
const int kNumOutput = 5;
static const char* Labels[kNumOutput] = {
        "alpha",  // equivalent plastic strain
	 "VM_Kirch",  // Kirchhoff Von Mises stress
	    "press",  // pressure
	    "s_max",  // max in-plane principal stress
	    "s_min"}; // min in-plane principal stress

/* constructor */
J2QL2DLinHardT::J2QL2DLinHardT(ifstreamT& in, const FSMatSupportT& support):
	QuadLog2D(in, support),
	J2PrimitiveT(in),
	fb_elastic(kNSD),
	fEPModuli(kNSD),

	/* work space */
	fa_inverse(kNSD),
	fMatrixTemp1(kNSD),
	fMatrixTemp2(kNSD),
	fMatrixTemp3(kNSD),
	fdev_beta(kNSD),
	fFtot(3),
	ffrel(3),
	
	/* 2D translation */
	fF_temp(2),
	fFtot_2D(2),
	ffrel_2D(2)	
{
	/* for intermediate config update */
	fa_inverse.Inverse(fEigMod);
}

/* update internal variables */
void J2QL2DLinHardT::UpdateHistory(void)
{
	/* get flags */
	ElementCardT& element = CurrentElement();	
	iArrayT& Flags = element.IntegerData();
	
	/* update plastic variables */
	for (int ip = 0; ip < NumIP(); ip++)
		if (Flags[ip] != kReset && Flags[ip] != kNotInit)
		//if (Flags[ip] == kIsPlastic || Flags[ip] == kP0)
		{
			/* do not repeat if called again. Also, needed
			 * so ComputeOutput does not repeat the update */
			Flags[ip] = kReset;
						
			/* fetch element data */
			LoadData(element,ip);
	
			/* plastic increment */
			fInternal[kalpha] += sqrt23*fInternal[kdgamma];
			
			/* principal values - plane strain */
			fb_2D.ReduceFrom3D(fb_tr);
			fb_2D.PrincipalValues(fEigs);
			fEigs[2] = fb_tr(2,2); //only out-of-plane value
	
			/* set spectral decomposition (don't perturb) */
			fSpectral.SpectralDecomp(fb_tr, fEigs, false);
					
			/* updated log stretches */
			fa_inverse.Multx(fbeta_tr, floge);
					
			/* compute intermediate config */
			fb_n.SetToCombination(exp(2.0*floge[0]), fSpectral.Rank1_Principal(0),
			                      exp(2.0*floge[1]), fSpectral.Rank1_Principal(1),
			                      exp(2.0*floge[2]), fSpectral.Rank1_Principal(2));
		}
}

/* reset internal variables to last converged solution */
void J2QL2DLinHardT::ResetHistory(void)
{
	/* flag not to update again */
	//(element->IntegerData()) = kReset;
	//need to check if initialized!
	
	/* get flags */
	ElementCardT& element = CurrentElement();
	iArrayT& Flags = element.IntegerData();

	for (int i = 0; i < Flags.Length(); i++)
		if (Flags[i] == kIsElastic || Flags[i] == kIsPlastic)
			Flags[i] = kReset;
}

/* print parameters */
void J2QL2DLinHardT::Print(ostream& out) const
{
	/* inherited */
	QuadLog2D::Print(out);
	J2PrimitiveT::Print(out);
}

/* modulus */
const dMatrixT& J2QL2DLinHardT::c_ijkl(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();

	/* get trial stretch */
	int ip = CurrIP();
	const dSymMatrixT& b_tr = TrialStretch(fFtot, ffrel, ip);

	/* principal values - plane strain */
	fb_2D.ReduceFrom3D(b_tr);
	fb_2D.PrincipalValues(fEigs);
	fEigs[2] = b_tr(2,2); //only out-of-plane value

	/* treat undeformed state separately */
	if ( fabs(fEigs[0] - 1.0) < kSmall &&
	     fabs(fEigs[1] - 1.0) < kSmall &&
	     fabs(fEigs[2] - 1.0) < kSmall ) //now explicitly check 3rd dim
	{
		IsotropicT::ComputeModuli2D(fModulus2D, fConstraintOption);
	}
	/* compute moduli */
	else
	{
		/* full spectral decomposition (and perturb) */
		fSpectral.DecompAndModPrep(b_tr, true);

		/* logarithmic stretches */
		fEigs = fSpectral.Eigenvalues();
		LogStretches(fEigs);

		/* principal stresses */
		fEigMod.Multx(floge,fBeta);
	
		/* Compute elastoplastic moduli and fetch principal stresses */
		fEPModuli = fEigMod;
		ElastoPlasticCorrection(fEPModuli, fBeta, ip);

		/* stiffness part */
		fModulus = fSpectral.EigsToRank4(fEPModuli);
		
		/* stress part */
		fModulus.AddScaled(2.0*fBeta[0], fSpectral.SpatialTensor(b_tr, 0));
		fModulus.AddScaled(2.0*fBeta[1], fSpectral.SpatialTensor(b_tr, 1));
		fModulus.AddScaled(2.0*fBeta[2], fSpectral.SpatialTensor(b_tr, 2));

		/* factor of J */
		fModulus /= sqrt( fEigs[0]*fEigs[1]*fEigs[2] );

		/* 3D -> 2D */
		fModulus2D.Rank4ReduceFrom3D(fModulus);
	}

	fModulus2D *= fThickness;
	return fModulus2D;
}
	
/* stress */
const dSymMatrixT& J2QL2DLinHardT::s_ij(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();

	/* get trial stretch */
	int ip = CurrIP();
	const dSymMatrixT& b_tr = TrialStretch(fFtot, ffrel, ip);

	/* spectral decomposition (don't perturb) */
	fSpectral.SpectralDecomp_new(b_tr, false);

	/* logarithmic stretches */
	fEigs = fSpectral.Eigenvalues();
	LogStretches(fEigs);
	
	/* principal stresses */
	fEigMod.Multx(floge, fBeta);

	/* plastic correction */
	ReturnMapping(b_tr, fBeta, ip);

	/* Kirchhoff -> Cauchy */
	fBeta /= sqrt(fEigs[0]*fEigs[1]*fEigs[2]);

	/* QuadLog3D stress */
	fStress = fSpectral.EigsToRank2(fBeta);

	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(fStress);
	fStress2D *= fThickness;

	return fStress2D;
}

/* strain energy density */
double J2QL2DLinHardT::StrainEnergyDensity(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();
	
	const dSymMatrixT& b_tr = TrialStretch(fFtot, ffrel, CurrIP());
	
	/* principal values - plane strain */
	fb_2D.ReduceFrom3D(b_tr);
	fb_2D.PrincipalValues(fEigs);
	fEigs[2] = b_tr(2,2); //only out-of-plane value
	
	/* logarithmic stretches */
	LogStretches(fEigs);

	return fThickness*ComputeEnergy(floge);
}

/*
* Returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default.
*/
int J2QL2DLinHardT::NumOutputVariables(void) const { return kNumOutput; }
void J2QL2DLinHardT::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void J2QL2DLinHardT::ComputeOutput(dArrayT& output)
{
	/* compute Cauchy stress (sets fBeta) */
	s_ij();
	output[2] = fStress.Trace()/3.0;

	/* in-plane principal values (Cauchy) */
	if (fBeta[0] > fBeta[1])
	{
		output[3] = fBeta[0];
		output[4] = fBeta[1];
	}
	else
	{
		output[3] = fBeta[1];
		output[4] = fBeta[0];
	}

	/* Cauchy -> Kirchhoff */
	fBeta *= sqrt( fEigs[0]*fEigs[1]*fEigs[2] );
	
	/* Von Mises equivalent (Kirchhoff) stress */
	fBeta    -= fBeta.Average(); //deviatoric part
	output[1] = fBeta.Magnitude()/sqrt23;

	/* plastic evolution parameter */
	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		output[0] = fInternal[kalpha];

		/* status flags */
		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == kIsPlastic) // output with update
			output[0] += sqrt23*fInternal[kdgamma];
	}
	else
		output[0] = 0.0;
}

/***********************************************************************
* Protected
***********************************************************************/

/* returns the elastic stretch */
const dSymMatrixT& J2QL2DLinHardT::TrialStretch(const dMatrixT& F_total,
	const dMatrixT& f_relative, int ip)
{
	/* element pointer */
	ElementCardT& element = CurrentElement();
	
	/* compute left Cauchy-Green */
	if (element.IsAllocated()) /* element is plastic */
	{
		/* load internal variables */
		LoadData(element, ip);
	
		/* check intermediate state */
		iArrayT& Flags = element.IntegerData();
		if (Flags[ip] <= kNotInit)
		{
			/* initialize intermediate config */
			InitIntermediate(F_total, f_relative);
			
			/* reset EP flag */
			Flags[ip] = (Flags[ip] == kP0) ? kIsPlastic : kIsElastic;
		}
	
		/* compute (and store) trial stretch */
		fb_tr.MultQBQT(f_relative, fb_n);

		return fb_tr;
	}
	else /* element is elastic */
	{
		fb_elastic.MultAAT(F_total);
	
		return fb_elastic;
	}
}

/*
* Return the correction to principal stress vector computed by the
* mapping the stress back to the yield surface, if needed, returning
* 1 of a plastic correction is needed (even possible and valid)
*/
void J2QL2DLinHardT::ReturnMapping(const dSymMatrixT& b_tr, dArrayT& beta, int ip)
{
	/* deviatoric part of principal stresses */
	fdev_beta = beta;
	fdev_beta -= fdev_beta.Average();

	/* element pointer */
	ElementCardT& element = CurrentElement();

	/* check consistency and initialize plastic element */
	if (PlasticLoading(fdev_beta, element, ip) &&
	    !element.IsAllocated())
	{
		/* new element */
		AllocateElement(element);

		/* set element data */
		PlasticLoading(fdev_beta, element, ip);

		/* initialize trial stretch */
		fb_tr = b_tr;
	}

	/* for plastic elements */
	if (element.IsAllocated())
	{				
		/* fetch data */
		double  ftrial = fInternal[kftrial];
		double& dgamma = fInternal[kdgamma];
		
		/* return mapping (single step) - set plastic increment */
		if (ftrial > kYieldTol)
		{	
			/* update internal variables */
			dgamma = (3.0*ftrial)/(6.0*Mu() + 2.0*fH_bar);
		
			/* corrected principal stresses */
			beta.AddScaled(-2.0*Mu()*dgamma,fUnitNorm);
		}
		else
			dgamma = 0.0;
			
		/* store */
		fbeta_tr = beta;
	}
}	

/*
* Return the correction to moduli due to plasticity (if any)
*
* Note: Return mapping occurs during the call to StressCorrection.
*       The element passed in is already assumed to carry current
*       internal variable values.
*/
void J2QL2DLinHardT::ElastoPlasticCorrection(dSymMatrixT& a_ep, dArrayT& beta,
	int ip)
{
	/* element pointer */
	const ElementCardT& element = CurrentElement();

	if (element.IsAllocated() &&
	    (element.IntegerData())[ip] == kIsPlastic)
	{
		/* load internal variables */
		LoadData(element, ip);
		
		/* copy beta */
		beta = fbeta_tr;
		
		/* correct elastoplastic moduli */
		if ((element.IntegerData())[ip] == kIsPlastic)
		{
			/* compute constants */
			double alpha    = fInternal[kalpha];
			double thetahat = 2.0*Mu()*fInternal[kdgamma]/
			                          fInternal[kstressnorm];
			double thetabar = (3.0/(3.0 + (dK(alpha) + dH(alpha))/Mu())) - thetahat;
			
			/* moduli corrections */
			fMatrixTemp3.Outer(fUnitNorm);
			a_ep.AddScaled(-2.0*Mu()*thetabar, fMatrixTemp3);
			a_ep.AddScaled(-2.0*Mu()*thetahat, fDevOp3);
		}					
	}
}	
	 	
/*
* Return a pointer to a new plastic element object constructed with
* the data from element
*/
void J2QL2DLinHardT::AllocateElement(ElementCardT& element)
{
	/* number of integration points */
	int numint = NumIP();

	/* determine storage */
	int i_size = 0;
	i_size += numint; //flags

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*numint; //fb_n
	d_size += dSymMatrixT::NumValues(kNSD)*numint; //fb_tr
	d_size += kNSD*numint;                  //fbeta_tr
	d_size += kNSD*numint;                  //fUnitNorm - principal space
	d_size += kNumInternal*numint;          //fInternal

	/* construct new plastic element */
	element.Dimension(i_size, d_size);

	/* initialize values */
	element.IntegerData() = kNotInit;
	element.DoubleData()  = 0.0;
}

/***********************************************************************
* Private
***********************************************************************/

/* compute F_total and f_relative */
void J2QL2DLinHardT::ComputeGradients(void)
{
	/* compute relative displacement */
	fFtot_2D = F();
	fF_temp.Inverse(F_total_last());
	ffrel_2D.MultAB(fFtot_2D,fF_temp);

	/* 2D -> 3D */
	fFtot.Rank2ExpandFrom2D(fFtot_2D);
	fFtot(2,2) = 1.0;

	ffrel.Rank2ExpandFrom2D(ffrel_2D);
	ffrel(2,2) = 1.0;	
}

/* initialize intermediate state from F_n (for ) */
void J2QL2DLinHardT::InitIntermediate(const dMatrixT& F_total,
	const dMatrixT& f_relative)
{
	/* compute F_n */
	fMatrixTemp1.Inverse(f_relative);
	fMatrixTemp2.MultAB(fMatrixTemp1, F_total);
	
	/* compute (and store) b */
	fb_n.MultAAT(fMatrixTemp2);
}

/* load element data for the specified integration point */
void J2QL2DLinHardT::LoadData(const ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	const dArrayT& d_array = element.DoubleData();

	/* decode */
	int stressdim = dSymMatrixT::NumValues(kNSD);
	dSymMatrixT::DimensionT dim = dSymMatrixT::int2DimensionT(kNSD);
	int blocksize = stressdim + stressdim + kNSD + kNSD + kNumInternal;
	int dex = ip*blocksize;
	
	     fb_n.Alias(        dim, &d_array[dex             ]);
	    fb_tr.Alias(        dim, &d_array[dex += stressdim]);
	 fbeta_tr.Alias(        kNSD, &d_array[dex += stressdim]);
	fUnitNorm.Alias(        kNSD, &d_array[dex += kNSD]);
	fInternal.Alias(kNumInternal, &d_array[dex += kNSD     ]);     	
}

/*
* Returns 1 if the trial elastic strain state lies outside of the
* yield surface.
*
* NOTE: pass (element = NULL) if element is not yet plastic, ie. has
*       no stored internal variables.
*/
int J2QL2DLinHardT::PlasticLoading(const dArrayT& dev_beta,
	ElementCardT& element, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated())
		return YieldCondition(dev_beta, 0.0) > kYieldTol;
	/* already plastic */
	else
	{
		/* get flags */
		iArrayT& Flags = element.IntegerData();

		/* load internal variables */
		LoadData(element,ip); // load in TrialStretch only ??
		
		/* set internal variables */		
		fInternal[kftrial] = YieldCondition(dev_beta, fInternal[kalpha]);
		
		/* plastic */
		if (fInternal[kftrial] > kYieldTol)
		{
			/* compute unit normal */
			double& norm = fInternal[kstressnorm];

			norm = dev_beta.Magnitude();
			fUnitNorm.SetToScaled(1.0/norm, dev_beta);
		
			/* set flag */
			Flags[ip] = (Flags[ip] == kNotInit) ? kP0 : kIsPlastic;
	
			return 1;
		}
		/* elastic */
		else
		{
			/* set flag */
			Flags[ip] = (Flags[ip] == kNotInit) ? kE0 : kIsElastic;
			
			return 0;
		}
	}
}

double J2QL2DLinHardT::YieldCondition(const dArrayT& devpstress,
	double alpha) const
{
	return devpstress.Magnitude() - sqrt23*K(alpha);
}
