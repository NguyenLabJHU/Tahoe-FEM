/* $Id: J2SSKStV1D.cpp,v 1.9 2004-07-15 08:28:12 paklein Exp $ */ 
#include "J2SSKStV1D.h"
#include "SSMatSupportT.h"
#include "ElementCardT.h"
#include "StringT.h"

#include <iostream.h>
#include <math.h>

#include "iArrayT.h"
#include "ifstreamT.h"

/* hardening functions */
#include "CubicSplineT.h"
#include "LinearExponentialT.h"

using namespace Tahoe;

/* parameters */
const double sqrt23 = sqrt(2.0/3.0);
const int    kNumInternal = 4; // number of internal variables
const double kYieldTol    = 1.0e-10;

const int kNSD = 1;

/* element output data */
const int kNumOutput = 4;
static const char* Labels[kNumOutput] = {
	"alpha",   // equivalent plastic strain
	"VM",      // Von Mises stress
	"press",   // pressure
	"dgamma"}; // dgamma

/* constructor */
J2SSKStV1D::J2SSKStV1D(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("J2SSKStV1D"),
	HookeanMatT(kNSD),

	fNumIP(NumIP()),
	fYoung(Young()),
	fK(NULL),
	
	fElasticStrain(kNSD),
	fStressCorr(kNSD),
	fModuliCorr(dSymMatrixT::NumValues(kNSD)),
	fDevStrain(kNSD),
	fTensorTemp(dSymMatrixT::NumValues(kNSD)),
	
	fStress_3D(3),
	
	fStress(kNSD),
	fModulus(dSymMatrixT::NumValues(kNSD))
{


	/* construct hardening function from stream */
	ConstructHardeningFunction(in);
}

J2SSKStV1D::J2SSKStV1D(void):
	ParameterInterfaceT("J2SSKStV1D")
{

}

/* initialization */
void J2SSKStV1D::Initialize(void)
{
ExceptionT::GeneralFail("J2SSKStV1D::Initialize", "out of date");
#if 0
	/* inherited */
	HookeanMatT::Initialize();
#endif
}

/* update internal variables */
void J2SSKStV1D::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Update(element);
}

/* reset internal variables to last converged solution */
void J2SSKStV1D::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Reset(element);
}

#if 0
/* print parameters */
void J2SSKStV1D::Print(ostream& out) const
{
	/* inherited */
	SSSolidMatT::Print(out);
	IsotropicT::Print(out);

	/* hardening function parameters */
	out << " Hardening function:\n";
	fK->Print(out);
}

/* print name */
void J2SSKStV1D::PrintName(ostream& out) const
{
	/* inherited */
	SSSolidMatT::PrintName(out);

	out << "    J2 Isotropic/Kinematic\n";
	out << "    Hardening with Radial Return\n";
	out << "    Small Strain\n";
	out << "    Kirchhoff-St.Venant\n";
}
#endif

/* modulus */
const dMatrixT& J2SSKStV1D::c_ijkl(void)
{
	/* elastoplastic correction */
	fModulus.SumOf(HookeanMatT::Modulus(), ModuliCorrection(CurrentElement(), CurrIP()));	
	return fModulus;
}

/* stress */
const dSymMatrixT& J2SSKStV1D::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);

	/* elastic stress */
	HookeanStress(e_els, fStress);

	/* modify Cauchy stress (return mapping) */
	int iteration = fSSMatSupport->IterationNumber();
	if (iteration > -1) /* elastic iteration */
		fStress += StressCorrection(e_els, element, ip);
	return fStress;	
}

/* returns the strain energy density for the specified strain */
double J2SSKStV1D::StrainEnergyDensity(void)
{
	return HookeanEnergy(ElasticStrain(e(), CurrentElement(), CurrIP()));		
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default. */
int J2SSKStV1D::NumOutputVariables(void) const  { return kNumOutput; }
void J2SSKStV1D::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void J2SSKStV1D::ComputeOutput(dArrayT& output)
{
	/* stress tensor (loads element data and sets fStress) */
	s_ij();

	/* 1D -> 3D */
	fStress_3D = 0.;
	fStress_3D[0] = fStress[0];
	
	/* pressure */
	output[2] = fStress_3D.Trace()/3.0;
	
	/* deviatoric Von Mises stress */
	fStress_3D.Deviatoric();
	double J2 = fStress_3D.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[1] = sqrt(3.0*J2);
	
	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		/* isotropic hardening */
		output[0] = fInternal[kalpha];

		/* status flags */
		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == kIsPlastic) // output with update
			output[0] += fInternal[kdgamma];

		output[3] = fInternal[kdgamma];
	}
	else
	{
		output[0] = 0.0;

		output[3] = 0.0;
	}
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void J2SSKStV1D::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli1D(modulus);
}

/*************************************************************************
* Private
*************************************************************************/

/* returns elastic strain */
const dSymMatrixT& J2SSKStV1D::ElasticStrain(const dSymMatrixT& totalstrain,
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
const dSymMatrixT& J2SSKStV1D::StressCorrection(const dSymMatrixT& trialstrain,
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
		double  alpha   = fInternal[kalpha];

		/* return mapping (single step) */
		if (ftrial > kYieldTol)
		{
			if (fType == kLinear)
			{
				/* plastic increment */
				dgamma = ftrial/(fYoung+dK(alpha));
			}
			else
			{	
				cout << "\n J2SSKStV1D::StressCorrection: unknown hardening function type: " 
				     << fType << endl;
				throw ExceptionT::kBadInputValue;
			}

			/* plastic increment stress correction */
			fStressCorr.SetToScaled(-dgamma*fYoung/fabs(fStress[0]), fStress);
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
const dMatrixT& J2SSKStV1D::ModuliCorrection(const ElementCardT& element,
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
		double alpha = fInternal[kalpha];

		/* moduli corrections */
		fTensorTemp = 1.0;
		fModuliCorr.AddScaled(-fYoung*fYoung/(fYoung+dK(alpha)), fTensorTemp);
	}

	return fModuliCorr;
}	
	 	
/* return a pointer to a new plastic element object constructed with
* the data from element */
void J2SSKStV1D::AllocateElement(ElementCardT& element)
{
	/* determine storage */
	int i_size = 0;
	i_size += fNumIP; //fFlags

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fPlasticStrain
	d_size += kNumInternal*fNumIP;	  //fInternal

	/* construct new plastic element */
	element.Dimension(i_size, d_size);
	
	/* initialize values */
	element.IntegerData() = kIsElastic;
	element.DoubleData()  = 0.0;
}

double J2SSKStV1D::YieldCondition(double alpha) const
{
	 return fabs(fStress[0]) - K(alpha);
}

/* element level data */
void J2SSKStV1D::Update(ElementCardT& element)
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
			fInternal[kalpha] += dgamma;
	
			/* dev plastic strain increment	*/
			fPlasticStrain.AddScaled(dgamma/fabs(fStress[0]), fStress);			
		}
}

/* resets to the last converged solution */
void J2SSKStV1D::Reset(ElementCardT& element)
{
	/* flag not to update again */
	(element.IntegerData()) = kReset;
}

/* load element data for the specified integration point */
void J2SSKStV1D::LoadData(const ElementCardT& element, int ip)
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

	fPlasticStrain.Alias(    dim, &d_array[                     dex]);
	fInternal.Alias(kNumInternal, &d_array[offset + ip*kNumInternal]);     	
}

/* returns 1 if the trial elastic strain state lies outside of the
* yield surface */
int J2SSKStV1D::PlasticLoading(const dSymMatrixT& trialstrain,
	ElementCardT& element, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated())
		return( YieldCondition(0.0) > kYieldTol );
	/* already plastic */
	else
	{
		/* get flags */
		iArrayT& Flags = element.IntegerData();

		/* load internal variables */
		LoadData(element, ip);

		fInternal[kftrial] = YieldCondition(fInternal[kalpha]);

		/* plastic */
		if (fInternal[kftrial] > kYieldTol)
		{		
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

/* construct isotropic hardening function */
void J2SSKStV1D::ConstructHardeningFunction(ifstreamT& in)
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
		default:
		{
			cout << "\n J2SSKStV1D::ConstructHardeningFunction: unknown hardening function type: "
			     << type << endl;
			throw ExceptionT::kBadInputValue;
		}
	}
}
