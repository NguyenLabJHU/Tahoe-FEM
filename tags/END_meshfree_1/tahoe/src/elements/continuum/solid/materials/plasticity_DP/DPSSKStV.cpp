/* $Id: DPSSKStV.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: myip (06/01/1999)                                             */

#include "DPSSKStV.h"
#include "ElementCardT.h"
#include "StringT.h"

/* parameters */
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 4;
static const char* Labels[kNumOutput] = {
	"alpha_dev",  // deviatoric part of equivalent plastic strain
	"alpha_vol",  // volumetric part of equivalent plastic strain
	       "VM",  // Von Mises stress
	    "press"}; // pressure

/* constructor */
DPSSKStV::DPSSKStV(ifstreamT& in, const ElasticT& element):
	SSStructMatT(in, element),
	KStV(in),	
	DPSSLinHardT(in, NumIP(), Mu(), Lambda()),
	fStress(3),
	fModulus(dSymMatrixT::NumValues(3)),
	fElasticModulus(dSymMatrixT::NumValues(3))
{
	/* inherited */
	KStV::SetModulus(fElasticModulus);
}

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT DPSSKStV::TangentType(void) const { return GlobalT::kNonSymmetric; }

/* update internal variables */
void DPSSKStV::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Update(element);
}

/* reset internal variables to last converged solution */
void DPSSKStV::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Reset(element);
}

/* print parameters */
void DPSSKStV::Print(ostream& out) const
{
	/* inherited */
	SSStructMatT::Print(out);
	KStV::Print(out);
	DPSSLinHardT::Print(out);
}

/* print name */
void DPSSKStV::PrintName(ostream& out) const
{
	/* inherited */
	SSStructMatT::PrintName(out);
	KStV::PrintName(out);
	DPSSLinHardT::PrintName(out);
}

/* modulus */
const dMatrixT& DPSSKStV::c_ijkl(void)
{
	/* elastoplastic correction */
	fModulus.SumOf(fElasticModulus,	ModuliCorrection(CurrentElement(), CurrIP()));
	return fModulus;
}

/* stress */
const dSymMatrixT& DPSSKStV::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);

	/* elastic stress */
	HookeanStress(fElasticModulus, e_els, fStress);

	/* modify Cauchy stress (return mapping) */
	fStress += StressCorrection(e_els, element, ip);
	return fStress;	
}

/* returns the strain energy density for the specified strain */
double DPSSKStV::StrainEnergyDensity(void)
{
	return HookeanEnergy(fElasticModulus,
		ElasticStrain(e(), CurrentElement(), CurrIP()));
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default. */
int DPSSKStV::NumOutputVariables(void) const  { return kNumOutput; }
void DPSSKStV::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Allocate(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void DPSSKStV::ComputeOutput(dArrayT& output)
{
	/* stress tensor (loads element data) */
	s_ij();

	/* pressure */
	output[3] = fStress.Trace()/3.0;
	
	/* deviatoric Von Mises stress */
	fStress.Deviatoric();
	double J2 = fStress.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[2] = sqrt(3.0*J2);
	
	/* equivalent plastic strains */
	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		output[0] = fInternal[kalpha_dev];
	  	output[1] = fInternal[kalpha_vol];

		/* status flags */
		iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == kIsPlastic) // output with update
		{
			output[0] -= fH_prime*fInternal[kdgamma];
		  	output[1] -= sqrt(3.0)*fdilation*fK_prime*fInternal[kdgamma];
		}
		// alpha not incremented until Update(), which
		// hasn't occurred yet
	}
	else
	{
		output[0] = 0.0;
		output[1] = 0.0;
	}
}

/***********************************************************************
* Protected
***********************************************************************/
