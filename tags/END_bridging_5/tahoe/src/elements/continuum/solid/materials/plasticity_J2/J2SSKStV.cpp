/* $Id: J2SSKStV.cpp,v 1.9 2004-01-10 04:41:20 paklein Exp $ */
/* created: paklein (06/18/1997) */
#include "J2SSKStV.h"
#include "SSMatSupportT.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* parameters */
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {
	"alpha",  // equivalent plastic strain
	   "VM",  // Von Mises stress
	"press"}; // pressure

/* constructor */
J2SSKStV::J2SSKStV(ifstreamT& in, const SSMatSupportT& support):
	SSSolidMatT(in, support),
	IsotropicT(in),
	HookeanMatT(3),
//	J2SSLinHardT(in, NumIP(), Mu()),
	J2SSC0HardeningT(in, NumIP(), Mu()),
	fStress(3),
	fModulus(dSymMatrixT::NumValues(3))
{

}

/* initialization */
void J2SSKStV::Initialize(void)
{
	/* inherited */
	HookeanMatT::Initialize();
}

/* update internal variables */
void J2SSKStV::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Update(element);
}

/* reset internal variables to last converged solution */
void J2SSKStV::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Reset(element);
}

/* print parameters */
void J2SSKStV::Print(ostream& out) const
{
	/* inherited */
	SSSolidMatT::Print(out);
	IsotropicT::Print(out);
//	J2SSLinHardT::Print(out);
	J2SSC0HardeningT::Print(out);
}

/* print name */
void J2SSKStV::PrintName(ostream& out) const
{
	/* inherited */
	SSSolidMatT::PrintName(out);
//	J2SSLinHardT::PrintName(out);
	J2SSC0HardeningT::PrintName(out);
	out << "    Kirchhoff-St.Venant\n";
}

/* modulus */
const dMatrixT& J2SSKStV::c_ijkl(void)
{
	/* elastoplastic correction */
	fModulus.SumOf(HookeanMatT::Modulus(), ModuliCorrection(CurrentElement(), CurrIP()));	
	return fModulus;
}

/* stress */
const dSymMatrixT& J2SSKStV::s_ij(void)
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
double J2SSKStV::StrainEnergyDensity(void)
{
	return HookeanEnergy(ElasticStrain(e(), CurrentElement(), CurrIP()));		
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default. */
int J2SSKStV::NumOutputVariables(void) const  { return kNumOutput; }
void J2SSKStV::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void J2SSKStV::ComputeOutput(dArrayT& output)
{
	/* stress tensor (loads element data and sets fStress) */
	s_ij();

	/* pressure */
	output[2] = fStress.Trace()/3.0;
	
	/* deviatoric Von Mises stress */
	fStress.Deviatoric();
	double J2 = fStress.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[1] = sqrt(3.0*J2);

	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		/* plastic strain */
		output[0] = fInternal[kalpha];
		
		/* status flags */
		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == kIsPlastic) // output with update
			output[0] += sqrt23*fInternal[kdgamma];
	}
	else
		output[0] = 0.0;
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void J2SSKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}
