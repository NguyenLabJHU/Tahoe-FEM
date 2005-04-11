/* $Id: J2SSKStV.cpp,v 1.9.28.3 2005-04-11 19:40:37 thao Exp $ */
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

	/*	if (CurrElementNumber() == 0 ) {
	  cout << "\nIP: "<<CurrIP();
	  cout << "\nStress: "<<fStress;
	  if (element.IsAllocated()) {
	    cout << "\n: Kalpha: "<<K(fInternal[kalpha]);
	    cout << "\nYield: "<<YieldCondition(fStress, fInternal[kalpha]);
	  }
	  }*/
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

const iArrayT& J2SSKStV::InternalDOF(void) const
{
  return(fInternalDOF);
}

const dArrayT& J2SSKStV::InternalStrainVars(void)
{
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) {
	        s_ij();
		double* p = fInternalStrainVars.Pointer();

		/* plastic increment */
      		double& dgamma = fInternal[kdgamma];

		*p++ = fInternal[kalpha] + sqrt23*dgamma;

		*p++ = fPlasticStrain[0] + dgamma*fUnitNorm[0];
		*p++ = fPlasticStrain[1] + dgamma*fUnitNorm[1];
		*p++ = fPlasticStrain[2] + dgamma*fUnitNorm[2];
		*p++ = fPlasticStrain[3] + dgamma*fUnitNorm[3];
		*p++ = fPlasticStrain[4] + dgamma*fUnitNorm[4];
		*p++ = fPlasticStrain[5] + dgamma*fUnitNorm[5];

		*p++ = fPlasticStrain[0] + dgamma*fUnitNorm[0];
		*p++ = fPlasticStrain[1] + dgamma*fUnitNorm[1];
		*p++ = fPlasticStrain[2] + dgamma*fUnitNorm[2];
		*p++ = fPlasticStrain[3] + dgamma*fUnitNorm[3];
		*p++ = fPlasticStrain[4] + dgamma*fUnitNorm[4];
		*p++ = fPlasticStrain[5] + dgamma*fUnitNorm[5];

		/*               const iArrayT& flags = element.IntegerData();
      		if (flags[CurrIP()] == kIsPlastic && CurrElementNumber() == 0) {
		  cout << "\nElement: "<<CurrElementNumber()<< "\t IP: "<<CurrIP();
		cout << "\n:Internal Strains: "<<fInternalStrainVars;
		}*/
	}
        else fInternalStrainVars = 0.0;
	return(fInternalStrainVars);
}

const dArrayT& J2SSKStV::InternalStressVars(void)
{
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) {
	        s_ij();
	        double* p = fInternalStressVars.Pointer();
		double& dgamma = fInternal[kdgamma];
		double alpha = fInternal[kalpha];
		alpha += sqrt23*dgamma;

		*p++ = -(K(alpha)-fYield);

		*p++ = -fBeta[0]; //-(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
		*p++ = -fBeta[1]; //-(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
		*p++ = -fBeta[2]; //-(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
		*p++ = -fBeta[3]; //-(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
		*p++ = -fBeta[4]; //-(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
		*p++ = -fBeta[5]; //-(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);

		*p++ = fStress[0];
		*p++ = fStress[1];
		*p++ = fStress[2];
		*p++ = fStress[3];
		*p++ = fStress[4];
		*p++ = fStress[5];

		/*		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == kIsPlastic && CurrElementNumber() == 0) {
		  cout << "\nElement: "<<CurrElementNumber()<< "\t IP: "<<CurrIP();
		  cout << "\n:Internal Stress: "<<fInternalStressVars;
		  }*/
	}
        else fInternalStressVars = 0.0;
	return(fInternalStressVars);
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
