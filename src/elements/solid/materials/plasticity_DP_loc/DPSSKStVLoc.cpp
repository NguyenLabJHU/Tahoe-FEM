/* $Id: DPSSKStVLoc.cpp,v 1.1 2004-03-20 23:35:32 raregue Exp $ */
/* created: myip (06/01/1999) */
#include "DPSSKStVLoc.h"
#include "SSMatSupportT.h"

#include "ElementCardT.h"
#include "StringT.h"
#include "DetCheckT.h"
#include <iostream.h>

using namespace Tahoe;

/* parameters */
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 11;
static const char* Labels[kNumOutput] = {
	    "alpha",  // stress-like internal state variable (isotropic linear hardening)
	    "VM",  // Von Mises stress
	    "press", // pressure
	    "loccheck",
	    "loccheckd", // localization check
	    "n1", // x1 component of normal n for contbif
	    "n2", // x2 component of normal n for contbif
	    "n3", // x3 component of normal n for contbif
	    "nd1", // x1 component of normal n for discbif	
	    "nd2", // x2 component of normal n for discbif	
	    "nd3"}; // x3 component of normal n for discbif	    

/* constructor */
DPSSKStVLoc::DPSSKStVLoc(ifstreamT& in, const SSMatSupportT& support):
	SSSolidMatT(in, support),
	IsotropicT(in),
	HookeanMatT(3),
	DPSSLinHardLocT(in, NumIP(), Mu(), Lambda()),
	fStress(3),
	fModulus(dSymMatrixT::NumValues(3)),
	fModulusdisc(dSymMatrixT::NumValues(3))
{
 
}

/* initialization */
void DPSSKStVLoc::Initialize(void)
{
	/* inherited */
	HookeanMatT::Initialize();
}

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT DPSSKStVLoc::TangentType(void) const { return GlobalT::kNonSymmetric; }

/* update internal variables */
void DPSSKStVLoc::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Update(element);
}

/* reset internal variables to last converged solution */
void DPSSKStVLoc::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Reset(element);
}

/* print parameters */
void DPSSKStVLoc::Print(ostream& out) const
{
	/* inherited */
	SSSolidMatT::Print(out);
	IsotropicT::Print(out);
	DPSSLinHardLocT::Print(out);
}

/* print name */
void DPSSKStVLoc::PrintName(ostream& out) const
{
	/* inherited */
	SSSolidMatT::PrintName(out);
	DPSSLinHardLocT::PrintName(out);
	out << "    Kirchhoff-St.Venant\n";
}

/* modulus */
const dMatrixT& DPSSKStVLoc::c_ijkl(void)
{

	fModulus.SumOf(HookeanMatT::Modulus(),
	ModuliCorrection(CurrentElement(), CurrIP()));
	
	return fModulus;
}

/*discontinuous modulus */
const dMatrixT& DPSSKStVLoc::cdisc_ijkl(void)
{
	/* elastoplastic correction */
	fModulusdisc.SumOf(HookeanMatT::Modulus(),
	ModuliCorrDisc(CurrentElement(), CurrIP()));
	return fModulusdisc;
}

/* stress */
const dSymMatrixT& DPSSKStVLoc::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);

	/* elastic stress */
	HookeanStress(e_els, fStress);

	/* modify Cauchy stress (return mapping) */
	fStress += StressCorrection(e_els, element, ip);
	return fStress;	
}


/*
* Test for localization using "current" values for Cauchy
* stress and the spatial tangent moduli. Returns 1 if the
* determinant of the acoustic tensor is negative and returns
* the normal for which the determinant is minimum. Returns 0
* of the determinant is positive.
*/
int DPSSKStVLoc::IsLocalized(dArrayT& normal)
{
        DetCheckT checker(fStress, fModulus);
        checker.SetfStructuralMatSupport(*fSSMatSupport);

        int loccheck= checker.IsLocalized(normal);
        return loccheck;
}


/* returns the strain energy density for the specified strain */
double DPSSKStVLoc::StrainEnergyDensity(void)
{
	return HookeanEnergy(ElasticStrain(e(), CurrentElement(), CurrIP()));
}

/* returns the number of variables computed for nodal extrapolation
 * during for element output, ie. internal variables. Returns 0
 * by default. */
int DPSSKStVLoc::NumOutputVariables(void) const  { return kNumOutput; } 
void DPSSKStVLoc::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void DPSSKStVLoc::ComputeOutput(dArrayT& output)
{
	
	/* stress tensor (load state) */
	const dSymMatrixT& stress = s_ij();

	/* pressure */
	output[2] = fStress.Trace()/3.0;

	/* deviatoric Von Mises stress */
	fStress.Deviatoric();
	double J2 = fStress.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[1] = sqrt(3.0*J2);
	
	/* stress-like internal variable alpha */
	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		output[0] = fInternal[kalpha];
		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == kIsPlastic)
		  {
			output[0] -= fH_prime*fInternal[kdgamma];
			
			// check for localization
			// compute modulus 
			const dMatrixT& modulus = c_ijkl();

			// continuous localization condition checker
			/*DetCheckT checker(stress, modulus);
			dArrayT normal(stress.Rows());
			output[3] = checker.IsLocalized_SS(normal);
			output[5] = normal[0];
			output[6] = normal[1];
			if (normal.Length() == 3)
				output[7] = normal[2];
			else
				output[7] = 0.0;
			*/
		    output[3] = 0.0;
		    output[5] = 0.0;
			output[6] = 0.0;
			output[7] = 0.0;
	
			/* compute discontinuous bifurcation modulus */
			const dMatrixT& modulusdisc = cdisc_ijkl();

			/* discontinuous localization condition checker */
			DetCheckT checkerdisc(stress, modulusdisc);
			dArrayT normaldisc(stress.Rows());
			output[4] = checkerdisc.IsLocalized_SS(normaldisc);
			output[8] = normaldisc[0];
			output[9] = normaldisc[1];
			if (normaldisc.Length() == 3)
				output[10] = normaldisc[2];
			else
				output[10] = 0.0;
			/*
			output[4] = 0.0;
		    output[8] = 0.0;
			output[9] = 0.0;
			output[10] = 0.0;*/
			
		  }
	}
	else
	{
		output[0] = 0.0;
		output[3] = 0.0;
		output[4] = 0.0;
		output[5] = 0.0;
		output[6] = 0.0;
		output[7] = 0.0;
		output[8] = 0.0;
		output[9] = 0.0;
		output[10] = 0.0;
	}

	
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void DPSSKStVLoc::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}

