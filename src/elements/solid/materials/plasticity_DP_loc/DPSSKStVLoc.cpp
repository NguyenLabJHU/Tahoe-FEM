/* $Id: DPSSKStVLoc.cpp,v 1.4 2004-07-15 08:28:56 paklein Exp $ */
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
const int kNumOutput = 4;
static const char* Labels[kNumOutput] = {
	"alpha",  // stress-like internal state variable (isotropic linear hardening)
	"VM",  // Von Mises stress
	"press", // pressure
	"loccheck"}; // localization check   
	
// need to store the normals somewhere.  as ISVs?

/* constructor */
DPSSKStVLoc::DPSSKStVLoc(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("DPSSKStVLoc"),
//	SSSolidMatT(in, support),
	IsotropicT(in),
	HookeanMatT(3),
	DPSSLinHardLocT(in, NumIP(), Mu(), Lambda()),
	fStress(3),
	fModulus(dSymMatrixT::NumValues(3)),
	fModulusPerfPlas(dSymMatrixT::NumValues(3))
{
 
}

/* initialization */
void DPSSKStVLoc::Initialize(void)
{
ExceptionT::GeneralFail("DPSSKStVLoc::Initialize", "out of date");
#if 0
	/* inherited */
	HookeanMatT::Initialize();
#endif
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

/* modulus */
const dMatrixT& DPSSKStVLoc::c_ijkl(void)
{
	fModulus.SumOf(HookeanMatT::Modulus(),
	ModuliCorrection(CurrentElement(), CurrIP()));

	return fModulus;
}

/*discontinuous modulus */
const dMatrixT& DPSSKStVLoc::c_perfplas_ijkl(void)
{
	/* elastoplastic correction */
	fModulusPerfPlas.SumOf(HookeanMatT::Modulus(),
	ModuliCorrPerfPlas(CurrentElement(), CurrIP()));
	return fModulusPerfPlas;
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
// not used
// see ComputeOutput
int DPSSKStVLoc::IsLocalized(dArrayT& normal)
{
	DetCheckT checker(fStress, fModulus, fModulusCe);
	checker.SetfStructuralMatSupport(*fSSMatSupport);

	int loccheck = checker.IsLocalized(normal);
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
	for (int i = 0; i < kNumOutput; i++) labels[i] = Labels[i];
}

void DPSSKStVLoc::ComputeOutput(dArrayT& output)
{
	dMatrixT Ce = HookeanMatT::Modulus();
	
	/* stress tensor (load state) */
	const dSymMatrixT& stress = s_ij();

	/* pressure */
	output[2] = fStress.Trace()/3.0;

	/* deviatoric Von Mises stress */
	fStress.Deviatoric();
	double J2 = fStress.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[1] = sqrt(3.0*J2);
	
	/* output stress-like internal variable alpha, and check for bifurcation */
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
			//const dMatrixT& modulus = c_ijkl();
			const dMatrixT& modulus = c_perfplas_ijkl();

			/* localization condition checker */
			DetCheckT checker(stress, modulus, Ce);
			AutoArrayT <dArrayT> normals;
			AutoArrayT <dArrayT> slipdirs;
			normals.Dimension(3);
			slipdirs.Dimension(3);
			output[3] = checker.IsLocalized_SS(normals,slipdirs);
		  }
	}
	else
	{
		output[0] = 0.0;
		output[3] = 0.0;
	}

}

/* information about subordinate parameter lists */
void DPSSKStVLoc::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	IsotropicT::DefineSubs(sub_list);
	SSSolidMatT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* DPSSKStVLoc::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = IsotropicT::NewSub(name);
	if (sub)
		return sub;
	else
		return SSSolidMatT::NewSub(name);
}

/* accept parameter list */
void DPSSKStVLoc::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	IsotropicT::TakeParameterList(list);
	SSSolidMatT::TakeParameterList(list);
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void DPSSKStVLoc::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}

