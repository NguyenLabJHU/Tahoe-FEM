/* $Id: DPSSKStVLoc.cpp,v 1.5 2004-07-21 20:52:46 raregue Exp $ */
/* created: myip (06/01/1999) */

#include "DPSSKStVLoc.h"
#include "SSMatSupportT.h"
#include "DPSSLinHardLocT.h"

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
DPSSKStVLoc::DPSSKStVLoc(void):
	ParameterInterfaceT("small_strain_StVenant_DP_Loc"),
	HookeanMatT(3),
	fDP(NULL)
{
 
}

/* destructor */
DPSSKStVLoc::~DPSSKStVLoc(void) { delete fDP; }

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT DPSSKStVLoc::TangentType(void) const { return GlobalT::kNonSymmetric; }

/* update internal variables */
void DPSSKStVLoc::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) fDP->Update(element);
}

/* reset internal variables to last converged solution */
void DPSSKStVLoc::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) fDP->Reset(element);
}

const dSymMatrixT& DPSSKStVLoc::ElasticStrain(const dSymMatrixT& totalstrain, const ElementCardT& element, int ip) {
	return fDP->ElasticStrain(totalstrain, element, ip);
}

/* modulus */
const dMatrixT& DPSSKStVLoc::c_ijkl(void)
{
	fModulus.SumOf(HookeanMatT::Modulus(),
	fDP->ModuliCorrection(CurrentElement(), CurrIP()));

	return fModulus;
}

/*discontinuous modulus */
const dMatrixT& DPSSKStVLoc::c_perfplas_ijkl(void)
{
	/* elastoplastic correction */
	fModulusPerfPlas.SumOf(HookeanMatT::Modulus(),
	fDP->ModuliCorrPerfPlas(CurrentElement(), CurrIP()));
	
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
	fStress += fDP->StressCorrection(e_els, element, ip);
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
	return HookeanEnergy(fDP->ElasticStrain(e(), CurrentElement(), CurrIP()));
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
		dArrayT& internal = fDP->Internal();
		output[0] = internal[DPSSLinHardLocT::kalpha];
		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == DPSSLinHardLocT::kIsPlastic)
		{
			output[0] -= fDP->H_prime()*internal[DPSSLinHardLocT::kdgamma];
			
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

/* describe the parameters needed by the interface */
void DPSSKStVLoc::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSSolidMatT::DefineParameters(list);
	IsotropicT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void DPSSKStVLoc::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	IsotropicT::DefineSubs(sub_list);
	SSSolidMatT::DefineSubs(sub_list);
	
	/* parameters for Drucker-Prager plasticity with localization */
	sub_list.AddSub("DP_Loc_SS_linear_hardening");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* DPSSKStVLoc::NewSub(const StringT& name) const
{
	if (name == "DP_Loc_SS_linear_hardening")
		return new DPSSLinHardLocT(0, 0.0, 0.0);
	else
	{
		/* inherited */
		ParameterInterfaceT* params = SSSolidMatT::NewSub(name);
		if (params) 
			return params;
		else
			return IsotropicT::NewSub(name);
	}
}

/* accept parameter list */
void DPSSKStVLoc::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	IsotropicT::TakeParameterList(list);
	SSSolidMatT::TakeParameterList(list);
	
	fStress.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(3));
	fModulusCe.Dimension(dSymMatrixT::NumValues(3));
	fModulusPerfPlas.Dimension(dSymMatrixT::NumValues(3));

	/* set modulus */
	HookeanMatT::Initialize();

	/* construct Drucker-Prager solver */
	fDP = new DPSSLinHardLocT(NumIP(), Mu(), Lambda());
	fDP->TakeParameterList(list.GetList("DP_Loc_SS_linear_hardening"));
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void DPSSKStVLoc::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}

