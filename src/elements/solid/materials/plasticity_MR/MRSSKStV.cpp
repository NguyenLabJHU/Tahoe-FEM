/* created: Majid T. Manzari (04/16/2001) */
#include "MRSSKStV.h"
#include "SSMatSupportT.h"
#include "MRSSNLHardT.h"

#include "ElementCardT.h"
#include "StringT.h"
#include "DetCheckT.h"
#include <iostream.h>

using namespace Tahoe;

/* parameters */
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 7;
static const char* Labels[kNumOutput] = {
	    "chi",
	    "cohesion",
	    "Friction Angle",
	    "Dilation Angle",
	    "VM",  // Von Mises stress
	    "press", // pressurefmo
	    "loccheck"}; // localization check	    

/* constructor */
MRSSKStV::MRSSKStV(void):
	ParameterInterfaceT("small_strain_StVenant_MR"),
	HookeanMatT(3),
	fMR(NULL)
{

}
	
/* destructor */
MRSSKStV::~MRSSKStV(void) { delete fMR; }

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT MRSSKStV::TangentType(void) const { return GlobalT::kNonSymmetric; }

/* update internal variables */
void MRSSKStV::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) fMR->Update(element);
}

/* reset internal variables to last converged solution */
void MRSSKStV::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) fMR->Reset(element);
}

const dSymMatrixT& MRSSKStV::ElasticStrain(const dSymMatrixT& totalstrain, const ElementCardT& element, int ip) 
{
	return fMR->ElasticStrain(totalstrain, element, ip);
}

/* modulus */
const dMatrixT& MRSSKStV::c_ijkl(void)
{
	fModulus = fMR->Moduli(CurrentElement(), CurrIP());
	return fModulus;
}

/*perfectly plastic modulus */
const dMatrixT& MRSSKStV::c_perfplas_ijkl(void)
{
	fModulusPerfPlas = fMR->ModuliPerfPlas(CurrentElement(), CurrIP());
	return fModulusPerfPlas;
}

/* stress */
const dSymMatrixT& MRSSKStV::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);

	/* Updated Cauchy stress (return mapping) */
	fStress = fMR->StressCorrection(e_els, element, ip);
	return fStress;	
}

/* returns the strain energy density for the specified strain */
double MRSSKStV::StrainEnergyDensity(void)
{
	return 0.;
}

/* returns the number of variables computed for nodal extrapolation
 * during for element output, ie. internal variables. Returns 0
 * by default. */
int MRSSKStV::NumOutputVariables(void) const  { return kNumOutput; } 
void MRSSKStV::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void MRSSKStV::ComputeOutput(dArrayT& output)
{
	dMatrixT Ce = HookeanMatT::Modulus();
	
	/* stress tensor (load state) */
	const dSymMatrixT& stress = s_ij();

	/* pressure */
	output[5] = fStress.Trace()/3.0;

	/* deviatoric Von Mises stress */
	fStress.Deviatoric();
	double J2 = fStress.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[4] = sqrt(3.0*J2);
	
	/* stress-like internal variable chi */
	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		dArrayT& internal = fMR->Internal();
		output[0] = internal[MRSSNLHardT::kchi];
		output[1] = internal[MRSSNLHardT::kc];
		output[2] = internal[MRSSNLHardT::ktanphi];
		output[3] = internal[MRSSNLHardT::ktanpsi];
		
		// check for localization
		// compute modulus 
		const dMatrixT& modulus = c_ijkl();
		// perfectly plastic modulus not implemented yet
		//const dMatrixT& modulus = c_perfplas_ijkl();

		/* localization condition checker */
		DetCheckT checker(stress, modulus, Ce);
		AutoArrayT <dArrayT> normals;
		AutoArrayT <dArrayT> slipdirs;
		normals.Dimension(3);
		slipdirs.Dimension(3);
		output[6] = checker.IsLocalized_SS(normals,slipdirs);
	}	
	else
	{
		output[6] = 0.0;
	}
}

/* describe the parameters needed by the interface */
void MRSSKStV::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSIsotropicMatT::DefineParameters(list);
	HookeanMatT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void MRSSKStV::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSIsotropicMatT::DefineSubs(sub_list);
	HookeanMatT::DefineSubs(sub_list);
	
	/* parameters for pressure sensitive plasticity with localization */
	sub_list.AddSub("MR_SS_nonlinear_hardening");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MRSSKStV::NewSub(const StringT& name) const
{
	if (name == "MR_SS_nonlinear_hardening")
		return new MRSSNLHardT(0, 0.0, 0.0);
	else
	{
		/* inherited */
		ParameterInterfaceT* params = SSIsotropicMatT::NewSub(name);
		if (params) 
			return params;
		else
			return HookeanMatT::NewSub(name);
	}
}

/* accept parameter list */
void MRSSKStV::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSIsotropicMatT::TakeParameterList(list);
	HookeanMatT::TakeParameterList(list);
	
	fStress.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(3));
	fModulusCe.Dimension(dSymMatrixT::NumValues(3));
	fModulusPerfPlas.Dimension(dSymMatrixT::NumValues(3));

	/* construct MR solver */
	fMR = new MRSSNLHardT(NumIP(), Mu(), Lambda());
	fMR->TakeParameterList(list.GetList("MR_SS_nonlinear_hardening"));
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void MRSSKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}
