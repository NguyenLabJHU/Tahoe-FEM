/* $Id: SMRSSKStV.cpp,v 1.1 2006-07-27 13:20:08 kyonten Exp $ */
/* created: Majid T. Manzari (04/16/2001) */
#include "SMRSSKStV.h"
#include "SSMatSupportT.h"
#include "SMRSSNLHardT.h"

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
SMRSSKStV::SMRSSKStV(void):
	ParameterInterfaceT("small_strain_StVenant_SMR"),
	HookeanMatT(3),
	fSMR(NULL)
{

}
	
/* destructor */
SMRSSKStV::~SMRSSKStV(void) { delete fSMR; }

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT SMRSSKStV::TangentType(void) const { return GlobalT::kNonSymmetric; }

/* update internal variables */
void SMRSSKStV::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) fSMR->Update(element);
}

/* reset internal variables to last converged solution */
void SMRSSKStV::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) fSMR->Reset(element);
}

const dSymMatrixT& SMRSSKStV::ElasticStrain(const dSymMatrixT& totalstrain, const ElementCardT& element, int ip) 
{
	return fSMR->ElasticStrain(totalstrain, element, ip);
}

/* modulus */
const dMatrixT& SMRSSKStV::c_ijkl(void)
{
	fModulus = fSMR->Moduli(CurrentElement(), CurrIP());
	return fModulus;
}

/*perfectly plastic modulus */
const dMatrixT& SMRSSKStV::c_perfplas_ijkl(void)
{
	fModulusPerfPlas = fSMR->ModuliPerfPlas(CurrentElement(), CurrIP());
	return fModulusPerfPlas;
}

/* stress */
const dSymMatrixT& SMRSSKStV::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);

	/* Updated Cauchy stress (return mapping) */
	//Stress0 = // initial stress
	fStress = fSMR->StressCorrection(e_els, element, ip);
	return fStress;	
}

/* returns the strain energy density for the specified strain */
double SMRSSKStV::StrainEnergyDensity(void)
{
	return 0.;
}

/* returns the number of variables computed for nodal extrapolation
 * during for element output, ie. internal variables. Returns 0
 * by default. */
int SMRSSKStV::NumOutputVariables(void) const  { return kNumOutput; } 
void SMRSSKStV::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void SMRSSKStV::ComputeOutput(dArrayT& output)
{
	dMatrixT Ce = HookeanMatT::Modulus();
	
	/* stress tensor (load state) */
	s_ij();

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
		dArrayT& internal = fSMR->Internal();
		output[0] = internal[SMRSSNLHardT::kchi];
		output[1] = internal[SMRSSNLHardT::kc];
		output[2] = internal[SMRSSNLHardT::ktanphi];
		output[3] = internal[SMRSSNLHardT::ktanpsi];
		
		// check for localization
		// compute modulus 
		//const dMatrixT& modulus = c_ijkl();
		// perfectly plastic modulus not implemented yet
		//const dMatrixT& modulus = c_perfplas_ijkl();

		/* localization condition checker */
		/*
		DetCheckT checker(stress, modulus, Ce);
		AutoArrayT <dArrayT> normals;
		AutoArrayT <dArrayT> slipdirs;
		normals.Dimension(3);
		slipdirs.Dimension(3);
		bool checkloc;
		double detA;
		checkloc = checker.IsLocalized_SS(normals,slipdirs);
		if (checkloc) output[6] = 1.0;
		else output[6] = 0.0;
		*/
		output[6] = 0.0;
	}	
	else
	{
		output[6] = 0.0;
	}
}

/* describe the parameters needed by the interface */
void SMRSSKStV::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSIsotropicMatT::DefineParameters(list);
	HookeanMatT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void SMRSSKStV::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSIsotropicMatT::DefineSubs(sub_list);
	HookeanMatT::DefineSubs(sub_list);
	
	/* parameters for pressure sensitive plasticity with localization */
	sub_list.AddSub("SMR_SS_nonlinear_hardening");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SMRSSKStV::NewSub(const StringT& name) const
{
	if (name == "SMR_SS_nonlinear_hardening")
		return new SMRSSNLHardT(0, 0.0, 0.0);
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
void SMRSSKStV::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSIsotropicMatT::TakeParameterList(list);
	HookeanMatT::TakeParameterList(list);
	
	fStress.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(3));
	fModulusCe.Dimension(dSymMatrixT::NumValues(3));
	fModulusPerfPlas.Dimension(dSymMatrixT::NumValues(3));

	/* construct MR solver */
	fSMR = new SMRSSNLHardT(NumIP(), Mu(), Lambda());
	fSMR->TakeParameterList(list.GetList("SMR_SS_nonlinear_hardening"));
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void SMRSSKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}
