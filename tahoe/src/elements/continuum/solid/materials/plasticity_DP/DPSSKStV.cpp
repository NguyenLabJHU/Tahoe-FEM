/* $Id: DPSSKStV.cpp,v 1.29 2006-07-02 18:44:31 cfoster01 Exp $ */
/* created: myip (06/01/1999) */
#include "DPSSKStV.h"
#include "SSMatSupportT.h"
#include "DPSSLinHardT.h"

#include "ElementCardT.h"
#include "StringT.h"
#include <iostream.h>

using namespace Tahoe;

/* parameters */
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {
	    "alpha",  // stress-like internal state variable (isotropic linear hardening)
	    "VM",  	// Von Mises stress
	    "press"}; // pressure	    

/* constructor */
DPSSKStV::DPSSKStV(void):
	ParameterInterfaceT("small_strain_StVenant_DP"),
	HookeanMatT(3),
	fDP(NULL)
{
 
}

DPSSKStV::~DPSSKStV(void) { delete fDP; }

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT DPSSKStV::TangentType(void) const { return GlobalT::kNonSymmetric; }

/* update internal variables */
void DPSSKStV::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) fDP->Update(element);
}

/* reset internal variables to last converged solution */
void DPSSKStV::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) fDP->Reset(element);
}

const dSymMatrixT& DPSSKStV::ElasticStrain(const dSymMatrixT& totalstrain, const ElementCardT& element, int ip) {
	return fDP->ElasticStrain(totalstrain, element, ip);
}

/* modulus */
const dMatrixT& DPSSKStV::c_ijkl(void)
{
	fModulus.SumOf(HookeanMatT::Modulus(),
	fDP->ModuliCorrection(CurrentElement(), CurrIP()));	
	//if (fModulus != HookeanMatT::Modulus())
	//	cout << "plastic modulus, element number " <<  CurrElementNumber() << endl;
	return fModulus;
}

/* elastic modulus */
const dMatrixT& DPSSKStV::ce_ijkl(void)
{
  return HookeanMatT::Modulus();
}


/* stress */
const dSymMatrixT& DPSSKStV::s_ij(void)
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


/* returns the strain energy density for the specified strain */
double DPSSKStV::StrainEnergyDensity(void)
{
	return HookeanEnergy(fDP->ElasticStrain(e(), CurrentElement(), CurrIP()));
}

/* returns the number of variables computed for nodal extrapolation
 * during for element output, ie. internal variables. Returns 0
 * by default. */
int DPSSKStV::NumOutputVariables(void) const  { return kNumOutput; } 
void DPSSKStV::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void DPSSKStV::ComputeOutput(dArrayT& output)
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
		dArrayT& internal = fDP->Internal();
		output[0] = internal[DPSSLinHardT::kalpha];
		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == DPSSLinHardT::kIsPlastic)
		  {
			output[0] -= fDP->H_prime()*internal[DPSSLinHardT::kdgamma];
		  }
	}
	else
	{
		output[0] = 0.0;
	}
}

/* describe the parameters needed by the interface */
void DPSSKStV::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSIsotropicMatT::DefineParameters(list);
	HookeanMatT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void DPSSKStV::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSIsotropicMatT::DefineSubs(sub_list);
	HookeanMatT::DefineSubs(sub_list);

	/* parameters for Drucker-Prager plasticity */
	sub_list.AddSub("DP_SS_linear_hardening");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* DPSSKStV::NewSub(const StringT& name) const
{
	if (name == "DP_SS_linear_hardening")
		return new DPSSLinHardT(0, 0.0, 0.0);
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
void DPSSKStV::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSIsotropicMatT::TakeParameterList(list);
	HookeanMatT::TakeParameterList(list);

	fStress.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(3));

	/* construct Drucker-Prager solver */
	fDP = new DPSSLinHardT(NumIP(), Mu(), Lambda());
	fDP->TakeParameterList(list.GetList("DP_SS_linear_hardening"));
}

bool DPSSKStV::IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, 
                             AutoArrayT <double> &detAs, AutoArrayT <double> &dissipations_fact)
 {
 #pragma unused(dissipations_fact)
 
     /* elastic modulus */
     /* this uses same space as c_ijkl(), so save separatley first */
     const dMatrixT modulus_e = ce_ijkl();
 
     /* localization condition checker */
     DetCheckT checker(s_ij(), c_ijkl(), modulus_e);
     normals.Dimension(NumSD());
     slipdirs.Dimension(NumSD());
     return checker.IsLocalized_SS(normals, slipdirs, detAs);
}
 
bool DPSSKStV::IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, double &detA)
 {


	 
	 /* elastic modulus */
     /* this uses same space as c_ijkl(), so save separatley first */
     dMatrixT modulus_e = ce_ijkl();
	 dMatrixT modulus = c_ijkl();

	 if ( modulus == modulus_e)
	    return false;

     /* localization condition checker */
     DetCheckT checker(s_ij(), c_ijkl(), modulus_e);
     normals.Dimension(NumSD());
     slipdirs.Dimension(NumSD());
     bool dummy = checker.IsLocalized_SS(normals, slipdirs, detA);
	 
	/*plasticity = loc */
	return true;
	
	/*traditional */
	//return dummy;
}
 
bool DPSSKStV::IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs)
{
     double dummyDetA = 0.0;
     return IsLocalized(normals, slipdirs, dummyDetA);
}


/*************************************************************************
 * Protected
 *************************************************************************/

/* set modulus */
void DPSSKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}

