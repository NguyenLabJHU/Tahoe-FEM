/* created: Karma Yonten (03/04/2004)                   
   MR version modified to incorporate gradient plasticity 
   theory.
*/
#include "GRAD_MRSSKStV.h"
#include "SSMatSupportT.h"
#include "GRAD_MRSSNLHardT.h"  

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
	    "loccheck"}; //	localization check    

/* constructor */
GRAD_MRSSKStV::GRAD_MRSSKStV(void):
	ParameterInterfaceT("small_strain_StVenant_GRAD_MR"),
	HookeanMatT(3),
	fGRAD_MR(NULL),
	fYieldFunction(0.0)
{

}

/* destructor */
GRAD_MRSSKStV::~GRAD_MRSSKStV(void) { delete fGRAD_MR; }

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT GRAD_MRSSKStV::TangentType(void) const { return GlobalT::kNonSymmetric; }

/* update internal variables */
void GRAD_MRSSKStV::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) fGRAD_MR->Update(element);
}

/* reset internal variables to last converged solution */
void GRAD_MRSSKStV::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) fGRAD_MR->Reset(element);
}

const dSymMatrixT& GRAD_MRSSKStV::ElasticStrain(const dSymMatrixT& totalstrain, const ElementCardT& element, int ip) 
{
	return fGRAD_MR->ElasticStrain(totalstrain, element, ip);
}

const dSymMatrixT& GRAD_MRSSKStV::LapElasticStrain(const dSymMatrixT& lap_totalstrain, const ElementCardT& element, int ip) 
{
	return fGRAD_MR->LapElasticStrain(lap_totalstrain, element, ip);
}


/* modulus */
const dMatrixT& GRAD_MRSSKStV::c_ijkl(void)
{
	fModulus =	fGRAD_MR->Moduli(CurrentElement(), CurrIP());
	return fModulus;
}

/*perfectly plastic modulus */
const dMatrixT& GRAD_MRSSKStV::c_perfplas_ijkl(void)
{
	/* elastoplastic correction */
	fModulusPerfPlas =	fGRAD_MR->ModuliPerfPlas(CurrentElement(), CurrIP());
	return fModulusPerfPlas;
}

/* yield function */
const double& GRAD_MRSSKStV::YieldF(void)
{
	fYieldFunction = fGRAD_MR->YieldFunction(CurrentElement(), CurrIP());
	return fYieldFunction;
}

/* stress */
const dSymMatrixT& GRAD_MRSSKStV::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	const double& dlam = CurrentDLam();           // passed from the global level
	const double& lap_dlam = CurrentLap_DLam(); // passed from the global level
	const dSymMatrixT& e_tot = e();          //remove thermal strain ??
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);
	const dSymMatrixT& lap_e_tot = lap_e(); // lap_e??
	const dSymMatrixT& lap_e_els = LapElasticStrain(lap_etot, element, ip);
// Note:  This part of the code is incomplete.  Basically, the
// laplacians of the strain and the plastic multiplier need
// to be passed from the higher level to the consititutive 
// routine
	/* Updated Cauchy stress (return mapping) */
	fStress = fGRAD_MR->StressCorrection(e_els, lap_e_els, dlam, lap_dlam, element, ip);
	return fStress;	
}


/* returns the strain energy density for the specified strain */
double GRAD_MRSSKStV::StrainEnergyDensity(void)
{
	return 0.;
}

/* returns the number of variables computed for nodal extrapolation
 * during for element output, ie. internal variables. Returns 0
 * by default. */
int GRAD_MRSSKStV::NumOutputVariables(void) const  { return kNumOutput; } 
void GRAD_MRSSKStV::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void GRAD_MRSSKStV::ComputeOutput(dArrayT& output)
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
	
	/* stress-like internal variable Chi */
	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		
		dArrayT& internal = fGRAD_MR->Internal();
		output[0] = internal[GRAD_MRSSNLHardT::kchi];
		output[1] = internal[GRAD_MRSSNLHardT::kc];
		output[2] = internal[GRAD_MRSSNLHardT::ktanphi];
		output[3] = internal[GRAD_MRSSNLHardT::ktanpsi];
		
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
void GRAD_MRSSKStV::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSIsotropicMatT::DefineParameters(list);
	HookeanMatT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void GRAD_MRSSKStV::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSIsotropicMatT::DefineSubs(sub_list);
	HookeanMatT::DefineSubs(sub_list);
	
	/* parameters for pressure sensitive plasticity with localization */
	sub_list.AddSub("GRAD_MR_SS_nonlinear_hardening");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* GRAD_MRSSKStV::NewSub(const StringT& name) const
{
	if (name == "GRAD_MR_SS_nonlinear_hardening")
		return new GRAD_MRSSNLHardT(0, 0.0, 0.0);
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
void GRAD_MRSSKStV::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSIsotropicMatT::TakeParameterList(list);
	HookeanMatT::TakeParameterList(list);
	
	fStress.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(3));
	fModulusCe.Dimension(dSymMatrixT::NumValues(3));
	fModulusPerfPlas.Dimension(dSymMatrixT::NumValues(3));
	//fYieldFunction(0.0); // kyonten

	/* construct GRAD_MR solver */
	fGRAD_MR = new GRAD_MRSSNLHardT(NumIP(), Mu(), Lambda());
	fGRAD_MR->TakeParameterList(list.GetList("GRAD_MR_SS_nonlinear_hardening"));
}


/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void GRAD_MRSSKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}
