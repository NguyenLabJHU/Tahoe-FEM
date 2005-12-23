/* $Id: GRAD_MRSSKStV.cpp,v 1.30 2005-12-23 14:49:15 kyonten Exp $ */
/* created: Karma Yonten (03/04/2004)                   
   Gradient Enhanced MR Model
*/
#include "GRAD_MRSSKStV.h"
#include "GRAD_MRSSNLHardT.h"  
#include "MFGPMatSupportT.h" 
#include "ElementCardT.h"
#include "StringT.h"
#include "DetCheckT.h"
#include <iostream.h>

using namespace Tahoe;

/* parameters */
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 8;
static const char* Labels[kNumOutput] = {
	    "chi",
	    "cohesion",
	    "Friction Angle",
	    "Dilation Angle",
	    "deviator stress",  // deviator stress, tau or q
	    "press", // pressure (mean stress, p)
	    "VM",  // Von Mises stress
	    "loccheck"}; //	localization check    

/* constructor */
GRAD_MRSSKStV::GRAD_MRSSKStV(void):
	ParameterInterfaceT("small_strain_StVenant_MR_grad_3D"),
	HookeanMatT(3),
	fGRAD_MR(NULL)
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

/* perfectly plastic modulus */
const dMatrixT& GRAD_MRSSKStV::c_perfplas_ijkl(void)
{
	/* elastoplastic correction */
	fModulusPerfPlas =	fGRAD_MR->ModuliPerfPlas(CurrentElement(), CurrIP());
	return fModulusPerfPlas;
}

const dMatrixT& GRAD_MRSSKStV::c_UU1_ijkl(void)
{ 
	fGRAD_MR->Moduli(CurrentElement(), CurrIP()); // call Moduli first
	fModulusUU1 = fGRAD_MR->Moduli_UU1();
	return fModulusUU1;
}

const dMatrixT& GRAD_MRSSKStV::c_UU2_ijkl(void)
{
	fGRAD_MR->Moduli(CurrentElement(), CurrIP()); // call Moduli first
	fModulusUU2 = fGRAD_MR->Moduli_UU2();
	return fModulusUU2;
}

const dMatrixT& GRAD_MRSSKStV::c_ULam1_ij(void)
{
	fGRAD_MR->Moduli(CurrentElement(), CurrIP()); // call Moduli first
	fModulusULam1 =	fGRAD_MR->Moduli_ULam1();
	return fModulusULam1;
}

const dMatrixT& GRAD_MRSSKStV::c_ULam2_ij(void)
{
	fGRAD_MR->Moduli(CurrentElement(), CurrIP()); // call Moduli first
	fModulusULam2 =	fGRAD_MR->Moduli_ULam2();
	return fModulusULam2;
}

const dMatrixT& GRAD_MRSSKStV::c_LamU1_ij(void)
{
	fGRAD_MR->Moduli(CurrentElement(), CurrIP()); // call Moduli first
	fModulusLamU1 =	fGRAD_MR->Moduli_LamU1();
	return fModulusLamU1;
}

const dMatrixT& GRAD_MRSSKStV::c_LamU2_ij(void)
{
	fGRAD_MR->Moduli(CurrentElement(), CurrIP()); // call Moduli first
	fModulusLamU2 =	fGRAD_MR->Moduli_LamU2();
	return fModulusLamU2;
}

const dMatrixT& GRAD_MRSSKStV::c_LamLam1(void)
{
	fGRAD_MR->Moduli(CurrentElement(), CurrIP()); // call Moduli first
	fModulusLamLam1 = fGRAD_MR->Moduli_LamLam1();
	return fModulusLamLam1;
}

const dMatrixT& GRAD_MRSSKStV::c_LamLam2(void)
{
	fGRAD_MR->Moduli(CurrentElement(), CurrIP()); // call Moduli first
	fModulusLamLam2 = fGRAD_MR->Moduli_LamLam2();
	return fModulusLamLam2;
}

/* yield function */
const double& GRAD_MRSSKStV::YieldF()
{
	s_ij(); // call s_ij first
	fYieldFunction = fGRAD_MR->YieldFunction();
	return fYieldFunction;
}

/* stress */
const dSymMatrixT& GRAD_MRSSKStV::s_ij(void)
{
	const char caller[] = "GRAD_MRSSKStV::s_ij";
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	const dSymMatrixT& eps = e();
	const dSymMatrixT& lap_eps = lap_e();
	const dArrayT& lam = pm();
	const dArrayT& lap_lam = lap_pm();
	const dSymMatrixT& e_els = ElasticStrain(eps, element, ip); 
	const dSymMatrixT& lap_e_els = LapElasticStrain(lap_eps, element, ip);
	
	/* check for correct lambda and it's laplacian */
    if (lam[0] < 0.) 
    	ExceptionT::GeneralFail(caller, "negative lambda %e", lam[0]);
    
	/* updated Cauchy stress (return mapping) */
	fStress = fGRAD_MR->StressCorrection(e_els, lap_e_els, lam, lap_lam, element, ip);
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
	output.Dimension(kNumOutput);
	dMatrixT Ce = HookeanMatT::Modulus();
	
	/* stress tensor (load state) */
	const dSymMatrixT& stress = s_ij();

	/* pressure, p */
	output[5] = fStress.Trace()/3.0;

	/* deviatoric Von Mises stress */
	fStress.Deviatoric();
	double J2 = fStress.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[4] = sqrt(J2); /* deviator stress, tau or q */
	output[6] = sqrt(3.0*J2);
	
	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		dArrayT& internal = fGRAD_MR->Internal();
		
		/* stress-like internal variable Chi */
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
		/*
		DetCheckT checker(stress, modulus, Ce);
		AutoArrayT <dArrayT> normals;
		AutoArrayT <dArrayT> slipdirs;
		normals.Dimension(3);
		slipdirs.Dimension(3);
		output[6] = checker.IsLocalized_SS(normals,slipdirs);
		bool checkloc;
		double detA;
		//checkloc = checker.IsLocalized_SS(normals,slipdirs,detA);
		checkloc = checker.IsLocalized_SS(normals,slipdirs);
		if (checkloc) output[18] = 1.0;
		else output[18] = 0.0;
		*/
		output[7] = 0.0;
	}
	else 
	{
		/* initial values of the variables */
		dArrayT& internal = fGRAD_MR->IniInternal();
		
		/* stress-like internal variable Chi */
		output.CopyIn(0, internal);
	}		
}

/* describe the parameters needed by the interface */
void GRAD_MRSSKStV::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	MFGPSSSolidMatT::DefineParameters(list);
	IsotropicT::DefineParameters(list);
	HookeanMatT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void GRAD_MRSSKStV::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	MFGPSSSolidMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
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
		ParameterInterfaceT* params1 = MFGPSSSolidMatT::NewSub(name);
		ParameterInterfaceT* params2 = IsotropicT::NewSub(name);
		if (params1) 
			return params1;
		else if (params2)
			return params2;
		else
			return HookeanMatT::NewSub(name);
	}
}

/* accept parameter list */
void GRAD_MRSSKStV::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	MFGPSSSolidMatT::TakeParameterList(list);
	IsotropicT::TakeParameterList(list);
	HookeanMatT::TakeParameterList(list);
	
	fStress.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(3));
	fModulusCe.Dimension(dSymMatrixT::NumValues(3));
	fModulusPerfPlas.Dimension(dSymMatrixT::NumValues(3));
	fModulusUU1.Dimension(dSymMatrixT::NumValues(3));
	fModulusUU2.Dimension(dSymMatrixT::NumValues(3));
	fModulusULam1.Dimension(dSymMatrixT::NumValues(3),1);
	fModulusULam2.Dimension(dSymMatrixT::NumValues(3),1);
	fModulusLamU1.Dimension(1,dSymMatrixT::NumValues(3));
	fModulusLamU2.Dimension(1,dSymMatrixT::NumValues(3));
	fModulusLamLam1.Dimension(1,1);
	fModulusLamLam2.Dimension(1,1);
	
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
