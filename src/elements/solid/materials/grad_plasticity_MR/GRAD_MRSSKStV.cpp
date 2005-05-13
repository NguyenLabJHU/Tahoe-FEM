/* $Id: GRAD_MRSSKStV.cpp,v 1.15 2005-05-13 22:03:22 kyonten Exp $ */
/* created: Karma Yonten (03/04/2004)                   
   MR version modified to incorporate gradient plasticity 
   theory.
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
const int kNumOutput = 19;
/* Note: For the time being stresses and strains are
   passed to MFGP_AssemblyT as internal state variables
   of the material
*/
static const char* Labels[kNumOutput] = {
	    "eps11","eps22","eps33","eps23","eps13","eps12",
	    "s11","s22","s33","s23","s13","s12",
	    "chi",
	    "cohesion",
	    "Friction Angle",
	    "Dilation Angle",
	    "VM",  // Von Mises stress
	    "press", // pressurefmo
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

/* initialize laplacian of strain and lambda, and lambda, all at ip */
void GRAD_MRSSKStV::Initialize(ElementCardT element, int ip, int n_ip,
                    dSymMatrixT strain_ip, dSymMatrixT strain_lap_ip, 
					dArrayT lambda_ip, dArrayT lambda_lap_ip)
{
	/*fStrain_IP = strain_ip;
	fLapStrain_IP = strain_lap_ip;
    fLambdaPM_IP = lambda_ip;
    fLapLambdaPM_IP = lambda_lap_ip;
    curr_element = element;
    curr_ip = ip;
    num_ip = n_ip;*/
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

const dMatrixT& GRAD_MRSSKStV::c_UU1_ijkl(void)
{
	fModulusUU1 = fGRAD_MR->Moduli_UU1();
	return fModulusUU1;
}

const dMatrixT& GRAD_MRSSKStV::c_UU2_ijkl(void)
{
	fModulusUU2 = fGRAD_MR->Moduli_UU2();
	return fModulusUU2;
}

const dMatrixT& GRAD_MRSSKStV::c_ULam1_ij(void)
{
	fModulusULam1 =	fGRAD_MR->Moduli_ULam1();
	return fModulusULam1;
}

const dMatrixT& GRAD_MRSSKStV::c_ULam2_ij(void)
{
	fModulusULam2 =	fGRAD_MR->Moduli_ULam2();
	return fModulusULam2;
}

const dMatrixT& GRAD_MRSSKStV::c_LamU1_ij(void)
{
	fModulusLamU1 =	fGRAD_MR->Moduli_LamU1();
	return fModulusLamU1;
}

const dMatrixT& GRAD_MRSSKStV::c_LamU2_ij(void)
{
	fModulusLamU2 =	fGRAD_MR->Moduli_LamU2();
	return fModulusLamU2;
}

const dMatrixT& GRAD_MRSSKStV::c_LamLam1(void)
{
	fModulusLamLam1 = fGRAD_MR->Moduli_LamLam1();
	return fModulusLamLam1;
}

const dMatrixT& GRAD_MRSSKStV::c_LamLam2(void)
{
	fModulusLamLam2 = fGRAD_MR->Moduli_LamLam2();
	return fModulusLamLam2;
}

/* yield function */
const double& GRAD_MRSSKStV::YieldF(void)
{
	fYieldFunction = fGRAD_MR->YieldFunction();
	return fYieldFunction;
}

/* stress */
const dSymMatrixT& GRAD_MRSSKStV::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	const dSymMatrixT& eps = e();
	const dSymMatrixT& lap_eps = lap_e();
	const dArrayT& lam = pm();
	const dArrayT& lap_lam = lap_pm();
	const dSymMatrixT& e_els = ElasticStrain(eps, element, ip); 
	const dSymMatrixT& lap_e_els = LapElasticStrain(lap_eps, element, ip);
	
	/* Updated Cauchy stress (return mapping) */
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

	/* stress components */
	output[6] = stress(0,0);
	output[7] = stress(1,1);
	output[8] = stress(2,2);
	output[9] = stress(1,2);
	output[10] = stress(0,2);
	output[11] = stress(0,1);
	/* pressure */
	output[17] = fStress.Trace()/3.0;

	/* deviatoric Von Mises stress */
	fStress.Deviatoric();
	double J2 = fStress.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[16] = sqrt(3.0*J2);
	
	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		
		dArrayT& internal = fGRAD_MR->Internal();
		/* strain components */
		output[0] = internal[GRAD_MRSSNLHardT::keps11];
		output[1] = internal[GRAD_MRSSNLHardT::keps22];
		output[2] = internal[GRAD_MRSSNLHardT::keps33];
		output[3] = internal[GRAD_MRSSNLHardT::keps23];
		output[4] = internal[GRAD_MRSSNLHardT::keps13];
		output[5] = internal[GRAD_MRSSNLHardT::keps12];
		
		/* stress-like internal variable Chi */
		output[12] = internal[GRAD_MRSSNLHardT::kchi];
		output[13] = internal[GRAD_MRSSNLHardT::kc];
		output[14] = internal[GRAD_MRSSNLHardT::ktanphi];
		output[15] = internal[GRAD_MRSSNLHardT::ktanpsi];
		
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
		output[18] = 0.0;
	}	
	else
	{
		output[18] = 0.0;
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
