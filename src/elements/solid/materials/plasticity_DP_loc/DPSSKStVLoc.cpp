/* $Id: DPSSKStVLoc.cpp,v 1.21 2005-04-25 05:05:57 raregue Exp $ */
/* created: myip (06/01/1999) */
#include "DPSSKStVLoc.h"
#include "SSEnhLocMatSupportT.h"
#include "DPSSLinHardLocT.h"

#include "ElementCardT.h"
#include "StringT.h"
#include "DetCheckT.h"
#include <iostream.h>

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#endif

using namespace Tahoe;

/* parameters */
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {
	"kappa",	// stress-like internal state variable (isotropic linear hardening)
	"VM",		// Von Mises stress
	"press"};	// pressure

/* constructor */
DPSSKStVLoc::DPSSKStVLoc(void):
	ParameterInterfaceT("small_strain_StVenant_DP_Loc"),
	HookeanMatT(3),
	fDP(NULL),
	fSSEnhLocMatSupport(NULL)
{
 
}

/* destructor */
DPSSKStVLoc::~DPSSKStVLoc(void) 
{ 
	delete fDP;
}

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT DPSSKStVLoc::TangentType(void) const { return GlobalT::kNonSymmetric; }

/* update internal variables */
void DPSSKStVLoc::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) fDP->Update(element, fSSMatSupport->TimeStep());
}

/* reset internal variables to last converged solution */
void DPSSKStVLoc::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) fDP->Reset(element);
}

const dSymMatrixT& DPSSKStVLoc::ElasticStrain(const dSymMatrixT& totalstrain, 
											const ElementCardT& element, int ip) 
{
	//cout << "totalstrain= \n" << totalstrain <<endl << endl;

	return fDP->ElasticStrain(totalstrain, element, ip);
}

/* modulus */
const dMatrixT& DPSSKStVLoc::c_ijkl(void)
{
	fModulus.SumOf(HookeanMatT::Modulus(), 
		fDP->ModuliCorrection(CurrentElement(), CurrIP(), fSSMatSupport->TimeStep()));

	return fModulus;
}

/* elastic modulus */
const dMatrixT& DPSSKStVLoc::ce_ijkl(void)
{
	fModulusCe = HookeanMatT::Modulus();
	return fModulusCe;
}

/* elastoplastic modulus */
const dMatrixT& DPSSKStVLoc::c_ep_ijkl(void)
{
	fModulusEP.SumOf(HookeanMatT::Modulus(),
		fDP->ModuliCorrectionEP(CurrentElement(), CurrIP()));

	return fModulusEP;
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
	int elem = CurrElementNumber();
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);

	//cout << "e_tot= \n" << e_tot <<endl << endl;
	//cout << "e_els= \n" << e_els <<endl << endl; 
	
#ifdef ENHANCED_STRAIN_LOC_DEV	
	int element_locflag = 0;
	if (element.IsAllocated()) 
	{
		element_locflag = fSSEnhLocMatSupport->ElementLocflag();
	}
	if ( element_locflag == 2 )
	{
		fStress = fSSEnhLocMatSupport->ElementStress(ip);
	}
	else
	{
		/* elastic stress */
		HookeanStress(e_els, fStress);

		/* modify Cauchy stress (return mapping) */
		fStress += fDP->StressCorrection(e_els, element, ip, fSSMatSupport->TimeStep());
	}
#else
	/* elastic stress */
	HookeanStress(e_els, fStress);

	/* modify Cauchy stress (return mapping) */
	fStress += fDP->StressCorrection(e_els, element, ip, fSSMatSupport->TimeStep());
#endif

	return fStress;
}


/*
* Test for localization using "current" values for Cauchy
* stress and the spatial tangent moduli. Returns true if the
* determinant of the acoustic tensor is negative and returns
* the normals and slipdirs. Returns false if the determinant is positive.
*/

//#if 0
bool DPSSKStVLoc::IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, 
							AutoArrayT <double> &detAs, AutoArrayT <double> &dissipations_fact)
{
	/* stress tensor */
	const dSymMatrixT& stress = s_ij();
			
	/* elasto-plastic tangent moduli */
	const dMatrixT& modulus = c_perfplas_ijkl();
	//const dMatrixT& modulus = c_ep_ijkl();
	//const dMatrixT& modulus = c_ijkl();
	
	/* elastic modulus */
	const dMatrixT& modulus_e = ce_ijkl();

	/* localization condition checker */
	DetCheckT checker(stress, modulus, modulus_e);
	normals.Dimension(NumSD());
	slipdirs.Dimension(NumSD());
	normals.Free();
	slipdirs.Free();
	detAs.Free();
	bool checkloc = checker.IsLocalized_SS(normals,slipdirs,detAs);
	
	if (checkloc)
	{
		/* calculate dissipation for each normal and slipdir */
		normals.Top();
		slipdirs.Top();
		dArrayT normal_tmp, slipdir_tmp;
		normal_tmp.Dimension(NumSD());
		slipdir_tmp.Dimension(NumSD());
		
		dissipations_fact.Free();
		
		double sigmn_scalar, nm, psi, cospsi;
		
		//dSymMatrixT devsig(NumSD());
		//devsig.Deviatoric(stress);
		
		dArrayT& internal = fDP->Internal();
		double kappaISV = internal[DPSSLinHardLocT::kkappa];
		const ElementCardT& element = CurrentElement();
		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == DPSSLinHardLocT::kIsPlastic)
			kappaISV -= fDP->H()*internal[DPSSLinHardLocT::kdgamma];
		
		while (normals.Next())
		{
			normal_tmp = normals.Current();
			slipdirs.Next();
			slipdir_tmp = slipdirs.Current();
			sigmn_scalar = stress.MultmBn(normal_tmp, slipdir_tmp);
			nm = dArrayT::Dot(normal_tmp, slipdir_tmp);
			psi = asin(nm);
			cospsi = cos(psi);
			double dissip = sigmn_scalar;
			dissip -= kappaISV*cospsi;
			dissipations_fact.Append(dissip);
		}
	}
	
	return checkloc;
}
//#endif

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
	
	/* output stress-like internal variable kappa, and check for bifurcation */
	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		dArrayT& internal = fDP->Internal();
		output[0] = internal[DPSSLinHardLocT::kkappa];
		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == DPSSLinHardLocT::kIsPlastic)
		{
			output[0] -= fDP->H()*internal[DPSSLinHardLocT::kdgamma];
			
			// check for localization
			// compute modulus 
			//const dMatrixT& modulus = c_ijkl();
			//const dMatrixT& modulus = c_perfplas_ijkl();
			//const dMatrixT& modulus = c_ep_ijkl();

			/* localization condition checker */
			/*
			DetCheckT checker(stress, modulus, Ce);
			AutoArrayT <dArrayT> normals;
			AutoArrayT <dArrayT> slipdirs;
			normals.Dimension(3);
			slipdirs.Dimension(3);
			output[3] = 0.0;
			if(checker.IsLocalized_SS(normals,slipdirs)) output[3] = 1.0;
			*/
			//output[3] = 0.0;
		  }
	}
	else
	{
		output[0] = 0.0;
		//output[3] = 0.0;
	}

}

/* describe the parameters needed by the interface */
void DPSSKStVLoc::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSIsotropicMatT::DefineParameters(list);
	HookeanMatT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void DPSSKStVLoc::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSIsotropicMatT::DefineSubs(sub_list);
	HookeanMatT::DefineSubs(sub_list);
	
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
		ParameterInterfaceT* params = SSIsotropicMatT::NewSub(name);
		if (params) 
			return params;
		else
			return HookeanMatT::NewSub(name);
	}
}

/* accept parameter list */
void DPSSKStVLoc::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSIsotropicMatT::TakeParameterList(list);
	HookeanMatT::TakeParameterList(list);
	
	fStress.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(3));
	fModulusEP.Dimension(dSymMatrixT::NumValues(3));
	fModulusCe.Dimension(dSymMatrixT::NumValues(3));
	fModulusPerfPlas.Dimension(dSymMatrixT::NumValues(3));

	/* construct Drucker-Prager solver */
	fDP = new DPSSLinHardLocT(NumIP(), Mu(), Lambda());
	fDP->TakeParameterList(list.GetList("DP_Loc_SS_linear_hardening"));
		
	/* cast to small strain embedded discontinuity material pointer */
	fSSEnhLocMatSupport = TB_DYNAMIC_CAST(const SSEnhLocMatSupportT*, fSSMatSupport);
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void DPSSKStVLoc::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}

