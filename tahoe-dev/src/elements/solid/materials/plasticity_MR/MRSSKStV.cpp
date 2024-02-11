/* $Id: MRSSKStV.cpp,v 1.14 2011-12-01 20:38:10 beichuan Exp $ */
/* created: Majid T. Manzari (04/16/2001) */
#include "MRSSKStV.h"
#include "SSEnhLocMatSupportT.h"
//#include "SSMatSupportT.h"
#include "MRSSNLHardT.h"

#include "ElementCardT.h"
#include "StringT.h"
#include "DetCheckT.h"
#include <iostream>

#include "DevelopmentElementsConfig.h"

using namespace Tahoe;

/* parameters */
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 10;
static const char* Labels[kNumOutput] = {
	    "chi",
	    "cohesion",
	    "Friction Angle",
	    "Dilation Angle",
	    "VM",  // Von Mises stress
	    "press", // pressure
	    "fyield", // yield function value
	    "dlambda", // plastic multiplier increment
	    "ip_loccheck", // localization check
	    "el_locflag"}; // element localization flag

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

/* elastic modulus */
const dMatrixT& MRSSKStV::ce_ijkl(void)
{
	fModulusCe = HookeanMatT::Modulus();
	return fModulusCe;
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
	int elem = CurrElementNumber();
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);

#ifdef ENHANCED_STRAIN_LOC_DEV	
	element_locflag = 0;
	if (element.IsAllocated()) 
	{
		element_locflag = fSSEnhLocMatSupport->ElementLocflag(elem);
	}
	if ( element_locflag == 2 )
	{
		fStress = fSSEnhLocMatSupport->ElementStress(elem,ip);
	}
	else
	{
		/* modify Cauchy stress (return mapping) */
		fStress = fMR->StressCorrection(e_els, element, ip, fSSMatSupport->GroupIterationNumber());
	}
#else
	/* modify Cauchy stress (return mapping) */
	fStress = fMR->StressCorrection(e_els, element, ip, fSSMatSupport->GroupIterationNumber());
#endif

	/* Updated Cauchy stress (return mapping) */
	//fStress = fMR->StressCorrection(e_els, element, ip);
	//fStress = fMR->StressCorrection(e_els, element, ip, fSSMatSupport->GroupIterationNumber());
	return fStress;	
}

/*
* Test for localization using "current" values for Cauchy
* stress and the spatial tangent moduli. Returns true if the
* determinant of the acoustic tensor is negative and returns
* the normals and slipdirs. Returns false if the determinant is positive.
*/
//#if 0
bool MRSSKStV::IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, 
							AutoArrayT <double> &detAs, AutoArrayT <double> &dissipations_fact)
{
	/* stress tensor */
	const dSymMatrixT& stress = s_ij();
	
	ElementCardT& element = CurrentElement();
	iArrayT& Flags = element.IntegerData();
	dArrayT& internal = fMR->Internal();
	int ip = CurrIP();
	if (internal[MRSSNLHardT::kdlambda] > 0.0) Flags[ip] = 0; //kIsPlastic=0
	/* elasto-plastic tangent moduli */
	//const dMatrixT& modulus = c_perfplas_ijkl(); //not implemented
	const dMatrixT& modulus = c_ijkl();
	
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
		/* calculate dissipation for each normal and slipdir;
		need to update this code
		*/
		double dissip = 1.0;
		normals.Top();
		dissipations_fact.Free();
		while (normals.Next())
		{
			//this is just temporary until implemented properly
			dissipations_fact.Append(dissip);
		}
	}
	
	return checkloc;
}
//#endif

/* read initial stress from file */
void MRSSKStV::InitialStress(dSymMatrixT& Stress0)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	
	/* read data from file */
	const StringT& input_file = fSSMatSupport->InputFile();
	StringT InitialStressData;
	InitialStressData.Root(input_file);
	// read element by element, ip by ip
	// if data is generated by a program
	Stress0=0.;	
}

/* collect global and local iterations info for each time step */
void MRSSKStV::GetIterationInfo(bool get_iters, int loc_iters)
{
	if(get_iters) {
		
		/* write data */
		const StringT& input_file = fSSMatSupport->InputFile();
	
		/* output filenames */
		StringT iterdata;
		iterdata.Root(input_file);
	
		iterdata.Append(".iterdata.", fSSMatSupport->StepNumber());
		
		/* open output streams */
		ofstreamT out_idata(iterdata);

		/* write */
		out_idata << "Time step number:       " << setw(kIntWidth) 
		          << fSSMatSupport->StepNumber() << endl;
		out_idata << "Global iteration number:" << setw(kIntWidth) 
		          << fSSMatSupport->GroupIterationNumber() << endl;
		out_idata << "Number of local iterations to converge:" << setw(kIntWidth) 
		          << loc_iters << endl;
		out_idata.close(); /* close */
		
	} // if(get_iters)
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
	//const ElementCardT& element = CurrentElement();
	ElementCardT& element = CurrentElement();
	int elem = CurrElementNumber();
	int ip = CurrIP();
	/* get flags */
	iArrayT& Flags = element.IntegerData();
	if (element.IsAllocated())
	//if(element.IsAllocated() && (element.IntegerData())[ip] == 0)
	{
		dArrayT& internal = fMR->Internal();
		output[0] = internal[MRSSNLHardT::kchi];
		output[1] = internal[MRSSNLHardT::kc];
		output[2] = internal[MRSSNLHardT::ktanphi];
		output[3] = internal[MRSSNLHardT::ktanpsi];
		output[6] = internal[MRSSNLHardT::kftrial];
		output[7] = internal[MRSSNLHardT::kdlambda];
		
		// check for localization
		// set flag first for Moduli to be elastoplastic for localization check
		if (internal[MRSSNLHardT::kdlambda] > 0.0) Flags[ip] = 0; //kIsPlastic=0
		const dMatrixT& modulus = c_ijkl();
		//const dMatrixT& modulus = c_perfplas_ijkl(); //not implemented
		DetCheckT checker(stress, modulus, Ce);
		AutoArrayT <dArrayT> normals;
		AutoArrayT <dArrayT> slipdirs;
		normals.Dimension(3);
		slipdirs.Dimension(3);
		output[8] = 0.0;
		double detA=1.0;
		if(checker.IsLocalized_SS(normals,slipdirs,detA)) output[8] = 1.0;
		
		// element localization flag
		output[9] = 0;
		#ifdef ENHANCED_STRAIN_LOC_DEV	
		element_locflag = fSSEnhLocMatSupport->ElementLocflag(elem);
		if (element_locflag > 0) output[9] = element_locflag;
		#endif
	}	
	else
	{
		output[8] = 0.0;
		output[9] = 0.0;
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
	
	/* cast to small strain embedded discontinuity material pointer */
	fSSEnhLocMatSupport = TB_DYNAMIC_CAST(const SSEnhLocMatSupportT*, fSSMatSupport);
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void MRSSKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}
