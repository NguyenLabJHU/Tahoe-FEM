/* $Id: SolidMatList1DT.cpp,v 1.21 2004-06-28 22:41:30 hspark Exp $ */

#include "SolidMatList1DT.h"

#include "SolidMaterialsConfig.h"
#include "ifstreamT.h"
#include "SolidMatSupportT.h"

/* 1D material types codes */
/* Add small strain linear elastic material here */
#include "SSHookean1D.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#endif

#ifdef GRAD_SMALL_STRAIN_DEV
#include "GradJ2SS1D.h"
#include "J2SSKStV1D.h"
#endif

#ifdef CAUCHY_BORN_MATERIAL
#include "Chain1D.h"
#endif

using namespace Tahoe;

/* constructor */
SolidMatList1DT::SolidMatList1DT(int length, const SolidMatSupportT& support):
	SolidMatListT(length, support)
{
	SetName("solid_materials_1D");
}

SolidMatList1DT::SolidMatList1DT(void)
{
	SetName("solid_materials_1D");
}

/* read material data from the input stream */
void SolidMatList1DT::ReadMaterialData(ifstreamT& in)
{
	const char caller[] = "SolidMatList1DT::ReadMaterialData";

	int i, matnum;
	SolidT::TypeT matcode;
	try {

	/* read material data */
	for (i = 0; i < fLength; i++)
 	{
		in >> matnum; matnum--;
		in >> matcode;

		/* checks */
		if (matnum < 0 || matnum >= fLength) 
			ExceptionT::BadInputValue(caller, "material number %d is out of range [1, %d]",
				matnum+1, fLength);
		
		/* repeated material number */
		if (fArray[matnum] != NULL)
			ExceptionT::BadInputValue(caller, "repeated material number %d", matnum+1);
		
		/* add to the list of matxxerials */
		switch (matcode)
		{
			case kSSKStV:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
				fArray[matnum] = new SSHookean1D(in, *fSSMatSupport);
				break;
		  	}
			case kGradJ2SS:
			{
#ifdef GRAD_SMALL_STRAIN_DEV
				/* check */
				if (!fGradSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new GradJ2SS1D(in, *fGradSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue("SolidMatList1DT::ReadMaterialData", "GRAD_SMALL_STRAIN_DEV not enabled: %d", matcode);
#endif
			}
			case kJ2SSKStV1D:
			{
#ifdef GRAD_SMALL_STRAIN_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new J2SSKStV1D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
#else
				ExceptionT::BadInputValue("SolidMatList1DT::ReadMaterialData", "GRAD_SMALL_STRAIN_DEV not enabled: %d", matcode);
#endif
			}
			case kChain1D:
			{
#ifdef CAUCHY_BORN_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new Chain1D(in, *fFSMatSupport);
				break;
#else
				ExceptionT::BadInputValue("SolidMatList1DT::ReadMaterialData", "CAUCHY_BORN_MATERIAL not enabled: %d", matcode);
#endif
			}
			default:
			{
				cout << "\n SolidMatList1DT::ReadMaterialData: unknown material code: ";
				cout << matcode << '\n' << endl;
				throw ExceptionT::kBadInputValue;
			}
		}

		/* safe cast since all structural */
		SolidMaterialT* pmat = (SolidMaterialT*) fArray[matnum];
		/* verify construction */
		if (!pmat) throw ExceptionT::kOutOfMemory;
		
		/* set thermal LTf pointer */
		int LTfnum = pmat->ThermalStrainSchedule();
		if (LTfnum > -1)
		{
			pmat->SetThermalSchedule(fSolidMatSupport->Schedule(LTfnum));
			
			/* set flag */
			fHasThermal = true;
		}
		
		/* perform initialization */
		pmat->Initialize();
	}  } /* end try */
	
	catch (ExceptionT::CodeT error)
	{
		cout << "\n SolidMatList1DT::ReadMaterialData: exception constructing material " << i+1
		     << '\n' << "     index " << matnum+1 << ", code " << matcode << endl;
		throw error;
	}
}

/* information about subordinate parameter lists */
void SolidMatList1DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);

	/* list of materials an array of choices */
	sub_list.AddSub("solid_material_list_1D", ParameterListT::OnePlus, true);
}

/* return the description of the given inline subordinate parameter list */
void SolidMatList1DT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	if (sub == "solid_material_list_1D")
	{
		order = ParameterListT::Choice;
	
		sub_sub_list.AddSub("small_strain_Hookean_1D");

#ifdef GRAD_SMALL_STRAIN_DEV
		sub_sub_list.AddSub("GradJ2SS1D");
		sub_sub_list.AddSub("J2SSKStV1D");
#endif		
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SolidMatList1DT::NewSub(const StringT& list_name) const
{
	if (list_name == "small_strain_Hookean_1D")
		return new SSHookean1D;
#ifdef GRAD_SMALL_STRAIN_DEV
	else if (list_name == "GradJ2SS1D")
		return new GradJ2SS1D;
	else if (list_name == "J2SSKStV1D")
		return new J2SSKStV1D;
#endif		
	else /* inherited */
		return SolidMatListT::NewSub(list_name);
}

/* error messages */

void SolidMatList1DT::Error_no_small_strain(ostream& out, int matcode) const
{
	out << "\n SolidMatList1DT: material " << matcode
		<< " requires a small strain element" << endl;
	throw ExceptionT::kBadInputValue;
}

void SolidMatList1DT::Error_no_finite_strain(ostream& out, int matcode) const
{
	out << "\n SolidMatList1DT: material " << matcode
		<< " requires a finite strain element" << endl;
	throw ExceptionT::kBadInputValue;
}

void SolidMatList1DT::Error_no_multi_scale(ostream& out, int matcode) const
{
	out << "\n SolidMatList1DT: material " << matcode
		<< " requires a variational multi-scale element" << endl;
	throw ExceptionT::kBadInputValue;
}

