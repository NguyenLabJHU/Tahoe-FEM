/* $Id: SSSolidMatList1DT.cpp,v 1.1.4.2 2004-06-19 23:28:09 paklein Exp $ */
#include "SSSolidMatList1DT.h"
#include "SSMatSupportT.h"
#include "fstreamT.h"

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

using namespace Tahoe;

/* constructor */
SSSolidMatList1DT::SSSolidMatList1DT(int length, const SSMatSupportT& support):
	SolidMatListT(length, support),
	fSSMatSupport(&support),
	fGradSSMatSupport(NULL)
{
#pragma message("check spatial dimension of material support")
	SetName("small_strain_material_1D");

#ifdef GRAD_SMALL_STRAIN_DEV
	/* cast to gradient enhanced small strain support */
	fGradSSMatSupport = TB_DYNAMIC_CAST(const GradSSMatSupportT*, fSSMatSupport);
#endif
}

SSSolidMatList1DT::SSSolidMatList1DT(void):
	fSSMatSupport(NULL),
	fGradSSMatSupport(NULL)	
{
	SetName("small_strain_material_1D");
}

/* read material data from the input stream */
void SSSolidMatList1DT::ReadMaterialData(ifstreamT& in)
{

ExceptionT::GeneralFail("SSSolidMaterialList1DT::ReadMaterialData");
#if 0
	int i, matnum;
	SolidT::TypeT matcode;
	try {

	/* read material data */
	for (i = 0; i < fLength; i++)
 	{
		in >> matnum; matnum--;
		in >> matcode;
		/* checks */
		if (matnum < 0  || matnum >= fLength) throw ExceptionT::kBadInputValue;
		
		/* repeated material number */
		if (fArray[matnum] != NULL)
		{
			cout << "\n SSSolidMatList1DT::ReadMaterialData: repeated material number: ";
			cout << matnum + 1 << endl;
			throw ExceptionT::kBadInputValue;
		}
		
		/* add to the list of matxxerials */
		switch (matcode)
		{
			case kSSKStV:
			{
				fArray[matnum] = new SSHookean1D(in, *fSSMatSupport);
				break;
		  	}
			case kGradJ2SS:
			{
#ifdef GRAD_SMALL_STRAIN_DEV
				fArray[matnum] = new GradJ2SS1D(in, *fGradSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue("SSSolidMatList1DT::ReadMaterialData", "GRAD_SMALL_STRAIN_DEV not enabled: %d", matcode);
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
				ExceptionT::BadInputValue("SSSolidMatList1DT::ReadMaterialData", "GRAD_SMALL_STRAIN_DEV not enabled: %d", matcode);
#endif
			}
			default:
			{
				cout << "\n SSSolidMatList1DT::ReadMaterialData: unknown material code: ";
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
		cout << "\n SSSolidMatList1DT::ReadMaterialData: exception constructing material " << i+1
		     << '\n' << "     index " << matnum+1 << ", code " << matcode << endl;
		throw error;
	}
#endif
}

/* information about subordinate parameter lists */
void SSSolidMatList1DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("ss_material_list_1D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void SSSolidMatList1DT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
		SubListT& sub_sub_list) const
{
	if (sub == "ss_material_list_1D")
	{
		order = ParameterListT::Choice;
	
		sub_sub_list.AddSub("linear_material_1D");
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSSolidMatList1DT::NewSub(const StringT& list_name) const
{
	/* try to construct material */
	SSSolidMatT* ss_solid_mat = NewSSSolidMat(list_name);
	if (ss_solid_mat)
		return ss_solid_mat;
	else /* inherited */
		return SolidMatListT::NewSub(list_name);
}

/* accept parameter list */
void SSSolidMatList1DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	AutoArrayT<SSSolidMatT*> materials;
	const ArrayT<ParameterListT>& subs = list.Lists();
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		SSSolidMatT* mat = NewSSSolidMat(sub.Name());
		if (mat) {
			materials.Append(mat);
			mat->TakeParameterList(sub);

			/* set flags */
			if (mat->HasHistory()) fHasHistory = true;	
			if (mat->HasThermalStrain()) fHasThermal = true;
			if (mat->HasLocalization()) fHasLocalizers = true;
		}
	}

	/* transfer */
	Dimension(materials.Length());
	for (int i = 0; i < materials.Length(); i++)
		fArray[i] = materials[i];
}

/* construct the specified material or NULL if the request cannot be completed */
SSSolidMatT* SSSolidMatList1DT::NewSSSolidMat(const StringT& list_name) const
{
	SSSolidMatT* mat = NULL;

	if (list_name == "linear_material_1D")
		mat = new SSHookean1D;

	/* set support */
	if (mat) mat->SetSSMatSupport(fSSMatSupport);

	return mat;
}

