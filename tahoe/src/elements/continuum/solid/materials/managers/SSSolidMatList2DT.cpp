/* $Id: SSSolidMatList2DT.cpp,v 1.1.4.1 2004-04-08 07:33:04 paklein Exp $ */
#include "SSSolidMatList2DT.h"
#include "SSMatSupportT.h"

#include "fstreamT.h"
#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#endif

#include "SSKStV2D.h"
#include "SSCubic2DT.h"

#ifdef PLASTICITY_J2_MATERIAL
#include "J2SSKStV2D.h"
#include "LocalJ2SSNonlinHard2D.h"
#include "GradJ2SSNonlinHard2D.h"
#endif

#ifdef PLASTICITY_DP_MATERIAL
#include "DPSSKStV2D.h"
#endif

#ifdef PLASTICITY_MR_MATERIAL_DEV
#include "MRSSKStV2D.h"
#endif

#ifdef VISCOELASTIC_MATERIALS_DEV
#include "SSSV_KStV2D.h"
#endif

#ifdef VISCOELASTICITY
#include "SSLinearVE2D.h"
#endif

#ifdef J2PLASTICITY_MATERIALS_DEV
#include "SSJ2LinHard2D.h"
#include "SSJ2LinHard3Dplane.h"
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
#include "ABAQUS_SS_BCJ_ISO.h"
#endif
#endif

#ifdef FOSSUM_MATERIAL_DEV
#include "FossumSSIso2DT.h"
#endif

using namespace Tahoe;

/* constructor */
SSSolidMatList2DT::SSSolidMatList2DT(int length, const SSMatSupportT& support):
	SolidMatListT(length, support),
	fSSMatSupport(&support),
	fGradSSMatSupport(NULL)
{
	SetName("small_strain_material_2D");
	if (fSSMatSupport->NumSD() != 2)
		ExceptionT::GeneralFail("SSSolidMatList2DT::SSSolidMatList2DT");

#ifdef __NO_RTTI__
	cout << "\n SSSolidMatList2DT::SSSolidMatList2DT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;
#endif

#ifdef GRAD_SMALL_STRAIN_DEV
	/* cast to gradient enhanced small strain support */
	fGradSSMatSupport = TB_DYNAMIC_CAST(const GradSSMatSupportT*, fSSMatSupport);
#endif
}

SSSolidMatList2DT::SSSolidMatList2DT(void):
	fSSMatSupport(NULL),
	fGradSSMatSupport(NULL)	
{
	SetName("small_strain_material_2D");

#ifdef __NO_RTTI__
	cout << "\n SSSolidMatList2DT::SSSolidMatList2DT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;
#endif
}

/* read material data from the input stream */
void SSSolidMatList2DT::ReadMaterialData(ifstreamT& in)
{
	const char caller[] = "SSSolidMatList2DT::ReadMaterialData";

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
			ExceptionT::BadInputValue(caller, "repeated material number: %d", matnum + 1);
		
		/* add to the list of materials */
		switch (matcode)
		{
			case kSSKStV:
			{
				fArray[matnum] = new SSKStV2D(in, *fSSMatSupport);
				break;
			}
			case kSSCubic:
			{ 
				fArray[matnum] = new SSCubic2DT(in, *fSSMatSupport);
				break;
			}
			case kJ2SSKStV:
			{
#ifdef PLASTICITY_J2_MATERIAL
				fArray[matnum] = new J2SSKStV2D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kDPSSKStV:
			{
#ifdef PLASTICITY_DP_MATERIAL
				fArray[matnum] = new DPSSKStV2D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_DP_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kMRSSKStV:
			{
#ifdef PLASTICITY_MR_MATERIAL_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new MRSSKStV2D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_MR_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kFossumSSIso:
			{
#ifdef FOSSUM_MATERIAL_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new FossumSSIso2DT(in, *fSSMatSupport);
				fHasHistory = true;		       		      				
				break;
#else
				ExceptionT::BadInputValue(caller, "FOSSUM_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kSSLinearVE:
			{
#ifdef VISCOELASTICITY
				fArray[matnum] = new SSLinearVE2D(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTICITY not enabled: %d", matcode);
#endif
			}
			case kLocJ2SSNlHard:
			{
#ifdef PLASTICITY_J2_MATERIAL
				fArray[matnum] = new LocalJ2SSNonlinHard2D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kGrdJ2SSNlHard:
			{
#ifdef PLASTICITY_J2_MATERIAL
				fArray[matnum] = new GradJ2SSNonlinHard2D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kABAQUS_BCJ_ISO:
			{
#ifdef __F2C__
#if defined(ABAQUS_MATERIAL) && defined(ABAQUS_BCJ_MATERIAL_DEV)
	
				/* small strain version */
				fArray[matnum] = new ABAQUS_SS_BCJ_ISO(in, *fSSMatSupport);

				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "ABAQUS_MATERIAL or ABAQUS_BCJ_MATERIAL_DEV not enabled: %d", matcode);
#endif
#else
				ExceptionT::BadInputValue(caller, "model requires f2c support: %d", kABAQUS_BCJ_ISO);
#endif /* __F2C__ */	
			}			
			case kSSSVKStV:
			{
#ifdef VISCOELASTIC_MATERIALS_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new SSSV_KStV2D(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTIC_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			case kSSJ2LinHard:
			{
#ifdef J2PLASTICITY_MATERIALS_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new SSJ2LinHard2D(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "J2PLASITICITY_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			case kSSJ2LinHardplane:
			{
#ifdef J2PLASTICITY_MATERIALS_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new SSJ2LinHard3Dplane(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "J2PLASITICITY_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			default:
				ExceptionT::BadInputValue(caller, "unknown material code: %d", matcode);
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
		ExceptionT::Throw(error, caller, "exception constructing material %d, index %d, code %d",
			i+1, matnum+1, matcode);
	}
}

/* return true if the list contains plane stress models */
bool SSSolidMatList2DT::HasPlaneStress(void) const
{
	/* check materials */
	for (int i = 0; i < Length(); i++)
	{
		/* get pointer to Material2DT */
		const ContinuumMaterialT* cont_mat = fArray[i];
		const SolidMaterialT* sol_mat = TB_DYNAMIC_CAST(const SolidMaterialT*, cont_mat);
		
		/* assume materials that don't have Material2DT are plane strain */
		if (sol_mat && sol_mat->Constraint() == SolidMaterialT::kPlaneStress) 
			return true;
	}
	return false;
}

/* information about subordinate parameter lists */
void SSSolidMatList2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("ss_material_list_2D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void SSSolidMatList2DT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
		SubListT& sub_sub_list) const
{
	if (sub == "ss_material_list_2D")
	{
		order = ParameterListT::Choice;
	
		sub_sub_list.AddSub("small_strain_cubic_2D");
		sub_sub_list.AddSub("small_strain_StVenant_2D");
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSSolidMatList2DT::NewSub(const StringT& list_name) const
{
	/* try to construct material */
	SSSolidMatT* ss_solid_mat = NewSSSolidMat(list_name);
	if (ss_solid_mat)
		return ss_solid_mat;
	else /* inherited */
		return SolidMatListT::NewSub(list_name);
}

/* accept parameter list */
void SSSolidMatList2DT::TakeParameterList(const ParameterListT& list)
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
		}
	}

	/* transfer */
	Dimension(materials.Length());
	for (int i = 0; i < materials.Length(); i++)
		fArray[i] = materials[i];
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct the specified material or NULL if the request cannot be completed */
SSSolidMatT* SSSolidMatList2DT::NewSSSolidMat(const StringT& list_name) const
{
	SSSolidMatT* mat = NULL;

	if (list_name == "small_strain_cubic_2D")
		mat = new SSCubic2DT;
	else if (list_name == "small_strain_StVenant_2D")
		mat = new SSKStV2D;

	/* set support */
	if (mat) mat->SetSSMatSupport(fSSMatSupport);

	return mat;
}
