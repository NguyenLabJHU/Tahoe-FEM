/* $Id: SSSolidMatList2DT.cpp,v 1.13 2006-07-21 20:04:31 tdnguye Exp $ */
#include "SSSolidMatList2DT.h"
#include "SSMatSupportT.h"

#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#include "DevelopmentElementsConfig.h"
#endif

#include "SSHookeanMat2DT.h"
#include "SSKStV2D.h"
#include "SSCubic2DT.h"

#ifdef PLASTICITY_J2_MATERIAL
#include "J2SSKStV2D.h"
#include "J2SSKStV2DPlaneStress.h"
#include "LocalJ2SSNonlinHard2D.h"
#include "GradJ2SSNonlinHard2D.h"
#endif

#ifdef PLASTICITY_DP_MATERIAL
#include "DPSSKStV2D.h"
#endif

#ifdef VISCOELASTICITY
#include "SSLinearVE2D.h"
#endif

/* development module materials require solid element development to be enabled */
#ifdef SOLID_ELEMENT_DEV

#ifdef PLASTICITY_DP_LOC_MATERIAL_DEV
#include "DPSSKStVLoc2D.h"
#endif

#ifdef PLASTICITY_MR_MATERIAL_DEV
#include "MRSSKStV2D.h"
#endif

#ifdef VISCOELASTIC_MATERIALS_DEV
#include "SSSV_KStV2D.h"
#endif

#ifdef J2PLASTICITY_MATERIALS_DEV
#include "SSJ2LinHard2D.h"
#include "SSJ2LinHard3Dplane.h"
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
#include "ABAQUS_SS_BCJ_ISO.h"
#endif
#ifdef ABAQUS_DEV
#include "ABAQUS_UMAT_IsoPlast.h"
#endif
#endif /* ABAQUS_MATERIAL */

#ifdef FOSSUM_MATERIAL_DEV
#include "FossumSSIso2DT.h"
#endif

#endif /* SOLID_ELEMENT_DEV */

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
void SSSolidMatList2DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "ss_material_list_2D")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("small_strain_Hookean_2D");
		sub_lists.AddSub("small_strain_cubic_2D");
		sub_lists.AddSub("small_strain_StVenant_2D");

#ifdef PLASTICITY_J2_MATERIAL
		sub_lists.AddSub("small_strain_StVenant_J2_PlaneStrain");
		sub_lists.AddSub("small_strain_StVenant_J2_PlaneStress");
#endif

#ifdef PLASTICITY_DP_MATERIAL
		sub_lists.AddSub("small_strain_StVenant_DP_2D");
#endif

#ifdef PLASTICITY_MR_MATERIAL_DEV
		sub_lists.AddSub("small_strain_StVenant_MR_2D");
#endif

#ifdef PLASTICITY_DP_LOC_MATERIAL_DEV
		sub_lists.AddSub("small_strain_StVenant_DP_Loc_2D");
#endif

#ifdef VISCOELASTICITY
		sub_lists.AddSub("linear_viscoelastic_2D");
#endif
	
#ifdef FOSSUM_MATERIAL_DEV
		sub_lists.AddSub("Fossum_small_strain_2D");
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
		sub_lists.AddSub("ABAQUS_UMAT_SS_BCJ_iso-damage");
#endif
#ifdef ABAQUS_DEV
		sub_lists.AddSub("ABAQUS_UMAT_IsoPlast");
#endif
#endif

	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSSolidMatList2DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	SSSolidMatT* ss_solid_mat = NewSSSolidMat(name);
	if (ss_solid_mat)
		return ss_solid_mat;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void SSSolidMatList2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	const ArrayT<ParameterListT>& subs = list.Lists();
	int count = 0;
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		SSSolidMatT* mat = NewSSSolidMat(sub.Name());
		if (mat) {
			
			/* store pointer */
			(*this)[count++] = mat;

			/* initialize material */
			mat->TakeParameterList(sub);

			/* set flags */
			if (mat->HasHistory()) fHasHistory = true;	
			if (mat->HasThermalStrain()) fHasThermal = true;
			if (mat->HasLocalization()) fHasLocalizers = true;
		}
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct the specified material or NULL if the request cannot be completed */
SSSolidMatT* SSSolidMatList2DT::NewSSSolidMat(const StringT& name) const
{
	SSSolidMatT* mat = NULL;

	if (name == "small_strain_Hookean_2D")
		mat = new SSHookeanMat2DT;
	else if (name == "small_strain_cubic_2D")
		mat = new SSCubic2DT;
	else if (name == "small_strain_StVenant_2D")
		mat = new SSKStV2D;

#ifdef PLASTICITY_J2_MATERIAL
	else if (name == "small_strain_StVenant_J2_PlaneStrain")
		mat = new J2SSKStV2D;
	else if (name == "small_strain_StVenant_J2_PlaneStress")
		mat = new J2SSKStV2DPlaneStress;
#endif

#ifdef PLASTICITY_DP_MATERIAL
	else if (name == "small_strain_StVenant_DP_2D")
		mat = new DPSSKStV2D;
#endif

#ifdef PLASTICITY_MR_MATERIAL_DEV
	else if (name == "small_strain_StVenant_MR_2D")
		mat = new MRSSKStV2D;
#endif

#ifdef PLASTICITY_DP_LOC_MATERIAL_DEV
	else if (name == "small_strain_StVenant_DP_Loc_2D")
		mat = new DPSSKStVLoc2D;
#endif

#ifdef VISCOELASTICITY
	else if (name == "linear_viscoelastic_2D")
		mat = new SSLinearVE2D;
#endif
	
#ifdef FOSSUM_MATERIAL_DEV
	else if (name == "Fossum_small_strain_2D")
	  mat = new FossumSSIso2DT;
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
	else if (name == "ABAQUS_UMAT_SS_BCJ_iso-damage")
		mat = new ABAQUS_SS_BCJ_ISO;
#endif
#ifdef ABAQUS_DEV
	else if (name == "ABAQUS_UMAT_IsoPlast")
		mat = new ABAQUS_UMAT_IsoPlast;
#endif
#endif /* ABAQUS_MATERIAL */

	/* set support */
	if (mat) mat->SetSSMatSupport(fSSMatSupport);

	return mat;
}
