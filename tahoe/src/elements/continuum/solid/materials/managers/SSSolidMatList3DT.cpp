/* $Id: SSSolidMatList3DT.cpp,v 1.13 2006-07-27 13:37:42 kyonten Exp $ */
#include "SSSolidMatList3DT.h"
#include "SSMatSupportT.h"
#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#include "DevelopmentElementsConfig.h"
#endif

#include "SSKStV.h"
#include "SSCubicT.h"

#ifdef VISCOELASTICITY
#include "SSLinearVE3D.h"
#endif

#ifdef PLASTICITY_J2_MATERIAL
#include "J2SSKStV.h"
#include "LocalJ2SSNonlinHard.h"
#include "GradJ2SSNonlinHard.h"
#endif

#ifdef PLASTICITY_DP_MATERIAL
#include "DPSSKStV.h"
#endif

/* development module materials require solid element development to be enabled */
#ifdef SOLID_ELEMENT_DEV

#ifdef VISCOELASTIC_MATERIALS_DEV
#include "SSSV_KStV3D.h"
#endif

#ifdef J2PLASTICITY_MATERIALS_DEV
#include "SSJ2LinHardT.h"
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
#include "ABAQUS_SS_BCJ_ISO.h"
#endif
#ifdef ABAQUS_DEV
#include "ABAQUS_UMAT_IsoPlast.h"
#endif
#endif /* ABAQUS_MATERIAL */

#ifdef PLASTICITY_DP_LOC_MATERIAL_DEV
#include "DPSSKStVLoc.h"
#endif

#ifdef FOSSUM_MATERIAL_DEV
#include "FossumSSIsoT.h"
#endif

#ifdef PLASTICITY_MR_MATERIAL_DEV
#include "MRSSKStV.h"
#endif

#ifdef PLASTICITY_SMR_MATERIAL_DEV
#include "SMRSSKStV.h"
#endif

#endif /* SOLID_ELEMENT_DEV */

using namespace Tahoe;

/* constructors */
SSSolidMatList3DT::SSSolidMatList3DT(int length, const SSMatSupportT& support):
	SolidMatListT(length, support),
	fSSMatSupport(&support)
{
	SetName("small_strain_material_3D");
}

SSSolidMatList3DT::SSSolidMatList3DT(void):
	fSSMatSupport(NULL)
{
	SetName("small_strain_material_3D");
}

/* information about subordinate parameter lists */
void SSSolidMatList3DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("ss_material_list_3D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void SSSolidMatList3DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "ss_material_list_3D")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("small_strain_Hookean");
		sub_lists.AddSub("small_strain_cubic");
		sub_lists.AddSub("small_strain_StVenant");

#ifdef PLASTICITY_J2_MATERIAL
		sub_lists.AddSub("small_strain_StVenant_J2");
#endif

#ifdef PLASTICITY_DP_MATERIAL
		sub_lists.AddSub("small_strain_StVenant_DP");
#endif

#ifdef PLASTICITY_MR_MATERIAL_DEV
		sub_lists.AddSub("small_strain_StVenant_MR");
#endif

#ifdef PLASTICITY_SMR_MATERIAL_DEV
		sub_lists.AddSub("small_strain_StVenant_SMR");
#endif

#ifdef PLASTICITY_DP_LOC_MATERIAL_DEV
		sub_lists.AddSub("small_strain_StVenant_DP_Loc");
#endif

#ifdef VISCOELASTICITY
		sub_lists.AddSub("linear_viscoelastic");
#endif

#ifdef FOSSUM_MATERIAL_DEV
		sub_lists.AddSub("Fossum_small_strain");
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
ParameterInterfaceT* SSSolidMatList3DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	SSSolidMatT* ss_solid_mat = NewSSSolidMat(name);
	if (ss_solid_mat)
		return ss_solid_mat;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void SSSolidMatList3DT::TakeParameterList(const ParameterListT& list)
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

/* construct the specified material or NULL if the request cannot be completed */
SSSolidMatT* SSSolidMatList3DT::NewSSSolidMat(const StringT& name) const
{
	SSSolidMatT* mat = NULL;

	if (name == "small_strain_Hookean")
		mat = new SSHookeanMatT;
	else if (name == "small_strain_cubic")
		mat = new SSCubicT;
	else if (name == "small_strain_StVenant")
		mat = new SSKStV;

#ifdef PLASTICITY_J2_MATERIAL
	else if (name == "small_strain_StVenant_J2")
		mat = new J2SSKStV;
#endif

#ifdef PLASTICITY_DP_MATERIAL
	else if (name == "small_strain_StVenant_DP")
		mat = new DPSSKStV;
#endif

#ifdef PLASTICITY_MR_MATERIAL_DEV
	else if (name == "small_strain_StVenant_MR")
		mat = new MRSSKStV;
#endif

#ifdef PLASTICITY_SMR_MATERIAL_DEV
	else if (name == "small_strain_StVenant_SMR")
		mat = new SMRSSKStV;
#endif

#ifdef PLASTICITY_DP_LOC_MATERIAL_DEV
	else if (name == "small_strain_StVenant_DP_Loc")
		mat = new DPSSKStVLoc;
#endif

#ifdef VISCOELASTICITY
	else if (name == "linear_viscoelastic")
		mat = new SSLinearVE3D;
#endif

#ifdef FOSSUM_MATERIAL_DEV
	else if (name == "Fossum_small_strain")
	  mat = new FossumSSIsoT;
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
