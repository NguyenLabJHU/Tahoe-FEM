/* $Id: SSSolidMatList3DT.cpp,v 1.3 2004-07-21 20:51:14 raregue Exp $ */
#include "SSSolidMatList3DT.h"
#include "SSMatSupportT.h"


#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#endif

#include "SSKStV.h"
#include "SSCubicT.h"

#ifdef VISCOELASTICITY
#include "SSLinearVE3D.h"
#endif

#ifdef VISCOELASTIC_MATERIALS_DEV
#include "SSSV_KStV3D.h"
#endif

#ifdef J2PLASTICITY_MATERIALS_DEV
#include "SSJ2LinHardT.h"
#endif

#ifdef PLASTICITY_J2_MATERIAL
#include "J2SSKStV.h"
#include "LocalJ2SSNonlinHard.h"
#include "GradJ2SSNonlinHard.h"
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
#include "ABAQUS_SS_BCJ_ISO.h"
#endif
#endif

#ifdef PLASTICITY_DP_MATERIAL
#include "DPSSKStV.h"
#endif

#ifdef PLASTICITY_DP_LOC_MATERIAL_DEV
#include "DPSSKStVLoc.h"
#endif

#ifdef FOSSUM_MATERIAL_DEV
#include "FossumSSIsoT.h"
#endif

#ifdef PLASTICITY_MR_MATERIAL_DEV
#include "MRSSKStV.h"
#endif

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
	
		sub_lists.AddSub("small_strain_cubic");
		sub_lists.AddSub("small_strain_StVenant");

#ifdef PLASTICITY_J2_MATERIAL
		sub_lists.AddSub("small_strain_StVenant_J2");
#endif

#ifdef PLASTICITY_DP_MATERIAL
		sub_lists.AddSub("small_strain_StVenant_DP");
#endif

#ifdef PLASTICITY_DP_LOC_MATERIAL_DEV
		sub_lists.AddSub("small_strain_StVenant_DP_Loc");
#endif

#ifdef VISCOELASTICITY
		sub_lists.AddSub("linear_viscoelastic");
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
SSSolidMatT* SSSolidMatList3DT::NewSSSolidMat(const StringT& name) const
{
	SSSolidMatT* mat = NULL;

	if (name == "small_strain_cubic")
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

#ifdef PLASTICITY_DP_LOC_MATERIAL_DEV
	else if (name == "small_strain_StVenant_DP_Loc")
		mat = new DPSSKStVLoc;
#endif

#ifdef VISCOELASTICITY
	else if (name == "linear_viscoelastic")
		mat = new SSLinearVE3D;
#endif

	/* set support */
	if (mat) mat->SetSSMatSupport(fSSMatSupport);

	return mat;
}
