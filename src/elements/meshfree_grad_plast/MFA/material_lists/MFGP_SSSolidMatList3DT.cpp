/* $Id: MFGP_SSSolidMatList3DT.cpp,v 1.5 2005-01-06 22:54:25 kyonten Exp $ */
#include "MFGP_SSSolidMatList3DT.h"
//#include "SSMatSupportT.h"
#include "DevelopmentMaterialsConfig.h"
#include "MFGP_MaterialSupportT.h"
#include "GRAD_MRSSKStV.h"

using namespace Tahoe;

/* constructors */
MFGP_SSSolidMatList3DT::MFGP_SSSolidMatList3DT(int length, const MFGP_MaterialSupportT& support):
	MFGP_SolidMatListT(length, support),
	fMFGPMaterialSupport(&support)
{
	SetName("mfgp_small_strain_material_3D");
}

MFGP_SSSolidMatList3DT::MFGP_SSSolidMatList3DT(void):
	fMFGPMaterialSupport(NULL)
{
	SetName("mfgp_small_strain_material_3D");
}

/* information about subordinate parameter lists */
void MFGP_SSSolidMatList3DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	MFGP_SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("mfgp_ss_material_list_3D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void MFGP_SSSolidMatList3DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "mfgp_ss_material_list_3D")
	{
		sub_lists.AddSub("small_strain_StVenant_MR_grad");
	}
	else /* inherited */
		MFGP_SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MFGP_SSSolidMatList3DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	MFGP_SSSolidMatT* ss_solid_mat = NewSSSolidMat(name);
	if (ss_solid_mat)
		return ss_solid_mat;
	else /* inherited */
		return MFGP_SolidMatListT::NewSub(name);
}

/* accept parameter list */
void MFGP_SSSolidMatList3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	MFGP_SolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	const ArrayT<ParameterListT>& subs = list.Lists();
	int count = 0;
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		MFGP_SSSolidMatT* mat = NewSSSolidMat(sub.Name());
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
MFGP_SSSolidMatT* MFGP_SSSolidMatList3DT::NewSSSolidMat(const StringT& name) const
{
	MFGP_SSSolidMatT* mat = NULL;

	if (name == "small_strain_StVenant_MR_grad")
		mat = new GRAD_MRSSKStV;

	/* set support */
	if (mat) mat->SetMFGPMaterialSupport(fMFGPMaterialSupport);

	return mat;
}
