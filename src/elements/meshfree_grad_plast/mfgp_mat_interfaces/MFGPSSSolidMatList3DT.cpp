/* $Id: MFGPSSSolidMatList3DT.cpp,v 1.1 2005-04-26 22:27:01 kyonten Exp $ */
#include "MFGPSSSolidMatList3DT.h"
#include "MFGPMatSupportT.h"

#include "DevelopmentMaterialsConfig.h"

#include "SSKStV.h"
#include "SSCubicT.h"

#ifdef GRAD_PLASTICITY_MR_MATERIAL_DEV
#include "GRAD_MRSSKStV.h"
#endif

using namespace Tahoe;

/* constructors */
MFGPSSSolidMatList3DT::MFGPSSSolidMatList3DT(int length, const MFGPMatSupportT& support):
	MFGPSolidMatListT(length, support),
	fMFGPMatSupport(&support)
{
	SetName("mfgp_material_3D");
}

MFGPSSSolidMatList3DT::MFGPSSSolidMatList3DT(void):
	fMFGPMatSupport(NULL)
{
	SetName("mfgp_material_3D");
}

/* information about subordinate parameter lists */
void MFGPSSSolidMatList3DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	MFGPSolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("mfgp_material_list_3D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void MFGPSSSolidMatList3DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "mfgp_material_list_3D")
	{

#ifdef GRAD_PLASTICITY_MR_MATERIAL_DEV
		sub_lists.AddSub("small_strain_StVenant_MR_grad_3D");
#endif

	}
	else /* inherited */
		MFGPSolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MFGPSSSolidMatList3DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	MFGPSSSolidMatT* ss_solid_mat = NewMFGPSSSolidMat(name);
	if (ss_solid_mat)
		return ss_solid_mat;
	else /* inherited */
		return MFGPSolidMatListT::NewSub(name);
}

/* accept parameter list */
void MFGPSSSolidMatList3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	MFGPSolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	const ArrayT<ParameterListT>& subs = list.Lists();
	int count = 0;
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		MFGPSSSolidMatT* mat = NewMFGPSSSolidMat(sub.Name());
		if (mat) {

			/* store pointer */
			(*this)[count++] = mat;

			/* initialize material */
			mat->TakeParameterList(sub);

			/* set flags */
			if (mat->HasHistory()) fHasHistory = true;	
		}
	}
}

/* construct the specified material or NULL if the request cannot be completed */
MFGPSSSolidMatT* MFGPSSSolidMatList3DT::NewMFGPSSSolidMat(const StringT& name) const
{
	MFGPSSSolidMatT* mat = NULL;

#ifdef GRAD_PLASTICITY_MR_MATERIAL_DEV
	if (name == "small_strain_StVenant_MR_grad_3D")
		mat = new GRAD_MRSSKStV;
#endif

	/* set support */
	if (mat) mat->SetMFGPMatSupport(fMFGPMatSupport);

	return mat;
}
