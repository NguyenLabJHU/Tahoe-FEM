/* $Id: GradSSSolidMatList3DT.cpp,v 1.1 2004-09-02 18:25:04 rdorgan Exp $ */
#include "GradSSSolidMatList3DT.h"
#include "GradSSMatSupportT.h"
#include "SolidMaterialsConfig.h"

/* gradient enhanced material type codes */
#include "GradJ2SSKStV.h"

using namespace Tahoe;

/* constructor */
GradSSSolidMatList3DT::GradSSSolidMatList3DT(int length, const GradSSMatSupportT& support):
	SSSolidMatList3DT(length, support),
	fGradSSMatSupport(&support)
{
	SetName("grad_small_strain_material_3D");
}

GradSSSolidMatList3DT::GradSSSolidMatList3DT(void):
	fGradSSMatSupport(NULL)	
{
	SetName("grad_small_strain_material_3D");
}

/* information about subordinate parameter lists */
void GradSSSolidMatList3DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("grad_ss_material_list_3D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void GradSSSolidMatList3DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "grad_ss_material_list_3D")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("grad_small_strain_StVenant_J2");
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* GradSSSolidMatList3DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	GradSSSolidMatT* grad_ss_solid_mat = NewGradSSSolidMat(name);
	if (grad_ss_solid_mat)
		return grad_ss_solid_mat;
	else /* inherited */
		return SSSolidMatList3DT::NewSub(name);
}

/* accept parameter list */
void GradSSSolidMatList3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSSolidMatList3DT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	AutoArrayT<GradSSSolidMatT*> materials;
	const ArrayT<ParameterListT>& subs = list.Lists();
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		GradSSSolidMatT* mat = NewGradSSSolidMat(sub.Name());
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
GradSSSolidMatT* GradSSSolidMatList3DT::NewGradSSSolidMat(const StringT& name) const
{
	GradSSSolidMatT* mat = NULL;

	if (name == "grad_small_strain_StVenant_J2")
		mat = new GradJ2SSKStV;

	/* set support */
	if (mat) mat->SetGradSSMatSupport(fGradSSMatSupport);

	return mat;
}

