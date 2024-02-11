/* $Id: GradSSSolidMatList2DT.cpp,v 1.1 2004-09-02 18:25:04 rdorgan Exp $ */
#include "GradSSSolidMatList2DT.h"
#include "GradSSMatSupportT.h"
#include "SolidMaterialsConfig.h"

#include "GradJ2SSKStV2D.h"

using namespace Tahoe;

/* constructor */
GradSSSolidMatList2DT::GradSSSolidMatList2DT(int length, const GradSSMatSupportT& support):
	SSSolidMatList2DT(length, support),
	fGradSSMatSupport(&support)
{
	SetName("grad_small_strain_material_2D");

	if (fGradSSMatSupport->NumSD() != 2)
		ExceptionT::GeneralFail("GradSSSolidMatList2DT::GradSSSolidMatList2DT");
}

GradSSSolidMatList2DT::GradSSSolidMatList2DT(void):
	fGradSSMatSupport(NULL)
{
	SetName("grad_small_strain_material_2D");
}

/* return true if the list contains plane stress models */
bool GradSSSolidMatList2DT::HasPlaneStress(void) const
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
void GradSSSolidMatList2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("grad_ss_material_list_2D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void GradSSSolidMatList2DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "grad_ss_material_list_2D")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("grad_small_strain_StVenant_J2_2D");
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* GradSSSolidMatList2DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	GradSSSolidMatT* grad_ss_solid_mat = NewGradSSSolidMat(name);
	if (grad_ss_solid_mat)
		return grad_ss_solid_mat;
	else /* inherited */
		return SSSolidMatList2DT::NewSub(name);
}

/* accept parameter list */
void GradSSSolidMatList2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSSolidMatList2DT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	const ArrayT<ParameterListT>& subs = list.Lists();
	int count = 0;
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		GradSSSolidMatT* mat = NewGradSSSolidMat(sub.Name());
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
GradSSSolidMatT* GradSSSolidMatList2DT::NewGradSSSolidMat(const StringT& name) const
{
	GradSSSolidMatT* mat = NULL;

	if (name == "grad_small_strain_StVenant_J2_2D")
		mat = new GradJ2SSKStV2D;

	/* set support */
	if (mat) mat->SetGradSSMatSupport(fGradSSMatSupport);

	return mat;
}
