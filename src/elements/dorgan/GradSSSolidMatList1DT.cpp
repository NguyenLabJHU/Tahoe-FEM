/* $Id: GradSSSolidMatList1DT.cpp,v 1.1 2004-07-20 23:41:03 rdorgan Exp $ */
#include "GradSSSolidMatList1DT.h"
#include "GradSSMatSupportT.h"

/* 1D gradient enhanced material type codes */
#include "GradJ2SSKStV1D.h"

using namespace Tahoe;

/* constructor */
GradSSSolidMatList1DT::GradSSSolidMatList1DT(int length, const GradSSMatSupportT& support):
	SSSolidMatList1DT(length, support),
	fGradSSMatSupport(&support)
{
	SetName("grad_small_strain_material_1D");
}

GradSSSolidMatList1DT::GradSSSolidMatList1DT(void):
	fGradSSMatSupport(NULL)	
{
	SetName("grad_small_strain_material_1D");
}

/* information about subordinate parameter lists */
void GradSSSolidMatList1DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("grad_ss_material_list_1D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void GradSSSolidMatList1DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "grad_ss_material_list_1D")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("grad_small_strain_StVenant_J2_1D");
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* GradSSSolidMatList1DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	GradSSSolidMatT* grad_ss_solid_mat = NewGradSSSolidMat(name);
	if (grad_ss_solid_mat)
		return grad_ss_solid_mat;
	else /* inherited */
		return SSSolidMatList1DT::NewSub(name);
}

/* accept parameter list */
void GradSSSolidMatList1DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSSolidMatList1DT::TakeParameterList(list);

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
GradSSSolidMatT* GradSSSolidMatList1DT::NewGradSSSolidMat(const StringT& name) const
{
	GradSSSolidMatT* mat = NULL;

	if (name == "grad_small_strain_StVenant_J2_1D")
		mat = new GradJ2SSKStV1D;

	/* set support */
	if (mat) mat->SetGradSSMatSupport(fGradSSMatSupport);

	return mat;
}

