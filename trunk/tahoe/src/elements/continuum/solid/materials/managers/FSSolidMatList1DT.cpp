/* $Id: FSSolidMatList1DT.cpp,v 1.3 2004-07-16 16:22:51 paklein Exp $ */
#include "FSSolidMatList1DT.h"
#include "FSMatSupportT.h"

#include "SolidMaterialsConfig.h"
#include "FSSolidMatT.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#endif

#ifdef CAUCHY_BORN_MATERIAL
#include "Chain1D.h"
#endif

using namespace Tahoe;

/* constructor */
FSSolidMatList1DT::FSSolidMatList1DT(int length, const FSMatSupportT& support):
	SolidMatListT(length, support),
	fFSMatSupport(&support)
{
	SetName("large_strain_material_1D");
}

FSSolidMatList1DT::FSSolidMatList1DT(void):
	fFSMatSupport(NULL)
{
	SetName("large_strain_material_1D");
}

/* information about subordinate parameter lists */
void FSSolidMatList1DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 1D materials */
	sub_list.AddSub("fs_material_list_1D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FSSolidMatList1DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "fs_material_list_1D")
	{
		order = ParameterListT::Choice;

#ifdef CAUCHY_BORN_MATERIAL
		sub_lists.AddSub("chain_1D");
#endif
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSSolidMatList1DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	FSSolidMatT* fs_solid_mat = NewFSSolidMat(name);
	if (fs_solid_mat)
		return fs_solid_mat;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void FSSolidMatList1DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	AutoArrayT<FSSolidMatT*> materials;
	const ArrayT<ParameterListT>& subs = list.Lists();
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		FSSolidMatT* mat = NewFSSolidMat(sub.Name());
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
FSSolidMatT* FSSolidMatList1DT::NewFSSolidMat(const StringT& name) const
{
	FSSolidMatT* mat = NULL;

	if (false) /* dummy */
		mat = NULL;

#ifdef CAUCHY_BORN_MATERIAL
	else if (name == "chain_1D")
		mat = new Chain1D;
#endif

	/* set support */
	if (mat) mat->SetFSMatSupport(fFSMatSupport);

	return mat;

}
