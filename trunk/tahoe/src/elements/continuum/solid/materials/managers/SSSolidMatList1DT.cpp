/* $Id: SSSolidMatList1DT.cpp,v 1.4 2004-08-08 02:06:32 paklein Exp $ */
#include "SSSolidMatList1DT.h"
#include "SSMatSupportT.h"


/* 1D material types codes */
/* Add small strain linear elastic material here */
#include "SSHookean1D.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#endif

#ifdef GRAD_SMALL_STRAIN_DEV
#include "J2SSKStV1D.h"
#endif

using namespace Tahoe;

/* constructor */
SSSolidMatList1DT::SSSolidMatList1DT(int length, const SSMatSupportT& support):
	SolidMatListT(length, support),
	fSSMatSupport(&support)
{
	SetName("small_strain_material_1D");
}

SSSolidMatList1DT::SSSolidMatList1DT(void):
	fSSMatSupport(NULL)
{
	SetName("small_strain_material_1D");
}

/* information about subordinate parameter lists */
void SSSolidMatList1DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("ss_material_list_1D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void SSSolidMatList1DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "ss_material_list_1D")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("linear_material_1D");

#ifdef GRAD_SMALL_STRAIN_DEV
		sub_lists.AddSub("small_strain_StVenant_J2_1D");
#endif
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSSolidMatList1DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	SSSolidMatT* ss_solid_mat = NewSSSolidMat(name);
	if (ss_solid_mat)
		return ss_solid_mat;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void SSSolidMatList1DT::TakeParameterList(const ParameterListT& list)
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
SSSolidMatT* SSSolidMatList1DT::NewSSSolidMat(const StringT& name) const
{
	SSSolidMatT* mat = NULL;

	if (name == "linear_material_1D")
		mat = new SSHookean1D;

#ifdef GRAD_SMALL_STRAIN_DEV
	else if (name == "small_strain_StVenant_J2_1D")
		mat = new J2SSKStV1D;
#endif

	/* set support */
	if (mat) mat->SetSSMatSupport(fSSMatSupport);

	return mat;
}

