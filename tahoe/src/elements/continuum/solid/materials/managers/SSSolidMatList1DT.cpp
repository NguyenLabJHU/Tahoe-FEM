/* $Id: SSSolidMatList1DT.cpp,v 1.1.6.4 2004-07-15 06:25:36 paklein Exp $ */
#include "SSSolidMatList1DT.h"
#include "SSMatSupportT.h"


/* 1D material types codes */
/* Add small strain linear elastic material here */
#include "SSHookean1D.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#endif

#ifdef GRAD_SMALL_STRAIN_DEV
#include "GradJ2SS1D.h"
#include "J2SSKStV1D.h"
#include "GradSSMatSupportT.h"
#endif

using namespace Tahoe;

/* constructor */
SSSolidMatList1DT::SSSolidMatList1DT(int length, const SSMatSupportT& support):
	SolidMatListT(length, support),
	fSSMatSupport(&support),
	fGradSSMatSupport(NULL)
{
	SetName("small_strain_material_1D");

#ifdef GRAD_SMALL_STRAIN_DEV
	/* cast to gradient enhanced small strain support */
	fGradSSMatSupport = TB_DYNAMIC_CAST(const GradSSMatSupportT*, fSSMatSupport);
#endif
}

SSSolidMatList1DT::SSSolidMatList1DT(void):
	fSSMatSupport(NULL),
	fGradSSMatSupport(NULL)	
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
SSSolidMatT* SSSolidMatList1DT::NewSSSolidMat(const StringT& name) const
{
	SSSolidMatT* mat = NULL;

	if (name == "linear_material_1D")
		mat = new SSHookean1D;

	/* set support */
	if (mat) mat->SetSSMatSupport(fSSMatSupport);

	return mat;
}

