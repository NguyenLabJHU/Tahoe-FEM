/* $Id: FSFiberOptimize_MatListT.cpp,v 1.1 2009-04-23 03:03:50 thao Exp $ */
/* created: paklein (02/14/1997) */
#include "FSFiberOptimize_MatListT.h"
#include "FSFiberMatSupportT.h"
#include "FSFiberMatT.h"

#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#include "DevelopmentElementsConfig.h"
#endif

#include "AnisoCorneaVisco_Opt.h"

using namespace Tahoe;

/* constructors */
FSFiberOptimize_MatListT::FSFiberOptimize_MatListT(int length,  FSFiberMatSupportT& support):
	FSFiberMatListT(length, support),
	fSupport(&support)
{
	SetName("fiber_opto_material");
}

FSFiberOptimize_MatListT::FSFiberOptimize_MatListT(void):
	fSupport(NULL)
{
	SetName("fiber_opto_material");
}

/* information about subordinate parameter lists */
void FSFiberOptimize_MatListT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("fiber_opto_material_list", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FSFiberOptimize_MatListT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "fiber_opto_material_list")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("aniso_cornea_visco_opto");
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSFiberOptimize_MatListT::NewSub(const StringT& name) const
{
	/* try to construct material */
	FSFiberMatT* fs_fiber_mat = NewFSFiberMat(name);
	if (fs_fiber_mat)
		return fs_fiber_mat;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void FSFiberOptimize_MatListT::TakeParameterList(const ParameterListT& list)
{		
	/* inherited */
	SolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	const ArrayT<ParameterListT>& subs = list.Lists();
	int count = 0;
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		FSFiberMatT* mat = NewFSFiberMat(sub.Name());
		if (mat) {
			
			/* store pointer */
			(*this)[count++] = mat;

			/* initialize material */
			mat->TakeParameterList(sub);

			/* set flags */
			if (mat->HasHistory()) fHasHistory = true;	
			if (mat->HasThermalStrain()) fHasThermal = false;
			if (mat->HasLocalization()) fHasLocalizers = false;
		}
	}
}

/* construct the specified material or NULL if the request cannot be completed */
FSFiberMatT* FSFiberOptimize_MatListT::NewFSFiberMat(const StringT& name) const
{
	FSFiberMatT* mat = NULL;
	FSFiberOptimize_MatT* mat_opto = NULL;

	if  (name == "aniso_cornea_visco_opto")
	{
		mat = new AnisoCorneaVisco_Opt;
		mat_opto = TB_DYNAMIC_CAST(FSFiberOptimize_MatT*, mat);
	}
	
	/* set support */
	
	
	if (mat_opto) {
		mat_opto->SetSupport(fSupport);
	}

	if (mat) mat->SetFSFiberMatSupport(fSupport);
	
	return mat;

}
