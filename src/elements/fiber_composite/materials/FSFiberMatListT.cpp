/* $Id: FSFiberMatListT.cpp,v 1.6 2007-12-19 23:35:54 thao Exp $ */
/* created: paklein (02/14/1997) */
#include "FSFiberMatListT.h"
#include "FSFiberMatSupportT.h"

#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#include "DevelopmentElementsConfig.h"
#endif

#include "FSFiberMatT.h"
#include "AnisoFiber3D.h"
#include "AnisoCornea.h"
#include "AnisoCorneaVisco.h"
#include "AnisoCorneaIVisco.h"
#include "CFiberViscT.h"
#include "C2FiberViscT.h"

//#include "AnisoCorneaVisco2.h"

using namespace Tahoe;

/* constructors */
FSFiberMatListT::FSFiberMatListT(int length, const FSFiberMatSupportT& support):
	FSSolidMatList3DT(length, support),
	fFSFiberMatSupport(&support)
{
	SetName("fiber_comp_material");
}

FSFiberMatListT::FSFiberMatListT(void):
	fFSFiberMatSupport(NULL)
{
	SetName("fiber_comp_material");
}

/* information about subordinate parameter lists */
void FSFiberMatListT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("fiber_comp_material_list", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FSFiberMatListT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "fiber_comp_material_list")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("aniso_cornea");
		sub_lists.AddSub("aniso_fiber_3D");
		sub_lists.AddSub("aniso_viscoelastic_cornea");
		sub_lists.AddSub("aniso_scalar_visco_cornea");
/*		sub_lists.AddSub("aniso_viscoelastic2_cornea");*/
		sub_lists.AddSub("fiberI4_viscoelastic");
		sub_lists.AddSub("fiberI5_viscoelastic");
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSFiberMatListT::NewSub(const StringT& name) const
{
	/* try to construct material */
	FSFiberMatT* fs_fiber_mat = NewFSFiberMat(name);
	if (fs_fiber_mat)
		return fs_fiber_mat;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void FSFiberMatListT::TakeParameterList(const ParameterListT& list)
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
FSFiberMatT* FSFiberMatListT::NewFSFiberMat(const StringT& name) const
{
	FSFiberMatT* mat = NULL;

	if (name == "aniso_cornea")
		mat = new AnisoCornea;
	else if (name == "aniso_fiber_3D")
		mat = new AnisoFiber3D;
	else if (name == "aniso_viscoelastic_cornea")
		mat = new AnisoCorneaVisco;
	else if (name == "aniso_scalar_visco_cornea")
		mat = new AnisoCorneaIVisco;
	else if (name == "fiberI5_viscoelastic")
		mat = new C2FiberViscT;
	else if (name == "fiberI4_viscoelastic")
		mat = new CFiberViscT;

	/* set support */
	if (mat) mat->SetFSFiberMatSupport(fFSFiberMatSupport);

	return mat;

}
