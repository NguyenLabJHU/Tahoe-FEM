/* $Id: MFGP_SSSolidMatList2DT.cpp,v 1.5 2005-01-06 22:54:25 kyonten Exp $ */
#include "MFGP_SSSolidMatList2DT.h"
//#include "SSMatSupportT.h"
#include "MFGP_MaterialSupportT.h"
#include "DevelopmentMaterialsConfig.h"
#include "GRAD_MRSSKStV2D.h"

using namespace Tahoe;

/* constructor */
//MFGP_SSSolidMatList2DT::MFGP_SSSolidMatList2DT(int length, const SSMatSupportT& support):
MFGP_SSSolidMatList2DT::MFGP_SSSolidMatList2DT(int length, const MFGP_MaterialSupportT& support):
	MFGP_SolidMatListT(length,support),
	fMFGPMaterialSupport(&support)
{
	SetName("mfgp_small_strain_material_2D");
	if (fMFGPMaterialSupport->NumSD() != 2)
		ExceptionT::GeneralFail("MFGP_SSSolidMatList2DT::MFGP_SSSolidMatList2DT");

#ifdef __NO_RTTI__
	cout << "\n MFGP_SSSolidMatList2DT::MFGP_SSSolidMatList2DT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;
#endif

}

MFGP_SSSolidMatList2DT::MFGP_SSSolidMatList2DT(void):
	fMFGPMaterialSupport(NULL)	
{
	SetName("mfgp_small_strain_material_2D");

#ifdef __NO_RTTI__
	cout << "\n MFGP_SSSolidMatList2DT::MFGP_SSSolidMatList2DT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;
#endif
}

/* return true if the list contains plane stress models */
bool MFGP_SSSolidMatList2DT::HasPlaneStress(void) const
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
void MFGP_SSSolidMatList2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	MFGP_SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("mfgp_ss_material_list_2D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void MFGP_SSSolidMatList2DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "mfgp_ss_material_list_2D")
	{
		order = ParameterListT::Choice;
		sub_lists.AddSub("small_strain_StVenant_MR_grad_2D");
	}
	else /* inherited */
		MFGP_SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MFGP_SSSolidMatList2DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	MFGP_SSSolidMatT* ss_solid_mat = NewSSSolidMat(name);
	if (ss_solid_mat)
		return ss_solid_mat;
	else /* inherited */
		return MFGP_SolidMatListT::NewSub(name);
}

/* accept parameter list */
void MFGP_SSSolidMatList2DT::TakeParameterList(const ParameterListT& list)
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

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct the specified material or NULL if the request cannot be completed */
MFGP_SSSolidMatT* MFGP_SSSolidMatList2DT::NewSSSolidMat(const StringT& name) const
{
	MFGP_SSSolidMatT* mat = NULL;

	if (name == "small_strain_StVenant_MR_grad_2D")
		mat = new GRAD_MRSSKStV2D;

	/* set support */
	if (mat) mat->SetMFGPMaterialSupport(fMFGPMaterialSupport);

	return mat;
}
