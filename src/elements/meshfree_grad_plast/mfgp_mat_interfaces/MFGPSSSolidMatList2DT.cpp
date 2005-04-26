/* $Id: MFGPSSSolidMatList2DT.cpp,v 1.1 2005-04-26 22:27:01 kyonten Exp $ */
#include "MFGPSSSolidMatList2DT.h"
#include "MFGPMatSupportT.h"

#include "DevelopmentMaterialsConfig.h"

#include "SSHookeanMat2DT.h"
#include "SSKStV2D.h"
#include "SSCubic2DT.h"

#include "GRAD_MRSSKStV2D.h"


using namespace Tahoe;

/* constructor */
MFGPSSSolidMatList2DT::MFGPSSSolidMatList2DT(int length, const MFGPMatSupportT& support):
	MFGPSolidMatListT(length, support),
	fMFGPMatSupport(&support)
{
	SetName("mfgp_material_2D");
	if (fMFGPMatSupport->NumSD() != 2)
		ExceptionT::GeneralFail("MFGPSSSolidMatList2DT::MFGPSSSolidMatList2DT");

#ifdef __NO_RTTI__
	cout << "\n MFGPSSSolidMatList2DT::MFGPSSSolidMatList2DT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;
#endif

}

MFGPSSSolidMatList2DT::MFGPSSSolidMatList2DT(void):
	fMFGPMatSupport(NULL)	
{
	SetName("mfgp_material_2D");

#ifdef __NO_RTTI__
	cout << "\n MFGPSSSolidMatList2DT::MFGPSSSolidMatList2DT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;
#endif
}

/* return true if the list contains plane stress models */
bool MFGPSSSolidMatList2DT::HasPlaneStress(void) const
{
	/* check materials */
	for (int i = 0; i < Length(); i++)
	{
		/* get pointer to Material2DT */
		const MFGPMaterialT* cont_mat = fArray[i];
		
		/* assume materials that don't have Material2DT are plane strain */
		if (cont_mat && cont_mat->Constraint() == MFGPMaterialT::kPlaneStress) 
			return true;
	}
	return false;
}

/* information about subordinate parameter lists */
void MFGPSSSolidMatList2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	MFGPSolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("mfgp_material_list_2D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void MFGPSSSolidMatList2DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "mfgp_material_list_2D")
	{
		order = ParameterListT::Choice;
		/*sub_lists.AddSub("small_strain_Hookean_2D");
		sub_lists.AddSub("small_strain_cubic_2D");
		sub_lists.AddSub("small_strain_StVenant_2D");*/
		
		sub_lists.AddSub("small_strain_StVenant_MR_grad_2D");

	}
	else /* inherited */
		MFGPSolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MFGPSSSolidMatList2DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	MFGPSSSolidMatT* ss_solid_mat = NewMFGPSSSolidMat(name);
	if (ss_solid_mat)
		return ss_solid_mat;
	else /* inherited */
		return MFGPSolidMatListT::NewSub(name);
}

/* accept parameter list */
void MFGPSSSolidMatList2DT::TakeParameterList(const ParameterListT& list)
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

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct the specified material or NULL if the request cannot be completed */
MFGPSSSolidMatT* MFGPSSSolidMatList2DT::NewMFGPSSSolidMat(const StringT& name) const
{
	MFGPSSSolidMatT* mat = NULL;
		
	/*if (name == "small_strain_Hookean_2D")
		mat = new SSHookeanMat2DT;
	else if (name == "small_strain_cubic_2D")
		mat = new SSCubic2DT;
	else if (name == "small_strain_StVenant_2D")
		mat = new SSKStV2D;*/

	if (name == "small_strain_StVenant_MR_grad_2D")
		mat = new GRAD_MRSSKStV2D;

	/* set support */
	if (mat) mat->SetMFGPMatSupport(fMFGPMatSupport);

	return mat;
}
