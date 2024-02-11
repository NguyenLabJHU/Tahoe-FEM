/* $Id: SSOptimize_MatList2DT.cpp,v 1.1 2009-04-23 03:03:51 thao Exp $ */
#include "SSOptimize_MatList2DT.h"
#include "SSOptimize_MatSupportT.h"
#include "SSOptimize_MatT.h"
#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#include "DevelopmentElementsConfig.h"
#endif

/* development module materials require solid element development to be enabled */
//#ifdef SOLID_ELEMENT_DEV
//#include "SSKStV2D_OptAdj.h"

//#endif /* SOLID_ELEMENT_DEV */

using namespace Tahoe;

/* constructor */
SSOptimize_MatList2DT::SSOptimize_MatList2DT(int length, const SSOptimize_MatSupportT& support):
	SSSolidMatList2DT(length, support),
	fSSOptimize_MatSupport(&support)
{
	SetName("small_strain_optimize_material_2D");
	if (fSSOptimize_MatSupport->NumSD() != 2)
		ExceptionT::GeneralFail("SSOptimize_MatList2DT::SSOptimize_MatList2DT");

}

SSOptimize_MatList2DT::SSOptimize_MatList2DT(void):
	fSSOptimize_MatSupport(NULL)
{
	SetName("small_strain_optimize_material_2D");
}

/* return the description of the given inline subordinate parameter list */
void SSOptimize_MatList2DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "ss_material_list_2D")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("ss_optimize_StVenant_2D");
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSOptimize_MatList2DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	SSOptimize_MatT* ss_adj_mat = NewMat(name);
	if (ss_adj_mat)
		return ss_adj_mat;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void SSOptimize_MatList2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	const ArrayT<ParameterListT>& subs = list.Lists();
	int count = 0;
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		SSOptimize_MatT* mat = NewMat(sub.Name());
		if (mat) {
			
			/* store pointer */
			(*this)[count++] = mat;

			/* initialize material */
			mat->TakeParameterList(sub);

			/* set flags */
			if (mat->HasThermalStrain()) fHasThermal = true;
			if (mat->HasLocalization()) fHasLocalizers = true;
			if (mat->HasHistory()) fHasHistory = true;	
		}
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct the specified material or NULL if the request cannot be completed */
SSOptimize_MatT* SSOptimize_MatList2DT::NewMat(const StringT& name) const
{
	SSOptimize_MatT* mat = NULL;

/*	if (name == "ss_optimize_StVenant_2D")
		mat = new SSKStV2D_OptAdj;
*/
	/* set support */
	if (mat) mat->SetSupport(fSSOptimize_MatSupport);

	return mat;
}
