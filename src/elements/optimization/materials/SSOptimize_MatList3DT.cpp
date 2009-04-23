/* $Id: SSOptimize_MatList3DT.cpp,v 1.1 2009-04-23 03:03:51 thao Exp $ */
#include "SSOptimize_MatList3DT.h"
#include "SSMatSupportT.h"
#include "SSOptimize_MatT.h"

#include "SolidMaterialsConfig.h"
#include "DevelopmentMaterialsConfig.h"
#include "DevelopmentElementsConfig.h"


/* development module materials require solid element development to be enabled */
#include "SSKStV_Optimize.h";

#ifdef VISCOELASTICITY
#include "SSVisco_Optimize.h";
#include "SSVE_Opto_test.h";
#endif

using namespace Tahoe;

/* constructor */
SSOptimize_MatList3DT::SSOptimize_MatList3DT(int length, SSMatSupportT& support):
	SSSolidMatList3DT(length, support),
	fSupport(&support)
{
	SetName("ss_opto_material");
	if (fSupport->NumSD() != 3)
		ExceptionT::GeneralFail("SSOptimize_MatList3DT::SSOptimize_MatList3DT");

}

SSOptimize_MatList3DT::SSOptimize_MatList3DT(void):
	fSupport(NULL)
{
	SetName("ss_opto_material");
}


void SSOptimize_MatList3DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("ss_opto_mat_list_3D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void SSOptimize_MatList3DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "ss_opto_mat_list_3D")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("ss_optimize_StVenant");
#ifdef VISCOELASTICITY
		sub_lists.AddSub("linear_prony_optimize_ve");
		sub_lists.AddSub("ssve_opto_test");
#endif		
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSOptimize_MatList3DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	SSSolidMatT* ss_adj_mat = NewMat(name);
	if (ss_adj_mat)
		return ss_adj_mat;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void SSOptimize_MatList3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	const ArrayT<ParameterListT>& subs = list.Lists();
	int count = 0;
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		SSSolidMatT* mat = NewMat(sub.Name());
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
SSSolidMatT* SSOptimize_MatList3DT::NewMat(const StringT& name) const
{
	const char caller[] = "SSOptimize_MatList3DT::NewMat";
	
	SSSolidMatT* mat = NULL;
	SSOptimize_MatT* mat_opto = NULL;

	if (name == "ss_optimize_StVenant")
	{
		mat = new SSKStV_Optimize;
		mat_opto = TB_DYNAMIC_CAST(SSOptimize_MatT*, mat);
	}
#ifdef VISCOELASTICITY
	else if (name == "linear_prony_optimize_ve")
	{
		mat = new SSVisco_Optimize;
		mat_opto = TB_DYNAMIC_CAST(SSOptimize_MatT*, mat);		
	}
	else if (name == "ssve_opto_test")
	{
		mat = new SSVE_Opto_test;
		mat_opto = TB_DYNAMIC_CAST(SSOptimize_MatT*, mat);		
	}
#endif
	/* set support */
	if (mat_opto) {
		mat_opto->SetSupport(fSupport);
	}
//	else ExceptionT::GeneralFail(caller, "could not set support for SSOptimize_MatT");
	
	if (mat){
		mat->SetSSMatSupport(fSupport);
	}
	return mat;
}
