/* $Id: MixtureSpeciesT.cpp,v 1.1 2004-11-05 22:53:49 paklein Exp $ */
#include "MixtureSpeciesT.h"
#include "UpdatedLagMixtureT.h"
#include "ShapeFunctionT.h"

using namespace Tahoe;

/* constructor */
MixtureSpeciesT::MixtureSpeciesT(const ElementSupportT& support):
	NLDiffusionElementT(support),
	fUpdatedLagMixture(NULL),
	fIndex(-1)
{
	SetName("mixture_species");
}

/* describe the parameters needed by the interface */
void MixtureSpeciesT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NLDiffusionElementT::DefineParameters(list);

	/* associated solid element group */
	list.AddParameter(ParameterT::Integer, "solid_element_group");

#if 0
	/* velocity of species is calculated wrt this reference frame */
	ParameterT frame(ParameterT::Word, "reference_frame");
	frame.SetDefault("global");
	species->AddParameter(frame);
#endif
}

/* accept parameter list */
void MixtureSpeciesT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MixtureSpeciesT::TakeParameterList";

	/* inherited */
	NLDiffusionElementT::TakeParameterList(list);

	/* resolve background solid element group */
	int solid_element_group = list.GetParameter("solid_element_group");
	solid_element_group--;
	ElementBaseT& element = ElementSupport().ElementGroup(solid_element_group);	
	fUpdatedLagMixture = TB_DYNAMIC_CAST(UpdatedLagMixtureT*, &element);
	if (!fUpdatedLagMixture)
		ExceptionT::GeneralFail(caller, "group %d is not a mixture", solid_element_group+1);
	
	/* resolve species index */
	const ArrayT<StringT>& labels = Field().Labels();
	if (labels.Length() != 1) 
		ExceptionT::GeneralFail(caller, "field dim 1 not %d", labels.Length());
	fIndex = fUpdatedLagMixture->SpeciesIndex(labels[0]);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form the residual force vector */
void MixtureSpeciesT::RHSDriver(void)
{
	/* project partial stresses to the nodes */
	fUpdatedLagMixture->ProjectPartialStress(fIndex);

	/* inherited */
	NLDiffusionElementT::RHSDriver();
}

/* calculate the internal force contribution ("-k*d") */
void MixtureSpeciesT::FormKd(double constK)
{
/* include effects from:
 * (1) source/sink terms
 * (2) divergence of flux
 */

	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	int nsd = NumSD();
	dMatrixT grad;
	dArrayT field;
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* set field gradient */
		grad.Set(1, nsd, fGradient_list[CurrIP()].Pointer());
		IP_ComputeGradient(fLocDisp, grad);
		
		/* interpolate field */
		field.Set(1, fField_list.Pointer(CurrIP()));
		IP_Interpolate(fLocDisp, field);

		/* get strain-displacement matrix */
		B(fShapes->CurrIP(), fB);

		/* compute heat flow */
//		fB.MultTx(fCurrMaterial->q_i(), fNEEvec);

		/* accumulate */
		fRHS.AddScaled(-constK*(*Weight++)*(*Det++), fNEEvec);
	}	
}
