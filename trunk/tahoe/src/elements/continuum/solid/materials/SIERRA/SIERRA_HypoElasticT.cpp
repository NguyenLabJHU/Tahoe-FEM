/* $Id: SIERRA_HypoElasticT.cpp,v 1.3 2004-07-27 03:16:04 paklein Exp $ */
#include "SIERRA_HypoElasticT.h"

using namespace Tahoe;

/* registration */
extern "C" {
void SIERRA_HypoElastic_reg(void);
}

/* constructor */
SIERRA_HypoElasticT::SIERRA_HypoElasticT(void):
	ParameterInterfaceT("SIERRA_hypoelastic")
{

}

/* returns the strain energy density for the specified strain */
double SIERRA_HypoElasticT::StrainEnergyDensity(void)
{
	/* load stored data */
	Load(CurrentElement(), CurrIP());

	return fstate_new[kStrainEnergyDensity];
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* call the SIERRA registration function */
void SIERRA_HypoElasticT::Register_SIERRA_Material(void) const
{
	/* call C function */
	SIERRA_HypoElastic_reg();
}

void SIERRA_HypoElasticT::SetOutputVariables(iArrayT& variable_index, 
	ArrayT<StringT>& output_labels) const
{
	/* no additional material outputs */
	variable_index.Dimension(0);
	output_labels.Dimension(0);
}
