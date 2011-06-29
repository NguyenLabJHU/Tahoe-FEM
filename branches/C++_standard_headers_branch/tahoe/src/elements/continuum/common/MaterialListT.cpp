/* $Id: MaterialListT.cpp,v 1.6 2004-07-15 08:26:13 paklein Exp $ */
/* created: paklein (02/16/1997) */
#include "MaterialListT.h"
#include "ContinuumMaterialT.h"

using namespace Tahoe;

/* constructors */
MaterialListT::MaterialListT(int length):
	pArrayT<ContinuumMaterialT*>(length),
	ParameterInterfaceT("material_list"),
	fHasHistory(false)
{

}

MaterialListT::MaterialListT(void):
	ParameterInterfaceT("material_list"),
	fHasHistory(false)
{

}

/* apply pre-conditions at the current time step */
void MaterialListT::InitStep(void)
{
	/* loop over list */
	for (int i = 0; i < fLength; i++)
		fArray[i]->InitStep();
}

/* finalize the current time step */
void MaterialListT::CloseStep(void)
{
	/* loop over list */
	for (int i = 0; i < fLength; i++)
		fArray[i]->CloseStep();
}
