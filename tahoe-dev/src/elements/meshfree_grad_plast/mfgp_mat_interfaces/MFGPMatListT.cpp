/* $Id: MFGPMatListT.cpp,v 1.1 2005-04-26 22:27:01 kyonten Exp $ */
#include "MFGPMatListT.h"
#include "MFGPMaterialT.h"

using namespace Tahoe;

/* constructors */
MFGPMatListT::MFGPMatListT(int length):
	pArrayT<MFGPMaterialT*>(length),
	ParameterInterfaceT("mfgp_material_list"),
	fHasHistory(false)
{

}

MFGPMatListT::MFGPMatListT(void):
	ParameterInterfaceT("mfgp_material_list"),
	fHasHistory(false)
{

}

/* apply pre-conditions at the current time step */
void MFGPMatListT::InitStep(void)
{
	/* loop over list */
	for (int i = 0; i < fLength; i++)
		fArray[i]->InitStep();
}

/* finalize the current time step */
void MFGPMatListT::CloseStep(void)
{
	/* loop over list */
	for (int i = 0; i < fLength; i++)
		fArray[i]->CloseStep();
}
