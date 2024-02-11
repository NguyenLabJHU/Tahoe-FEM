/* $Id: MFGPSolidMatListT.cpp,v 1.1 2005-04-26 22:27:01 kyonten Exp $ */
#include "MFGPSolidMatListT.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#include "DevelopmentElementsConfig.h"
#endif

#include "MFGPMaterialT.h"
#include "MFGPMatSupportT.h"

using namespace Tahoe;

/* constructors */
MFGPSolidMatListT::MFGPSolidMatListT(int length, const MFGPMatSupportT& support):
	MFGPMatListT(length),
	fHasLocalizers(false),
	fMFGPMatSupport(&support)
{
#ifdef __NO_RTTI__
	cout << "\n MFGPSolidMatListT::MFGPSolidMatListT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;
#endif
}

MFGPSolidMatListT::MFGPSolidMatListT(void):
	fHasLocalizers(false),
	fMFGPMatSupport(NULL)
{

}