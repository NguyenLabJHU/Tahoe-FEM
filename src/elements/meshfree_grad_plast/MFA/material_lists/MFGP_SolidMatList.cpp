/* $Id: MFGP_SolidMatList.cpp,v 1.1 2005-01-06 22:49:49 kyonten Exp $ */
#include "MFGP_SolidMatListT.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#include "DevelopmentElementsConfig.h"
#endif

#include "SolidMaterialT.h"
#include "MFGP_MaterialSupportT.h"

using namespace Tahoe;

/* constructors */
MFGP_SolidMatListT::MFGP_SolidMatListT(int length, const MFGP_MaterialSupportT& support):
	MaterialListT(length),
	fHasLocalizers(false),
	fHasThermal(false),
	fMFGPMaterialSupport(&support)
{
#ifdef __NO_RTTI__
	cout << "\n MFGP_SolidMatListT::MFGP_SolidMatListT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;
#endif
}

MFGP_SolidMatListT::MFGP_SolidMatListT(void):
	fHasLocalizers(false),
	fHasThermal(false),
	fMFGPMaterialSupport(NULL)
{

}

/* return true if the contains materials that generate heat */
bool MFGP_SolidMatListT::HasHeatSources(void) const
{
	/* check materials */
	bool has_heat = false;
	for (int i = 0; !has_heat && i < Length(); i++)
	{
		const ContinuumMaterialT* cont_mat = fArray[i];
		const SolidMaterialT* struct_mat = TB_DYNAMIC_CAST(const SolidMaterialT*, cont_mat);
		if (!struct_mat) throw ExceptionT::kGeneralFail;

		/* test */
		has_heat = struct_mat->HasIncrementalHeat();
	}

	return has_heat;
}
