/* $Id: SolidMatListT.cpp,v 1.14 2004-07-15 08:28:29 paklein Exp $ */
#include "SolidMatListT.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#include "DevelopmentElementsConfig.h"
#endif

#include "SolidMaterialT.h"
#include "SSMatSupportT.h"
#include "FSMatSupportT.h"
#ifdef GRAD_SMALL_STRAIN_DEV
#include "GradSSMatSupportT.h"
#endif

using namespace Tahoe;

/* constructors */
SolidMatListT::SolidMatListT(int length, const SolidMatSupportT& support):
	MaterialListT(length),
	fHasLocalizers(false),
	fHasThermal(false),
	fSolidMatSupport(&support)
{
#ifdef __NO_RTTI__
	cout << "\n SolidMatListT::SolidMatListT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;
#endif
}

SolidMatListT::SolidMatListT(void):
	fHasLocalizers(false),
	fHasThermal(false),
	fSolidMatSupport(NULL)
{

}

/* return true if the contains materials that generate heat */
bool SolidMatListT::HasHeatSources(void) const
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
