/* $Id: SolidMatListT.cpp,v 1.4 2002-10-05 20:04:16 paklein Exp $ */
#include "SolidMatListT.h"
#include "SolidMaterialT.h"

using namespace Tahoe;

/* constructors */
SolidMatListT::SolidMatListT(int length):
	MaterialListT(length),
	fHasLocalizers(false),
	fHasThermal(false)
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
		const SolidMaterialT* struct_mat = 
			dynamic_cast<const SolidMaterialT*>(cont_mat);
		if (!struct_mat) throw eGeneralFail;

		/* test */
		has_heat = struct_mat->HasIncrementalHeat();
	}

	return has_heat;
}
