/* $Id: SolidMatListT.cpp,v 1.6 2002-11-14 17:06:22 paklein Exp $ */
#include "SolidMatListT.h"
#include "SolidMaterialT.h"
#include "SSMatSupportT.h"
#include "FDMatSupportT.h"

using namespace Tahoe;

/* constructors */
SolidMatListT::SolidMatListT(int length, const SolidMatSupportT& support):
	MaterialListT(length),
	fHasLocalizers(false),
	fHasThermal(false),
	fSolidMatSupport(support),
	fSSMatSupport(NULL),
	fFDMatSupport(NULL)
{
#ifdef __NO_RTTI__
	cout << "\n SolidMatListT::SolidMatListT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;

	/* cast and hope for the best */
	fSSMatSupport = (const SSMatSupportT*) &fSolidMatSupport;
	fFDMatSupport = (const FDMatSupportT*) &fSolidMatSupport;
#else

	/* cast to small strain support */
	fSSMatSupport = dynamic_cast<const SSMatSupportT*>(&fSolidMatSupport);

	/* cast to finite strain support */
	fFDMatSupport = dynamic_cast<const FDMatSupportT*>(&fSolidMatSupport);
	
	/* must have at least one */
	if (!fSSMatSupport && !fFDMatSupport)
	{
		cout << "\n SolidMatListT::SolidMatListT: could not cast element group to\n" 
		     <<   "     neither SSMatSupportT nor SSMatSupportT" << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif
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
		if (!struct_mat) throw ExceptionT::kGeneralFail;

		/* test */
		has_heat = struct_mat->HasIncrementalHeat();
	}

	return has_heat;
}
