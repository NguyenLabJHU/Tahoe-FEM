/* $Id: SolidMatListT.cpp,v 1.8 2003-07-29 21:17:49 rdorgan Exp $ */
#include "SolidMatListT.h"
#include "SolidMaterialT.h"
#include "SSMatSupportT.h"
#include "FSMatSupportT.h"
#include "GradSSMatSupportT.h"

using namespace Tahoe;

/* constructors */
SolidMatListT::SolidMatListT(int length, const SolidMatSupportT& support):
	MaterialListT(length),
	fHasLocalizers(false),
	fHasThermal(false),
	fSolidMatSupport(support),
	fSSMatSupport(NULL),
	fFSMatSupport(NULL),
	fGradSSMatSupport(NULL)
{
#ifdef __NO_RTTI__
	cout << "\n SolidMatListT::SolidMatListT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;

	/* cast and hope for the best */
	fSSMatSupport = (const SSMatSupportT*) &fSolidMatSupport;
	fFSMatSupport = (const FSMatSupportT*) &fSolidMatSupport;
	fGradSSMatSupport = (const GradSSMatSupportT*) &fSolidMatSupport;
#else

	/* cast to small strain support */
	fSSMatSupport = dynamic_cast<const SSMatSupportT*>(&fSolidMatSupport);

	/* cast to finite strain support */
	fFSMatSupport = dynamic_cast<const FSMatSupportT*>(&fSolidMatSupport);
	
	/* cast to gradient enhanced small strain support */
	fGradSSMatSupport = dynamic_cast<const GradSSMatSupportT*>(&fSolidMatSupport);

	/* must have at least one */
	if (!fSSMatSupport && !fFSMatSupport && !fGradSSMatSupport)
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
