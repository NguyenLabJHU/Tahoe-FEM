/* $Id: SolidMatListT.cpp,v 1.13 2003-12-28 08:23:33 paklein Exp $ */
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
	fSolidMatSupport(&support),
	fSSMatSupport(NULL),
	fFSMatSupport(NULL),
	fGradSSMatSupport(NULL)
{
#ifdef __NO_RTTI__
	cout << "\n SolidMatListT::SolidMatListT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;
#endif

	/* cast to small strain support */
	fSSMatSupport = TB_DYNAMIC_CAST(const SSMatSupportT*, fSolidMatSupport);

	/* cast to finite strain support */
	fFSMatSupport = TB_DYNAMIC_CAST(const FSMatSupportT*, fSolidMatSupport);
#ifdef GRAD_SMALL_STRAIN_DEV
	/* cast to gradient enhanced small strain support */
	fGradSSMatSupport = TB_DYNAMIC_CAST(const GradSSMatSupportT*, fSolidMatSupport);
#endif

	/* must have at least one */
	if (!fSSMatSupport && !fFSMatSupport && !fGradSSMatSupport)
	{
		cout << "\n SolidMatListT::SolidMatListT: could not cast element group to\n" 
		     <<   "     neither SSMatSupportT nor SSMatSupportT" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

SolidMatListT::SolidMatListT(void):
	fHasLocalizers(false),
	fHasThermal(false),
	fSolidMatSupport(NULL),
	fSSMatSupport(NULL),
	fFSMatSupport(NULL),
	fGradSSMatSupport(NULL)
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
