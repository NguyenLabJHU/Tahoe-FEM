/* $Id: SolidMatListT.h,v 1.2 2002-06-08 20:20:40 paklein Exp $ */

#ifndef _STRUCT_MAT_LIST_T_H_
#define _STRUCT_MAT_LIST_T_H_

/* base class */
#include "MaterialListT.h"

class SolidMatListT: public MaterialListT
{
public:

	/** constructor */
	SolidMatListT(int length);

	/** returns true if any of the materials in the list can undergo
	 * strain localization */
	bool HasLocalizingMaterials(void) const;

	/** returns true if any of the materials in the list is going to
	 * be subject to thermal loading */
	bool HasThermalStrains(void) const;

	/** return true if the contains materials that generate heat */
	bool HasHeatSources(void) const;

protected:

	/* flags for material properties */
	bool fHasLocalizers; /* materials that localize        */
	bool fHasThermal;    /* materials with thermal loading */
};

inline bool SolidMatListT::HasLocalizingMaterials(void) const { return fHasLocalizers; }
inline bool SolidMatListT::HasThermalStrains(void) const { return fHasThermal; }

#endif /* _STRUCT_MAT_LIST_T_H_ */