/* $Id: SolidMatListT.h,v 1.4 2002-10-05 20:04:16 paklein Exp $ */

#ifndef _STRUCT_MAT_LIST_T_H_
#define _STRUCT_MAT_LIST_T_H_

/* base class */
#include "MaterialListT.h"


namespace Tahoe {

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
	
	/** return true if the list contains plane stress models */
	virtual bool HasPlaneStress(void) const { return false; };

protected:

	/* flags for material properties */
	bool fHasLocalizers; /* materials that localize        */
	bool fHasThermal;    /* materials with thermal loading */
};

inline bool SolidMatListT::HasLocalizingMaterials(void) const { return fHasLocalizers; }
inline bool SolidMatListT::HasThermalStrains(void) const { return fHasThermal; }

} // namespace Tahoe 
#endif /* _STRUCT_MAT_LIST_T_H_ */
