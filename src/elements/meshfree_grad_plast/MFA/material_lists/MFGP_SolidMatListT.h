/* $Id: MFGP_SolidMatListT.h,v 1.1 2005-01-06 22:49:49 kyonten Exp $ */
#ifndef _MFGP_STRUCT_MAT_LIST_T_H_
#define _MFGP_STRUCT_MAT_LIST_T_H_

/* base class */
#include "MaterialListT.h"

namespace Tahoe {

/* forward declarations */
class MFGP_MaterialSupportT;

/** list of materials for structural analysis */
class MFGP_SolidMatListT: public MaterialListT
{
public:

	/** constructor */
	MFGP_SolidMatListT(int length, const MFGP_MaterialSupportT& support);
	MFGP_SolidMatListT(void);

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

	/** \name flags for material properties */
	/*@{*/
	bool fHasLocalizers; /**< materials that localize */
	bool fHasThermal;    /**< materials with thermal loading */
	/*@}*/

	/** base class for structural material support */
	const MFGP_MaterialSupportT* fMFGPMaterialSupport; 
};

inline bool MFGP_SolidMatListT::HasLocalizingMaterials(void) const { return fHasLocalizers; }
inline bool MFGP_SolidMatListT::HasThermalStrains(void) const { return fHasThermal; }

} // namespace Tahoe 
#endif /* _MFGP_STRUCT_MAT_LIST_T_H_ */
