/* $Id: SolidMatListT.h,v 1.4.6.1 2002-11-13 08:44:20 paklein Exp $ */
#ifndef _STRUCT_MAT_LIST_T_H_
#define _STRUCT_MAT_LIST_T_H_

/* base class */
#include "MaterialListT.h"

namespace Tahoe {

/* forward declarations */
class SolidMatSupportT;
class SSMatSupportT;
class FDMatSupportT;

/** list of materials for structural analysis */
class SolidMatListT: public MaterialListT
{
public:

	/** constructor */
	SolidMatListT(int length, const SolidMatSupportT& support);

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

	/** \name material support classes 
	 * These are dynamically cast from the MaterialSupportT passed
	 * into the constructor */
	/*@{*/
	/** base class for structural material support */
	const SolidMatSupportT& fSolidMatSupport; 
	
	/** support for small strain materials */
	const SSMatSupportT* fSSMatSupport;

	/** support for finite strain materials */
	const FDMatSupportT* fFDMatSupport;
	/*@}*/
};

inline bool SolidMatListT::HasLocalizingMaterials(void) const { return fHasLocalizers; }
inline bool SolidMatListT::HasThermalStrains(void) const { return fHasThermal; }

} // namespace Tahoe 
#endif /* _STRUCT_MAT_LIST_T_H_ */
