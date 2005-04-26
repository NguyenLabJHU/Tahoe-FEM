/* $Id: MFGPSolidMatListT.h,v 1.1 2005-04-26 22:27:01 kyonten Exp $ */
#ifndef _MFGP_STRUCT_MAT_LIST_T_H_
#define _MFGP_STRUCT_MAT_LIST_T_H_

/* base class */
#include "MFGPMatListT.h"

namespace Tahoe {

/* forward declarations */
class MFGPMatSupportT;

/** list of materials for structural analysis */
class MFGPSolidMatListT: public MFGPMatListT
{
public:

	/** constructor */
	MFGPSolidMatListT(int length, const MFGPMatSupportT& support);
	MFGPSolidMatListT(void);

	/** returns true if any of the materials in the list can undergo
	 * strain localization */
	bool HasLocalizingMaterials(void) const;
	
	/** return true if the list contains plane stress models */
	virtual bool HasPlaneStress(void) const { return false; };

protected:

	/** \name flags for material properties */
	/*@{*/
	bool fHasLocalizers; /**< materials that localize */
	/*@}*/

	/** base class for structural material support */
	const MFGPMatSupportT* fMFGPMatSupport; 
};

inline bool MFGPSolidMatListT::HasLocalizingMaterials(void) const { return fHasLocalizers; }
} // namespace Tahoe 
#endif /* _MFGP_STRUCT_MAT_LIST_T_H_ */
