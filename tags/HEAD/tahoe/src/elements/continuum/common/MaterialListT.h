/* $Id: MaterialListT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (02/16/1997)                                          */

#ifndef _MATERIAL_LIST_T_H_
#define _MATERIAL_LIST_T_H_

/* base class */
#include "pArrayT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ContinuumMaterialT;
class ifstreamT;

class MaterialListT: public pArrayT<ContinuumMaterialT*>
{
public:

	/* constructors */
	MaterialListT(int length);

	/* read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in) = 0;

	/* write material data to the output stream */
	void WriteMaterialData(ostream& out) const;

	/* apply pre-conditions at the current time step */
	void InitStep(void);

	/* returns true if any of the materials in the list makes
	 * use of history variables. If it does, Update/Reset
	 * of these variables needs to be taken care of */
	bool HasHistoryMaterials(void) const;

	/* returns true if any of the materials in the list can undergo
	 * strain localization */
	bool HasLocalizingMaterials(void) const;

	/* returns true if any of the materials in the list is going to
	 * be subject to thermal loading */
	bool HasThermalStrains(void) const;

protected:

	/* flags for material properties */
	bool fHasHistory;    /* internal history variables     */
	bool fHasLocalizers; /* materials that localize        */
	bool fHasThermal;    /* materials with thermal loading */
	
};

/* inlines */
inline bool MaterialListT::HasHistoryMaterials(void) const { return fHasHistory;  }
inline bool MaterialListT::HasLocalizingMaterials(void) const { return fHasLocalizers; }
inline bool MaterialListT::HasThermalStrains(void) const { return fHasThermal; }

#endif /* _MATERIAL_LIST_T_H_ */
