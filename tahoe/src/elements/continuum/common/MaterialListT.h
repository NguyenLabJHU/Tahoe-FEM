/* $Id: MaterialListT.h,v 1.2.2.1 2002-05-17 01:23:27 paklein Exp $ */
/* created: paklein (02/16/1997) */

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

	/** constructor */
	MaterialListT(int length);

	/** read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in) = 0;

	/** write material data to the output stream */
	void WriteMaterialData(ostream& out) const;

	/** apply pre-conditions at the current time step */
	void InitStep(void);

	/** finalize the current time step */
	void CloseStep(void);
	
	/** returns true if any of the materials in the list makes
	 * use of history variables. If it does, Update/Reset
	 * of these variables needs to be taken care of */
	bool HasHistoryMaterials(void) const;

protected:

	/** true if list contains materials with history variables */
	bool fHasHistory;    

};

/* inlines */
inline bool MaterialListT::HasHistoryMaterials(void) const { return fHasHistory;  }

#endif /* _MATERIAL_LIST_T_H_ */
