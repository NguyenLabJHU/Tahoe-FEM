/* $Id: MaterialListT.h,v 1.9 2004-07-15 08:26:13 paklein Exp $ */
/* created: paklein (02/16/1997) */
#ifndef _MATERIAL_LIST_T_H_
#define _MATERIAL_LIST_T_H_

/* base classes */
#include "pArrayT.h"
#include "ParameterInterfaceT.h"
#include "ContinuumMaterialT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** base class for materials lists */
class MaterialListT: public pArrayT<ContinuumMaterialT*>, public ParameterInterfaceT
{
public:

	/** constructor */
	MaterialListT(int length);
	MaterialListT(void);

	/* destructor */
	virtual ~MaterialListT(void) { };

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

} // namespace Tahoe 
#endif /* _MATERIAL_LIST_T_H_ */
