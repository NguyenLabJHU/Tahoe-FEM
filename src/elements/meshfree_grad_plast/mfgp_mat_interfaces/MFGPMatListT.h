/* $Id: MFGPMatListT.h,v 1.1 2005-04-26 22:27:01 kyonten Exp $ */
#ifndef _MFGP_MAT_LIST_T_H_
#define _MFGP_MAT_LIST_T_H_

/* base classes */
#include "pArrayT.h"
#include "ParameterInterfaceT.h"
#include "MFGPMaterialT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** base class for materials lists */
class MFGPMatListT: public pArrayT<MFGPMaterialT*>, public ParameterInterfaceT
{
public:

	/** constructor */
	MFGPMatListT(int length);
	MFGPMatListT(void);

	/* destructor */
	virtual ~MFGPMatListT(void) { };

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
inline bool MFGPMatListT::HasHistoryMaterials(void) const { return fHasHistory;  }

} // namespace Tahoe 
#endif /* _MFGP_MAT_LIST_T_H_ */
