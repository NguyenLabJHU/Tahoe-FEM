/* $Id: MFGP_SSSolidMatList3DT.h,v 1.4 2005-01-06 22:54:25 kyonten Exp $ */
/* created: paklein (02/14/1997) */
#ifndef _MFGP_SS_MATLIST_3D_T_H_
#define _MFGP_SS_MATLIST_3D_T_H_

/* base class */
#include "MFGP_SolidMatListT.h"
#include "SolidT.h"

//#include "SSMatSupportT.h"
#include "MFGP_MaterialSupportT.h"

namespace Tahoe {

/* forward declarations */
class MFGP_SSSolidMatT;

/** materials list for 3D structural analysis */
class MFGP_SSSolidMatList3DT: public MFGP_SolidMatListT, public SolidT
{
public:

	/** constructors */
	MFGP_SSSolidMatList3DT(int length, const MFGP_MaterialSupportT& support);
	MFGP_SSSolidMatList3DT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** construct the specified material or NULL if the request cannot be completed */
	MFGP_SSSolidMatT* NewSSSolidMat(const StringT& name) const;

private:

	/** support for finite strain materials */
	//const SSMatSupportT* fSSMatSupport;
	const MFGP_MaterialSupportT* fMFGPMaterialSupport;
	
};

} // namespace Tahoe 
#endif /* _MFGP_SS_MATLIST_3D_T_H_ */
