/* $Id: MFGP_SSSolidMatList3DT.h,v 1.3 2004-10-30 00:50:12 raregue Exp $ */
/* created: paklein (02/14/1997) */
#ifndef _MFGP_SS_MATLIST_3D_T_H_
#define _MFGP_SS_MATLIST_3D_T_H_

/* base class */
#include "SolidMatListT.h"
#include "SolidT.h"

#include "SSMatSupportT.h"

namespace Tahoe {

/* forward declarations */
class SSSolidMatT;

/** materials list for 3D structural analysis */
class MFGP_SSSolidMatList3DT: public SolidMatListT, public SolidT
{
public:

	/** constructors */
	MFGP_SSSolidMatList3DT(int length, const SSMatSupportT& support);
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
	SSSolidMatT* NewSSSolidMat(const StringT& name) const;

private:

	/** support for finite strain materials */
	const SSMatSupportT* fSSMatSupport;
	
};

} // namespace Tahoe 
#endif /* _MFGP_SS_MATLIST_3D_T_H_ */
