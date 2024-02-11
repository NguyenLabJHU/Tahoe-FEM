/* $Id: SSOptimize_MatList3DT.h,v 1.1 2009-04-23 03:03:52 thao Exp $ */
/* created: paklein (02/14/1997) */
#ifndef _SSOptimize_MatList_3D_T_H_
#define _SSOptimize_MatList_3D_T_H_

/* base classes */
#include "SSSolidMatList3DT.h"

namespace Tahoe {

/* forward declaration */
//class SSSolidMatT;
class SSSolidMatT;
class SSOptimize_MatSupportT;

/** small strain materials list for 3D structural analysis */
class SSOptimize_MatList3DT: public SSSolidMatList3DT
{
public:

	/** constructor */
	SSOptimize_MatList3DT(int length,  SSMatSupportT& support);
	SSOptimize_MatList3DT(void);
	
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
	SSSolidMatT* NewMat(const StringT& name) const;

private:

	/** support for small strain materials */
	SSMatSupportT* fSupport;
};

} /* namespace Tahoe  */

#endif /* _MATLIST_3D_T_H_ */
