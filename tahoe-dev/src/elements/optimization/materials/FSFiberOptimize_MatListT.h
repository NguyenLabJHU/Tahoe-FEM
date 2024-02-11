/* $Id: FSFiberOptimize_MatListT.h,v 1.1 2009-04-23 03:03:50 thao Exp $ */
/* created: paklein (02/14/1997) */
#ifndef _FSFiberOptimize_MatListT_
#define _FSFiberOptimize_MatListT_

/* base classes */
#include "FSFiberMatListT.h"

namespace Tahoe {

/* forward declaration */
//class SSSolidMatT;
class FSFiberMatT;
class FSFiberMatSupportT;

/** small strain materials list for 3D structural analysis */
class FSFiberOptimize_MatListT: public FSFiberMatListT
{
public:

	/** constructor */
	FSFiberOptimize_MatListT(int length,  FSFiberMatSupportT& support);
	FSFiberOptimize_MatListT(void);
	
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
	FSFiberMatT* NewFSFiberMat(const StringT& name) const;

private:

	/** support for small strain materials */
	FSFiberMatSupportT* fSupport;
};

} /* namespace Tahoe  */

#endif /* _MATLIST_3D_T_H_ */
