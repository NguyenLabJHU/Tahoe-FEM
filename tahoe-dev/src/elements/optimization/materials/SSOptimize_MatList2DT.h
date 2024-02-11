/* $Id: SSOptimize_MatList2DT.h,v 1.1 2009-04-23 03:03:51 thao Exp $ */
/* created: paklein (02/14/1997) */
#ifndef _SS_OPT_LIST_2D_T_H_
#define _SS_OPT_LIST_2D_T_H_

/* base classes */
#include "SSSolidMatList2DT.h"
namespace Tahoe {

/* forward declaration */
//class SSSolidMatT;
class SSOptimize_MatT;
class SSOptimize_MatSupportT;

/** small strain materials list for 2D structural analysis */
class SSOptimize_MatList2DT: public SSSolidMatList2DT
{
public:

	/** constructor */
	SSOptimize_MatList2DT(int length, const SSOptimize_MatSupportT& support);
	SSOptimize_MatList2DT(void);

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** construct the specified material or NULL if the request cannot be completed */
	SSOptimize_MatT* NewMat(const StringT& name) const;

private:

	/** support for small strain materials */
	const SSOptimize_MatSupportT* fSSOptimize_MatSupport;

};

} /* namespace Tahoe  */

#endif /* _MATLIST_2D_T_H_ */
