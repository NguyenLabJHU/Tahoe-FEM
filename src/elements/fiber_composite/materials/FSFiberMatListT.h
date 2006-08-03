/* $Id: FSFiberMatListT.h,v 1.1 2006-08-03 01:10:41 thao Exp $ */
#ifndef _FIBER_MATLIST_3D_T_H_
#define _FIBER_MATLIST_3D_T_H_

/* base class */
#include "FSSolidMatList3DT.h"

namespace Tahoe {

/* forward declarations */
class FSFiberMatT;
class FSFiberMatSupportT;

/** materials list for 3D structural analysis */
class FSFiberMatListT: public FSSolidMatList3DT
{
public:

	/** constructors */
	FSFiberMatListT(int length, const FSFiberMatSupportT& support);
	FSFiberMatListT(void);

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
	FSFiberMatT* NewFSFiberMat(const StringT& name) const;

private:

	/** support for finite strain fiber composite materials */
	const FSFiberMatSupportT* fFSFiberMatSupport;
};

} /* namespace Tahoe */

#endif /* _FIBER_MATLIST_3D_T_H_ */
