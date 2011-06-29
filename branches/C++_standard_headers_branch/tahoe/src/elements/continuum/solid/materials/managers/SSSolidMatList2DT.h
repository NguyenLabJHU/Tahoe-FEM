/* $Id: SSSolidMatList2DT.h,v 1.2 2004-07-15 08:28:28 paklein Exp $ */
/* created: paklein (02/14/1997) */
#ifndef _SS_MATLIST_2D_T_H_
#define _SS_MATLIST_2D_T_H_

/* base classes */
#include "SolidMatListT.h"
#include "SolidT.h"

namespace Tahoe {

/* forward declaration */
class SSSolidMatT;

/** small strain materials list for 2D structural analysis */
class SSSolidMatList2DT: public SolidMatListT, public SolidT
{
public:

	/** constructor */
	SSSolidMatList2DT(int length, const SSMatSupportT& support);
	SSSolidMatList2DT(void);

	/** return true if the list contains plane stress models */
	virtual bool HasPlaneStress(void) const;

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

	/** support for small strain materials */
	const SSMatSupportT* fSSMatSupport;

	/** support for gradient enhanced small strain materials */
	const GradSSMatSupportT* fGradSSMatSupport;
};

} /* namespace Tahoe  */

#endif /* _MATLIST_2D_T_H_ */
