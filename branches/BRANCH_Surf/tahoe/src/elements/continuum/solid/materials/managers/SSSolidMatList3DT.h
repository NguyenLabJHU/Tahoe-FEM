/* $Id: SSSolidMatList3DT.h,v 1.2 2004-07-15 08:28:28 paklein Exp $ */
/* created: paklein (02/14/1997) */
#ifndef _MATLIST_3D_T_H_
#define _MATLIST_3D_T_H_

/* base class */
#include "SolidMatListT.h"
#include "SolidT.h"

namespace Tahoe {

/* forward declarations */
class SSSolidMatT;

/** materials list for 3D structural analysis */
class SSSolidMatList3DT: public SolidMatListT, public SolidT
{
public:

	/** constructors */
	SSSolidMatList3DT(int length, const SSMatSupportT& support);
	SSSolidMatList3DT(void);

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
#endif /* _MATLIST_3D_T_H_ */
