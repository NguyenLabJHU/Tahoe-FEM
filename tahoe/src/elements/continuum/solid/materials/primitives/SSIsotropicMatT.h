/* $Id: SSIsotropicMatT.h,v 1.1.6.1 2004-07-06 06:54:06 paklein Exp $ */
#ifndef _SS_ISOTROPIC_MAT_T_H_
#define _SS_ISOTROPIC_MAT_T_H_

/* base classes */
#include "SSSolidMatT.h"
#include "IsotropicT.h"

namespace Tahoe {

/** defines the interface for small strain isotropic materials */
class SSIsotropicMatT: public SSSolidMatT, public IsotropicT
{
public:

	/** constructor */
	SSIsotropicMatT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _SS_ISOTROPIC_MAT_T_H_ */
