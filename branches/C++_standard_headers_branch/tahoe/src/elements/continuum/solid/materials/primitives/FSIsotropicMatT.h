/* $Id: FSIsotropicMatT.h,v 1.2 2004-07-15 08:29:19 paklein Exp $ */
#ifndef _FS_ISOTROPIC_MAT_T_H_
#define _FS_ISOTROPIC_MAT_T_H_

/* base classes */
#include "FSSolidMatT.h"
#include "IsotropicT.h"

namespace Tahoe {

/** defines the interface for small strain isotropic materials */
class FSIsotropicMatT: public FSSolidMatT, public IsotropicT
{
public:

	/** constructor */
	FSIsotropicMatT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _FS_ISOTROPIC_MAT_T_H_ */
