/* $Id: MFGP_SSIsotropicMatT.h,v 1.1 2005-01-06 22:51:27 kyonten Exp $ */
#ifndef _MFGP_SS_ISOTROPIC_MAT_T_H_
#define _MFGP_SS_ISOTROPIC_MAT_T_H_

/* base classes */
#include "MFGP_SSSolidMatT.h"
#include "IsotropicT.h"

namespace Tahoe {

/** defines the interface for small strain isotropic materials */
class MFGP_SSIsotropicMatT: public MFGP_SSSolidMatT, public IsotropicT
{
public:

	/** constructor */
	MFGP_SSIsotropicMatT(void);

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

#endif /* _MFGP_SS_ISOTROPIC_MAT_T_H_ */
