/* $Id: FSFiberMatT.h,v 1.1 2006-08-03 01:10:41 thao Exp $ */
/* created: paklein (06/09/1997) */
#ifndef _FD_FIB_MAT_T_H_
#define _FD_FIB_MAT_T_H_

/* base class */
#include "FSSolidMatT.h"

/* direct members */
#include "FSFiberMatSupportT.h"

namespace Tahoe {

/* forward declarations */
class UpLagFiberCompT;

/** base class for finite deformation fiber composite constitutive models. The interface *
 * provides access to the element-computed fiber orientation vectors in the global (lab) *
 * cartesian coordinates.                                                                */
class FSFiberMatT: public FSSolidMatT
{
public:

	/** constructor */
	FSFiberMatT(void);

	/** set the material support or pass NULL to clear */
	virtual void SetFSFiberMatSupport(const FSFiberMatSupportT* support);

	/** fiber materials support */
	const FSFiberMatSupportT& FiberMatSupportT(void) const;
	/*@}*/

protected:

	/** support for finite strain materials */
	const FSFiberMatSupportT* fFSFiberMatSupport;
};

/* fiber element materials support */
inline const FSFiberMatSupportT& FSFiberMatT::FiberMatSupportT(void) const
{ 
#if __option(extended_errorcheck)
	if (!fFSFiberMatSupport) 
		ExceptionT::GeneralFail("FSFiberMatT::FSFiberMatSupport", "pointer not set");
#endif

	return *fFSFiberMatSupport; 
}

} /* namespace Tahoe */

#endif /* _FD_STRUCT_MAT_T_H_ */
