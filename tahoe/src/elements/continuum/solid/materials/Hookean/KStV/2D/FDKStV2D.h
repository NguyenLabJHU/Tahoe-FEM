/* $Id: FDKStV2D.h,v 1.5.30.1 2004-03-02 17:46:15 paklein Exp $ */
/* created: paklein (06/10/97) */
#ifndef _FD_KSTV_2D_H_
#define _FD_KSTV_2D_H_

/* base classes */
#include "FDKStV.h"
#include "IsotropicT.h"

namespace Tahoe {

class FDKStV2D: public FDKStV
{
public:

	/* constructor */
	FDKStV2D(ifstreamT& in, const FSMatSupportT& support);

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

private:

	/* set inverse of thermal transformation - return true if active */
	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);  			
};

} // namespace Tahoe 
#endif /* _FD_KSTV_2D_H_ */
