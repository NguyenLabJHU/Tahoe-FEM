/* $Id: FDKStV2D.h,v 1.3 2002-07-02 19:55:46 cjkimme Exp $ */
/* created: paklein (06/10/97)                                            */

#ifndef _FD_KSTV_2D_H_
#define _FD_KSTV_2D_H_

/* base classes */
#include "FDKStV.h"
#include "IsotropicT.h"


namespace Tahoe {

class FDKStV2D: public FDKStV, public Material2DT
{
public:

	/* constructor */
	FDKStV2D(ifstreamT& in, const FiniteStrainT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

private:

	/* set inverse of thermal transformation - return true if active */
	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);  			
};

} // namespace Tahoe 
#endif /* _FD_KSTV_2D_H_ */
