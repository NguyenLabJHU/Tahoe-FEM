/* $Id: FDKStV2D.h,v 1.1.1.1.2.1 2001-06-06 16:20:45 paklein Exp $ */
/* created: paklein (06/10/97)                                            */

#ifndef _FD_KSTV_2D_H_
#define _FD_KSTV_2D_H_

/* base classes */
#include "FDKStV.h"
#include "IsotropicT.h"

class FDKStV2D: public FDKStV, public Material2DT
{
public:

	/* constructor */
	FDKStV2D(ifstreamT& in, const ElasticT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

private:

	/* set inverse of thermal transformation - return true if active */
	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);  			
};

#endif /* _FD_KSTV_2D_H_ */
