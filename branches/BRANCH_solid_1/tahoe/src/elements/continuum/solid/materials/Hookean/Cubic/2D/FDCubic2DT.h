/* $Id: FDCubic2DT.h,v 1.1.1.1.2.1 2001-06-06 16:22:03 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#ifndef _FD_CUBIC_2D_T_H_
#define _FD_CUBIC_2D_T_H_

/* base classes */
#include "FDCubicT.h"
#include "Material2DT.h"

class FDCubic2DT: public FDCubicT, public Material2DT
{
public:

	/* constructor */
	FDCubic2DT(ifstreamT& in, const ElasticT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);
	
private:

	/* set inverse of thermal transformation - return true if active */
	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);  			
};

#endif /* _FD_CUBIC_2D_T_H_ */
