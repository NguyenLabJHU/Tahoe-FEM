/* $Id: FDCubic2DT.h,v 1.2.6.1 2002-06-27 18:03:11 cjkimme Exp $ */
/* created: paklein (06/11/1997)                                          */

#ifndef _FD_CUBIC_2D_T_H_
#define _FD_CUBIC_2D_T_H_

/* base classes */
#include "FDCubicT.h"
#include "Anisotropic2DT.h"
#include "Material2DT.h"


namespace Tahoe {

class FDCubic2DT: public FDCubicT, public Anisotropic2DT, public Material2DT
{
public:

	/* constructor */
	FDCubic2DT(ifstreamT& in, const FiniteStrainT& element);

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
#endif /* _FD_CUBIC_2D_T_H_ */
