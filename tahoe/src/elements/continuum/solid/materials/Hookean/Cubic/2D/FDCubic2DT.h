/* $Id: FDCubic2DT.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#ifndef _FD_CUBIC_2D_T_H_
#define _FD_CUBIC_2D_T_H_

/* base classes */
#include "FDHookeanMatT.h"
#include "Cubic2DT.h"

class FDCubic2DT: public FDHookeanMatT, public Cubic2DT
{
public:

	/* constructor */
	FDCubic2DT(ifstreamT& in, const ElasticT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
private:

	/* set inverse of thermal transformation - return true if active */
	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);  			
};

#endif /* _FD_CUBIC_2D_T_H_ */
