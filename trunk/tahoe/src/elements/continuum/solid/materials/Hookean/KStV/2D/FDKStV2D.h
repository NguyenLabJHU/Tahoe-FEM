/* $Id: FDKStV2D.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/10/97)                                            */

#ifndef _FD_KSTV_2D_H_
#define _FD_KSTV_2D_H_

/* base classes */
#include "FDHookeanMatT.h"
#include "KStV2D.h"

class FDKStV2D: public FDHookeanMatT, public KStV2D
{
public:

	/* constructor */
	FDKStV2D(ifstreamT& in, const ElasticT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

private:

	/* set inverse of thermal transformation - return true if active */
	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);  			
};

#endif /* _FD_KSTV_2D_H_ */
