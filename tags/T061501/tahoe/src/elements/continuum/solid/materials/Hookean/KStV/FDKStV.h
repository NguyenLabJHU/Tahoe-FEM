/* $Id: FDKStV.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#ifndef _FD_KSTV_H_
#define _FD_KSTV_H_

/* base classes */
#include "FDHookeanMatT.h"
#include "KStV.h"

class FDKStV: public FDHookeanMatT, public KStV
{
public:

	/* constructor */
	FDKStV(ifstreamT& in, const ElasticT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
};

#endif /* _FD_KSTV_H_ */
