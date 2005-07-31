/* $Id: SSKStV.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#ifndef _SS_KSTV_H_
#define _SS_KSTV_H_

/* base classes */
#include "SSHookeanMatT.h"
#include "KStV.h"

class SSKStV: public SSHookeanMatT, public KStV
{
public:

	/* constructor */
	SSKStV(ifstreamT& in, const ElasticT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
};

#endif /* _SS_KSTV_H_ */
