/* $Id: EAMFCC3DSym.h,v 1.3 2002-07-05 22:28:10 paklein Exp $ */
/* created: paklein (12/06/1996)                                          */
/* EAMFCC3DSym.h                                                          */

#ifndef _EAMFCC3DSYM_H_
#define _EAMFCC3DSYM_H_

/* base class */
#include "EAMFCC3D.h"

namespace Tahoe {

/* forward declarations */
class dMatrixT;

/* bond parameters */
const int kEAMFCC3DSymNumBonds = 27;

class EAMFCC3DSym: public EAMFCC3D
{
public:

	/* constructor */
	EAMFCC3DSym(ifstreamT& in, int EAMcode, int numspatialdim,
		int numbonds = kEAMFCC3DSymNumBonds);
	EAMFCC3DSym(ifstreamT& in, const dMatrixT& Q, int EAMcode, int numspatialdim,
		int numbonds = kEAMFCC3DSymNumBonds);
	    	
protected:

	/* initialize bond table values */
	virtual void LoadBondTable(void);
};

} // namespace Tahoe 
#endif /* _EAMFCC3DSYM_H_ */
