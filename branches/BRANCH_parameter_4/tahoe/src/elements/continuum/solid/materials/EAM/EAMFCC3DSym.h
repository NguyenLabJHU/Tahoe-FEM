/* $Id: EAMFCC3DSym.h,v 1.3.64.1 2004-07-06 06:53:26 paklein Exp $ */
/* created: paklein (12/06/1996) */
#ifndef _EAMFCC3DSYM_H_
#define _EAMFCC3DSYM_H_

/* base class */
#include "EAMFCC3D.h"

namespace Tahoe {

/* forward declarations */
class dMatrixT;

class EAMFCC3DSym: public EAMFCC3D
{
public:

	/** constructor */
	EAMFCC3DSym(void);
//	EAMFCC3DSym(ifstreamT& in, int EAMcode, int nsd);

protected:

	/** initialize bond table values */
	virtual void LoadBondTable(void);
};

} /* namespace Tahoe */

#endif /* _EAMFCC3DSYM_H_ */
