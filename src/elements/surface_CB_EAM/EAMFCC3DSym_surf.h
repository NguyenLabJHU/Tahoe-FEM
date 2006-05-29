/* $Id: EAMFCC3DSym_surf.h,v 1.2 2006-05-29 17:22:56 paklein Exp $ */
/* created: paklein (12/06/1996) */
#ifndef _EAMFCC3DSYM_SURF_H_
#define _EAMFCC3DSYM_SURF_H_

/* base class */
#include "EAMFCC3D_surf.h"

namespace Tahoe {

/* forward declarations */
class dMatrixT;

class EAMFCC3DSym_surf: public EAMFCC3D_surf
{
public:

	/** constructor */
	EAMFCC3DSym_surf(int nshells, int normal);

protected:

	/** initialize bond table values */
	virtual void LoadBondTable(void);
};

} /* namespace Tahoe */

#endif /* _EAMFCC3DSYM_SURF_H_ */
