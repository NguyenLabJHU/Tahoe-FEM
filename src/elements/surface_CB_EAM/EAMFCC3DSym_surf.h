/* $Id: EAMFCC3DSym_surf.h,v 1.1 2006-05-21 15:55:19 hspark Exp $ */
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
	EAMFCC3DSym_surf(void);

protected:

	/** initialize bond table values */
	virtual void LoadBondTable(void);
};

} /* namespace Tahoe */

#endif /* _EAMFCC3DSYM_SURF_H_ */
