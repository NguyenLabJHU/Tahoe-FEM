/* $Id: EAMFCC3DSym_110surf.h,v 1.1 2007-06-12 01:50:02 hspark Exp $ */
/* created: paklein (12/06/1996) */
#ifndef _EAMFCC3DSYM_110SURF_H_
#define _EAMFCC3DSYM_110SURF_H_

/* base class */
/* This class defines {110} surfaces for <100> oriented nanostructures
#include "EAMFCC3D_surf.h"

namespace Tahoe {

/* forward declarations */
class dMatrixT;

class EAMFCC3DSym_110surf: public EAMFCC3D_surf
{
public:

	/** constructor */
	EAMFCC3DSym_110surf(int nshells, int normal);

protected:

	/** initialize bond table values */
	virtual void LoadBondTable(void);
	
private:

	/** normal code to do bond table rotation */
	int fNormalCode;
	
	/** Return rotation matrix for bond table */
	dMatrixT RotationMatrixA(const double angle);
	
	/** Reeturn other rotation matrix for bond table */
	dMatrixT RotationMatrixB(const double angle);
	
};

} /* namespace Tahoe */

#endif /* _EAMFCC3DSYM_110SURF_H_ */
