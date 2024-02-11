/* $Id: EAMFCC3DSym_edge.h,v 1.2 2009-06-04 16:25:46 hspark Exp $ */
/* created: hspark (6/2/2009) */
#ifndef _EAMFCC3DSYM_EDGE_H_
#define _EAMFCC3DSYM_EDGE_H_

/* base class */
#include "EAMFCC3D_edge.h"

namespace Tahoe {

/* forward declarations */
class dMatrixT;

class EAMFCC3DSym_edge: public EAMFCC3D_edge
{
public:

	/** constructor */
	EAMFCC3DSym_edge(int nshells, int normal);

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

#endif /* _EAMFCC3DSYM_EDGE_H_ */
