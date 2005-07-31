/* $Id: SWMaterial2D.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (08/25/1996)                                          */

#ifndef _SWMATERIAL2D_H_
#define _SWMATERIAL2D_H_

/* base classes */
#include "NL_E_RotMat2DT.h"
#include "SWDataT.h"

class SWMaterial2D: public NL_E_RotMat2DT, public SWDataT
{
public:

	/* constructor */
	SWMaterial2D(ifstreamT& in, const ElasticT& element);
	
	/* print parameters */
	virtual void Print(ostream& out) const;	
	
};

#endif /* _SWMATERIAL2D_H_ */
