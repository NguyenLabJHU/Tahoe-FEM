/* $Id: SWMaterial2D.h,v 1.5 2003-01-29 07:34:37 paklein Exp $ */
/* created: paklein (08/25/1996) */
#ifndef _SWMATERIAL2D_H_
#define _SWMATERIAL2D_H_

/* base classes */
#include "NL_E_RotMat2DT.h"
#include "SWDataT.h"

namespace Tahoe {

class SWMaterial2D: public NL_E_RotMat2DT, public SWDataT
{
public:

	/* constructor */
	SWMaterial2D(ifstreamT& in, const FSMatSupportT& support);
	
	/* print parameters */
	virtual void Print(ostream& out) const;	
	
};

} // namespace Tahoe 
#endif /* _SWMATERIAL2D_H_ */
