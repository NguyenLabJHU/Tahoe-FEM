/* $Id: SolidMatList2DT.h,v 1.12 2003-12-02 17:12:22 paklein Exp $ */
/* created: paklein (02/14/1997) */
#ifndef _MATLIST_2D_T_H_
#define _MATLIST_2D_T_H_

/* base classes */
#include "SolidMatListT.h"
#include "SolidT.h"

namespace Tahoe {

/** materials list for 2D structural analysis */
class SolidMatList2DT: public SolidMatListT, public SolidT
{
public:

	/** constructor */
	SolidMatList2DT(int length, const SolidMatSupportT& support);
	SolidMatList2DT(void);

	/** read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

	/** return true if the list contains plane stress models */
	virtual bool HasPlaneStress(void) const;

private:
	
	/** \name errror messages */
	/*@{*/
	void Error_no_small_strain(ostream& out, int matcode) const;
	void Error_no_finite_strain(ostream& out, int matcode) const;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _MATLIST_2D_T_H_ */
