/* $Id: QuadLog2D.h,v 1.4 2002-11-14 17:06:11 paklein Exp $ */
/* created: paklein (06/28/1997) */
#ifndef _QUAD_LOG_2D_
#define _QUAD_LOG_2D_

/* base classes */
#include "QuadLog3D.h"
#include "Material2DT.h"

namespace Tahoe {

/** (2D <-> 3D) translator for the QuadLog3D */
class QuadLog2D: public QuadLog3D, public Material2DT
{
public:

	/* constructor */
	QuadLog2D(ifstreamT& in, const FDMatSupportT& support);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

protected:

	/* return values */
	dSymMatrixT fStress2D;
	dMatrixT    fModulus2D;

	/* workspace */
	dSymMatrixT fb_2D;
};

} // namespace Tahoe 
#endif /* _QUAD_LOG_2D_ */
