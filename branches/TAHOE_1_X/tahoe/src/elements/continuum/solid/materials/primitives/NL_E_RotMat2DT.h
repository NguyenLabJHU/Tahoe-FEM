/* $Id: NL_E_RotMat2DT.h,v 1.5 2003-01-29 07:35:08 paklein Exp $ */
/* created: paklein (06/13/1997) */
#ifndef _NL_E_ROTMAT_2D_T_H_
#define _NL_E_ROTMAT_2D_T_H_

/* base classes */
#include "NL_E_Mat2DT.h"
#include "Anisotropic2DT.h"

namespace Tahoe {

/** base class for materials with 2D nonlinear elastic behavior
 * with in-plane orientation with respect to global coordinate
 * axes, ie. the moduli, stress, and strain energy density functions
 * are formulated in the material's natural coordinates.
 * (See notes in NL_E_Mat2DT)
 */
class NL_E_RotMat2DT: public NL_E_Mat2DT, public Anisotropic2DT
{
public:

	/* constructor */
	NL_E_RotMat2DT(ifstreamT& in, const FSMatSupportT& support, ConstraintOptionT constraint);

	/* print parameters */
	virtual void Print(ostream& out) const;
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/* strain energy density */
	virtual double StrainEnergyDensity(void);
	
};

} // namespace Tahoe 
#endif /* _NL_E_ROTMAT_2D_T_H_ */
