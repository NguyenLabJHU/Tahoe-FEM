/* $Id: QuadLog2D.h,v 1.6 2004-07-15 08:27:35 paklein Exp $ */
/* created: paklein (06/28/1997) */
#ifndef _QUAD_LOG_2D_
#define _QUAD_LOG_2D_

/* base classes */
#include "QuadLog3D.h"

namespace Tahoe {

/** (2D <-> 3D) translator for the QuadLog3D */
class QuadLog2D: public QuadLog3D
{
public:

	/* constructor */
	QuadLog2D(void);

	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* return values */
	dSymMatrixT fStress2D;
	dMatrixT    fModulus2D;

	/* workspace */
	dSymMatrixT fb_2D;
};

} // namespace Tahoe 
#endif /* _QUAD_LOG_2D_ */
