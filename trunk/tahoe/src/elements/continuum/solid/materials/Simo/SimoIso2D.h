/* $Id: SimoIso2D.h,v 1.9 2004-09-10 22:39:12 paklein Exp $ */
/* created: paklein (03/04/1997) */
#ifndef _SIMO_ISO_2D_H_
#define _SIMO_ISO_2D_H_

/* base classes */
#include "SimoIso3D.h"

namespace Tahoe {

/** (2D <-> 3D) translator for the SimoIso3D */
class SimoIso2D: public SimoIso3D
{
public:

	/** constructor */
	SimoIso2D(void);

	/** initialize step. Verify that the thermal dilatation deformation
	 * gradient is equibiaxial. This state is assumed when computing
	 * the state of plain strain deformation. */
	virtual void InitStep(void);

	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** strain energy density */
	virtual double StrainEnergyDensity(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** compute 3D stretch tensor \b b from the 2D deformation state. 
	 * \todo Make this a FSSolidMatT function? */
	void Compute_b_3D(dSymMatrixT& b_3D);

protected:

	/* return values */
	dSymMatrixT fStress2D;  /**< return value for the Cauchy stress */
	dMatrixT    fModulus2D; /**< return value for the spatial tangent modulus */
	
	/** workspace */
	dSymMatrixT fb_2D;		 	 	
};

} // namespace Tahoe 
#endif /* _SIMO_ISO_2D_H_ */
