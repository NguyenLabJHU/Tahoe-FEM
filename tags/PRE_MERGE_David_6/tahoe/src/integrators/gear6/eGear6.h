/* $Id: eGear6.h,v 1.4 2004-07-15 08:30:43 paklein Exp $ */
#ifndef _E_GEAR_06_H_
#define _E_GEAR_06_H_

/* base classes */
#include "Gear6.h"
#include "eIntegratorT.h"

namespace Tahoe {

/** Element controller for an explicit 4th order velocity verlet
  * algorithm */
class eGear6: public virtual Gear6, public eIntegratorT
{
public:

	/** constructor */
	eGear6(void);

	/** \name elements of the effective mass matrix
	 * returns 1 if the algorithm requires M, C, or K and sets const equal
	 * to the coefficient for the linear combination of components in the
	 * element effective mass matrix */
	/*@{*/
	virtual int FormM(double& constM) const;
	virtual int FormC(double& constC) const;
	virtual int FormK(double& constK) const;
	/*@}*/

	/** \name elements of the residual
	 * components of the internal force vector */
	/*@{*/
	virtual int FormMa(double& constMa) const;
	virtual int FormCv(double& constCv) const;
	virtual int FormKd(double& constKd) const;
	/*@}*/

protected:  	
	
	/** recalculate constants */
	virtual void eComputeParameters(void);

private:
	
	/** effective mass coefficients */
	double fconstC;
	
};

} // namespace Tahoe

#endif /* _E_Gear_06_H_ */
