/* $Id: eVerlet.h,v 1.4 2004-07-15 08:30:57 paklein Exp $ */
#ifndef _E_VERLET_H_
#define _E_VERLET_H_

/* base classes */
#include "Verlet.h"
#include "eIntegratorT.h"

namespace Tahoe {

/** Element controller for an explicit 4th order velocity verlet
  * algorithm */
class eVerlet: public virtual Verlet, public eIntegratorT
{
public:

	/** constructor */
	eVerlet(void);

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

#endif /* _E_VERLET_H_ */
