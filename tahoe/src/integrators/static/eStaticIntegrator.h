/* $Id: eStaticIntegrator.h,v 1.5 2002-07-05 22:27:55 paklein Exp $ */
/* created: paklein (10/14/1996) */

#ifndef _E_STATIC_CONTROLLER_H_
#define _E_STATIC_CONTROLLER_H_

/* base classes */
#include "StaticT.h"
#include "eIntegratorT.h"

namespace Tahoe {

/** element component of the time integrator for quasi-static problems */ 
class eStaticIntegrator: public virtual StaticT, public eIntegratorT
{
public:

	/** constructor */
	eStaticIntegrator(void);

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

};

} // namespace Tahoe 
#endif /* _E_STATICCONTROLLER_H_ */
