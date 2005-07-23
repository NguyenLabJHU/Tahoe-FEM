/* $Id: eExplicitCD.h,v 1.5 2003-01-29 07:35:15 paklein Exp $ */
/* created: paklein (03/23/1997) */

#ifndef _E_EXP_CD_H_
#define _E_EXP_CD_H_

/* base classes */
#include "ExplicitCD.h"
#include "eIntegratorT.h"

namespace Tahoe {

/** Element controller for an explicit 2nd order
 * accurate, central difference time-stepping algorithm */
class eExplicitCD: public virtual ExplicitCD, public eIntegratorT
{
public:

	/** constructor */
	eExplicitCD(void);

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
#endif /* _E_EXP_CD_H_ */
