/* $Id: eLinearHHTalpha.h,v 1.6 2004-07-15 08:30:28 paklein Exp $ */
/* created: paklein (10/14/1996) */
#ifndef _E_LINEARHHT_A_H_
#define _E_LINEARHHT_A_H_

/* base classes */
#include "HHTalpha.h"
#include "eIntegratorT.h"

namespace Tahoe {

/** element component of the HHT-\f$\alpha\f$ time integrator */
class eLinearHHTalpha: virtual public HHTalpha, public eIntegratorT
{
public:

	/** constructor */
	eLinearHHTalpha(double alpha);

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
	
	/** \name effective mass coefficients */
	/*@{*/
	double	fconstM;
	double	fconstC;
	double	fconstK;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _E_LINEARHHT_A_H_ */
