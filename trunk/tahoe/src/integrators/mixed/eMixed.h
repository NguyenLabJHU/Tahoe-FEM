/* $Header: /home/cvs/t/ta/tahoe/tahoe/src/integrators/mixed/eMixed.h,v 1.2 2006-08-18 01:15:50 a-kopacz Exp $ */
/* created: a-kopacz (08/08/2006) */

#ifndef _E_MIXED_H_
#define _E_MIXED_H_

/* base classes */
#include "Mixed.h"
#include "eIntegratorT.h"


namespace Tahoe {

class eMixed: public virtual Mixed, public eIntegratorT
{
public:

	/** constructor */
	eMixed(void);

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
	double	fconstC;
	double	fconstK;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _E_MIXED_H_ */
