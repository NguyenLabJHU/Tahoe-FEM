/* $Id: eIntegratorT.h,v 1.2 2002-04-02 23:19:22 paklein Exp $ */
/* created: paklein (10/14/1996) */

#ifndef _E_CONTROLLERT_H_
#define _E_CONTROLLERT_H_

/* base class */
#include "IntegratorT.h"

/** Base class for a general (upto) second order element controller.  
 * This is the interface between the elements and the Controller 
 * class heirarchy. */
class eIntegratorT: virtual public IntegratorT
{
public:

	/** constructor */
	eIntegratorT(void);

	/** destructor */
	virtual ~eIntegratorT(void);

	/** \name elements of the effective mass matrix
	 * returns 1 if the algorithm requires M, C, or K and sets const equal
	 * to the coefficient for the linear combination of components in the
	 * element effective mass matrix */
	/*@{*/
	virtual int FormM(double& constM) const = 0;
	virtual int FormC(double& constC) const = 0;
	virtual int FormK(double& constK) const = 0;
	/*@}*/

	/** \name elements of the residual
	 * components of the internal force vector */
	/*@{*/
	virtual int FormMa(double& constMa) const = 0;
	virtual int FormCv(double& constCv) const = 0;
	virtual int FormKd(double& constKd) const = 0;
	/*@}*/

protected:
	
	/** recalculate element branch constants */
	virtual void eComputeParameters(void) = 0;
};

#endif /* _E_CONTROLLERT_H_ */
