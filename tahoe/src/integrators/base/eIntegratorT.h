/* $Id: eIntegratorT.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/14/1996)                                          */
/* Base class for a general (upto) second order element controller.  This */
/* is the interface between the elements and the Controller class heirarchy. */

#ifndef _E_CONTROLLERT_H_
#define _E_CONTROLLERT_H_

class eIntegratorT
{
public:

	/* time integration flags */
	enum StatDynFlagT {kStatic, kDynamic};
	enum ImpExpFlagT  {kImplicit, kExplicit};

	/* constructor */
	eIntegratorT(void);

	/* destructor */
	virtual ~eIntegratorT(void);

	/* time integration parameters */
	virtual StatDynFlagT StaticDynamic(void) const = 0;
	virtual ImpExpFlagT  ImplicitExplicit(void) const = 0;
	
	/* return order time discretization */
	virtual int Order(void) const = 0;

	/* returns 1 if the algorithm requires M, C, or K and sets const equal
	 * to the coefficient for the linear combination of components in the
	 * element effective mass matrix.  Returning 0 implies that the constX
	 * are exactly zero */
	virtual int FormM(double& constM) const = 0;
	virtual int FormC(double& constC) const = 0;
	virtual int FormK(double& constK) const = 0;

	/* components of the internal force vector.  Returning 0 implies that
	 * the constXx are exactly zero */
	virtual int FormMa(double& constMa) const = 0;
	virtual int FormCv(double& constCv) const = 0;
	virtual int FormKd(double& constKd) const = 0;

protected:
	
	/* recalculate element branch constants */
	virtual void eComputeParameters(void) = 0;
};

#endif /* _E_CONTROLLERT_H_ */
