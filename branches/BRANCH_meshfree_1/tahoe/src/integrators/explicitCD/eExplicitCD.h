/* $Id: eExplicitCD.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (03/23/1997)                                          */
/* Element controller for an explicit 2nd order                           */
/* accurate, central difference time-stepping algorithm.                  */

#ifndef _E_EXP_CD_H_
#define _E_EXP_CD_H_

/* base classes */
#include "ControllerT.h"
#include "eControllerT.h"

class eExplicitCD: public virtual ControllerT, public eControllerT
{
public:

	/* constructor */
	eExplicitCD(void);

	/* time integration parameters */
	virtual StatDynFlagT StaticDynamic(void) const;
	virtual ImpExpFlagT ImplicitExplicit(void) const;

	/* return order time discretization */
	virtual int Order(void) const;

	/* returns 1 if the algorithm requires M, C, or K and sets const equal
	 * to the coefficient for the linear combination of components in the
	 * element effective mass matrix */
	virtual int FormM(double& constM) const;
	virtual int FormC(double& constC) const;
	virtual int FormK(double& constK) const;

	/* components of the internal force vector */
	virtual int FormMa(double& constMa) const;
	virtual int FormCv(double& constCv) const;
	virtual int FormKd(double& constKd) const;

protected:  	
	
	/* recalculate constants */
	virtual void eComputeParameters(void);

private:
	
	/* effective mass coefficients */
	double	fconstC; //assumes only mass proportionate Rayleigh damping
	
};

#endif /* _E_EXP_CD_H_ */
