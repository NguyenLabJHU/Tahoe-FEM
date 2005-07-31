/* $Id: IntegratorT.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/14/1996)                                          */
/* base class for controller roots (manages controller data)              */

#ifndef _CONTROLLER_T_H_
#define _CONTROLLER_T_H_

#include "Environment.h"

/* forward declarations */
class dArrayT;
class NodeManagerT;

class IntegratorT
{
public:

	/* constructor */
	IntegratorT(void);

	/* destructor */
	virtual ~IntegratorT(void);

	/* set the time step size */
	void SetTimeStep(double dt);

	/* take control of the time at which the external force vector
	 * is formed.  Default action to simply call the NodeManagerT's
	 * FormRHS function */
	virtual void FormNodalForce(NodeManagerT* nodeboss) const;

protected:

	/* called by SetParameters to compute the specific time stepping
	 * coefficients as needed */
	virtual void ComputeParameters(void) = 0;

protected:

	double	fdt; /* current time step value */		
};

#endif /* _CONTROLLER_T_H_ */
