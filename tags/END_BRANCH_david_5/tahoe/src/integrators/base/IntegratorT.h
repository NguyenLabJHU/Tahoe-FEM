/* $Id: IntegratorT.h,v 1.6 2004-07-15 08:30:33 paklein Exp $ */
/* created: paklein (10/14/1996) */
#ifndef _CONTROLLER_T_H_
#define _CONTROLLER_T_H_

#include "Environment.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;
class NodeManagerT;
class eIntegratorT;
class nIntegratorT;

/** abstract class for time integrators */
class IntegratorT
{
public:

	/** enum of integrator types */
	enum CodeT {
    kLinearStatic = 0,
          kStatic = 1,
       kTrapezoid = 2,
       kLinearHHT = 3,
    kNonlinearHHT = 4,
      kExplicitCD = 5,
          kVerlet = 6,
           kGear6 = 7
	};

	/** enum for implicit/explicit attribute */
	enum ImpExpFlagT {kImplicit, kExplicit};

	/** constructor */
	IntegratorT(void);

	/** destructor */
	virtual ~IntegratorT(void);

	/** factory method */
	static IntegratorT* New(int type, bool exception_on_fail);

	/** \name integrator parameters */
	/*@{*/
	/** return flag indicating whether integrator is implicit or explicit */
	virtual ImpExpFlagT ImplicitExplicit(void) const = 0;

	/** return order time discretization */
	virtual int Order(void) const = 0;

	/** return order field derivative which is treated as the primary 
	 * unknown value */
	virtual int OrderOfUnknown(void) const = 0;
	/*@}*/

	/** set the time step size */
	void SetTimeStep(double dt);

	/** take control of the time at which the external force vector
	 * is formed.  Default action to simply call the NodeManagerT's
	 * FormRHS function */
	virtual void FormNodalForce(NodeManagerT* nodeboss) const;

	/** \name casting up */
	/*@{*/
	virtual const eIntegratorT& eIntegrator(void) const = 0;
	virtual const nIntegratorT& nIntegrator(void) const = 0;
	/*@}*/

protected:

	/** called by SetParameters to compute the specific time stepping
	 * coefficients as needed */
	virtual void ComputeParameters(void) = 0;

protected:

	/** current time step value */
	double	fdt; 
};

} // namespace Tahoe 
#endif /* _CONTROLLER_T_H_ */
