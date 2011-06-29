/* $Id: HHTalpha.h,v 1.6 2004-07-15 08:30:27 paklein Exp $ */
/* created: paklein (10/14/1996) */

#ifndef _HHT_ALPHA_H_
#define _HHT_ALPHA_H_

#include "Environment.h"

/* base class */
#include "IntegratorT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dArrayT;

/** HHT-\f$\alpha\f$ time integrator */
class HHTalpha: virtual public IntegratorT
{
public:

	/** constructor 
	 * \param alpha damping parameters 
	 */
	HHTalpha(double alpha);

	/** \name integrator parameters */
	/*@{*/
	/** return flag indicating whether integrator is implicit or explicit */
	virtual ImpExpFlagT ImplicitExplicit(void) const { return kImplicit; };

	/** return order time discretization */
	virtual int Order(void) const { return 2; };

	/** return order field derivative which is treated as the primary 
	 * unknown value */
	virtual int OrderOfUnknown(void) const { return 2; };
	/*@}*/

protected:

	/** set time integration to single parameter 2nd order */
	void Set2ndOrder(double alpha);
	
protected:

	/** autoset parameters */
	bool fAuto2ndOrder;

	/** \name time integration parameters */
	/*@{*/
	double fgamma;
	double fbeta;
	double falpha;
	/*@}*/		
};

} // namespace Tahoe 
#endif /* _HHT_ALPHA_H_ */
