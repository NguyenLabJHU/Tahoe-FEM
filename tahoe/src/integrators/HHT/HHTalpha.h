/* $Id: HHTalpha.h,v 1.5 2003-01-29 07:35:14 paklein Exp $ */
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
	 * \param in stream to extract parameters
	 * \param out stream to write messages
	 * \param auto2ndorder if true, sets the integrator parameters
	 *        to produce second order accuracy and does not read
	 *        any parameters from in. */
	HHTalpha(ifstreamT& in, ostream& out, bool auto2ndorder);

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
