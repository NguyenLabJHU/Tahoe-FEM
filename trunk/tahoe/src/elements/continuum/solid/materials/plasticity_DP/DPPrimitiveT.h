/* $Id: DPPrimitiveT.h,v 1.9 2004-07-15 08:28:48 paklein Exp $ */
/* created: myip (06/01/1999) */
#ifndef _DP_PRIMITIVET_H_
#define _DP_PRIMITIVET_H_

/* base class */
#include "ParameterInterfaceT.h"

namespace Tahoe {

/* forward declarations */
class dSymMatrixT;

/** base class for Drucker-Prager, nonassociative, small-strain,
 * pressure dependent plastic model with linear isotropic hardening
 */
class DPPrimitiveT: public ParameterInterfaceT
{
  public:

	/** constructor */
	DPPrimitiveT(void);

	/** destructor */
	virtual ~DPPrimitiveT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** \name accessors to parameters */
	/*@{*/
	double H_prime(void) const { return fH_prime; };
	/*@}*/

  protected:

  	/** returns the value of the yield function given the
  	 * Cauchy stress vector and state variables, where alpha is
  	 * the deviatoric stress-like internal state variable
  	 */
  	double YieldCondition(const dSymMatrixT& devstress, 
			      const double meanstress, double alpha) const;

  protected:

	/** \name parameters */
	/*@{*/
  	double falpha_bar; /**< cohesion-like strength parameter (falpha_bar >= 0.0) */
  	double ffriction;  /**< friction-like parameter (ffriction >= 0.0) */
  	double fdilation;  /**< dilation-like parameter (fdilation >= 0.0) */
  	double fH_prime;   /**< Deviatoric hardening parameter */
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _DP_PRIMITIVET_H_ */
