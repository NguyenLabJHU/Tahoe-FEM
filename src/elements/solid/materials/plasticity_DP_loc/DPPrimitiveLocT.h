/* $Id: DPPrimitiveLocT.h,v 1.4 2004-09-10 01:07:59 cfoster Exp $ */
/* created: myip (06/01/1999)                                      */

/*
 * Base class for Drucker-Prager, nonassociative, small-strain,
 * pressure dependent plastic model with linear isotropic hardening
 * with localization
 */

#ifndef _DP_PRIMITIVELOCT_H_
#define _DP_PRIMITIVELOCT_H_

/* base class */
#include "ParameterInterfaceT.h"

namespace Tahoe {

/* forward declarations */
class dSymMatrixT;

class DPPrimitiveLocT: public ParameterInterfaceT
{
public:

	/* constructor */
	DPPrimitiveLocT(void);

	/* destructor */
	virtual ~DPPrimitiveLocT(void);

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
	double H_delta(void) const { return fH_delta; };
	/*@}*/
	
protected:

  	/*
  	 * Returns the value of the yield function given the
  	 * Cauchy stress vector and state variables, where alpha is
  	 * the deviatoric stress-like internal state variable
  	 */
	double YieldCondition(const dSymMatrixT& devstress, 
				const double meanstress, double alpha) const;

protected:
	
	/** \name parameters */
	/*@{*/	
	double falpha_bar; /* cohesion-like strength parameter (falpha_bar >= 0.0) */
	double ffriction;  /* friction-like parameter (ffriction >= 0.0) */
	double fdilation;  /* dilation-like parameter (fdilation >= 0.0) */
	double fH_prime;   /* Deviatoric hardening parameter */
	double fH_delta;   /* Localized deviatoric hardening parameter (fH_delta < 0.0) */
	double fEta; /*fluidity parameter eta */
	/*@}*/	
};

} // namespace Tahoe 
#endif /* _DP_PRIMITIVELOCT_H_ */
