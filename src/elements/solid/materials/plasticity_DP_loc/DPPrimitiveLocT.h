/* $Id: DPPrimitiveLocT.h,v 1.2 2004-06-09 17:27:39 raregue Exp $ */
/* created: myip (06/01/1999)                                      */

/*
 * Base class for Drucker-Prager, nonassociative, small-strain,
 * pressure dependent plastic model with linear isotropic hardening
 * with localization
 */

#ifndef _DP_PRIMITIVELOCT_H_
#define _DP_PRIMITIVELOCT_H_

/* project headers */
#include "Environment.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dSymMatrixT;

class DPPrimitiveLocT
{
public:

	/* constructor */
	DPPrimitiveLocT(ifstreamT& in);

	/* destructor */
	virtual ~DPPrimitiveLocT(void);

	/* write parameters to stream */
   	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
protected:

  	/*
  	 * Returns the value of the yield function given the
  	 * Cauchy stress vector and state variables, where alpha is
  	 * the deviatoric stress-like internal state variable
  	 */
	double YieldCondition(const dSymMatrixT& devstress, 
				const double meanstress, double alpha) const;

protected:
	
	double falpha_bar; /* cohesion-like strength parameter (falpha_bar >= 0.0) */
	double ffriction;  /* friction-like parameter (ffriction >= 0.0) */
	double fdilation;  /* dilation-like parameter (fdilation >= 0.0) */
	double fH_prime;   /* Deviatoric hardening parameter */
	double fH_delta;   /* Localized deviatoric hardening parameter (fH_delta < 0.0) */
	
};

} // namespace Tahoe 
#endif /* _DP_PRIMITIVELOCT_H_ */
