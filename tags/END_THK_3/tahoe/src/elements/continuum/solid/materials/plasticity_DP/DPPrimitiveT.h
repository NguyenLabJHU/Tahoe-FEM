/* $Id: DPPrimitiveT.h,v 1.7 2002-07-05 22:28:24 paklein Exp $ */
/* created: myip (06/01/1999)                                      */
/*
 * Base class for Drucker-Prager, nonassociative, small-strain,
 * pressure dependent plastic model with linear isotropic hardening
 */

#ifndef _DP_PRIMITIVET_H_
#define _DP_PRIMITIVET_H_

/* project headers */
#include "Environment.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dSymMatrixT;

class DPPrimitiveT
{
  public:

	/* constructor */
	DPPrimitiveT(ifstreamT& in);

	/* destructor */
	virtual ~DPPrimitiveT(void);

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
#endif /* _DP_PRIMITIVET_H_ */
