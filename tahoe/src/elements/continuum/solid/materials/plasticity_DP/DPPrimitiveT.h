/* $Id: DPPrimitiveT.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: myip (06/01/1999)                                             */
/* Base class for a pressure dependent plastic material                   */

#ifndef _DP_PRIMITIVET_H_
#define _DP_PRIMITIVET_H_

/* project headers */
#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
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
	 * Returns the value value of the yield function given the
	 * Cauchy stress vector and state variables, where alpha_dev and
	 * alpha_vol are the deviatoric and volumetric part of internal
	 * variables, respectively.
	 */
	double YieldCondition(const dSymMatrixT& devstress, const double meanstress,
			      double alpha_dev, double alpha_vol) const; //**mien**//

protected:
	
	double falpha_bar;  /* cohesion-like strength parameter (falpha_bar >= 0.0) */
	double ffriction;   /* friction-like parameter (ffriction >= 0.0) */
	double fdilation;   /* dilation-like parameter (fdialtion >= 0.0) */
	double fH_prime;    /* Deviatoric hardening parameter */
	double fK_prime;    /* Volumetric hardening parameter */
	double fH_delta;    /* Localized deviatoric hardening parameter (fH_delta < 0.0) */
	double fK_delta;    /* Localized volumetric hardening parameter (fK_delta < 0.0) */
};

#endif /* _DP_PRIMITIVET_H_ */
