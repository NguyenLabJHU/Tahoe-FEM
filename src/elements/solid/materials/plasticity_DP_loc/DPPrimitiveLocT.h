/* $Id: DPPrimitiveLocT.h,v 1.5 2005-04-08 19:22:46 raregue Exp $ */
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
	double H(void) const { return fH; };
	/*@}*/
	
protected:

  	/*
  	 * Returns the value of the yield function given the
  	 * Cauchy stress vector and state variables, where kappa is
  	 * the deviatoric stress-like internal state variable
  	 */
	double YieldCondition(const dSymMatrixT& devstress, 
				const double meanstress, double kappa) const;

protected:
	
	/** \name parameters */
	/*@{*/	
	double fkappa; /* cohesion-like strength parameter (fkappa >= 0.0) */
	double ffriction;  /* friction-like parameter (ffriction >= 0.0) */
	double fdilation;  /* dilation-like parameter (fdilation >= 0.0) */
	double fH;   /* Deviatoric hardening parameter */
	double fEta; /*fluidity parameter eta */
	/*@}*/	
};

} // namespace Tahoe 
#endif /* _DP_PRIMITIVELOCT_H_ */
