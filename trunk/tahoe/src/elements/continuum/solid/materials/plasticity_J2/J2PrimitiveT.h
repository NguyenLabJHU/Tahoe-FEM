/* $Id: J2PrimitiveT.h,v 1.4 2004-07-15 08:28:54 paklein Exp $ */
/* created: paklein (02/17/1997) */
#ifndef _J2_PRIMITIVET_H_
#define _J2_PRIMITIVET_H_

/* base class */
#include "ParameterInterfaceT.h"

namespace Tahoe {

/* forward declarations */
class dSymMatrixT;

/** base class for a J2 plastic material with linear kinematic and
 * isotropic hardening laws defined by:
 	\f[
		H(\alpha) = (1 - \theta) \bar{H} \alpha
	\f]
	\f[
		K(\alpha) = Y + \theta \bar{H} \alpha
	\f]
 * where \f$ \alpha \f$ is the internal hardening variable
 */
class J2PrimitiveT: virtual public ParameterInterfaceT
{
public:

	/** constructor */
	J2PrimitiveT(void);

	/** destructor */
	virtual ~J2PrimitiveT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* returns the value value of the yield function given the
	 * Cauchy stress vector and state variables, where beta and alpha
	 * represent the kinematic and isotropic hardening, respectively */
	double YieldCondition(const dSymMatrixT& relstress, double alpha) const;

	/* hardening functions and their 1st derivatives.
	 *
	 *		H(a) = (1 - ftheta) fH_bar a
	 *      K(a)  = fYield + ftheta fH_bar a
	 */
	double  H(double a) const;
	double dH(double a) const;
	double  K(double a) const;
	double dK(double a) const;
	
protected:
	
	double fYield;	/* initial flow stress (fYield > 0) */
	double fH_bar;	/* hardening parameter (fH_bar > 0) */	
	double ftheta;	/* (0 < ftheta < 1) 				*/

};

/* hardening functions and their 1st derivatives */
inline double J2PrimitiveT::H(double a) const
{
	return ( (1.0 - ftheta)*fH_bar*a );
}

inline double J2PrimitiveT::dH(double a) const
{
#pragma unused(a)

	return ( (1.0 - ftheta)*fH_bar );
}

inline double J2PrimitiveT::K(double a) const
{
	return ( fYield + ftheta*fH_bar*a );
}

inline double J2PrimitiveT::dK(double a) const
{
#pragma unused(a)

	return ( ftheta*fH_bar );
}

} /* namespace Tahoe */

#endif /* _J2_PRIMITIVET_H_ */
