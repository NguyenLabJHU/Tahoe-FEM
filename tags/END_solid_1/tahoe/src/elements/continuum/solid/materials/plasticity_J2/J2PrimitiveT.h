/* $Id: J2PrimitiveT.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (02/17/1997)                                          */
/* Base class for a J2 plastic material with linear kinematic/            */
/* isotropic hardening laws defined by:                                   */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */

#ifndef _J2_PRIMITIVET_H_
#define _J2_PRIMITIVET_H_

/* project headers */
#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class dSymMatrixT;

class J2PrimitiveT
{
public:

	/* constructor */
	J2PrimitiveT(ifstreamT& in);

	/* destructor */
	virtual ~J2PrimitiveT(void);

	/* output parameters to stream */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
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

/*
* Hardening functions and their 1st derivatives.
*
*		H(a) = (1 - ftheta) fH_bar a
*      K(a) = fYield + ftheta fH_bar a
*/
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

#endif /* _J2_PRIMITIVET_H_ */
