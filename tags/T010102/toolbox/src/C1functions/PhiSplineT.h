/* $Id: PhiSplineT.h,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (01/30/2000)                                          */
/* PhiSplineT.h                                                           */
/* cubic spline form of the EAM pair potential from PRL v##, n##, 1986.   */
/* Pair potential is written as:                                          */
/* phi = Z_a(r) Z_b(r)/r                                                  */
/* Class only handles (a == b). Z(r) is given as a spline. For r > r_cut, */
/* f, f', f'' = 0 (discontinuously).                                      */

#ifndef _EAM_PHI_SPLINE_T_H_
#define _EAM_PHI_SPLINE_T_H_

/* base class */
#include "CubicSplineT.h"

class PhiSplineT: public CubicSplineT
{
public:

	/* constructor */
	PhiSplineT(const dArray2DT& points, FixityT fixity, double r_cut);

	/* returning values */
	virtual double Function(double r) const;
	virtual double DFunction(double r) const;
	virtual double DDFunction(double r) const;

	/* returning values in groups - returns refence to out to allow:
	 *
	 *	dArrayT& goodname = pfunc->MapFunction(in, tempspace);
	 */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;
	
	/*
	 * Return 0th, 1st, and 2nd derivative in the respective
	 * fields of the dArrayT.
	 */  	
	virtual void SetAll(double r, dArrayT& data) const;   	

private:

	/* cut off radius */
	double fr_cut;
};

#endif /* _EAM_PHI_SPLINE_T_H_ */
