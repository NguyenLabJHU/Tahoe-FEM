/* $Id: CubicSplineT.h,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (12/02/1996)                                          */
/* CubicSplineT.h -                                                       */
/* Interface for a piecewise cubic function.  The type of end             */
/* extensions for the spline are set by the coefficients in the           */
/* first and last rows of fCoefficients.  Note there are n+1              */
/* sets of coefficients for a spline with n knots.                        */
/* Note: CubicSplineT makes a deep copy of the data it's given.           */

#ifndef _CUBICSPLINE_T_H_
#define _CUBICSPLINE_T_H_

/* base class */
#include "C1FunctionT.h"

/* direct members */
#include "dRangeArrayT.h"
#include "dArray2DT.h"

/* forward declarations */
class dMatrixT;

class CubicSplineT: public C1FunctionT
{
public:

	enum FixityT {kParabolic = 0,  // extend with same curvature
	                kFreeRun = 1}; // extend with same slope

	/* constructor */
	CubicSplineT(const dArrayT& knots, const dArray2DT& coefficients);
	CubicSplineT(const dArray2DT& points, FixityT fixity);

	/* I/O */
	virtual void Print(ostream& out) const;     	    	   	
	virtual void PrintName(ostream& out) const;     	    	   	
	    	   	    	
	/* returning values */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;

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
	virtual void SetAll(double x, dArrayT& data) const;   	

	/* accessor to the spline information */
	const dArray2DT& Coefficients(void) const;

protected:

	/* non-virtual function calls */
	double function(double x) const;
	double Dfunction(double x) const;
	double DDfunction(double x) const;
	void all_functions(double x, double& f, double& Df, double& DDf) const;
	
private:

	/* compute spline coefficients */
	void SetSpline(const dArray2DT& points, FixityT fixity);
	
protected:   	

	dRangeArrayT fXPoints;	
	dArray2DT    fCoefficients;	
};

#endif /* _CUBICSPLINE_T_H_ */
