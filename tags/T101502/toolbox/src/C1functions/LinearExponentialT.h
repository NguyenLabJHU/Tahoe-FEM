/* $Id: LinearExponentialT.h,v 1.3 2002-07-05 22:26:16 paklein Exp $ */
/* created: paklein (05/04/2001)                                    */

#ifndef _LINEAR_EXPONENTIAL_T_H_
#define _LINEAR_EXPONENTIAL_T_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

/** implementation of the function:
 *
 *  f(x) = a + b x + c (1 - exp[-x/d])
 *
 * with parameters {a, b, c, d}
 */
class LinearExponentialT: public C1FunctionT
{
public:

	/** constructor */
	LinearExponentialT(double a, double b, double c, double d);

	/** print parameters */
	virtual void Print(ostream& out) const;

	/** print function name */
	virtual void PrintName(ostream& out) const;
	
	/** evaluate function */
	virtual double Function(double x) const;

	/** evaluate first derivative function */
	virtual double DFunction(double x) const;

	/** evaluate second derivative function */
	virtual double DDFunction(double x) const;

	/* Returning values in groups */

	/** multiple function evaluations */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;

	/** multiple first derivative evaluations */
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;

	/** multiple second derivative evaluations */
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;

private:

	/* parameters */
	double fa, fb, fc, fd;
};

} // namespace Tahoe 
#endif /* _LINEAR_EXPONENTIAL_T_H_ */
