/* $Id: LinearExponentialT.h,v 1.4 2004-06-09 06:24:13 paklein Exp $ */
/* created: paklein (05/04/2001)                                    */

#ifndef _LINEAR_EXPONENTIAL_T_H_
#define _LINEAR_EXPONENTIAL_T_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

/** implementation of the function:
	\f[
		f(x) = a + b x + c (1 - \exp[-x/d])
	\f]
 * with parameters {a, b, c, d}
 */
class LinearExponentialT: public C1FunctionT
{
public:

	/** \name constructors */
	/*@{*/
	LinearExponentialT(double a, double b, double c, double d);
	LinearExponentialT(void);
	/*@}*/

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

	/** \name returning values in groups */
	/*@{*/
	/** multiple function evaluations */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;

	/** multiple first derivative evaluations */
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;

	/** multiple second derivative evaluations */
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@{*/

private:

	/* parameters */
	double fa, fb, fc, fd;
};

} /* namespace Tahoe */

#endif /* _LINEAR_EXPONENTIAL_T_H_ */
