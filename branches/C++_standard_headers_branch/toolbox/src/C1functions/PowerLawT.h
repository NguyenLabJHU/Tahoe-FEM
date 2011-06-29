/* $Id: PowerLawT.h,v 1.2 2005-03-11 20:26:05 paklein Exp $ */
#ifndef _POWER_LAW_T_H_
#define _POWER_LAW_T_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

/** implementation of the function:
\f[
	f(x) = a (b + c x)^n
\f]
 *
 * with parameters {a, b, c, n}
 */
class PowerLawT: public C1FunctionT
{
public:

	/** constructor */
	PowerLawT(double a, double b, double c, double n);
	PowerLawT(void);
	
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

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface.
	 * \param list destination for the parameter descriptions. The list should have the
	 *        name corresponding to ParameterInterfaceT::Name. */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list.
	 * \param list input parameter list, which should be validated using ParameterInterfaceT::ValidateParameterList
	 *        to ensure the list conforms to the description defined by the interface. */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@{*/

private:

	/* parameters */
	double fa, fb, fc, fn;
};

} /* namespace Tahoe */

#endif /* _POWER_LAW_T_H_ */
