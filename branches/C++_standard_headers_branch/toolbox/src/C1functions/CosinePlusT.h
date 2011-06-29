/* $Id: CosinePlusT.h,v 1.2 2009-07-09 00:49:48 regueiro Exp $ */
#ifndef _COSINE_PLUS_T_H_
#define _COSINE_PLUS_T_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

/** implementation of the function:
\f[
	f(t) = H(t,t1-t0)*(a + b \cos (c t) + d \sin (e t) + f t \cos (g t) + p t \sin (q t))
\f]
 *
 * with parameters {t0, t1, a, b, c, d, e, f, g, p, q}
 */
class CosinePlusT: public C1FunctionT
{
public:

	/** constructor */
	CosinePlusT(double t0, double t1, double a, double b, double c, double d, double e, double f, double g,
		double p, double q);
	CosinePlusT(void);
	
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
	double ft0, ft1, fa, fb, fc, fd, fe, ff, fg, fp, fq; 
};

} /* namespace Tahoe */

#endif /* _COSINE_PLUS_T_H_ */
