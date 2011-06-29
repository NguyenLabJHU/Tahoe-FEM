/* $Id: LennardJones612.h,v 1.5 2004-06-07 21:37:33 paklein Exp $ */
/* created: paklein (10/30/1997) */
#ifndef _LJ_612_H_
#define _LJ_612_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

/** Lennard-Jones 6/12 function
 \f[
	f(x) = a (\frac{1}{2} (x/b)^{-12} - (x/b)^{-6})
 \f]
 *
 * with parameters {a, b} */
class LennardJones612: public C1FunctionT
{
public:

	/** \name constructors */
	/*@{*/
	LennardJones612(void);
	LennardJones612(double A, double B);
	/*@}*/

	/** \name I/O */
	/*@{*/
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	/*@}*/
	
	/** \name returning values */
	/*@{*/
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;
	/*@}*/

	/** returning values in groups. Derived classes should define
	 * their own non-virtual function called within this functon
	 * which maps in to out w/o requiring a virtual function call
	 * everytime. Default behavior is just to map the virtual functions
	 * above.
	 */
	/*@{*/
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list.
	 * \param list input parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@{*/

private:

	/** \name potential parameters */
	/*@{*/
	double fA;
	double fB;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _LJ_612_H_ */
