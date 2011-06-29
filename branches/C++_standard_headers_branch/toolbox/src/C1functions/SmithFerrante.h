/* $Id: SmithFerrante.h,v 1.5 2004-06-07 21:37:33 paklein Exp $ */
/* created: paklein (10/30/1997) */
#ifndef _SMITH_FERRANTE_H_
#define _SMITH_FERRANTE_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

/** the  cohesive force law:
 \f[
 	f(\Delta x) = a \Delta x \exp (-\Delta x/B)
 \f]
 * where \f$ \Delta x = x - x_0 \f$ 
 */
class SmithFerrante: public C1FunctionT
{
public:

	/** \name constructors */
	/*@{*/
	SmithFerrante(void);
	SmithFerrante(double A, double B, double l_0 = 1.0);
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

	/** \name Returning values in groups. derived classes should define
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

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** \name potential parameters */
	/*@{*/
	double fA;
	double fB;
	double fl_0; /**< equilibrium length */
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _SMITH_FERRANTE_H_ */
