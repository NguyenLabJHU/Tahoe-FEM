/* $Id: LinearT.h,v 1.6 2004-03-06 17:28:32 paklein Exp $ */
#ifndef _LINEAR_T_H_
#define _LINEAR_T_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

/** specialization of C1FunctionT to a linear function. The two
 * parameter function defined by:
   \f[
		f(x) = A x + B
   \f]
 */
class LinearT: public C1FunctionT
{
public:

	/** \name constructors */
	/*@{*/
	LinearT(double A, double B);
	LinearT(void);
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

	/** \name returning values in groups */
	/*@{*/
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;
	/*@}*/

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

	/** \name function parameters */
	/*@{*/
	double fA;
	double fB;
	/*@}*/
};

/* inlines */

/* returning values */
inline double LinearT::Function(double x) const { return (fA*x+fB); }
inline double LinearT::DFunction(double x) const 
{
#pragma unused(x)
	return fA; 
}
inline double LinearT::DDFunction(double x) const
{ 
#pragma unused(x)
	return (0.0); 
}

} /* namespace Tahoe */

#endif /* _LINEAR_T_H_ */
