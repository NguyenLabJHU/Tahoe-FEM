/* $Id: LennardJones612.h,v 1.3 2004-03-17 17:55:41 paklein Exp $ */
/* created: paklein (10/30/1997) */
#ifndef _LJ_612_H_
#define _LJ_612_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

/** Lennard-Jones 6/12 function */
class LennardJones612: public C1FunctionT
{
public:

	/** \name constructors */
	/*@{*/
	LennardJones612(void);
	LennardJones612(double A);
	/*@}*/

	/** \name I/O */
	/*@{*/
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	/*@}*/
	
	/** \name returning values
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

	/** potential parameters */
	double fA;
};

} // namespace Tahoe 
#endif /* _LJ_612_H_ */
