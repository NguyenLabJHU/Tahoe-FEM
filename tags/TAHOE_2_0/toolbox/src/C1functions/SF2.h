/* $Id: SF2.h,v 1.3 2004-06-19 23:27:18 paklein Exp $ */
#ifndef _SF2_H_
#define _SF2_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe{

/** cohesive force law
 * The force is given by
 \f[
	F(dr) = A dr \exp \left[ -\frac{dr^2}{B} \right]
 \f]
 * where \f$ dr = l - l_0 \f$. */
class SF2: public C1FunctionT
{
public:

	/* constructor */
	SF2(double A, double B, double l_0 = 1.0);
	SF2(void);	

	/* I/O */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* returning values */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;

	/** returning values in groups - derived classes should define
	 * their own non-virtual function called within this functon
	 * which maps in to out w/o requiring a virtual function call
	 * everytime. Default behavior is just to map the virtual functions
	 * above.
	 */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@{*/

private:

	/* potential parameters */
	double fA;
	double fB;
	double fl_0; //equilibrium length
};
}
#endif /* _SMITH_FERRANTE_H_ */
