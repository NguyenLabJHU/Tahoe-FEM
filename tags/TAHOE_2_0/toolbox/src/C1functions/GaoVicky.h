/* $Id: GaoVicky.h,v 1.3 2004-06-19 23:27:18 paklein Exp $ */
/* created: Ji (12/26/1998) */
#ifndef _GAO_VICKY_H_
#define _GAO_VICKY_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

/** cohesive force law:
	\f[
		F(dr) = A dr / (1 + \exp[(-B/C + dr)/D]
	\f]
 * where: \f$ dr = l - L \f$ and  \f$ D << 1 \f$.
 */
class GaoVicky: public C1FunctionT
{
public:

	/* constructor */
	GaoVicky(double A, double B, double C, double D, double L = 1.0);
	GaoVicky(void);	

	/* I/O */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* returning values */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;

	/* returning values in groups - derived classes should define
	 * their own non-virtual function called within this functon
	 * which maps in to out w/o requiring a virtual function call
	 * everytime. Default behavior is just to map the virtual functions
	 * above */
	 
	// virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
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
	double fC;
	double fD;
	double fL; // unstretched length
};

} // namespace Tahoe 
#endif /* _GAO_VICKY_H_ */
