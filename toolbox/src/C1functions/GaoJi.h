/* $Id: GaoJi.h,v 1.3 2004-06-19 23:27:18 paklein Exp $ */
/* created: Baohua Ji (25/02/2002) */
#ifndef _GAO_JI_H_
#define _GAO_JI_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

/** cohesive force law:
	\f[
		F(dr) = A B C \frac{\ln (1 + \frac{dr}{B C})}{(1 + \frac{x}{B C})^C}
	\f]
 * where: \f$ dr = l - L_0 \f$.
 */
class GaoJi: public C1FunctionT
{
public:

	/* constructor */
	GaoJi(double A, double B, double C, double L_0 = 1.0);
	GaoJi(void);	

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
	double fC;
	double fL_0; // unstretched length
};

} // namespace Tahoe 
#endif /* _GAO_JI_H_ */
