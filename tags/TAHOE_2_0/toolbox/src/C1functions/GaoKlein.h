/* $Id: GaoKlein.h,v 1.2 2002-07-02 19:56:31 cjkimme Exp $ */
/* created: paklein (12/26/1998)                                          */
/* Cohesive force law:                                                    */
/* F(dr) = A dr Exp[B/(dr - C)]                                           */
/* where: dr = l - L.                                                     */
/* 	                                                                      */

#ifndef _GAO_KLEIN_H_
#define _GAO_KLEIN_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class GaoKlein: public C1FunctionT
{
public:

	/* constructor */
	GaoKlein(double A, double B, double C, double L = 1.0);

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

private:

	/* potential parameters */
	double fA;
	double fB;
	double fC;
	double fL; // unstretched length
};

} // namespace Tahoe 
#endif /* _GAO_KLEIN_H_ */
