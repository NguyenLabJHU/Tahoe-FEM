/* $Id: GaoJi2.h,v 1.2 2002-07-02 19:56:31 cjkimme Exp $ */
/* created: Baohua Ji (25/02/2002)                                            */
/* Cohesive force law:                                                        */
/* F(dl) = A*x/(1 + (x/B1)^2)^N                                               */
/* where: dl = l - L_0.                                                       */
/* 	                                                                      */

#ifndef _GAO_JI2_H_
#define _GAO_JI2_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class GaoJi2: public C1FunctionT
{
public:

	/* constructor */
	GaoJi2(double A, double B, double C, double L_0 = 1.0);

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

private:

	/* potential parameters */
	double fA;
	double fB;
	double fN;
	double fL_0; // unstretched length
//        double B1;
};

} // namespace Tahoe 
#endif /* _GAO_JI2_H_ */
