/* $Id: ModBessel.h,v 1.3 2002-07-02 19:56:31 cjkimme Exp $ */

#ifndef _MOD_BESSEL_H_
#define _MOD_BESSEL_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class ModBessel: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	ModBessel(double p);
	
	/*
	 * Destructor
	 */
	~ModBessel();
	
	/*
	 * Evaluation function
	 */
	 double Eval(int d, double x) const;
	
	/*
	 * I/O
	 */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/*
	 * Returning values
	 */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;

	/*
	 * Returning values in groups - derived classes should define
	 * their own non-virtual function called within this functon
	 * which maps in to out w/o requiring a virtual function call
	 * everytime. Default behavior is just to map the virtual functions
	 * above.
	 */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;

private:

	/* potential parameters */
	double power;
    double chebev(double a, double b, double c[], int m, double x) const;
    void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi) const;
    void bessk(double x, double xnu, double *rk, double *rkp) const;

};

} // namespace Tahoe 
#endif /* _MOD_BESSEL_H_ */


