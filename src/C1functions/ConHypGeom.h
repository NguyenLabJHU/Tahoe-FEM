/* $Id: ConHypGeom.h,v 1.3 2003-02-03 04:40:16 paklein Exp $ */
#ifndef _CON_HYP_GEOM_H_
#define _CON_HYP_GEOM_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class ConHypGeom: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	ConHypGeom(double A, double B);
	
	/*
	 * Destructor
	 */
	~ConHypGeom();
	
	/*
	 * Methods
	 */
	double PochhammerRat(double numer, double denom, int index) const;
	double PochhammerProd(double lhs, double rhs, int index) const;

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
	double fA, fB;
	
	//could add more parameters

};

} // namespace Tahoe 
#endif /* _CON_HYP_GEOM_H_ */

