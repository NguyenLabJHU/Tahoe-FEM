/* $Id: LinearT.h,v 1.3 2002-09-12 16:00:20 paklein Exp $ */

#ifndef _LINEAR_T_H_
#define _LINEAR_T_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class LinearT: public C1FunctionT
{
public:

	/* constructor */
	LinearT(double A, double B);

	/* I/O */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* returning values */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;

	/* returning values in groups */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;

private:

	/* potential parameters */
	double fA;
	double fB;
};

/* inlines */

/* returning values */
inline double LinearT::Function(double x) const 
{ return (fA*x+fB); }

inline double LinearT::DFunction(double x) const 
{ return (fA*x); }
inline double LinearT::DDFunction(double x) const
{ 
#pragma unused(x)
	return (0.0); 
}

} // namespace Tahoe 
#endif /* _LINEAR_T_H_ */



