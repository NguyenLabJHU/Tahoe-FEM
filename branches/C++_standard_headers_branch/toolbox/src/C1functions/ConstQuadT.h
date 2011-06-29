/* $Id: ConstQuadT.h,v 1.2 2002-07-02 19:56:31 cjkimme Exp $ */
/* created: paklein (03/25/1999)                                          */

#ifndef _CONSTQUAD_T_H_
#define _CONSTQUAD_T_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class ConstQuadT: public C1FunctionT
{
public:

	/* constructor */
	ConstQuadT(double A,double B, double C);

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
	double fC;
};

/* inlines */

/* returning values */
inline double ConstQuadT::Function(double x) const 
{
	double fun; 

	if (x <= fB)
		fun = fC;
	else
		fun = fA*(x-fB)*(x-fB)+fC;
	return fun; 
}
inline double ConstQuadT::DFunction(double x) const 
{ 
	double dfun; 

	if (x <= fB)
		dfun = 0;
	else
		dfun = 2*fA*(x-fB);

	return dfun; 
}
inline double ConstQuadT::DDFunction(double x) const
{
	double ddfun; 

	if (x <= fB) 
		ddfun = 0;
	else 
		ddfun = 2*fA;

	return ddfun; 
}

} // namespace Tahoe 
#endif /* BIQUAD_T_H_ */
