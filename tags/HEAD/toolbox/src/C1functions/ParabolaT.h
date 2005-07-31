/* $Id: ParabolaT.h,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (03/25/1999)                                          */

#ifndef _PARABOLA_T_H_
#define _PARABOLA_T_H_

/* base class */
#include "C1FunctionT.h"

class ParabolaT: public C1FunctionT
{
public:

	/* constructor */
	ParabolaT(double k);

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
	double fk;
};

/* inlines */

/* returning values */
inline double ParabolaT::Function(double x) const { return fk*x*x; }
inline double ParabolaT::DFunction(double x) const { return fk*x; }
inline double ParabolaT::DDFunction(double x) const
{
#pragma unused(x)
	return fk;
}

#endif /* _PARABOLA_T_H_ */
