/* $Id: ParabolaPotT.h,v 1.1 2003-05-21 16:04:08 thao Exp $ */
/* created: paklein (03/25/1999)                                          */

#ifndef _PARABOLA_POT_T_H_
#define _PARABOLA_POT_T_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class ParabolaPotT: public C1FunctionT
{
public:

	/* constructor */
	ParabolaPotT(double k, double B, double l0=1.0);

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
	double fl0;
	double fB;
};

/* inlines */

/* returning values */
inline double ParabolaPotT::Function(double x) const { 
	return (0.5*fk*(x-fl0)*(x-fl0)-0.5*fk*fB); }
inline double ParabolaPotT::DFunction(double x) const { return fk*(x-fl0); }
inline double ParabolaPotT::DDFunction(double x) const
{
#pragma unused(x)
	return fk;
}

} // namespace Tahoe 
#endif /* _PARABOLA_POT_T_H_ */
