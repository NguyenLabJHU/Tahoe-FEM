/* $Id: WLCwRep.h,v 1.1 2005-04-20 23:47:27 thao Exp $ */

#ifndef _WLCwRep_H_
#define _WLCwRep_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class WLCwRep: public C1FunctionT
{
public:

	/* constructor */
	WLCwRep(double A, double B, double C, double n);

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
	/* parameters */
	double fA;
	double fB;
	double fC;
	double fn;
};


} // namespace Tahoe 
#endif /* _WORMLIKECHAIN_H_ */
