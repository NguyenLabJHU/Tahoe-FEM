/* $Id: PolyDistributionT.h,v 1.2 2003-06-30 22:07:25 rjones Exp $ */

#ifndef _POLYDISTRIBUTION_T_H_
#define _POLYDISTRIBUTION_T_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class PolyDistributionT: public C1FunctionT
{
public:

	/* constructor */
	PolyDistributionT(double p, double m, double w);

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
	double fPower;
	double fMean;
	double fWidth;
};


} // namespace Tahoe 
#endif /* _POLYDISTRIBUTION_T_H_ */
