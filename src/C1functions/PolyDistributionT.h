/* $Id: PolyDistributionT.h,v 1.1 2003-06-03 16:32:12 rjones Exp $ */

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
	double Mom0  (const double x, const double d) const;
	double Mom1  (const double x, const double d) const;
	double Mom1_5(const double x, const double d) const;
	double dMom0dx  (const double x, const double d) const;
	double dMom1dx  (const double x, const double d) const;
	double dMom1dd  (const double x, const double d) const;
	double dMom1_5dx(const double x, const double d) const;
	double dMom1_5dd(const double x, const double d) const;

	/* parameters */
	double fPower;
	double fMean;
	double fWidth;
};


} // namespace Tahoe 
#endif /* _POLYDISTRIBUTION_T_H_ */
