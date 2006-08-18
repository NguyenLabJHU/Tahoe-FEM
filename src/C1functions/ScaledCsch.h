/* $Id: ScaledCsch.h,v 1.1 2006-08-18 18:43:15 thao Exp $ */

#ifndef _ScaledCsch_H_
#define _ScaledCsch_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class ScaledCsch: public C1FunctionT
{
public:

	/* constructor */
	ScaledCsch(double A, double B);
	ScaledCsch(void);

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

	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);

private:
	/* parameters */
	double feta0;
	double ftau0;
};


} // namespace Tahoe 
#endif /* _ScaledCsch_H_ */
