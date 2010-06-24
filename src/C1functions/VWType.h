/* $Id: VWType.h,v 1.1 2010-06-24 13:32:38 thao Exp $ */

#ifndef _FungType_H_
#define _FungType_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class FungType: public C1FunctionT
{
public:

	/* constructor */
	FungType(double A, double B);
	FungType(void);

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
	double fA;
	double fB;
};


} // namespace Tahoe 
#endif /* _FungType_H_ */
