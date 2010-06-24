/* $Id: VWType.h,v 1.2 2010-06-24 13:48:30 thao Exp $ */

#ifndef _VWType_H_
#define _VWType_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class VWType: public C1FunctionT
{
public:

	/* constructor */
	VWType(double A, double B);
	VWType(void);

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
#endif /* _VWType_H_ */
