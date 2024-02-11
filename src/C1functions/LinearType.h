/* $Id: LinearType.h,v 1.1 2013-08-08 17:17:06 tahoe.vickynguyen Exp $ */

#ifndef _LinearType_H_
#define _LinearType_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class LinearType: public C1FunctionT
{
public:

	/* constructor */
	LinearType(double A);
	LinearType(void);

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
};


} // namespace Tahoe 
#endif /* _LinearType_H_ */
