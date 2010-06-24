/* $Id: FungType.h,v 1.4 2010-06-24 13:32:38 thao Exp $ */

#ifndef _FungType2_H_
#define _FungType2_H_

/* base class */
#include "C1FunctionT.h"

/*GA Holzapfel, TC Gasser, M Stadler, "A structural model for the viscoelastic behavior of arterial walls: continuum formulation and finite element simulations*/
/*Wf = alpha/beta ( Exp[beta( (l^2-1)/2)^2 ) ] -1 )*/
namespace Tahoe {

class FungType2: public C1FunctionT
{
public:

	/* constructor */
	FungType2(double A, double B);
	FungType2(void);

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
#endif /* _FungType2_H_ */
