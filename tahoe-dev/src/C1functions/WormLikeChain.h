/* $Id: WormLikeChain.h,v 1.1 2005-04-15 22:47:45 thao Exp $ */

#ifndef _WORMLIKECHAIN_H_
#define _WORMLIKECHAIN_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class WormLikeChain: public C1FunctionT
{
public:

	/* constructor */
	WormLikeChain(double A, double B);

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
};


} // namespace Tahoe 
#endif /* _WORMLIKECHAIN_H_ */
