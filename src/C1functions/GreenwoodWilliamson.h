/* $Id: GreenwoodWilliamson.h,v 1.3 2002-02-05 15:26:56 dzeigle Exp $ */

#ifndef _GREENWOOD_WILLIAMSON_H_
#define _GREENWOOD_WILLIAMSON_H_

/* base class */
#include "C1FunctionT.h"

class GreenwoodWilliamson: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	GreenwoodWilliamson(double MU, double SIGMA);

	/*
	 * I/O
	 */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/*
	 * Returning values
	 */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;
	virtual double ContactArea(double x) const;

	/*
	 * Returning values in groups - derived classes should define
	 * their own non-virtual function called within this functon
	 * which maps in to out w/o requiring a virtual function call
	 * everytime. Default behavior is just to map the virtual functions
	 * above.
	 */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapContactArea(const dArrayT& in, dArrayT& out) const;

private:

	/* potential parameters */
	double fM;
	double fS;
};

#endif /* _GREENWOOD_WILLIAMSON_H_ */



