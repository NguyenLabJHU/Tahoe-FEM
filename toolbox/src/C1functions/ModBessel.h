/* $Id: ModBessel.h,v 1.1 2002-01-28 18:39:47 dzeigle Exp $ */

#ifndef _MOD_BESS_H_
#define _MOD_BESS_H_

/* base class */
#include "C1FunctionT.h"

class ModBessel: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	ModBessel(double P);

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

private:

	/* potential parameters */
	double fP;
	
	//could add more parameters

};

#endif /* _MOD_BESS_H_ */

