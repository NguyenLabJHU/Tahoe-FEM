/* $Id: ErrorFunc.h,v 1.1 2002-02-05 15:26:56 dzeigle Exp $ */

#ifndef _ERR_FUN_H_
#define _ERR_FUN_H_

/* base class */
#include "C1FunctionT.h"

class ErrorFunc: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	ErrorFunc(double fS);

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
	double fS;

};

#endif /* _ERR_FUN_H_ */


