/* $Id: LennardJones612.h,v 1.2 2002-07-02 19:56:31 cjkimme Exp $ */
/* created: paklein (10/30/1997)                                          */

#ifndef _LJ_612_H_
#define _LJ_612_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class LennardJones612: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	LennardJones612(double A);

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
	double fA;
	
	//could add more parameters

};

} // namespace Tahoe 
#endif /* _LJ_612_H_ */
