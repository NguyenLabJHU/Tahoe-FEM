/* $Id $ */

#ifndef _MOD_SMITH_FERRANTE_H_
#define _MOD_SMITH_FERRANTE_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class ModSmithFerrante: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	ModSmithFerrante(double A, double B);

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
	double fB;
};

} // namespace Tahoe 
#endif /* _MOD_SMITH_FERRANTE_H_ */
