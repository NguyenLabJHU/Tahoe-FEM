/* $Id: GreenwoodWilliamson.h,v 1.8 2003-02-03 04:40:16 paklein Exp $ */
#ifndef _GREENWOOD_WILLIAMSON_H_
#define _GREENWOOD_WILLIAMSON_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

class GreenwoodWilliamson: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	GreenwoodWilliamson(double POWER, double MU, double SIGMA);
	
	/*
	 * Destructor
	 */
	~GreenwoodWilliamson();

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
	double fM;
	double fS;
};

} // namespace Tahoe 
#endif /* _GREENWOOD_WILLIAMSON_H_ */
