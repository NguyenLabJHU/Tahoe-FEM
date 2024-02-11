/* $Id: MajumdarBhushan.h,v 1.2 2003-05-12 22:01:28 dzeigle Exp $ */
#ifndef _MAJUMDAR_BHUSHAN_H_
#define _MAJUMDAR_BHUSHAN_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

class MajumdarBhushan: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	MajumdarBhushan(double FRACDIM, double SIGMA, double C);
	
	/*
	 * Destructor
	 */
	~MajumdarBhushan();

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
	double fD;
	double fS;
	double fC;
};

} // namespace Tahoe 
#endif /* _MAJUMDAR_BHUSHAN_H_ */
