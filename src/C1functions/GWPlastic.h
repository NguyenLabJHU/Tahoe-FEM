/* $Id: GWPlastic.h,v 1.1 2003-06-03 16:32:12 rjones Exp $ */
#ifndef _GW_PLASTIC_H_
#define _GW_PLASTIC_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

class GWPlastic: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	GWPlastic(double POWER, double MU, double SIGMA);
	
	/*
	 * Destructor
	 */
	~GWPlastic();

	/*
	 * Reset parameters
	 */
	void ResetParameters(double POWER, double MU, double SIGMA);

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
#endif /* _GW_PLASTIC_H_ */
