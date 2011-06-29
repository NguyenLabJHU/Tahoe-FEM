/* $Id: Gamma.h,v 1.2 2002-07-02 19:56:31 cjkimme Exp $ */

#ifndef _GAMMA_H_
#define _GAMMA_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class Gamma: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	Gamma();
	
	/*
	 * Destructor
	 */
	~Gamma();
	
	/*
	 * I/O
	 */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/*
	 * Returning values
	 */
	virtual double Function(double xx) const;
	virtual double DFunction(double xx) const;
	virtual double DDFunction(double xx) const;

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
	/* room for more		*/

};

} // namespace Tahoe 
#endif /* _GAMMA_H_ */



