/* $Id: ErrorFunc.h,v 1.3 2002-07-02 19:56:31 cjkimme Exp $ */

#ifndef _ERR_FUN_H_
#define _ERR_FUN_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class ErrorFunc: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	ErrorFunc();
	
	/*
	 * Destructor
	 */
	~ErrorFunc();
	
	/*
	* Methods
	*/
	double gammp(double a, double x) const;
	void gcf(double *gammcf, double a, double x, double *gln) const;
	void gser(double *gamser, double a, double x, double *gln) const;

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

};

} // namespace Tahoe 
#endif /* _ERR_FUN_H_ */



