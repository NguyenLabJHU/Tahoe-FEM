/* $Id: C1FunctionT.h,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (12/04/1996)                                          */
/* C1FunctionT.h                                                          */
/* Interface for a twice differentiable function.                         */

#ifndef _C2_FUNCTION_T_H_
#define _C2_FUNCTION_T_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class dArrayT;

class C1FunctionT
{
public:

	/* function codes of derived classes */
	enum TypesT {kLennardJones = 0,
	            kSmithFerrante = 1,
	                 kGaoKlein = 2,
	                kQuadratic = 3,
	              kCubicSpline = 4};

	/* constructor */
	C1FunctionT(void);

	/* destructor */
	virtual ~C1FunctionT(void);
	
	/* I/O */
	virtual void Print(ostream& out) const = 0;     	    	   	
	virtual void PrintName(ostream& out) const = 0;     	    	   	
	    	   	    	
	/* returning values */
	virtual double Function(double x) const = 0;
	virtual double DFunction(double x) const = 0;
	virtual double DDFunction(double x) const = 0;

	/* returning values in groups - derived classes should define
	 * their own non-virtual function called within this functon
	 * which maps in to out w/o requiring a virtual function call
	 * everytime. Default behavior is just to map the virtual functions
	 * above */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;
	
	/* Return 0th, 1st, and 2nd derivative in the respective
	 * fields of the dArrayT. Default behavior is just to call the
	 * virtual functions above */  	
	virtual void SetAll(double x, dArrayT& data) const;   	
	    	   	    	
};

#endif /* _C2_FUNCTION_T_H_ */
