/* $Id: dRangeArrayT.h,v 1.1.1.1 2001-01-25 20:56:24 paklein Exp $ */
/* created: paklein (12/02/1996)                                          */
/* dRangeArrayT.h                                                         */
/* Must set the knot values with the constructor.  This class             */
/* does not have all the behavior of a dArrayT.  The inheritance          */
/* is PRIVATE.                                                            */

#ifndef _DRANGEARRAY_T_H_
#define _DRANGEARRAY_T_H_

/* base class */
#include "dArrayT.h"

/* forward declarations */
class dArray2DT;

class dRangeArrayT: private dArrayT
{
public:

	/* constructor */
	dRangeArrayT(const dArrayT& values);
	dRangeArrayT(int colnum, const dArray2DT& values2D);

	/* I/O operators */
	friend ostream& operator<<(ostream& out, const dRangeArrayT& array);
	    	
	/*
	 * Return the range for the given value - the Range will is
	 * an integer { 0...Length() }, where 0 means the value is less
	 * than the first array value, Length() means it's greater than
	 * the last, and in general:
	 *
	 *		i(x): array[i-1] < x < array[i]
	 *
	 */
	int Range(double value) const;
	
	/* still allow length checks */
	int Length(void) const;
	
	/* make element accessor public */
	dArrayT::operator[];
	
private:
	
	/* returns 1 if the data is in ascending order */
	int IsSequential(void) const;
	
};

/* inlines */

/*  still allow length checks */
inline int dRangeArrayT::Length(void) const { return dArrayT::Length(); }	

#endif /* _DRANGEARRAY_T_H_ */
