/* $Id: dRangeArrayT.cpp,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (12/02/1996)                                          */
/* dRangeArrayT.cpp                                                       */

#include "dRangeArrayT.h"
#include "dArray2DT.h"

/* constructor */
dRangeArrayT::dRangeArrayT(const dArrayT& values)
{
	dArrayT::operator=(values);
	if ( !IsSequential()) throw eGeneralFail;
}

dRangeArrayT::dRangeArrayT(int colnum, const dArray2DT& values2D)
{
	/* allocate space */
	Allocate(values2D.MajorDim());
	
	/* copy values */
	values2D.ColumnCopy(colnum,*this);
}

/* I/O operators */
ostream& operator<<(ostream& out, const dRangeArrayT& array)
{
	/* inherited */
	const dArrayT& temp = array;
	return (out << temp);
}

/*
* Return the range for the given value - the Range will is
* an integer { 0...Length() }, where 0 means the value is less
* than the first array value, Length() means it's greater than
* the last, and in general:
*
*		i(x): array[i-1] < x < array[i]
*
*/
int dRangeArrayT::Range(double value) const
{
	/* less than smallest value */
	if (value < fArray[0])
		return 0;

	/* bisection */
	int lower = 0;
	int upper = Length();
	int dex;
	do {
		dex = (lower + upper)/2;
	
		if (value > fArray[dex])
			lower = dex;
		else
			upper = dex;
		
	} while (upper > lower + 1);
	
	return upper;
}
	
/* returns 1 if the data is in ascending order */
int dRangeArrayT::IsSequential(void) const
{
	/* must be in weak ascending order */
	for (int i = 1; i < Length(); i++)
		if (fArray[i] < fArray[i-1])
			return 0;
			
	return 1;
}
