/* $Id: dRangeArrayT.cpp,v 1.6 2002-12-05 08:23:03 paklein Exp $ */
/* created: paklein (12/02/1996) */

#include "dRangeArrayT.h"
#include "dArray2DT.h"

/* constructor */

using namespace Tahoe;

dRangeArrayT::dRangeArrayT(void) { }

dRangeArrayT::dRangeArrayT(const dArrayT& values)
{
	SetValues(values);
}

dRangeArrayT::dRangeArrayT(int colnum, const dArray2DT& values2D)
{
	/* allocate space */
	Dimension(values2D.MajorDim());
	
	/* copy values */
	values2D.ColumnCopy(colnum,*this);

	/* check */
	if (!IsSequential())
		ExceptionT::GeneralFail("dRangeArrayT::dRangeArrayT", 
			"array values must be sorted in ascending order"); 
}

namespace Tahoe {

/* I/O operators */
ostream& operator<<(ostream& out, const dRangeArrayT& array)
{
	/* inherited */
	const dArrayT& temp = array;
	return (out << temp);
}

} // namespace Tahoe

/* set values */
void dRangeArrayT::SetValues(const dArrayT& values)
{
	dArrayT::operator=(values);
	if (!IsSequential())
		ExceptionT::GeneralFail("dRangeArrayT::SetValues", 
			"array values must be sorted in ascending order"); 
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
