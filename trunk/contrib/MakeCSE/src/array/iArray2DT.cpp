/*
 * File: iArray2DT.h
 */

/*
 * created      : PAK (09/23/96)
 * last modified: PAK (04/21/98)
 */

#include <iostream.h>
#include <iomanip.h>

#include "Constants.h"

#include "iArray2DT.h"

/* constructor */
iArray2DT::iArray2DT(void) { }
iArray2DT::iArray2DT(int majordim, int minordim):
	nArray2DT<int>(majordim, minordim) { }
iArray2DT::iArray2DT(int majordim, int minordim, int* p):
	nArray2DT<int>(majordim, minordim, p) { }
iArray2DT::iArray2DT(const iArray2DT& source):
	nArray2DT<int>(source) { }

int iArray2DT::Count(int value) const
{
	int* p = Pointer();
	int  count = 0;

	for (int i = 0; i < Length(); i++)
		if (*p++ == value)
			count++;

	return count;			
}

int iArray2DT::HasValue(int value) const
{
	int* p = Pointer();

	for (int i = 0; i < Length(); i++)
		if (*p++ == value) return 1;

	return 0;			
}

int iArray2DT::RowHasValue(int row, int value, int& column) const
{
	int* prow = (*this)(row);
	for (int i = 0; i < fMinorDim; i++)
		if (*prow++ == value)
		{
			column = i;
			return 1;
		}

	return 0;
}

int iArray2DT::ColumnHasValue(int column, int value, int& row) const
{
	int* pcol = Pointer(column);
	for (int i = 0; i < fMajorDim; i++)
	{
		if (*pcol == value)
		{
			row = i;
			return 1;
		}
		
		pcol += fMinorDim;
	}

	return 0;
}
