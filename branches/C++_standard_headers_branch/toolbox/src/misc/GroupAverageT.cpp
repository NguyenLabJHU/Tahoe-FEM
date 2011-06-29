/* $Id: GroupAverageT.cpp,v 1.8 2003-11-21 22:41:57 paklein Exp $ */
/* created: paklein (10/03/1996) */
#include "GroupAverageT.h"

#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

#include "toolboxConstants.h"
#include "dArrayT.h"
#include "StringT.h"
#include "iArray2DT.h"

using namespace Tahoe;

/* constructor */
GroupAverageT::GroupAverageT(void):
	fNumRows(0),
	fIsAveraged(0),
	fCurrRow(0),
	fMemory(0)
{

}

void GroupAverageT::ResetAverage(int numcols)
{
	/* dimension workspace */
	fCounts.Dimension(fNumRows);
	fMemory.Dimension(fNumRows*numcols);
	fValues.Set(fNumRows, numcols, fMemory.Pointer());

	/* initialize data */
	fCounts = 0;
	fValues = 0.0;
	fIsAveraged = 0;
}

/* accumulate the values in the data */
void GroupAverageT::AssembleAverage(const iArrayT& nums,
	const dArray2DT& vals)
{
	/* increment count */
	for (int i = 0; i < nums.Length(); i++)
		fCounts[nums[i]]++;
		
	/* accumulate values */
	fValues.Accumulate(nums, vals);
}

const dArray2DT& GroupAverageT::OutputAverage(void)
{
	/* compute averages */
	Average();
	return fValues;
}

void GroupAverageT::OutputAverage(const iArrayT& rows,
	dArray2DT& average)
{
	/* compute averages */
	Average();

	/* collect rows */
	average.Dimension(rows.Length(), fValues.MinorDim());
	average.RowCollect(rows, fValues);
}

void GroupAverageT::OutputUsedAverage(dArray2DT& average)
{
	/* compute averages */
	Average();

	average.Dimension(NumRowsUsed(), fValues.MinorDim());
	int row = 0;
	int* pcount = fCounts.Pointer();
	for (int i = 0; i < fNumRows; i++)
		if (*pcount++ > 0)
		{
			average.SetRow(row, fValues(i));
			row++;
		}
}

/* return the number of rows in the current smoothing set */
int GroupAverageT::NumberOfAverageCols(void) const
{
	int nodecount = 0;
	const int* p = fCounts.Pointer();
	for (int i = 0; i < fNumRows; i++)
		if (*p++ > 0)
			nodecount++;
			
	return nodecount;
}

/* sequential smoothed value access */
void GroupAverageT::TopAverage(void)
{
	Average();
	fCurrRow = 0;
}

/* copy the next row into the values.  Returns -1 if there are
* no more rows, else returns row number */
int GroupAverageT::NextAverageRow(dArrayT& values)
{
	/* advance to next non-zero count */
	while (fCurrRow < fNumRows && fCounts[fCurrRow] < 1) fCurrRow++;
	
	if (fCurrRow == fNumRows)
		return -1;
	else
	{
		fValues.RowAlias(fCurrRow, values);
		return fCurrRow++;
	}
}

/* return value and position of largest (active) average value */
void GroupAverageT::MaxInColumn(int column, int& maxrow, double& maxval)
{
	/* check range */
	if (column < 0 || column >= fValues.MinorDim())
	{
		cout << "\n GroupAverageT::MaxInColumn: column out of range: "
		     << column << endl;
		throw ExceptionT::kOutOfRange;
	}

	/* compute averages */
	Average();

	/* initialize */
	maxrow = -1;
	maxval = 0.0;
	int    offset = fValues.MinorDim();
	int*    count = fCounts.Pointer();
	double* value = fValues.Pointer(column);
	for (int i = 0; i < fNumRows; i++)
	{
		if (*count > 0 && (maxrow == -1 || *value > maxval))
		{
			maxrow = i;
			maxval = *value;
		}
		
		count++;
		value += offset;
	}
}

/* returns the row numbers which where used by the current averaging
* calculation */
void GroupAverageT::RowsUsed(iArrayT& rowsused) const
{
	/* allocate */
	int count = NumRowsUsed();
	rowsused.Dimension(count);

	/* copy in used rows */
	const int* pcount = fCounts.Pointer();
	int* pused = rowsused.Pointer();
	for (int j = 0; j < fNumRows; j++)
		if (*pcount++ > 0) *pused++ = j;
}

int GroupAverageT::NumRowsUsed(void) const
{
	const int* pcount = fCounts.Pointer();
	int count = 0;
	for (int i = 0; i < fNumRows; i++)
		if (*pcount++ > 0) count++;

	return count;
}

/* clear the specified row numbers */
void GroupAverageT::ClearRows(const ArrayT<int>& rows)
{
	for (int i = 0; i < rows.Length(); i++)
		fCounts[rows[i]] = 0;
}

/**********************************************************************
* Protected
**********************************************************************/

/* divide the smoothing values by the counts */
void GroupAverageT::Average(void)
{
	if (!fIsAveraged)
	{
		for (int i = 0; i < fNumRows; i++)
		{
			int& count = fCounts[i];
		
			if (count > 0)
				fValues.ScaleRow(i, 1.0/count);
		}
		
		/* set flag */
		fIsAveraged = 1;
	}
}
