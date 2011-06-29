/* $Id: GroupAverageT.h,v 1.5 2003-01-27 06:42:48 paklein Exp $ */
/* created: paklein (10/03/1996) */
#ifndef _GROUPAVERAGE_T_H_
#define _GROUPAVERAGE_T_H_

/* direct members */
#include "iArrayT.h"
#include "dArray2DT.h"
#include "AutoArrayT.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;
class StringT;
class iArray2DT;

class GroupAverageT
{
public:

	/* constructor */
	GroupAverageT(void);
	
	/* set the number of rows in the smoothing set */
	void SetNumAverageRows(int numrows);
	
	/* set for averaging the given number of columns */
	void ResetAverage(int numcols);

	/* accumulate the values in the data */
	void AssembleAverage(const iArrayT& nums, const dArray2DT& vals);
	
	/* output the averaged values */
	const dArray2DT& OutputAverage(void);                        // for all nodes
	void OutputAverage(const iArrayT& rows, dArray2DT& average); // requested rows only
	void OutputUsedAverage(dArray2DT& average); // active rows only

	/* return the number of rows in the current set */
	int NumberOfAverageCols(void) const;
	
	/* copy the next row into the values.  Returns -1 if there are
	 * no more rows, else returns row number for the returned values */
	void TopAverage(void); //triggers average
	int NextAverageRow(dArrayT& values);
	
	/* return value and position of largest (active) average value */
	void MaxInColumn(int column, int& maxrow, double& maxval); //triggers averaging

	/* returns the row numbers which where used by the current averaging
	 * calculation */
	void RowsUsed(iArrayT& rowsused) const;
	int NumRowsUsed(void) const;

	/* clear the specified row numbers */
	void ClearRows(const ArrayT<int>& rows);

protected:

	/* return the values array */
	const dArray2DT& Values(void) const;

	/* perform the averaging */
	void Average(void);

private:

	int	fNumRows;
	int	fIsAveraged;
	int	fCurrRow;

	dArray2DT	       fValues;
	AutoArrayT<int>    fCounts;
	AutoArrayT<double> fMemory;
	 	
}; /* _GROUPAVERAGE_T_H_ */

/* inlines */
/* set the number of rows in the smoothing set */
inline void GroupAverageT::SetNumAverageRows(int numrows) { fNumRows = numrows; }

/* return the values array */
inline const dArray2DT& GroupAverageT::Values(void) const { return fValues; }

} // namespace Tahoe 
#endif
