/* $Id: iArray2DT.h,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (09/23/1996)                                          */

#ifndef _IARRAY2D_T_H_
#define _IARRAY2D_T_H_

/* base class */
#include "nArray2DT.h"

class iArray2DT: public nArray2DT<int>
{
public:

	/* constructors */
	iArray2DT(void);
	iArray2DT(int majordim, int minordim);
	iArray2DT(int majordim, int minordim, int* p);
	iArray2DT(const iArray2DT& source);

	/* assignment operators */
	iArray2DT& operator=(const iArray2DT& RHS);
	iArray2DT& operator=(const int value);

	/* flagging operations */
	int	Count(int value) const;
	int HasValue(int value) const;
	int RowHasValue(int row, int value, int& column) const;
	int ColumnHasValue(int column, int value, int& row) const;
};

/* inlines */

/* assigment operators */
inline iArray2DT& iArray2DT::operator=(const iArray2DT& RHS)
{
	nArray2DT<int>::operator=(RHS);
	return *this;
}

inline iArray2DT& iArray2DT::operator=(const int value)
{
	nArray2DT<int>::operator=(value);
	return *this;
}

#endif /* _IARRAY2D_T_H_ */
