/* $Id: iArray2DT.h,v 1.4 2003-11-21 22:41:30 paklein Exp $ */
/* created: paklein (09/23/1996) */
#ifndef _IARRAY2D_T_H_
#define _IARRAY2D_T_H_

/* base class */
#include "nArray2DT.h"

namespace Tahoe {

class iArray2DT: public nArray2DT<int>
{
public:

	/** \name constructors */
	/*@{*/
	iArray2DT(void);
	iArray2DT(int majordim, int minordim);
	iArray2DT(const iArray2DT& source);

	/** create an alias */
	iArray2DT(int majordim, int minordim, const int* p);
	/*@}*/

	/** \name assignment operators */
	/*@{*/
	iArray2DT& operator=(const iArray2DT& RHS);
	iArray2DT& operator=(const int value);
	/*@}*/

	/** \name flagging operations */
	/*@{*/
	int	Count(int value) const;
	int HasValue(int value) const;
	int RowHasValue(int row, int value, int& column) const;
	int ColumnHasValue(int column, int value, int& row) const;
	/*@}*/
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

} // namespace Tahoe 
#endif /* _IARRAY2D_T_H_ */
