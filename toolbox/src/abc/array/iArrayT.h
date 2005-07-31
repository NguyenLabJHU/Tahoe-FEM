/* $Id: iArrayT.h,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (08/10/1996)                                          */

#ifndef _IARRAY_T_H_
#define _IARRAY_T_H_

/* base class */
#include "nArrayT.h"

class iArrayT: public nArrayT<int>
{
public:

	/* constructors */
	iArrayT(void);
	iArrayT(int length);
	iArrayT(int length, int* p);
	iArrayT(const iArrayT& source);
	
	/* assigment operators */
	iArrayT& operator=(const iArrayT& RHS);
	iArrayT& operator=(const int value);
	
	/* flagging operations */
	int ChangeValue(int from, int to); // returns the number changed
	int	Count(int value) const;	
	int HasValue(int value) const;	           // returns 1 if present, 0 otherwise.
	int HasValue(int value, int& index) const; // first occurence at index.
	
	/* set array value to its position in the array */
	void SetValueToPosition(void);
	
	/* sort both arrays by the values in master */
	void SortAscending(void);
	void SortAscending(ArrayT<int>& master);
	void SortAscending(ArrayT<double>& master);
};

/* inlines */

/* assigment operators */
inline iArrayT& iArrayT::operator=(const iArrayT& RHS)
{
	nArrayT<int>::operator=(RHS);
	return *this;
}

inline iArrayT& iArrayT::operator=(const int value)
{
	nArrayT<int>::operator=(value);
	return *this;
}

inline void iArrayT::SortAscending(void)
{
	/* inherited */
	nArrayT<int>::SortAscending();
}

inline void iArrayT::SortAscending(ArrayT<int>& master)
{
	/* inherited */
	nArrayT<int>::SortAscending(master);
}

inline void iArrayT::SortAscending(ArrayT<double>& master)
{
	/* utility */
	::SortAscending(master, *this);
}

#endif /* _IARRAY_T_H_ */
