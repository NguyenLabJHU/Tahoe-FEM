/*
 * File: iArrayT.h
 *
 */

/*
 * created      : PAK (08/10/96)
 * last modified: PAK (02/12/98)
 */

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
	void ChangeValue(int from, int to);
	void PrintValued(ostream& out, int value, int wrapat = 4) const;
	int	Count(int value) const;
	
	int HasValue(int value) const;	           // returns 1 if present, 0 otherwise.
	int HasValue(int value, int& index) const; // first occurence at index.
	
	/* set array value to its position in the array */
	void SetValueToPosition(void);
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

#endif /* _IARRAY_T_H_ */
