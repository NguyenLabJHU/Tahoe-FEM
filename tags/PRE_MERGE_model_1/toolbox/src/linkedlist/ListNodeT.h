/* $Id: ListNodeT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (02/07/1996)                                          */
/* Container class for LinkedListT                                        */
/* Note: the TYPE stored in the list should have an appropriate           */
/* copy constructor                                                       */

#ifndef _LISTNODE_T_H_
#define _LISTNODE_T_H_

template <class TYPE> class LinkedListT;

template <class TYPE>
class ListNodeT
{
	friend class LinkedListT<TYPE>;

private:

	/* constructor */
	ListNodeT(const TYPE &value);

public:
	 	 	
	/* return data from the node */
	const TYPE& Value(void) const;

	/* pointer to the next node */
	ListNodeT<TYPE>* NextPtr(void) const;
	 	 	
private:

	TYPE		fValue;
	ListNodeT*	fNextPtr;
};

/*************************************************************************
*
* Implementation
*
*************************************************************************/

/* constructor */
template <class TYPE>
inline ListNodeT<TYPE>::ListNodeT(const TYPE &value):
	fValue(value),
	fNextPtr(NULL)
{
	
}

/* return data from the node */
template <class TYPE>
inline const TYPE& ListNodeT<TYPE>::Value(void) const
{
	return fValue;
}

/* pointer to the next node */
template <class TYPE>
inline ListNodeT<TYPE>* ListNodeT<TYPE>::NextPtr(void) const
{
	return fNextPtr;
}

#endif /* _LISTNODE_T_H_ */
