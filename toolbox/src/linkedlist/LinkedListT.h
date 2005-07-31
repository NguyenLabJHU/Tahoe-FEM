/* $Id: LinkedListT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (02/07/1996)                                          */
/* Basic linked list template                                             */
/* Note: the TYPE stored in the list should have an appropriate           */
/* copy constructor, an assignment ("=") operator, and a                  */
/* stream insertion ("<<"), and a test for equality ("==").               */

#ifndef _LINKEDLIST_T_H_
#define _LINKEDLIST_T_H_

/* ANSI headers */
#include <iostream.h>

#include "ExceptionCodes.h"

/* direct members */
#include "ListNodeT.h"

template <class TYPE>
class LinkedListT
{
public:

	/* constructor */
	LinkedListT(void);
	LinkedListT(const LinkedListT<TYPE>& source);

	/* destructor */
	~LinkedListT(void);
	 	
	/* append value to the end of the list */
	void Append(const TYPE& value);

	/* append value if not already in the list (using "==").
	 * returns 1 if the value was added, else returns 0 */
	int AppendUnique(const TYPE& value);

	/* delete/insert values at - error if out of range */
	void InsertAt(const TYPE& value, int position);
	void DeleteAt(int position);
		
	/* top/next loops */
	void Top(void);
	int Next(TYPE& value);
				
	/* clears list contents */
	void Clear(void);	

	/* return 1 if the list is empty else return 0 */
	int	IsEmpty(void) const;
		
	/* returns the number of values in the list */	
	int Length(void) const;
	
	/* assignment operator */
	LinkedListT<TYPE>& operator=(const LinkedListT<TYPE>& source);
	
private:

	ListNodeT<TYPE>* fCurrPtr;	/* pointer to the current node in the list */
	ListNodeT<TYPE>* fFirstPtr;	/* pointer to the first node in the list   */
	ListNodeT<TYPE>* fLastPtr;	/* pointer to the last node in the list   */
	
	int fAtTop; /* flag to indicate that the list has been reset */	
};

/* output operator */
template <class TYPE>
ostream& operator<<(ostream& out, const LinkedListT<TYPE>& list)
{
	ListNodeT<TYPE>* currPtr = list.fFirstPtr;
	while (currPtr != NULL)
	{
		out << currPtr->Value() << '\n';
		currPtr = currPtr->NextPtr();
	}
	return out;
};

/*************************************************************************
* Implementation
*************************************************************************/

template <class TYPE>
LinkedListT<TYPE>::LinkedListT(void):
	fCurrPtr(NULL),
	fFirstPtr(NULL),
	fLastPtr(NULL),
	fAtTop(0)
{

}

template <class TYPE>
LinkedListT<TYPE>::LinkedListT(const LinkedListT<TYPE>& source):
	fCurrPtr(NULL),
	fFirstPtr(NULL),
	fLastPtr(NULL),
	fAtTop(0)
{
	(*this) = source;
}
	
template <class TYPE>
LinkedListT<TYPE>::~LinkedListT(void)
{
	Clear();	
} 	
	

template <class TYPE>
void LinkedListT<TYPE>::Append(const TYPE &value)
{
	ListNodeT<TYPE>* newptr = new ListNodeT<TYPE>(value);
	if (!newptr) throw(eOutOfMemory);
	
	if (fFirstPtr == NULL)
	{
		/* empty list */
		fFirstPtr = newptr;
		fLastPtr  = newptr;
	}
	else 	
	{
		/* list not empty */
		fLastPtr->fNextPtr = newptr;
		fLastPtr = newptr;
	}
}

template <class TYPE>
int LinkedListT<TYPE>::AppendUnique(const TYPE &value)
{
	/* search for matching value */
	ListNodeT<TYPE>* currPtr = fFirstPtr;
	while (currPtr != NULL)
	{
		if (currPtr->fValue == value) return 0;
		currPtr = currPtr->fNextPtr;		
	}
		
	/* value not found */	
	Append(value);
	return 1;
}

/* delete/insert values at - error if out of range */
template <class TYPE>
void LinkedListT<TYPE>::InsertAt(const TYPE& value, int position)
{
	/* advance */
	int i = 0;
	ListNodeT<TYPE>* currPtr = fFirstPtr;
	ListNodeT<TYPE>* lastPtr = NULL;
	while (i++ < position && currPtr != NULL)
	{
		lastPtr = currPtr;
		currPtr = currPtr->fNextPtr;
	}
		
	/* check */	
	if (currPtr == NULL) throw eGeneralFail;

	/* new list node */
	ListNodeT<TYPE>* newptr = new ListNodeT<TYPE>(value);
	if (!newptr) throw(eOutOfMemory);
	
	if (lastPtr == NULL)
		fFirstPtr = newPtr;
	else
		lastPtr->fNextPtr = newPtr;

	newPtr->fNextPtr = currPtr;
	if (currPtr == fLastPtr) fLastPtr = newPtr;
}

template <class TYPE>
void LinkedListT<TYPE>::DeleteAt(int position)
{
	/* advance */
	int i = 0;
	ListNodeT<TYPE>* currPtr = fFirstPtr;
	ListNodeT<TYPE>* lastPtr = NULL;
	while (i++ < position && currPtr != NULL)
	{
		lastPtr = currPtr;
		currPtr = currPtr->fNextPtr;
	}
		
	/* check */	
	if (currPtr == NULL) throw eGeneralFail;
	
	if (lastPtr == NULL)
		fFirstPtr = currPtr->fNextPtr;
	else
		lastPtr->fNextPtr = currPtr->fNextPtr;

	if (currPtr == fLastPtr) fLastPtr = lastPtr;

	delete currPtr;	
}

template <class TYPE>
inline void LinkedListT<TYPE>::Top(void) { fAtTop = 1; }

template <class TYPE>
int LinkedListT<TYPE>::Next(TYPE &value)
{
	if (fFirstPtr == NULL)
		return (0);
	else
	{
		if (fAtTop)
		{
			fAtTop = 0;
			fCurrPtr = fFirstPtr;
		}
		else if (fCurrPtr == fLastPtr)
			return (0);
		else
			fCurrPtr = fCurrPtr->fNextPtr;
	
		value = fCurrPtr->fValue;	
		return (1);
	}
}

template <class TYPE>
void LinkedListT<TYPE>::Clear(void)
{
	ListNodeT<TYPE>* currPtr;
	ListNodeT<TYPE>* tempPtr;
	
	currPtr = fFirstPtr;
	while (currPtr != NULL)
	{
		tempPtr = currPtr;
		currPtr = currPtr->fNextPtr;
		delete tempPtr;
	}

	fFirstPtr = NULL;
	fLastPtr  = NULL;
	fCurrPtr  = NULL;
}

template <class TYPE>
inline int	LinkedListT<TYPE>::IsEmpty(void) const
{
	return (fFirstPtr == NULL);
}

template <class TYPE>
int LinkedListT<TYPE>::Length(void) const
{
	int length = 0;

	ListNodeT<TYPE>* ptr = fFirstPtr;
	while (ptr != NULL)
	{
		length++;
		ptr = ptr->fNextPtr;
	}
	
	return length;
}

/*
* Assignment operator
*/
template <class TYPE>
LinkedListT<TYPE>& LinkedListT<TYPE>::operator=(const LinkedListT<TYPE>& source)
{
	/* dispose of current list entries */
	Clear();

	ListNodeT<TYPE>* currPtr = source.fFirstPtr;
	while (currPtr != NULL)
	{
		Append(currPtr->fValue);
		currPtr = currPtr->fNextPtr;
	}

	return *this;
}

#endif /* _LINKEDLIST_T_H_ */
