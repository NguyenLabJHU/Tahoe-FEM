/* $Id: LinkedListT.h,v 1.8 2004-02-18 16:38:01 paklein Exp $ */
/* created: paklein (02/07/1996) */
#ifndef _LINKEDLIST_T_H_
#define _LINKEDLIST_T_H_

/* ANSI headers */
#include <iostream.h>

#include "ExceptionT.h"

/* direct members */
#include "ListNodeT.h"

namespace Tahoe {

/** basic templated linked list class
 * \note the TYPE stored in the list should have an appropriate
 * copy constructor, an assignment ("=") operator, and a
 * stream insertion ("<<"), and a test for equality ("==").
 */
template <class TYPE>
class LinkedListT
{
public:

	/** \name constructors */
	/*@{*/
	LinkedListT(void);
	LinkedListT(const LinkedListT<TYPE>& source);
	/*@}*/

	/** destructor */
	~LinkedListT(void);
	 	
	/** append value to the end of the list */
	void Append(const TYPE& value);
	
	/** append a given number of values */
	void AppendArray(int length, TYPE* values);
	
	/** append a given number of the same value */
	void AppendArray(int length, TYPE& value);

	/** append value if not already in the list (using "==").
	 * returns 1 if the value was added, else returns 0 */
	int AppendUnique(const TYPE& value);
	
	/** insert values at the current position in the list */
	void InsertAtCurrent(const TYPE& value);

	/** delete/insert values at - error if out of range */
	void InsertAt(const TYPE& value, int position);
	void DeleteAt(int position);
		
	/** \name top/next loops */
	/*@{*/
	void Top(void);
	
	int AtTop(void) { return fAtTop; };
	
	/** increment the list pointer and copy the contents of the next element into value */
	int Next(TYPE& value);
	
	/** return a pointer to the next value in the list. Returns NULL if there are no 
	 * more values */
	TYPE* Next(void);
	/*@}*/
				
	/** clears list contents */
	void Clear(void);	

	/** return 1 if the list is empty else return 0 */
	int	IsEmpty(void) const;
		
	/** returns the number of values in the list */	
	int Length(void) const;
	
	/** returns pointer to the current value in the list or NULL if empty or undefined */
	TYPE* CurrentValue(void) const;
	
	/** returns pointer to the next node's value or NULL */
	TYPE* PeekAhead(void) const;
	
	/** shallow copy */
	void Alias(const LinkedListT<TYPE>& RHS);
	
	/** deep or shallow? */
	bool IsAllocated(void) const;
	
	/** assignment operator */
	LinkedListT<TYPE>& operator=(const LinkedListT<TYPE>& source);
	
private:

	/** \name pointers into the list */
	/*@{*/
	ListNodeT<TYPE>* fCurrPtr;	/**< pointer to the current node in the list */
	ListNodeT<TYPE>* fFirstPtr;	/**< pointer to the first node in the list   */
	ListNodeT<TYPE>* fLastPtr;	/**< pointer to the last node in the list    */
	/*@}*/
	
	/** flag to indicate that the list has been reset */	
	int fAtTop;
	
	/** flag to indicate that this owns the memory */
	int fDelete;
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
}

/*************************************************************************
* Implementation
*************************************************************************/

template <class TYPE>
LinkedListT<TYPE>::LinkedListT(void):
	fCurrPtr(NULL),
	fFirstPtr(NULL),
	fLastPtr(NULL),
	fAtTop(0),
	fDelete(1)
{

}

template <class TYPE>
LinkedListT<TYPE>::LinkedListT(const LinkedListT<TYPE>& source):
	fCurrPtr(NULL),
	fFirstPtr(NULL),
	fLastPtr(NULL),
	fAtTop(0),
	fDelete(1)
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
	if (fDelete)
	{
		ListNodeT<TYPE>* newptr = new ListNodeT<TYPE>(value);
		if (!newptr) throw ExceptionT::kOutOfMemory;
	
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
	else
		ExceptionT::GeneralFail("LinkedListT::Append","Cannot Append to shallow list");
}

template <class TYPE>
void LinkedListT<TYPE>::AppendArray(int length, TYPE* values)
{
	if (fDelete)
	{
		if (length)
		{
			/* guarantee list is non-empty before looping */
			ListNodeT<TYPE>* newptr = new ListNodeT<TYPE>(*values++);
			if (!newptr) throw ExceptionT::kOutOfMemory;
	
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
			length--;
			
			for (int i = 0; i < length; i++)
			{
				ListNodeT<TYPE>* newptr = new ListNodeT<TYPE>(*values++);
				if (!newptr) throw ExceptionT::kOutOfMemory;	
			
				fLastPtr->fNextPtr = newptr;
				fLastPtr = newptr;
			}
		}
	}
	else
		ExceptionT::GeneralFail("LinkedListT::AppendArray","Cannot append to shallow Array \n");
}

template <class TYPE>
void LinkedListT<TYPE>::AppendArray(int length, TYPE& value)
{
	if (fDelete)
	{
		if (length)
		{
			/* guarantee list is non-empty before looping */
			ListNodeT<TYPE>* newptr = new ListNodeT<TYPE>(value);
			if (!newptr) throw ExceptionT::kOutOfMemory;
		
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
			length--;
			for (int i = 0; i < length; i++)
			{
				ListNodeT<TYPE>* newptr = new ListNodeT<TYPE>(value);
				if (!newptr) throw ExceptionT::kOutOfMemory;	
				
				fLastPtr->fNextPtr = newptr;
				fLastPtr = newptr;
			}
		}
	}
	else
		ExceptionT::GeneralFail("LinkedListT::AppendArray","Cannot append to shallow list");
}


template <class TYPE>
int LinkedListT<TYPE>::AppendUnique(const TYPE &value)
{
	if (fDelete)
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
	else
		ExceptionT::GeneralFail("LinkedListT::AppendUnique","Cannot append to shallow list");

	/* dummy */
	return 0;
}

/* insert value at the current node of the linked list */
template <class TYPE>
void LinkedListT<TYPE>::InsertAtCurrent(const TYPE& value)
{
	ListNodeT<TYPE>* newptr = new ListNodeT<TYPE>(value);
	if (!newptr) throw ExceptionT::kOutOfMemory;
	
	if (fFirstPtr == NULL)
	{
		/* empty list */
		fFirstPtr = newptr;
		fLastPtr  = newptr;
	}
	else 	
	{
		/* If insertion follows a call to Top(), nothing is set but this flag */
		if (fAtTop)
		{
			newptr->fNextPtr = fFirstPtr;
			fFirstPtr = newptr;
		}
		else
		{
			newptr->fNextPtr = fCurrPtr->fNextPtr;
			fCurrPtr->fNextPtr = newptr;
			fCurrPtr = fCurrPtr->fNextPtr;
			if (!fCurrPtr->fNextPtr)
				fLastPtr = fCurrPtr;
		}
	}
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
	if (currPtr == NULL) throw ExceptionT::kGeneralFail;

	/* new list node */
	ListNodeT<TYPE>* newptr = new ListNodeT<TYPE>(value);
	if (!newptr) throw ExceptionT::kOutOfMemory;
	
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
	if (fDelete)
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
		if (currPtr == NULL) throw ExceptionT::kGeneralFail;
	
		if (lastPtr == NULL)
			fFirstPtr = currPtr->fNextPtr;
		else
			lastPtr->fNextPtr = currPtr->fNextPtr;

		if (currPtr == fLastPtr) fLastPtr = lastPtr;

		delete currPtr;	
	}
	else
		ExceptionT::GeneralFail("LinkedListT::DeleteAt","Cannot delete from shallow list");
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
TYPE* LinkedListT<TYPE>::Next(void)
{
	if (fFirstPtr == NULL)
		return NULL;
	else
	{
		if (fAtTop)
		{
			fAtTop = 0;
			fCurrPtr = fFirstPtr;
		}
		else if (fCurrPtr == fLastPtr)
			return NULL;
		else
			fCurrPtr = fCurrPtr->fNextPtr;
	
		return &(fCurrPtr->fValue);
	}
}

template <class TYPE>
void LinkedListT<TYPE>::Clear(void)
{
	if (fDelete)
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
	}

	fFirstPtr = NULL;
	fLastPtr  = NULL;
	fCurrPtr  = NULL;
	
	fDelete = 1;
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

template <class TYPE>
TYPE* LinkedListT<TYPE>::CurrentValue(void) const
{
	if (!fFirstPtr) 
		return NULL;
	
	if (fAtTop)
		return &(fFirstPtr->fValue);
	else
	{
		if (!fCurrPtr) // this is actually an error, but ignore it for now
			return NULL;
	
		return &(fCurrPtr->fValue);
	}
}

template <class TYPE>
TYPE* LinkedListT<TYPE>::PeekAhead(void) const
{
	if (!fCurrPtr)
		return NULL;
	
	if (!fCurrPtr->fNextPtr)
		return NULL;
		
	return &(fCurrPtr->fNextPtr->fValue);
		
}

template <class TYPE>
void LinkedListT<TYPE>::Alias(const LinkedListT<TYPE>& RHS)
{
	if (this != &RHS)
	{
		Clear();
	
		fFirstPtr = RHS.fFirstPtr;
		fLastPtr = RHS.fLastPtr;
		fCurrPtr = RHS.fCurrPtr;
		fAtTop = RHS.fAtTop;
		fDelete = 0;
	}
}

/* deep or shallow ? */
template <class TYPE>
inline bool LinkedListT<TYPE>::IsAllocated(void) const { return fDelete != 0; }

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
	
} // namespace Tahoe 
#endif /* _LINKEDLIST_T_H_ */
