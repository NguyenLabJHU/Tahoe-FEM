/* $Id: AutoArrayT.h,v 1.1.1.1 2001-01-25 20:56:22 paklein Exp $ */
/* created: paklein (12/05/1997)                                          */
/* Array that automatically increases its dimensions when                 */
/* elements are inserted using Append() or AppendUnique.                  */
/* NOTE: currently the class only supports automatic expansion            */
/* of the memory space. a strategy for both expansion and                 */
/* contraction might be:                                                  */
/* at size b                                                              */
/* mem size is (1 + headroom) b                                           */
/* min size is b/(1 + headroom) = mem size/(1 + headroom)^2               */
/* when b > mem size -> allocate (1 + headroom) b                         */
/* when b < mem size -> allocate b/(1 + headroom)                         */
/* need to resolve when size reduction is active. working space           */
/* is often [0 b], and don't want memory operation at Reset()             */

#ifndef _AUTO_ARRAY_T_H_
#define _AUTO_ARRAY_T_H_

/* base class */
#include "ArrayT.h"

template <class TYPE>
class AutoArrayT: public ArrayT<TYPE>
{
public:

	/* constructors */
	AutoArrayT(void);
	AutoArrayT(int headroom);
	AutoArrayT(int length, int headroom);
	AutoArrayT(const ArrayT<TYPE>& source);
	AutoArrayT(const ArrayT<TYPE>& source, int headroom);

	/* set logical size using "smart" memory management */
	void Allocate(int length);
	void Resize(int new_length);
	void Resize(int new_length, const TYPE& fill);

	/* free memory (if allocated) and set size to zero */
	void Free(void);

	/* (re-)set size of the memory headroom - new headroom
	 * not used until the next memory allocation */
	void SetHeadRoom(int headroom);

	/* assignment operators - dimensions must be correct */
	AutoArrayT<TYPE>& operator=(const AutoArrayT<TYPE>& RHS); //CW was jumping to ArrayT<>::op=
	AutoArrayT<TYPE>& operator=(const ArrayT<TYPE>& RHS);
	AutoArrayT<TYPE>& operator=(const TYPE& value);
	
	/* add to the end of the logical size */
	void Append(const ArrayT<TYPE>& source);
	void Append(const TYPE& value);
	
	/* delete/insert values at - error if out of range */
	void InsertAt(const TYPE& value, int position);
	void DeleteAt(int position);

	/* returns 1 if the value is already present - NOTE: "==" must be defined
	 * for TYPE */
	bool HasValue(const TYPE& value) const;
	
	/* returns the index of the value, -1 if not present - NOTE: requires "==" */
	int PositionOf(const TYPE& value) const;
	
	/* add only if not already in the list - returns 1 if */
	/* value added 0 if not - NOTE: "==" must be defined */
	/* for TYPE */
	bool AppendUnique(const TYPE& value);
	int AppendUnique(const ArrayT<TYPE>& source); // returns the number of appended
	
	/* make *this into the union of the current contents and the values
	 * in the argument */
	void Union(const ArrayT<TYPE>& source);	
	
	/* copy logical size - OK as long as RHS >= *this in length */
	void CopyInto(ArrayT<TYPE>& RHS) const;

	/* Top/Next loop control */
	void Top(void);
	bool Next(TYPE** value);  // returns false when end of list encountered
	bool Next(void);          // just increment internal counter
	bool InRange(void) const; // returns list status without incrementing

	/* returns the current position in the list */
	const int& Position(void) const;
	TYPE& Current(void) const;

	/* stack-like operations */
	void Push(const TYPE& value);
	void Pop(void);

  private:
  
  	/* size plus headroom */
  	int WithHeadRoom(int length) { return (length*(100 + fHeadRoom))/100; };
  	
#if 0
private:

	/* safe memory copying - ie. like memcpy if fByteCopy is 1 */
	void MemCopy(TYPE* to, const TYPE* from, int length);
	void MemCopy(TYPE* to, const TYPE* from, int length) const;

	/* safe memory moving - ie. like memmove if fByteCopy is 1 */
	void MemMove(TYPE* to, const TYPE* from, int length);
#endif
	  	
private:
	
	// Size of allocated memory. fLength used to store
	// the logical size, ie. the number of initialized
	// elements in the array
	int fMemSize;

	// percent of overallocation to cut-down on calls
	// for memory allocation
	int fHeadRoom;
	
	// for Top/Next loops
	int	fCurrElement;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* default size */
const int kAutoDefSize     = 5;
const int kAutoDefHeadRoom = 20;

/* constructors */
template <class TYPE>
inline AutoArrayT<TYPE>::AutoArrayT(void):
	fMemSize(0),
	fHeadRoom(kAutoDefHeadRoom),
	fCurrElement(-1)
{
	/* initial memory */
	Allocate(kAutoDefSize);

	/* no logical size */
	fLength = 0;
}

template <class TYPE>
inline AutoArrayT<TYPE>::AutoArrayT(int headroom):
	fMemSize(0),
	fHeadRoom(headroom),
	fCurrElement(-1)
{
	/* check flags */
	if (fHeadRoom < 0) throw eGeneralFail;
}

template <class TYPE>
inline AutoArrayT<TYPE>::AutoArrayT(int length, int headroom):
	fMemSize(0),
	fHeadRoom(headroom),
	fCurrElement(-1)
{
	/* check flags */
	if (fHeadRoom < 0) throw eGeneralFail;
	Allocate(length);	
}

template <class TYPE>
inline AutoArrayT<TYPE>::AutoArrayT(const ArrayT<TYPE>& source):
	fMemSize(0),
	fHeadRoom(kAutoDefHeadRoom),
	fCurrElement(-1)
{
	operator=(source);	
}

template <class TYPE>
inline AutoArrayT<TYPE>::AutoArrayT(const ArrayT<TYPE>& source, int headroom):
	fMemSize(0),
	fHeadRoom(headroom),
	fCurrElement(-1)
{
	/* check flags */
	if (fHeadRoom < 0) throw eGeneralFail;

	operator=(source);	
}

/* set logical size using "smart" memory management */
template <class TYPE>
inline void AutoArrayT<TYPE>::Allocate(int length)
{
	/* only grows */
	if (length > fMemSize)
	{
		/* inherited */
		ArrayT<TYPE>::Allocate(WithHeadRoom(length));
	
		/* allocated size */
		fMemSize = fLength;
	}

	/* logical size */
	fLength = length;
}

template <class TYPE>
inline void AutoArrayT<TYPE>::Resize(int new_length)
{
	/* only grows */
	if (new_length > fMemSize)
	{
		/* inherited */
		ArrayT<TYPE>::Resize(WithHeadRoom(new_length));
		
		/* allocated size */
		fMemSize = fLength;
	}

	/* logical size */
	fLength = new_length;
}

template <class TYPE>
void AutoArrayT<TYPE>::Resize(int new_length, const TYPE& fill)
{
	/* only grows */
	if (new_length > fMemSize)
	{
		/* inherited */
		ArrayT<TYPE>::Resize(WithHeadRoom(new_length), fill);
		
		/* allocated size */
		fMemSize = fLength;
	}

	/* logical size */
	fLength = new_length;
}

/* free memory (if allocated) and set size to zero */
template <class TYPE>
inline void AutoArrayT<TYPE>::Free(void)
{
	/* inherited */
	ArrayT<TYPE>::Free();
	
	/* empty */
	fCurrElement = -1;
	fMemSize = 0;
}

/* (re-) set size of the memory headroom - new headroom
* not used until the next memory allocation */
template <class TYPE>
inline void AutoArrayT<TYPE>::SetHeadRoom(int headroom)
{
	fHeadRoom = headroom;
	
	/* check value */
	if (fHeadRoom < 0) throw eGeneralFail;
}

/*
* Assignment operators - dimensions must be correct
*/
template <class TYPE>
inline AutoArrayT<TYPE>& AutoArrayT<TYPE>::operator=(const TYPE& value)
{
	/* inherited */
	ArrayT<TYPE>::operator=(value);
	
	return *this;
}

template <class TYPE>
AutoArrayT<TYPE>& AutoArrayT<TYPE>::operator=(const AutoArrayT<TYPE>& RHS)
{
	/* no copies of self */
	if (fArray != RHS.Pointer())
	{
		/* set logical size */
		if (fMemSize < RHS.Length())
			Allocate(RHS.Length());
		else
			fLength = RHS.Length();
				
		/* copy data */
		MemCopy(fArray, RHS.Pointer(), fLength);
	}

	return *this;
}

template <class TYPE>
AutoArrayT<TYPE>& AutoArrayT<TYPE>::operator=(const ArrayT<TYPE>& RHS)
{
	/* no copies to self */
	if (fArray != RHS.Pointer())
	{
		/* set logical size */
		if (fMemSize < RHS.Length())
			Allocate(RHS.Length());
		else
			fLength = RHS.Length();
				
		/* copy data */
		MemCopy(fArray, RHS.Pointer(), fLength);
	}

	return *this;
}

/* resizing insertion functions */
template <class TYPE>
void AutoArrayT<TYPE>::Append(const ArrayT<TYPE>& source)
{
	/* increase memory if needed */
	if (fLength + source.Length() >= fMemSize)
	{
		/* old array */
		TYPE* olddata = fArray;
		
		/* allocate more memory (w/ headroom) */
		fMemSize += source.Length();
		fMemSize = (fMemSize*(100 + fHeadRoom))/100;
		fArray = New(fMemSize);		
		
		/* copy data into new space */
		MemCopy(fArray, olddata, fLength);
		
		/* free memory */
		delete[] olddata;
	}	

	/* append data from the list */
	MemCopy(fArray + fLength, source.Pointer(), source.Length());
	fLength += source.Length();

}

template <class TYPE>
void AutoArrayT<TYPE>::Append(const TYPE& value)
{
	if (fLength < fMemSize)
		fArray[fLength++] = value;
	else /* need more memory */
	{
		/* old array */
		TYPE* olddata = fArray;
		
		/* allocate larger block */
		fMemSize += (fMemSize*fHeadRoom)/100;
		if (fMemSize == fLength) fMemSize += 1; //insure growth
		fArray = New(fMemSize);		
		
		/* copy data into new space */
		MemCopy(fArray, olddata, fLength);
		
		/* append new value */
		fArray[fLength++] = value;
		
		/* free memory */
		delete[] olddata;
	}	
}

/* delete/insert values at - error if out of range */
template <class TYPE>
void AutoArrayT<TYPE>::InsertAt(const TYPE& value, int position)
{
	/* range check */
	if (position < 0 || position > fLength) throw eOutOfRange;

	/* empty or at end */
	if (position == fLength)
		Append(value);
	else
	{
		/* resize (by re-appending the last) */
		Append(fArray[fLength - 1]);
	
		/* move data */
		TYPE* from = fArray + position;
		TYPE*   to = from + 1;
		int length = fLength - (position + 2);
		MemMove(to, from, length);
		
		/* insert value */
		fArray[position] = value;
	}
}

template <class TYPE>
void AutoArrayT<TYPE>::DeleteAt(int position)
{
	/* range check */
	if (position < 0 || position >= fLength) throw eOutOfRange;

	/* move data */
	TYPE*   to = fArray + position;
	TYPE* from = to + 1;
	int length = fLength - (position + 1);
	MemMove(to, from, length);
	
	/* reset size */
	fLength--;
}

/* returns 1 if the value is already present - NOTE: "==" must be defined */
/* for TYPE */
template <class TYPE>
inline bool AutoArrayT<TYPE>::HasValue(const TYPE& value) const
{
	/* scan logical size for duplicates */
	TYPE* pthis = fArray;
	for (int i = 0; i < fLength; i++)
		if (*pthis++ == value)
			return true;

	return false;
}

/* returns the index of the value, -1 if not present - NOTE: requires "==" */
template <class TYPE>
inline int AutoArrayT<TYPE>::PositionOf(const TYPE& value) const
{
	/* scan logical size for duplicates */
	TYPE* pthis = fArray;
	for (int i = 0; i < fLength; i++)
		if (*pthis++ == value)
			return i;

	return -1;			
}

/* add only if no already in the list - returns 1 if */
/* value added 0 if not - NOTE: "==" must be defined */
/* for TYPE */
template <class TYPE>
bool AutoArrayT<TYPE>::AppendUnique(const TYPE& value)
{
	/* scan logical size for duplicates */
	TYPE* pthis = fArray;
	for (int i = 0; i < fLength; i++)
		if (*pthis++ == value)
			return false;
			
	/* append value on fall through */
	Append(value);			
	return true;
}

template <class TYPE>
int AutoArrayT<TYPE>::AppendUnique(const ArrayT<TYPE>& source)
{	
	TYPE* psrc = source.Pointer();
	int length = source.Length();
	int count = 0;
	for (int i = 0; i < length; i++)
		count += AppendUnique(*psrc++);
	return count;
}

/* copy logical size - OK as long as RHS >= *this in length */
template <class TYPE>
inline void AutoArrayT<TYPE>::CopyInto(ArrayT<TYPE>& RHS) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (fLength > RHS.Length()) throw eSizeMismatch;
#endif
	
	/* copy logical size */
	MemCopy(RHS.Pointer(), fArray, fLength);
}

/* Top/Next loop control */
template <class TYPE>
inline void AutoArrayT<TYPE>::Top(void) { fCurrElement = -1; }

template <class TYPE>
inline bool AutoArrayT<TYPE>::Next(TYPE** value)
{
	if (++fCurrElement < fLength)
	{
		*value = fArray + fCurrElement;
		return true;
	}
	else
		return false;
}

template <class TYPE>
inline bool AutoArrayT<TYPE>::Next(void) { return ++fCurrElement < fLength; }
template <class TYPE>
inline bool AutoArrayT<TYPE>::InRange(void) const
{
	return fCurrElement >= 0 && fCurrElement < fLength;
}

/* returns the current position in the list */
template <class TYPE>
inline const int& AutoArrayT<TYPE>::Position(void) const { return fCurrElement; }

template <class TYPE>
inline TYPE& AutoArrayT<TYPE>::Current(void) const
{
#if __option(extended_errorcheck)
	/* range check */
	if (fCurrElement < 0 || fCurrElement >= fLength) throw eOutOfRange;
#endif
	return *(fArray + fCurrElement);
}

/* stack-like operations */
template <class TYPE>
inline void AutoArrayT<TYPE>::Push(const TYPE& value)
{
	InsertAt(value, 0);
}

template <class TYPE>
inline void AutoArrayT<TYPE>::Pop(void)
{
	if (Length() > 0) DeleteAt(0);
}

/***********************************************************************
* Private
***********************************************************************/

#if 0
/* safe memory copying - ie. like memcpy if fByteCopy is 1 */
template <class TYPE>
inline void AutoArrayT<TYPE>::MemCopy(TYPE* to, const TYPE* from, int length)
{
	if (fByteCopy)	
		/* byte copy */
		memcpy(to, from, sizeof(TYPE)*length);
	else
		/* overloaded assigment operator */
		for (int i = 0; i < length; i++)
			*to++ = *from++;
}

template <class TYPE>
inline void AutoArrayT<TYPE>::MemCopy(TYPE* to, const TYPE* from, int length) const
{
	if (fByteCopy)	
		/* byte copy */
		memcpy(to, from, sizeof(TYPE)*length);
	else
		/* overloaded assigment operator */
		for (int i = 0; i < length; i++)
			*to++ = *from++;
}

/* safe memory moving - ie. like memmove if fByteCopy is 1 */
template <class TYPE>
void AutoArrayT<TYPE>::MemMove(TYPE* to, const TYPE* from, int length)
{
	if (fByteCopy)	
		memmove(to, from, sizeof(TYPE)*length);
	else
	{
		/* copy from top */
		if (from > to)
		{
			for (int i = 0; i < length; i++)
				*to++ = *from++; /* requires assigment operator */
		}
		/* copy from back */
		else
		{
			to   += length - 1;
			from += length - 1; 		
			for (int i = 0; i < length; i++)
				*to-- = *from--; /* requires assigment operator */
		}
	}
}
#endif

#endif /* _AUTO_ARRAY_T_H_ */
