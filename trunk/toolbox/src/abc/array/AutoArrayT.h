/* $Id: AutoArrayT.h,v 1.18 2004-04-23 20:21:53 paklein Exp $ */
/* created: paklein (12/05/1997) */
#ifndef _AUTO_ARRAY_T_H_
#define _AUTO_ARRAY_T_H_

/* base class */
#include "ArrayT.h"

namespace Tahoe {

/** array class with "smart" memory management. An allocation headroom
 * can be specified to make the allocated array larger than the logical
 * size of the array. This allows the array to change dimension without
 * always requiring calls to memory allocation routines.
 * \note currently the class only supports automatic expansion
 * of the memory space. a strategy for both expansion and
 * contraction might be:
 * at size b
 * mem size is (1 + headroom) b
 * min size is b/(1 + headroom) = mem size/(1 + headroom)^2
 * when b > mem size -> allocate (1 + headroom) b
 * when b < mem size -> allocate b/(1 + headroom)
 * need to resolve when size reduction is active. working space
 * is often [0 b], and don't want memory operation at Reset() */
template <class TYPE>
class AutoArrayT: public ArrayT<TYPE>
{
public:

	/** \name constructors */
	/*@{*/
	AutoArrayT(void);
	explicit AutoArrayT(int headroom);
	AutoArrayT(int length, int headroom);
	AutoArrayT(const ArrayT<TYPE>& source);
	AutoArrayT(const ArrayT<TYPE>& source, int headroom);
	/*@}*/

	/** \name dimensioning methods */
	/*@{*/
	/** set the array size to the given length. No change occurs if the array
	 * is already the specified length. The previous contents of the array is
	 * not preserved. To preserve the array contents while changing the dimension
	 * use AutoArrayT::Resize. */
	void Dimension(int length);
	
	/** dimensions this array to the same length as the source, but no data is copied */
	void Dimension(const ArrayT<TYPE>& source) { Dimension(source.Length()); };
	
	/** \deprecated replaced by AutoArrayT::Dimension on 02/13/2002 */
	void Allocate(int length) { Dimension(length); };

	/** dimension the array to the new length keeping as much of the previous
	 * data as fits in the new space. Added space is not initialized. */
	void Resize(int new_length);
	void Resize(int new_length, const TYPE& fill);

	/* free memory (if allocated) and set size to zero */
	void Free(void);

	/** (re-)set size of the memory headroom - new headroom
	 * not used until the next memory allocation */
	void SetHeadRoom(int headroom);
	int HeadRoom(void) const { return fHeadRoom; };
	/*@}*/

	/** \name assignment operators */
	/*@{*/
	AutoArrayT<TYPE>& operator=(const AutoArrayT<TYPE>& RHS);
	AutoArrayT<TYPE>& operator=(const ArrayT<TYPE>& RHS);
	AutoArrayT<TYPE>& operator=(const TYPE& value);
	/*@}*/
	
	/** \name add to the end of the logical size */
	/*@{*/
	void Append(const ArrayT<TYPE>& source);
	void Append(const TYPE& value);
	/*@}*/
	
	/** \name delete/insert values
	 * Methods to operate on specific locations in the logical size. 
     * Error if out of range. */
	/*@{*/
	void InsertAt(const TYPE& value, int position);
	void DeleteAt(int position);
	/*@}*/

	/** returns 1 if the value is already present - NOTE: "==" must be defined
	 * for TYPE */
	bool HasValue(const TYPE& value) const;
	
	/** returns the index of the value, -1 if not present - NOTE: requires "==" */
	int PositionOf(const TYPE& value) const;
	
	/** \name appending unique values
	 * add only if not already in the list - returns 1 if
     * value added 0 if not - NOTE: "==" must be defined
	 * for TYPE */
	/*@{*/
	/** append value if not matching any other entries with operator==(). */
	bool AppendUnique(const TYPE& value);
	
	/** append unique using the given comparator function.
	 * \param value potential new value in array
	 * \param comparator function that returns true if a and b are not unique */
	bool AppendUnique(const TYPE& value, bool (*comp)(const TYPE& a, const TYPE& b));

	/** append unique values from the list using operator==(). 
	 * \return the number of appended values */
	int AppendUnique(const ArrayT<TYPE>& source);
	/*@}*/
	
	/** make *this into the union of the current contents and the values
	 * in the argument */
	void Union(const ArrayT<TYPE>& source);	
	
	/** copy logical size - OK as long as RHS >= *this in length */
	void CopyInto(ArrayT<TYPE>& RHS) const;

	/** \name Top/Next loop control */
	/*@{*/
	void Top(void);
	
	/** returns false when end of list encountered */
	bool Next(TYPE** value);  

	/** just increment internal counter */
	bool Next(void);
	
	/** returns list status without incrementing */
	bool InRange(void) const;
	/*@}*/

	/** returns the current position in the list */
	const int& Position(void) const;

	/** returns a reference to the current element in the list */
	const TYPE& Current(void) const;

	/** returns a reference to the current element in the list */
	TYPE& Current(void);
	
	/** set the position of the current pointer in the list and return a reference
	 * to the element in that position */
	TYPE& Current(int position);
	
	/** \name stack-like operations */
	/*@{*/
	void Push(const TYPE& value);
	void Pop(void);
	/*@}*/

  private:
  
  	/** size plus headroom */
  	int WithHeadRoom(int length) { return (length*(100 + fHeadRoom))/100; };
  	
private:
	
	/** size of allocated memory. fLength used to store
     * the logical size, ie. the number of initialized
     * elements in the array. */
	int fMemSize;

	/** percent of overallocation to cut-down on calls
     * for memory allocation. */
	int fHeadRoom;
	
	/** for Top/Next loops */
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

}

template <class TYPE>
inline AutoArrayT<TYPE>::AutoArrayT(int headroom):
	fMemSize(0),
	fHeadRoom(headroom),
	fCurrElement(-1)
{
	/* check flags */
	if (fHeadRoom < 0) ExceptionT::GeneralFail();
}

template <class TYPE>
inline AutoArrayT<TYPE>::AutoArrayT(int length, int headroom):
	fMemSize(0),
	fHeadRoom(headroom),
	fCurrElement(-1)
{
	/* check flags */
	if (fHeadRoom < 0) ExceptionT::GeneralFail();
	Dimension(length);	
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
	if (fHeadRoom < 0) ExceptionT::GeneralFail();

	operator=(source);	
}

/* set logical size using "smart" memory management */
template <class TYPE>
inline void AutoArrayT<TYPE>::Dimension(int length)
{
	/* only grows */
	if (length > fMemSize)
	{
		/* inherited */
		ArrayT<TYPE>::Dimension(WithHeadRoom(length));
	
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
	if (fHeadRoom < 0) ExceptionT::GeneralFail();
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
			Dimension(RHS.Length());
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
			Dimension(RHS.Length());
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
		/* owns memory */
		bool was_allocated = IsAllocated();

		/* who owns the existing memory */
		int old_length = fLength;
		TYPE* olddata;
		if (was_allocated)
			ReleasePointer(&olddata);
		else
			olddata = Pointer();
		
		/* allocate more memory */
		Dimension(fMemSize + source.Length());
		
		/* reset logical size */
		fLength = old_length;
				
		/* copy data into new space */
		MemCopy(fArray, olddata, fLength);
		
		/* free memory */
		if (was_allocated) delete[] olddata;
	}	

	/* append data from the list */
	MemCopy(fArray + fLength, source.Pointer(), source.Length());

	/* reset logical size */
	fLength += source.Length();
}

template <class TYPE>
void AutoArrayT<TYPE>::Append(const TYPE& value)
{
	if (fLength < fMemSize)
		fArray[fLength++] = value;
	else /* need more memory */
	{
		/* owns memory */
		bool was_allocated = IsAllocated();

		/* who owns the existing memory */
		int old_length = fLength;
		TYPE* olddata;
		if (was_allocated)
			ReleasePointer(&olddata);
		else
			olddata = Pointer();

		/* allocate larger block */
		Dimension(fMemSize + 1);		

		/* reset logical size */
		fLength = old_length;
		
		/* copy data into new space */
		MemCopy(fArray, olddata, fLength);
		
		/* append new value */
		fArray[fLength++] = value;
		
		/* free memory */
		if (was_allocated) delete[] olddata;
	}	
}

/* delete/insert values at - error if out of range */
template <class TYPE>
void AutoArrayT<TYPE>::InsertAt(const TYPE& value, int position)
{
	/* range check */
	if (position < 0 || position > fLength) ExceptionT::OutOfRange();

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
	if (position < 0 || position >= fLength) ExceptionT::OutOfRange();

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

/* append unique using the given comparator function */
template <class TYPE>
bool AutoArrayT<TYPE>::AppendUnique(const TYPE& value, bool (*comp)(const TYPE& a, const TYPE& b))
{
	/* scan logical size for duplicates */
	TYPE* pthis = fArray;
	for (int i = 0; i < fLength; i++)
		if ((*comp)(*pthis++, value))
			return false;
			
	/* append value on fall through */
	Append(value);			
	return true;
}

template <class TYPE>
int AutoArrayT<TYPE>::AppendUnique(const ArrayT<TYPE>& source)
{	
	const TYPE* psrc = source.Pointer();
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
	if (fLength > RHS.Length()) ExceptionT::SizeMismatch();
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
	return fCurrElement > -1 && fCurrElement < fLength;
}

/* returns the current position in the list */
template <class TYPE>
inline const int& AutoArrayT<TYPE>::Position(void) const { return fCurrElement; }

template <class TYPE>
inline const TYPE& AutoArrayT<TYPE>::Current(void) const
{
#if __option(extended_errorcheck)
	/* range check */
	if (fCurrElement < 0 || fCurrElement >= fLength) 
		ExceptionT::OutOfRange("AutoArrayT<TYPE>::Current", 
			"position is out of range: %d", fCurrElement);
#endif
	return *(fArray + fCurrElement);
}

template <class TYPE>
inline TYPE& AutoArrayT<TYPE>::Current(void)
{
#if __option(extended_errorcheck)
	/* range check */
	if (fCurrElement < 0 || fCurrElement >= fLength) 
		ExceptionT::OutOfRange("AutoArrayT<TYPE>::Current", 
			"position is out of range: %d", fCurrElement);
#endif
	return *(fArray + fCurrElement);
}

template <class TYPE>
inline TYPE& AutoArrayT<TYPE>::Current(int position)
{
#if __option(extended_errorcheck)
	/* range check */
	if (position < 0 || position >= fLength) 
		ExceptionT::OutOfRange("AutoArrayT<TYPE>::Current", 
			"position is out of range: %d", position);
#endif
	fCurrElement = position;
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

} /* namespace Tahoe */

#endif /* _AUTO_ARRAY_T_H_ */
