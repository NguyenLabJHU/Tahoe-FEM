/* $Id: ArrayT.h,v 1.21 2005-06-05 06:18:28 paklein Exp $ */
/* created: paklein (06/19/1996) */
#ifndef _ARRAY_T_H_
#define _ARRAY_T_H_

/* Environmental */
#include "Environment.h"
#include "toolboxConstants.h"

/* ANSI headers */
#include <iostream.h>
#if defined(__JANUS__) && defined(__ROGUE_STL__)
#include <istream>
#include <ostream>
#endif
#include <string.h>
#ifdef __NEW_THROWS__
#include <new.h>
#endif

#include "ExceptionT.h"

namespace Tahoe {

/** templated class for arrays. Memory is dynamically allocated and freed.
 * Contents of the array is copied and moved using memcpy and memmove if
 * the static memory ArrayT::fByteCopy is true. Otherwise, the templated
 * type must define a propoer operator=. */
template <class TYPE>
class ArrayT
{
public:

	/** construct an empty array */
	ArrayT(void);

	/** construct a array of the specified length */
	explicit ArrayT(int length);

	/** construct an array of the specified length using the given memory.
	 * \param length dimension of the array
	 * \param TYPEPtr pointer to a block of memory with a dimension of at
	 *        least length. It is assumed that this memory will remain valid
	 *        for the lifetime of this array
	 * \note the pointer passed in to alias isn't modified by creating the
	 *       alias, but there is nothing preventing values in the data block
	 *       to be modified by this alias subsequently */
	ArrayT(int length, const TYPE* TYPEPtr);

	/** construct a deep copy of the source array */
	ArrayT(const ArrayT& source);

	/** destructor */
	~ArrayT(void);
	
	/** configure an array of the specified length using the given memory.
	 * Any memory previously owned by the array will be freed and previous
	 * content is lost.
	 * \param length dimension of the array
	 * \param TYPEPtr pointer to a block of memory with a dimension of at
	 *        least length. It is assumed that this memory will remain valid
	 *        for the lifetime of this array 
	 * \note the pointer passed in to alias isn't modified by creating the
	 *       alias, but there is nothing preventing values in the data block
	 *       to be modified by this alias subsequently */
	void Alias(int length, const TYPE* TYPEPtr);

	/** \deprecated replaced by ArrayT::Alias on 09/04/2003 */
	void Set(int length, TYPE* TYPEPtr);

	/** create shallow copy of the source array. Any memory previously owned 
	 * by the array will be freed and previous content is lost. The source array
	 * is assumed to stay valid for the lifetime of this array.
	 * \note the array passed in to alias isn't modified by creating the
	 *       alias, but there is nothing preventing values in its data
	 *       to be modified by this alias subsequently */
	void Alias(const ArrayT<TYPE>& RHS);

	/** set the array size to the given length. No change occurs if the array
	 * is already the specified length. The previous contents of the array is
	 * not preserved. To preserve the array contents while changing the dimension
	 * use ArrayT::Resize. */
	void Dimension(int length);
	
	/** dimensions this array to the same length as the source, but no data is copied */
	void Dimension(const ArrayT& source) { Dimension(source.Length()); };

	/** \deprecated replaced by ArrayT::Dimension on 02/13/2002 */
	void Allocate(int length) { Dimension(length); };

	/** return true if the array owns the memory it is using */
	bool IsAllocated(void) const;

	/** dimension the array to the new length keeping as much of the previous
	 * data as fits in the new space. Added space is not initialized. */
	void Resize(int new_length);

	/** dimension the array to the new length keeping as much of the previous
	 * data as fits in the new space 
	 * \param new_length dimension of the array following the function call
	 * \param value to initialize elements if the new length is greater than
	 *        the previous length */
	void Resize(int new_length, const TYPE& fill);
	
	/** length of the array */
	int Length(void) const;

	/** free memory (if allocated) and set size to zero */
	void Free(void);

	/** \name return a pointer specified element in the array
	 * offset must be 0 <= offset <= Length() <--- one beyond the end! */
	/*@{*/
	TYPE* Pointer(int offset = 0);
	const TYPE* Pointer(int offset = 0) const;
	/*@}*/

	/** \name element accessor */
	/*@{*/
	TYPE& operator[](int index);
	const TYPE& operator[](int index) const;
	/*@}*/

	/** reference to the first element in the array. Array must be at
	 * least length 1 */
	TYPE& First(void);
	const TYPE& First(void) const;

	/** reference to the last element in the array. Array must be at
	 * least length 1 */
	TYPE& Last(void);
	const TYPE& Last(void) const;

	/** \name assignment operators */
	/*@{*/		
	/** set all elements in the array to the given value */
	ArrayT<TYPE>& operator=(const TYPE& value);

	/** create a deep copy of the source array. */
	ArrayT<TYPE>& operator=(const ArrayT<TYPE>& RHS);
	/*@}*/		

	/** \name equality operators 
	 * Assumes TYPE has a suitably defined operator==. */
	/*@{*/
	/** element-by-element comparison */
	bool operator==(const ArrayT<TYPE>& RHS);

	/** element-by-element comparison assuming pRHS is as long as *this */
	bool operator==(const TYPE* pRHS);
	
	/** all values the same */
	bool operator==(const TYPE& value);
	/*@}*/

	/** \name inequality operators 
	 * Assumes TYPE has a suitably defined operator==. */
	/*@{*/
	/** element-by-element comparison */
	bool operator!=(const ArrayT<TYPE>& RHS);

	/** element-by-element comparison assuming pRHS is as long as *this */
	bool operator!=(const TYPE* pRHS);
	
	/** all values the same */
	bool operator!=(const TYPE& value);
	/*@}*/

	/** copy contents of the array. The valid portion of the source array is
	 * assumed to be at least as long as the length of this array */ 
	void Copy(const TYPE* pRHS);
	
	/** copy source into this array. The length of this array, accounting for
	 * the offset must be long enough to accept the contents of the source.
	 * \param offset destination in this array for the first element of source
	 * \param source source array */
	void CopyIn(int offset, const ArrayT<TYPE>& source);

	/** copy part of source into this array. The length of this array, accounting for
	 * the offsets must be long enough to accept the specified part of source.
	 * \param offset destination in this array for the first element copied 
	 *        from source
	 * \param source source array 
	 * \param source_offset offset to the first element copied from source
	 * \param length number of elements copied from source. */
	void CopyPart(int offset, const ArrayT<TYPE>& source, int source_offset, int length);

	/** \name copy selected values of source into this array
	 * \param keys indicies of the source array to copy. The number of
	 *        keys must be the same as the length of this array
	 * \param source source array. */
	/*@{*/
	void Collect(const ArrayT<int>& keys, const ArrayT<TYPE>& source);
	void Collect(const ArrayT<int>& keys, const TYPE* source);
	/*@}*/
	
	/** exchange data with the source array. Exchanges all fields of the arrays. */
	void Swap(ArrayT<TYPE>& source);

	/** transform this array into a shallow copy. The array gives up ownership of
	 * its memory. This memory must be released elsewhere. An exception is thrown
	 * if this array is already shallow 
	 * \param array returns with a pointer to the memory owned by the array */
	void ReleasePointer(TYPE** array);

	/** take ownership of the given memory. Any existing memory discarded. */
	void TakePointer(int length, TYPE* array);

	/** read array from stream as binary */
	void ReadBinary(istream& in);

	/** write array to stream as binary */
	void WriteBinary(ostream& out) const;

	/** \name useful functions for manipulating memory */
	/*@{*/
	/** return allocate memory and return a pointer (new C++ error handling)
	 * \return NULL on fail */
	static TYPE* New(int size);

	/** copying memory depending on the ArrayT:fByteCopy flag */
	static void MemCopy(TYPE* to, const TYPE* from, int length);

	/** moving memory depending on the ArrayT:fByteCopy flag */
	static void MemMove(TYPE* to, const TYPE* from, int length);
	/*@}*/
	  	
protected:
	
	int fLength;  /**< logical size of the array */
	TYPE* fArray; /**< the data */
	
private:
	
	int fDelete; /**< allocation flag */

//NOTE: could use fDelete as "memory protection" flag
//      0: not owned => no guarantees
//      1: owned     => will free on destruction
//      2: locked    => must free on destruction, i.e., no Swap or ReleasePointer

	/** must be uniquely defined at file scope if
	 * copying operations are instantiated */
	static const bool fByteCopy;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructor */
template <class TYPE>
inline ArrayT<TYPE>::ArrayT(void):
	fLength(0), 
	fArray(NULL),
	fDelete(0)
{

}

template <class TYPE>
inline ArrayT<TYPE>::ArrayT(int length):
	fLength(0), 
	fArray(NULL),
	fDelete(0)
{
	Dimension(length);	
}

template <class TYPE>
inline ArrayT<TYPE>::ArrayT(int length, const TYPE* TYPEPtr): 
	fLength(0), 
	fArray(NULL),
	fDelete(0)
{
	Alias(length, TYPEPtr);
}

template <class TYPE>
inline ArrayT<TYPE>::ArrayT(const ArrayT& source):
	fLength(0),
	fArray(NULL),
	fDelete(0)
{
	operator=(source);	
}

/* destructor */
template <class TYPE>
inline ArrayT<TYPE>::~ArrayT(void)
{
	if (fDelete)
	{
		delete[] fArray;
		fDelete = 0;	
	}
	fArray  = NULL;
	fLength = 0;
}

/* return allocate memory and return a pointer (new C++ error handling)
* throws ExceptionT::kOutOfMemory on fail */
template <class TYPE>
#if ! __option(extended_errorcheck)
inline
#endif
TYPE* ArrayT<TYPE>::New(int size)
{
	TYPE* p = NULL;
	if (size > 0)
	{
#ifdef __NEW_THROWS__
		try { p = new TYPE[size]; }
		catch (bad_alloc) { p = NULL; }
#else
		p = new TYPE[size];
#endif

		if (!p) ExceptionT::OutOfMemory("ArrayT<TYPE>::New");
	}
	else if (size < 0) 
		ExceptionT::GeneralFail("ArrayT<TYPE>::New", "invalid size: %d", size);
	
	return p;
}

/* safe memory copying - ie. like memcpy if fByteCopy is 1 */
template <class TYPE>
inline void ArrayT<TYPE>::MemCopy(TYPE* to, const TYPE* from, int length)
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
inline void ArrayT<TYPE>::MemMove(TYPE* to, const TYPE* from, int length)
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

/* set fields */
template <class TYPE>
inline void ArrayT<TYPE>::Alias(int length, const TYPE* TYPEPtr)
{
	/* release memory if allocated */
	if (fDelete) {
		delete[] fArray;
		fDelete = 0;
	}
	
#if __option(extended_errorcheck)
	const char caller[] = "ArrayT<TYPE>::Alias";
	if (length < 0) ExceptionT::GeneralFail(caller, "bad length %d", length);
	if (length > 0 && TYPEPtr == NULL) ExceptionT::GeneralFail(caller, "pointer is NULL");
#endif

	fLength = length;
	fArray = const_cast<TYPE*>(TYPEPtr);	
}

/* shallow copy/conversion */
template <class TYPE>
inline void ArrayT<TYPE>::Alias(const ArrayT<TYPE>& RHS)
{
	Alias(RHS.Length(), RHS.Pointer());
}

template <class TYPE>
inline void ArrayT<TYPE>::Set(int length, TYPE* TYPEPtr) { Alias(length, TYPEPtr); }

/* allocate an array of the specified size.*/
template <class TYPE>
void ArrayT<TYPE>::Dimension(int length)
{
	/* do nothing if the correct length already */
	if (length != fLength)
	{
		/* old data is lost if it has new length */
		if (fDelete) delete[] fArray;
	
		/* (try) allocate new memory */
		fArray = New(length);
		
		/* set dimensions */
		fLength = length;
		fDelete = (fLength > 0) ? 1 : 0;
	}
}

/* deep or shallow ? */
template <class TYPE>
inline bool ArrayT<TYPE>::IsAllocated(void) const { return fDelete != 0; }

/* if allocated, free memory and set size to zero */
template <class TYPE>
void ArrayT<TYPE>::Free(void)
{
	/* free memory */
	if (fDelete) delete[] fArray;

	/* "empty" parameters */
	fLength = 0;
	fArray  = NULL;
	fDelete = 0;
}

/* resize to new dimension, copying in at most what fits.
* extra space is initialized by specifying the fill. */
template <class TYPE>
void ArrayT<TYPE>::Resize(int new_length)
{
	/* quick return */
	if (new_length == fLength) return;

	/* abort on negative lengths */
	if (new_length < 0) 
		ExceptionT::GeneralFail("ArrayT<TYPE>::Resize", "bad size: %d", new_length);
	
	/* allocate new space */
	TYPE* new_data = New(new_length);
	int copy_size = (new_length < fLength) ? new_length : fLength;
	MemCopy(new_data, fArray, copy_size);

	/* free old memory */
	if (fDelete) delete[] fArray;
	fArray  = new_data;
	fLength = new_length;
	fDelete = 1;
}

/* with specified fill */
template <class TYPE>
void ArrayT<TYPE>::Resize(int new_length, const TYPE& fill)
{
	int old_length = fLength;
	
	/* resize the memory and copy in data */
	Resize(new_length);

	/* initialize added space */
	TYPE* pthis = fArray + old_length;		
	for (int i = old_length; i < fLength; i++)
		*pthis++ = fill;
}

/* returns the length of the array */
template <class TYPE>
inline int ArrayT<TYPE>::Length(void) const { return fLength; }

/* returns a pointer specified element in the array - offset
* must be 0 <= offset <= Length() <--- one passed the end! */
template <class TYPE>
inline TYPE* ArrayT<TYPE>::Pointer(int offset)
{
/* range checking */
#if __option (extended_errorcheck)
	if (offset < 0 || offset > fLength) 
		ExceptionT::OutOfRange("ArrayT<TYPE>::Pointer", "offset: !(0 <= %d <= %d)", offset, fLength);
#endif
	return fArray + offset;
}

template <class TYPE>
inline const TYPE* ArrayT<TYPE>::Pointer(int offset) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (offset < 0 || offset > fLength) 
		ExceptionT::OutOfRange("ArrayT<TYPE>::Pointer", "offset: !(0 <= %d <= %d)", offset, fLength);
#endif
	return fArray + offset;
}

/* element accessor */
template <class TYPE>
inline TYPE& ArrayT<TYPE>::operator[](int index)
{
/* range checking */
#if __option (extended_errorcheck)
	if (index < 0 || index >= fLength) 
		ExceptionT::OutOfRange("ArrayT<TYPE>::operator[]", "index !(0 <= %d < %d)", index, fLength);
#endif
	return fArray[index];
}

template <class TYPE>
inline const TYPE& ArrayT<TYPE>::operator[](int index) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (index < 0 || index >= fLength) 
		ExceptionT::OutOfRange("ArrayT<TYPE>::operator[]", "index !(0 <= %d < %d)", index, fLength);
#endif
	return fArray[index];
}

template <class TYPE>
inline TYPE& ArrayT<TYPE>::First(void)
{
#if __option(extended_errorcheck)
	if (fArray == NULL) ExceptionT::GeneralFail("ArrayT<TYPE>::First", "array is NULL");
#endif	
	return *fArray;
}

template <class TYPE>
inline const TYPE& ArrayT<TYPE>::First(void) const
{
#if __option(extended_errorcheck)
	if (fArray == NULL) ExceptionT::GeneralFail("ArrayT<TYPE>::First", "array is NULL");
#endif	
	return *fArray;
}

template <class TYPE>
inline TYPE& ArrayT<TYPE>::Last(void)
{
#if __option(extended_errorcheck)
	if (fArray == NULL) ExceptionT::GeneralFail("ArrayT<TYPE>::Last", "array is NULL");
#endif	
	return *(fArray + fLength - 1);
}

template <class TYPE>
inline const TYPE& ArrayT<TYPE>::Last(void) const
{
#if __option(extended_errorcheck)
	if (fArray == NULL) ExceptionT::GeneralFail("ArrayT<TYPE>::Last", "array is NULL");
#endif	
	return *(fArray + fLength - 1);
}

/* assignment operators */
template <class TYPE>
inline ArrayT<TYPE>& ArrayT<TYPE>::operator=(const TYPE& value)
{
	TYPE* p = fArray;
	for (int i = 0; i < fLength; i++)
		*p++ = value;
	return *this;
}

template <class TYPE>
inline ArrayT<TYPE>& ArrayT<TYPE>::operator=(const ArrayT<TYPE>& RHS)
{
	/* no copies to self */
	if (fArray != RHS.fArray)
	{
		/* allocate space if needed */
		if (fLength != RHS.fLength) Dimension(RHS.fLength); // old data is gone
				
		/* copy data */
		MemCopy(fArray, RHS.fArray, fLength);	
	}
	return *this;
}

/* element-by-element comparison */
template <class TYPE>
inline bool ArrayT<TYPE>::operator==(const ArrayT<TYPE>& RHS)
{
	if (fLength != RHS.fLength)
		return false;
	else
		return operator==(RHS.Pointer());
}

/* element-by-element comparison assuming pRHS is as long as *this */
template <class TYPE>
inline bool ArrayT<TYPE>::operator==(const TYPE* pRHS)
{
	if (fArray == pRHS) /* also catches this is empty && pRHS == NULL */
		return true;
	else if (fLength == 0)
		return false;
	else
	{
		TYPE* p = fArray;
		for (int i = 0; i < fLength; i++)
			if (*p++ != *pRHS++)
				return false;
	}
	return true;
}
	
/* all values the same */
template <class TYPE>
inline bool ArrayT<TYPE>::operator==(const TYPE& value)
{
	if (fLength == 0)
		return false;
	else
	{
		TYPE* p = fArray;
		for (int i = 0; i < fLength; i++)
			if (*p++ != value)
				return false;
	}
	return true;
}

/* element-by-element comparison */
template <class TYPE>
inline bool ArrayT<TYPE>::operator!=(const ArrayT<TYPE>& RHS)
{
	return !(operator==(RHS));
}

/* element-by-element comparison assuming pRHS is as long as *this */
template <class TYPE>
inline bool ArrayT<TYPE>::operator!=(const TYPE* pRHS)
{
	return !(operator==(pRHS));
}
	
/* all values the same */
template <class TYPE>
inline bool ArrayT<TYPE>::operator!=(const TYPE& value)
{
	return !(operator==(value));
}

template <class TYPE>
inline void ArrayT<TYPE>::Copy(const TYPE* pRHS)
{
	/* no copies to self */
	if (fArray != pRHS) MemCopy(fArray, pRHS, fLength);	
}

/* copy length elements of source beginning at start */
template <class TYPE>
inline void ArrayT<TYPE>::CopyIn(int offset, const ArrayT<TYPE>& source)
{
	CopyPart(offset, source, 0, source.Length());
}

template <class TYPE>
void ArrayT<TYPE>::CopyPart(int offset, const ArrayT<TYPE>& source,
	int source_offset, int length)
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if (offset + length > fLength) ExceptionT::SizeMismatch("ArrayT<TYPE>::CopyPart");
	if (source_offset + length > source.fLength) ExceptionT::OutOfRange("ArrayT<TYPE>::CopyPart");
#endif

	/* copy */
	MemCopy(fArray + offset, source.fArray + source_offset, length);
}

/* collect prescribed values from source */
template <class TYPE>
inline void ArrayT<TYPE>::Collect(const ArrayT<int>& keys, const ArrayT<TYPE>& source)
{
#if __option (extended_errorcheck)
	if (keys.Length() != Length()) ExceptionT::SizeMismatch("ArrayT<TYPE>::Collect");
#endif
	Collect(keys, source.Pointer());
}

template <class TYPE>
inline void ArrayT<TYPE>::Collect(const ArrayT<int>& keys, const TYPE* source)
{
	const int*  pkeys = keys.Pointer();
	TYPE* pthis = Pointer();
	for (int i = 0; i < Length(); i++)
		*pthis++ = source[*pkeys++];
}

/* exchange data */
template <class TYPE>
void ArrayT<TYPE>::Swap(ArrayT<TYPE>& source)
{
	int length = fLength;
	fLength = source.fLength;
	source.fLength = length;

	TYPE* array = fArray;
	fArray = source.fArray;
	source.fArray = array;
	
	int del = fDelete;
	fDelete = source.fDelete;	
	source.fDelete = del;
}

/* array no longer responsible for freeing memory. error if
* the array is already shallow. DANGER. */
template <class TYPE>
void ArrayT<TYPE>::ReleasePointer(TYPE** array)
{
	if (fArray != NULL && !fDelete)
		ExceptionT::GeneralFail("ArrayT<TYPE>::ReleasePointer", "array is shallow");
	else if (array == NULL)
		ExceptionT::GeneralFail("ArrayT<TYPE>::ReleasePointer", "array is NULL");
	else
		fDelete = 0;
	
	*array = fArray;
}

/* take ownership of the given memory. Any existing memory discarded. */
template <class TYPE>
void ArrayT<TYPE>::TakePointer(int length, TYPE* array)
{
	/* free existing */
	if (fDelete) delete[] fArray;

	/* take ownership */
	fArray  = array;
	fLength = length;
	fDelete = 1;
}

/* binary I/O */
template <class TYPE>
inline void ArrayT<TYPE>::ReadBinary(istream& in)
{
	in.read((char*) fArray, sizeof(TYPE)*fLength);
}

template <class TYPE>
inline void ArrayT<TYPE>::WriteBinary(ostream& out) const
{
	out.write((char*) fArray, sizeof(TYPE)*fLength);
}

} /* namespace Tahoe */

#endif /* _ARRAY_T_H_ */
