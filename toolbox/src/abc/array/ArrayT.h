/* $Id: ArrayT.h,v 1.1.1.1 2001-01-25 20:56:22 paklein Exp $ */
/* created: paklein (06/19/1996)                                          */
/* Base class for handling memory allocation for arrays of TYPE           */

#ifndef _ARRAY_T_H_
#define _ARRAY_T_H_

/* Environmental */
#include "Environment.h"
#include "Constants.h"

/* ANSI headers */
#include <iostream.h>
#include <string.h>
#ifdef __NEW_THROWS__
#include <new.h>
#endif

#include "ExceptionCodes.h"

template <class TYPE>
class ArrayT
{
public:

	/* constructors */
	ArrayT(void);
	ArrayT(int length);
	ArrayT(int length, TYPE* TYPEPtr);
	ArrayT(const ArrayT& source);

	/* destructor */
	~ArrayT(void);
	
	/* set fields */
	void Set(int length, TYPE* TYPEPtr);
	
	/* allocate an array of the specified size  */
	void Allocate(int length);
	bool IsAllocated(void) const;

	/* resize to new dimension, copying in at most what fits.
	 * extra space is initialized by specifying the fill. */
	void Resize(int new_length);
	void Resize(int new_length, const TYPE& fill);
	
	/* length of the array */
	int Length(void) const;

	/* free memory (if allocated) and set size to zero */
	void Free(void);

	/* returns a pointer specified element in the array - offset
	 * must be 0 <= offset <= Length() <--- one beyond the end! */
	TYPE* Pointer(int offset = 0) const;

	/* element accessor */
	TYPE& operator[](int index) const;
	TYPE& First(void) const;
	TYPE& Last(void) const;
		
	/* assignment operators */
	ArrayT<TYPE>& operator=(const TYPE& value);
	ArrayT<TYPE>& operator=(const ArrayT<TYPE>& RHS);

	/* shallow copy/conversion */
	void Alias(const ArrayT<TYPE>& RHS);
	
	/* copy length elements of source beginning at offset */
	void CopyIn(int offset, const ArrayT<TYPE>& source);
	void CopyPart(int offset, const ArrayT<TYPE>& source, int source_offset, int length);
	void Collect(const ArrayT<int>& keys, const ArrayT<TYPE>& source);
	
	/* exchange data */
	void Swap(ArrayT<TYPE>& source);

	/* array no longer responsible for freeing memory. error if
	 * the array is already shallow. DANGER. */
	void ReleasePointer(TYPE** array);

	/* binary I/O - entire length of fArray */
	void ReadBinary(istream& in);
	void WriteBinary(ostream& out) const;

protected:

	/* return allocate memory and return a pointer (new C++ error handling)
	 * returns NULL on fail */
	TYPE* New(int size) const;

	/* safe memory copying - ie. like memcpy if fByteCopy is true */
	static void MemCopy(TYPE* to, const TYPE* from, int length);
	static void MemMove(TYPE* to, const TYPE* from, int length);
	  	
protected:
	
	int fLength;  /* size */
	TYPE* fArray; /* data */
	
private:
	
	int fDelete; /* allocation flag */

//NOTE: could use fDelete as "memory protection" flag
//      0: not owned => no guarantees
//      1: owned     => will free on destruction
//      2: locked    => must free on destruction, i.e., no Swap or ReleasePointer

	/* must be uniquely defined at file scope if
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
	Allocate(length);	
}

template <class TYPE>
inline ArrayT<TYPE>::ArrayT(int length, TYPE* TYPEPtr): 
	fLength(length), 
	fArray(TYPEPtr),
	fDelete(0)
{

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
* throws eOutOfMemory on fail */
template <class TYPE>
#if ! __option(extended_errorcheck)
inline
#endif
TYPE* ArrayT<TYPE>::New(int size) const
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

		if (!p)
		{
			cout << "\n ArrayT<TYPE>::New: out of memory" << endl;
			throw eOutOfMemory;
		}
	}
	else if (size < 0)
	{
		cout << "\n ArrayT<TYPE>::New: invalid size: " << size << endl;
		throw eGeneralFail;
	}
	
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
inline void ArrayT<TYPE>::Set(int length, TYPE* TYPEPtr)
{
	/* release memory if allocated */
	if (fDelete)
	{
		delete[] fArray;
		fDelete = 0;
	}
	fLength = length;
	fArray  = TYPEPtr;	
}

/* allocate an array of the specified size.*/
template <class TYPE>
void ArrayT<TYPE>::Allocate(int length)
{
	/* abort on negative lengths */
	if (length < 0) throw eGeneralFail;

	if (length != fLength || !fDelete)
	{
		/* old data is lost if it has new length */
		if (fDelete) delete[] fArray;
	
		fDelete = 1;
		fLength = length;
		fArray = New(fLength);
	}
}

/* deep or shallow ? */
template <class TYPE>
inline bool ArrayT<TYPE>::IsAllocated(void) const { return fDelete != 0; }

/* if allocated, free memory and set size to zero */
template <class TYPE>
void ArrayT<TYPE>::Free(void)
{
	/* set size */
	fLength = 0;

	/* free memory */
	if (fDelete)
	{
		delete[] fArray;
		fArray  = NULL;
		fDelete = 0;
	}
}

/* resize to new dimension, copying in at most what fits.
* extra space is initialized by specifying the fill. */
template <class TYPE>
void ArrayT<TYPE>::Resize(int new_length)
{
	/* quick return */
	if (new_length == fLength) return;

	/* abort on negative lengths */
	if (new_length < 0) throw eGeneralFail;
	
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
inline TYPE* ArrayT<TYPE>::Pointer(int offset) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (offset < 0 || offset > fLength) throw eOutOfRange;
#endif
	return fArray + offset;
}

/* element accessor */
template <class TYPE>
inline TYPE& ArrayT<TYPE>::operator[](int index) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (index < 0 || index >= fLength) throw eOutOfRange;
#endif
	return fArray[index];
}

template <class TYPE>
TYPE& ArrayT<TYPE>::First(void) const
{
#if __option(extended_errorcheck)
	if (fArray == NULL) throw eGeneralFail;
#endif	
	return *fArray;
}

template <class TYPE>
inline TYPE& ArrayT<TYPE>::Last(void) const
{
#if __option(extended_errorcheck)
	if (fArray == NULL) throw eGeneralFail;
#endif	
	return *(fArray + fLength - 1);
}

/* assignment operators */
template <class TYPE>
ArrayT<TYPE>& ArrayT<TYPE>::operator=(const TYPE& value)
{
	TYPE* p = fArray;
	for (int i = 0; i < fLength; i++)
		*p++ = value;
	return *this;
}

template <class TYPE>
ArrayT<TYPE>& ArrayT<TYPE>::operator=(const ArrayT<TYPE>& RHS)
{
	/* no copies to self */
	if (fArray != RHS.fArray)
	{
		/* allocate space if needed */
		if (fLength != RHS.fLength) Allocate(RHS.fLength); // old data is gone
				
		/* copy data */
		MemCopy(fArray, RHS.fArray, fLength);	
	}
	return *this;
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
	if (offset + length > fLength) throw eSizeMismatch;
	if (source_offset + length > source.fLength) throw eOutOfRange;
#endif

	/* copy */
	MemCopy(fArray + offset, source.fArray + source_offset, length);
}

/* collect prescribed values from source */
template <class TYPE>
void ArrayT<TYPE>::Collect(const ArrayT<int>& keys, const ArrayT<TYPE>& source)
{
#if __option (extended_errorcheck)
	if (keys.Length() != Length()) throw eSizeMismatch;
#endif

	int*  pkeys = keys.Pointer();
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
	{
		cout << "\n ArrayT<TYPE>::ReleasePointer: array is shallow" << endl;
		throw eGeneralFail;
	}
	else if (array == NULL)
	{
		cout << "\n ArrayT<TYPE>::ReleasePointer: cannot release pointer to NULL" << endl;
		throw eGeneralFail;
	}
	else
		fDelete = 0;
	
	*array = fArray;
}

/* shallow copy/conversion */
template <class TYPE>
inline void ArrayT<TYPE>::Alias(const ArrayT<TYPE>& RHS)
{
	Set(RHS.Length(), RHS.Pointer());
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

#endif /* _ARRAY_T_H_ */
