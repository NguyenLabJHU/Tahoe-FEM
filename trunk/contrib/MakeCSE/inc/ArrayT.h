/*
 * File: ArrayT.h - TEMPLATE
 *
 * Base class for handling memory allocation for arrays of TYPE
 *
 */

/*
 * created      : PAK (06/19/96)
 * last modified: PAK (04/14/99)
 */

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
	int IsAllocated(void) const;

	/* free memory (if allocated) and set size to zero */
	void Free(void);
	
	/* resize to new dimension, copying in at most what fits. 
	 * extra space is initialized by specifying the fill. */
	void Resize(int newlength, bool copy_in, bool byte_copy); 
	void Resize(int newlength, const TYPE& fill, bool copy_in, bool byte_copy);

	/* length of the array */
	int Length(void) const;

	/* returns a pointer specified element in the array */
	TYPE* Pointer(int offset = 0) const;

	/* element accessor */
	TYPE& operator[](int index) const;
	TYPE& First(void) const;
	TYPE& Last(void) const;
		
	/* assignment operators */
 	ArrayT<TYPE>& operator=(const TYPE& value);
 	ArrayT<TYPE>& operator=(const ArrayT<TYPE>& RHS);
 	
 	/* copy length elements of source beginning */
 	void CopyPart(int offset, const ArrayT<TYPE>& source, int source_offset, int length);
	void Collect(const ArrayT<int>& keys, const ArrayT<TYPE>& source);
 	
 	/* exchange data */
 	void Swap(ArrayT<TYPE>& source);

	/* array no longer responsible for freeing memory. error if
	 * the array is already shallow. DANGER. */
	void ReleasePointer(TYPE** array);
 
  	/* shallow copy/conversion */
	void ShallowCopy(const ArrayT& RHS);

	/* binary I/O - entire length of fArray */
	void ReadBinary(istream& in);
	void WriteBinary(ostream& out) const;

  protected:
  
  	/* return allocate memory and return a pointer (new C++ error handling) 
  	 * returns NULL on fail */
	TYPE* New(int size) const;
          	  	
  protected:
	 
	int	  fLength; /* size */
	TYPE* fArray;  /* data */
	
  private:
	
  	int   fDelete; /* allocation flag */  
};

/*************************************************************************
 * Implementation
 *************************************************************************/

/* constructor */
template <class TYPE> 
inline ArrayT<TYPE>::ArrayT(void):fDelete(0), fLength(0), fArray(NULL)
{

}

template <class TYPE> 
inline ArrayT<TYPE>::ArrayT(int length):fDelete(0), fLength(0), fArray(NULL)
{
	Allocate(length);	
}

template <class TYPE>
inline ArrayT<TYPE>::ArrayT(int length, TYPE* TYPEPtr): fDelete(0),
	fLength(length), fArray(TYPEPtr)
{

}

template <class TYPE>
inline ArrayT<TYPE>::ArrayT(const ArrayT& source):fDelete(0), fLength(0), 
	fArray(NULL)
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
inline TYPE* ArrayT<TYPE>::New(int size) const
{
	TYPE* p;

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
	return p;
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
inline int ArrayT<TYPE>::IsAllocated(void) const { return fDelete; }

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
void ArrayT<TYPE>::Resize(int newlength, bool copy_in, bool byte_copy)
{
	/* quick return */
	if (newlength == fLength) return;

	/* abort on negative lengths */
	if (newlength < 0) throw eGeneralFail;

	/* allocate new space */
	TYPE* newdata = New(newlength);
	
	/* copy in old data */
	if (copy_in)
	{
		int copysize = (newlength < fLength) ? newlength : fLength;
		
		if (byte_copy)
			memcpy(newdata, fArray, sizeof(TYPE)*copysize);
		else
			for (int i = 0; i < copysize; i++)
				newdata[i] = fArray[i];
	}
	
	/* free old memory */
	if (fDelete) delete[] fArray;
	
	/* set parameters */
	fLength = newlength;
	fArray  = newdata;
	fDelete = 1;
}

/* with specified fill */
template <class TYPE> 
void ArrayT<TYPE>::Resize(int newlength, const TYPE& fill, bool copy_in, bool byte_copy)
{
	int oldlength = fLength;
	
	/* resize the memory and copy in data */
	Resize(newlength, copy_in, byte_copy);

	/* initialize added space */
	TYPE* pthis = fArray + oldlength;		
	for (int i = oldlength; i < fLength; i++)
		*pthis++ = fill; 
}

/* returns the length of the array */
template <class TYPE> 
inline int ArrayT<TYPE>::Length(void) const { return fLength; }

/* returns a pointer to the 1st element in the array */
template <class TYPE>
inline TYPE* ArrayT<TYPE>::Pointer(int offset) const { return fArray + offset; }

/* element accessor */
template <class TYPE>
inline TYPE& ArrayT<TYPE>::operator[](int index) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (index < 0 || index >= fLength) throw(eOutOfRange);
#endif

	return fArray[index];
}

template <class TYPE>
TYPE& ArrayT<TYPE>::First(void) const
{
	return *fArray;
}

template <class TYPE>
inline TYPE& ArrayT<TYPE>::Last(void) const
{
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
				
		/* copy bytes */	
		memcpy(fArray, RHS.fArray, sizeof(TYPE)*fLength);
		// not always safe	
	}

	return *this;
}

/* copy length elements of source beginning at start */
template <class TYPE>
void ArrayT<TYPE>::CopyPart(int offset, const ArrayT<TYPE>& source, 
	int source_offset, int length)
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if (offset + length > fLength) throw eSizeMismatch;
	if (source_offset + length > source.fLength) throw eOutOfRange;
#endif

	/* copy bytes */
	memcpy(fArray + offset, source.fArray + source_offset, sizeof(TYPE)*length);
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
inline void ArrayT<TYPE>::ShallowCopy(const ArrayT& RHS)
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
