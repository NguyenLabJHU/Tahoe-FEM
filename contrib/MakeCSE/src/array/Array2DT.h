/*
 * File: Array2DT.h - TEMPLATE
 *
 */

/*
 * created      : PAK (11/02/98)
 * last modified: PAK (11/02/98)
 */

#ifndef _ARRAY2D_T_H_
#define _ARRAY2D_T_H_

/* base class */
#include "ArrayT.h"

template <class TYPE>
class Array2DT: public ArrayT<TYPE>
{
  public:

	/* constructors */
	Array2DT(void);
	Array2DT(int majordim, int minordim);
	Array2DT(int majordim, int minordim, TYPE* TYPEPtr);
	Array2DT(const Array2DT& source); 

	/* set fields - convert to shallow object */
	void Set(int majordim, int minordim, TYPE* TYPEPtr);

	/* allocate an array of the specified size */
	void Allocate(int majordim, int minordim);

	/* free memory (if allocated) and set size to zero */
	void Free(void);

	/* resize to new major dimension, copying in at most what fits. 
	 * extra space is initialized by specifying the fill. */
	void Resize(int newmajordim, bool copy_in, bool byte_copy); 
	void Resize(int newmajordim, const TYPE& fill, bool copy_in, bool byte_copy); 

  	/* dimensions */
  	int MajorDim(void) const;
  	int MinorDim(void) const;

  	/* accessors */
  	TYPE& operator()(int majordim, int minordim) const;
  	TYPE* operator()(int majordim) const;

  	/* copy/assignment operators - by a scalar or element by element */
  	Array2DT<TYPE>& operator=(const Array2DT& RHS);
  	Array2DT<TYPE>& operator=(const TYPE& value);
  	void ShallowCopy(const Array2DT& RHS);

	/* by row - (length = -1) implies (length = fMinorDim) */
  	void WriteRow(int row, const TYPE* array, int length = -1);
  	void WriteRow(int row, const ArrayT<TYPE>& array);
  		//NOTE: WriteRow does not throw an exception for a size
  		//      mismatch. Copied length is the shorter of the
  		//      2 implied sizes.

  protected:
  
  	int fMajorDim;
  	int fMinorDim;
  	
};

/*************************************************************************
 *
 * Implementation
 *
 *************************************************************************/

/* constructor */
template <class TYPE> 
inline Array2DT<TYPE>::Array2DT(void): fMajorDim(0), fMinorDim(0) { }

template <class TYPE> 
inline Array2DT<TYPE>::Array2DT(int majordim, int minordim)
{
	Allocate(majordim, minordim);
}

template <class TYPE> 
inline Array2DT<TYPE>::Array2DT(int majordim, int minordim, TYPE* TYPEPtr): 
	ArrayT<TYPE>(majordim*minordim, TYPEPtr), 
	fMajorDim(majordim), 
	fMinorDim(minordim) 
{ 

}

template <class TYPE>
inline Array2DT<TYPE>::Array2DT(const Array2DT& source)
{
	operator=(source);	
}

/* set fields - convert to shallow object */
template <class TYPE>
inline void Array2DT<TYPE>::Set(int majordim, int minordim, 
	TYPE* TYPEPtr)
{
	/* inherited */
	ArrayT<TYPE>::Set(majordim*minordim,TYPEPtr);

	/* set dimensions */
	fMajorDim = majordim;
	fMinorDim = minordim;
}

/*
 * Allocate an array of the specified size - works only with
 * Array2DT's created by default construction.
 */
template <class TYPE> 
inline void Array2DT<TYPE>::Allocate(int majordim, int minordim)
{
	/* inherited */
	ArrayT<TYPE>::Allocate(majordim*minordim);

	/* set dimensions */
	fMajorDim = majordim;
	fMinorDim = minordim;
}

/* free memory (if allocated) and set size to zero */
template <class TYPE> 
inline void Array2DT<TYPE>::Free(void)
{
	/* inherited */
	ArrayT<TYPE>::Free();
	
	/* reset dimensions */
	fMajorDim = 0;
	fMinorDim = 0;
}

/* resize to new major dimension, copying in at most what fits. 
 * extra space is initialized by specifying the fill. */
template <class TYPE> 
void Array2DT<TYPE>::Resize(int newmajordim, bool copy_in, bool byte_copy)
{
	/* inherited */
	ArrayT<TYPE>::Resize(newmajordim*fMinorDim, copy_in, byte_copy);

	/* new dimension */
	fMajorDim = newmajordim;
}

template <class TYPE> 
void Array2DT<TYPE>::Resize(int newmajordim, const TYPE& fill, bool copy_in, 
	bool byte_copy)
{
	/* inherited */
	ArrayT<TYPE>::Resize(newmajordim*fMinorDim, fill, copy_in, byte_copy);

	/* new dimension */
	fMajorDim = newmajordim;
}

/* dimensions */
template <class TYPE>
inline int Array2DT<TYPE>::MajorDim(void) const { return fMajorDim; }

template <class TYPE>
inline int Array2DT<TYPE>::MinorDim(void) const { return fMinorDim; }

/* accessors */
template <class TYPE>
inline TYPE& Array2DT<TYPE>::operator()(int majordim, int minordim) const
{
/* range checking */
#if __option (extended_errorcheck)
if (majordim < 0 || majordim >= fMajorDim ||
	minordim < 0 || minordim >= fMinorDim) throw eOutOfRange;
#endif

	return fArray[majordim*fMinorDim + minordim];
}

template <class TYPE>
inline TYPE* Array2DT<TYPE>::operator()(int majordim) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (majordim < 0 || majordim >= fMajorDim) throw eOutOfRange;
#endif

	return fArray + majordim*fMinorDim ;
}

/* copy/assignment operators - by a scalar or element by element */
template <class TYPE>
inline Array2DT<TYPE>& Array2DT<TYPE>::operator=(const Array2DT& RHS)
{
	/* inherited */
	ArrayT<TYPE>::operator=(RHS);

	/* set dimensions */
	fMajorDim = RHS.fMajorDim;
	fMinorDim = RHS.fMinorDim;

	return *this;
}

template <class TYPE>
inline Array2DT<TYPE>& Array2DT<TYPE>::operator=(const TYPE& value)
{
	/* inherited */
	ArrayT<TYPE>::operator=(value);
	return *this;
}

template <class TYPE>
inline void Array2DT<TYPE>::ShallowCopy(const Array2DT& RHS)
{
	/* inherited */
	ArrayT<TYPE>::ShallowCopy(RHS);

	/* set dimensions */
	fMajorDim = RHS.fMajorDim;
	fMinorDim = RHS.fMinorDim;
}

template <class TYPE>  
inline void Array2DT<TYPE>::WriteRow(int row, const TYPE* array, int length)
{
	/* copy bytes */	
	length = (length > fMinorDim) ? fMinorDim : length;
	memcpy((*this)(row), array, sizeof(TYPE)*length);	
}

template <class TYPE>  
inline void Array2DT<TYPE>::WriteRow(int row, const ArrayT<TYPE>& array)
{	
	/* copy bytes */	
	int length = array.Length();
	length = (length > fMinorDim) ? fMinorDim : length;
	memcpy((*this)(row), array.Pointer(), sizeof(TYPE)*length);	
}

#endif /* _ARRAY2D_T_H_ */
