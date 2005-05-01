/* $Id: Array2DT.h,v 1.11 2005-05-01 19:27:45 paklein Exp $ */
/* created: paklein (11/02/1998) */
#ifndef _ARRAY2D_T_H_
#define _ARRAY2D_T_H_

/* base class */
#include "ArrayT.h"

namespace Tahoe {

/** 2D array */
template <class TYPE>
class Array2DT: public ArrayT<TYPE>
{
public:

	/** \name constructors */
	/*@{*/
	Array2DT(void);
	Array2DT(int majordim, int minordim);
	Array2DT(const Array2DT& source);

	/** construct an alias */
	Array2DT(int majordim, int minordim, const TYPE* TYPEPtr);
	/*@}*/

	/** destructor */
	~Array2DT(void);

	/** create alias */
	void Alias(int majordim, int minordim, const TYPE* TYPEPtr);

	/** \deprecated replaced with Array2DT::Alias */
	void Set(int majordim, int minordim, TYPE* TYPEPtr);

	/** set the array size to the given dimensions. No change occurs if the array
	 * is already the specified size. The previous contents of the array is
	 * not preserved. To preserve the array contents while changing the dimension
	 * use Array2DT::Resize. */
	void Dimension(int majordim, int minordim);

	/** dimensions this array to the same length as the source, but no data is copied */
	void Dimension(const Array2DT& source) { Dimension(source.MajorDim(), source.MinorDim()); };

	/** \deprecated replaced by Array2DT::Dimension on 02/13/2002 */
	void Allocate(int majordim, int minordim) { Dimension(majordim, minordim); };

	/** resize to new major dimension, copying in at most what fits.
	 * extra space is initialized by specifying the fill. */
	void Resize(int new_majordim);
	void Resize(int new_majordim, const TYPE& fill);

	/** exchange data */
	void Swap(Array2DT<TYPE>& source);

	/* free memory (if allocated) and set size to zero */
	void Free(void);

	/** \name dimensions */
	/*@{*/
	int MajorDim(void) const;
	int MinorDim(void) const;
	/*@}*/

	/** \name accessors */
	/*@{*/
	TYPE& operator()(int majordim, int minordim);
	const TYPE& operator()(int majordim, int minordim) const;

	TYPE* operator()(int majordim);
	const TYPE* operator()(int majordim) const;
	/*@}*/

	/** shallow copy of a row.
	 * \param row row number to alias
	 * \param array array to alias to the row data */
	void RowAlias(int row, ArrayT<TYPE>& array) const;

	/* copy/assignment operators - by a scalar or element by element */
	Array2DT<TYPE>& operator=(const Array2DT& RHS);
	Array2DT<TYPE>& operator=(const TYPE& value);
	void Alias(const Array2DT& RHS);

	/** \name set values in batch */
	/*@{*/
	void SetRow(int row, const TYPE& value);
	void SetRow(int row, const TYPE* array);
	void SetRow(int row, const ArrayT<TYPE>& array);
	/*@}*/

protected:

	int fMajorDim;
	int fMinorDim;  	
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructor */
template <class TYPE>
inline Array2DT<TYPE>::Array2DT(void):
	fMajorDim(0), 
	fMinorDim(0)
{ 

}

template <class TYPE>
inline Array2DT<TYPE>::Array2DT(int majordim, int minordim)
{
	Dimension(majordim, minordim);
}

template <class TYPE>
inline Array2DT<TYPE>::Array2DT(int majordim, int minordim, const TYPE* TYPEPtr):
	ArrayT<TYPE>(majordim*minordim, TYPEPtr),
	fMajorDim(majordim),
	fMinorDim(minordim)
{

}

/* destructor */
template <class TYPE>
inline Array2DT<TYPE>::~Array2DT(void)
{
	fMajorDim = 0;
	fMinorDim = 0;
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

template <class TYPE>
inline void Array2DT<TYPE>::Alias(int majordim, int minordim,
	const TYPE* TYPEPtr)
{
	/* inherited */
	ArrayT<TYPE>::Alias(majordim*minordim,TYPEPtr);

	/* set dimensions */
	fMajorDim = majordim;
	fMinorDim = minordim;
}

/*
* Allocate an array of the specified size - works only with
* Array2DT's created by default construction.
*/
template <class TYPE>
inline void Array2DT<TYPE>::Dimension(int majordim, int minordim)
{
	/* zero dimensions */
	fMajorDim = fMinorDim = 0;

	/* (try) inherited */
	ArrayT<TYPE>::Dimension(majordim*minordim);

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
inline void Array2DT<TYPE>::Resize(int new_majordim)
{
	/* inherited */
	ArrayT<TYPE>::Resize(new_majordim*fMinorDim);

	/* new dimension */
	fMajorDim = new_majordim;
}

template <class TYPE>
inline void Array2DT<TYPE>::Resize(int new_majordim, const TYPE& fill)
{
	/* inherited */
	ArrayT<TYPE>::Resize(new_majordim*fMinorDim, fill);

	/* new dimension */
	fMajorDim = new_majordim;
}

/* exchange data */
template <class TYPE>
inline void Array2DT<TYPE>::Swap(Array2DT<TYPE>& source)
{
	/* inherited */
	ArrayT<TYPE>::Swap(source);

	/* dimensions */
	int tmp = fMajorDim;
	fMajorDim = source.fMajorDim;
	source.fMajorDim = tmp;

	tmp = fMinorDim;
	fMinorDim = source.fMinorDim;
	source.fMinorDim = tmp;
}

/* dimensions */
template <class TYPE>
inline int Array2DT<TYPE>::MajorDim(void) const { return fMajorDim; }

template <class TYPE>
inline int Array2DT<TYPE>::MinorDim(void) const { return fMinorDim; }

/* accessors */
template <class TYPE>
inline TYPE& Array2DT<TYPE>::operator()(int majordim, int minordim)
{
/* range checking */
#if __option (extended_errorcheck)
if (majordim < 0 || majordim >= fMajorDim ||
	minordim < 0 || minordim >= fMinorDim) throw ExceptionT::kOutOfRange;
#endif

	return this->fArray[majordim*fMinorDim + minordim];
}

template <class TYPE>
inline const TYPE& Array2DT<TYPE>::operator()(int majordim, int minordim) const
{
/* range checking */
#if __option (extended_errorcheck)
if (majordim < 0 || majordim >= fMajorDim ||
	minordim < 0 || minordim >= fMinorDim) throw ExceptionT::kOutOfRange;
#endif

	return this->fArray[majordim*fMinorDim + minordim];
}

template <class TYPE>
inline TYPE* Array2DT<TYPE>::operator()(int majordim)
{
/* range checking */
#if __option (extended_errorcheck)
	if (majordim < 0 || majordim >= fMajorDim) throw ExceptionT::kOutOfRange;
#endif

	return this->fArray + majordim*fMinorDim ;
}

template <class TYPE>
inline const TYPE* Array2DT<TYPE>::operator()(int majordim) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (majordim < 0 || majordim >= fMajorDim) throw ExceptionT::kOutOfRange;
#endif

	return this->fArray + majordim*fMinorDim ;
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
inline void Array2DT<TYPE>::Alias(const Array2DT& RHS)
{
	/* inherited */
	ArrayT<TYPE>::Alias(RHS);

	/* set dimensions */
	fMajorDim = RHS.fMajorDim;
	fMinorDim = RHS.fMinorDim;
}

/* set values in batch */
template <class TYPE>
inline void Array2DT<TYPE>::SetRow(int row, const TYPE& value)
{
	TYPE* prow = (*this)(row);
	for (int i = 0; i < fMinorDim; i++)
		*prow++ = value;
}			

template <class TYPE>
inline void Array2DT<TYPE>::SetRow(int row, const TYPE* array)
{
	/* copy */	
	MemCopy((*this)(row), array, fMinorDim);	
}

template <class TYPE>
inline void Array2DT<TYPE>::SetRow(int row, const ArrayT<TYPE>& array)
{
/* range checking */
#if __option (extended_errorcheck)
	if (array.Length() != fMinorDim) throw ExceptionT::kSizeMismatch;
#endif
	
	/* copy */	
	MemCopy((*this)(row), array.Pointer(), fMinorDim);	
}

template <class TYPE>
inline void Array2DT<TYPE>::RowAlias(int row, ArrayT<TYPE>& array) const
{
	array.Alias(fMinorDim, (*this)(row));
}	

}//namespace Tahoe

#endif /* _ARRAY2D_T_H_ */
