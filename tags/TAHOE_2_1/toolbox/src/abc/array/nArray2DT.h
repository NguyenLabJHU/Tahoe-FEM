/* $Id: nArray2DT.h,v 1.21 2003-11-21 22:41:30 paklein Exp $ */
/* created: paklein (07/09/1996) */
#ifndef _NARRAY2D_T_H_
#define _NARRAY2D_T_H_

/* base class */
#include "nArrayT.h"

namespace Tahoe {

/** nArrayT with subdimension. Row major storage */
template <class nTYPE>
class nArray2DT: public nArrayT<nTYPE>
{
public:

	/** \name constructors */
	/*@{*/
	nArray2DT(void);
	nArray2DT(int majordim, int minordim);
	nArray2DT(const nArray2DT& source);

	/** construct an alias */
	nArray2DT(int majordim, int minordim, const nTYPE* MATHTYPEPtr);
	/*@}*/

	/** destructors */
	~nArray2DT(void);

	/** \name convert to a shallow object */
	/*@{*/
	void Alias(int majordim, int minordim, const nTYPE* MATHTYPEPtr);
	void Alias(const nArray2DT& RHS);
	/*@}*/
	
	/** \deprecated replaced by nArray2DT::Alias on 09/04/2003 */	
	void Set(int majordim, int minordim, nTYPE* MATHTYPEPtr);

	/** set the array size to the given dimensions. No change occurs if the array
	 * is already the specified size. The previous contents of the array is
	 * not preserved. To preserve the array contents while changing the dimension
	 * use nArray2DT::Resize. */
	void Dimension(int majordim, int minordim);
	
	/** dimensions this array to the same length as the source, but no data is copied */
	void Dimension(const nArray2DT& source) { Dimension(source.MajorDim(), source.MinorDim()); };
	
	/** \deprecated replaced by nArray2DT::Dimension on 02/13/2002 */
	void Allocate(int majordim, int minordim) { Dimension(majordim, minordim); };

	/* free memory (if allocated) and set size to zero */
	void Free(void);

	/* resize to new major dimension, copying in at most what fits.
	 * extra space is initialized by specifying the fill. */
	void Resize(int new_majordim);
	void Resize(int new_majordim, const nTYPE& fill);

	/** \name dimensions */
	/*@{*/
	int MajorDim(void) const; /**< major dimension */
	int MinorDim(void) const; /**< minor dimension */
	/*@}*/

	/** \name accessors */
	/*@{*/
	/** reference to an element in the array */
	nTYPE& operator()(int majordim, int minordim);

	/** const reference to an element in the array */
	const nTYPE& operator()(int majordim, int minordim) const;

	/** pointer to the row in the array */
	nTYPE* operator()(int majordim);

	/** pointer to the row in the array */
	const nTYPE* operator()(int majordim) const;
	/*@}*/

	/* copy/assignment operators - by a scalar or element by element */
	nArray2DT<nTYPE>& operator=(const nArray2DT& RHS);
	nArray2DT<nTYPE>& operator=(const nTYPE& value);
	
	/** exchange data */
	void Swap(nArray2DT<nTYPE>& source);

	/** return transposed list (re-dimensions) */
	void Transpose(const nArray2DT<nTYPE>& source);
	  	
	/** copy the specified row into array.
	 * \param array destination for row data. Returns with
	 *        length equation to the minor dimension of *this. */
	void RowCopy(int row, nArrayT<nTYPE>& array) const;

	/** copy the specified row into array without range checking. */
	void RowCopy(int row, nTYPE* array) const;

	/** copy the specified column into array. */
	void ColumnCopy(int col, nArrayT<nTYPE>& array) const;

	/** copy the specified column into array without range checking. */
	void ColumnCopy(int col, nTYPE* array) const;

	/** copy the specified column out of the source array. 
	 * \param col destination in this of the column
	 * \param a source array
	 * \param col_a column to copy out of the source array */
	void ColumnCopy(int col, const nArray2DT<nTYPE>& a, int col_a);

	/** shallow copy of a row.
	 * \param row row number to alias
	 * \param array array to alias to the row data */
	void RowAlias(int row, nArrayT<nTYPE>& array) const;
	
	/* set values in batch */
	void SetRow(int row, const nTYPE& value);
	void SetRow(int row, const nTYPE* array);
	void SetRow(int row, const nArrayT<nTYPE>& array);
	
	void SetColumn(int col, const nTYPE& value);
	void SetColumn(int col, const nTYPE* array);
	void SetColumn(int col, const nArrayT<nTYPE>& array);

	/** return the sum of the values in the given row */
	nTYPE RowSum(int row) const;

	/** return the sum of the values in the given column */
	nTYPE ColumnSum(int col) const;

	/** scale all values in the given row by a factor */	
	void ScaleRow(int row, const nTYPE& scale);

	/** scale all values in the given column by a factor */	
	void ScaleColumn(int col, const nTYPE& scale);

	/* dot the specified row number with the array */
	nTYPE DotRow(int row, const nTYPE* array) const;
	nTYPE DotRow(int row, const nArrayT<nTYPE>& array) const;
	nTYPE DotColumn(int col, const nArrayT<nTYPE>& array) const;

	/* row to row operations */
	void CopyRowFromRow(int copyrow, int fromrow);
	void SwapRows(int row1, int row2);

	/* copying selected rows of array */
	void RowCollect(const ArrayT<int>& rows, const nArray2DT& array);
	void RowCollect(const int* rows, const nArray2DT& array);

	/* create sublist with transposed indexing */
	void SetLocal(const ArrayT<int>& rows, nArrayT<nTYPE>& sublist) const;

	/** set the rows of vals to the existing rows of data */
	void Assemble(const ArrayT<int>& rows, const nArray2DT& vals);

	/** \name add the rows of vals to the existing rows of data */
	/*@{*/
	void Accumulate(const ArrayT<int>& rows, const nArray2DT& vals);
	void Accumulate(int col, const ArrayT<int>& rows, const nArrayT<nTYPE>& vals);
	void Accumulate(int col, const ArrayT<int>& rows, const nTYPE* vals);
	/*@}*/

	/* commonly used single row operations;
	 *
	 *      this_i =  scale*RHS_i;
	 *		this_i += scale*RHS_i;
	 *
	 *  Note: this and RHS only need to have same minor dimension, as
	 *        long as row is in range for both.
	 */		
	void SetToRowScaled(int row, const nTYPE& scale, const nArray2DT& RHS);
	void AddRowScaled(int row, const nTYPE& scale, const nArray2DT& RHS);

	/** \name add scaled values to rows */
	/*@{*/
	/** check dimensions and add */
	void AddToRowScaled(int row, const nTYPE& scale, const nArrayT<nTYPE>& array);

	/** add assuming array is length of a row */
	void AddToRowScaled(int row, const nTYPE& scale, const nTYPE* array);
	/*@}*/

	/** add the same array to all columns */
	void AddToRowsScaled(const nTYPE& scale, const nArrayT<nTYPE>& array);

	/** \name add scaled values to columns */
	/*@{*/
	/** check dimensions and add */
	void AddToColumnScaled(int col, const nTYPE& scale, const nArrayT<nTYPE>& array);

	/** add assuming array is length of a column */
	void AddToColumnScaled(int col, const nTYPE& scale, const nTYPE* array);
	/*@}*/

	/** copy all values from source into this array. This array must be
	 * large enough to copy the contents of the source array at the
	 * specified location
	 * \param source source of data. All values are copied
	 * \param start row in this array where data will start being written */
	void BlockRowCopyAt(const nArray2DT& source, int start);

	/** copy all values from source into this array. This array must be
	 * large enough to copy the contents of the source array at the
	 * specified location
	 * \param source source of data. All values are copied
	 * \param start column in this array where data will start being written */
	void BlockColumnCopyAt(const nArray2DT& source, int start);

	/* read/write with line numbers (1,...) */
	void ReadNumbered(istream& in);
	void WriteNumbered(ostream& out) const;

	/* output by row */
	void PrintRow(int row, ostream& out) const;
	void PrintRow(int row, int rowlength, ostream& out) const;
	 		
protected:

	int fMajorDim;
	int fMinorDim;
};

/* I/O operators */	
template <class nTYPE>
ostream& operator<<(ostream& out, const nArray2DT<nTYPE>& array)
{
	const nTYPE* p = array.Pointer();
	int width = OutputWidth(out, p);
	for (int i = 0; i < array.MajorDim(); i++)
	{
		if (i > 0) out << '\n';
		for (int j = 0; j < array.MinorDim(); j++)
			out << setw(width) << *p++;			
	}
		return out;
}

/*************************************************************************
* Implementation
*************************************************************************/

/* constructor */
template <class nTYPE>
inline nArray2DT<nTYPE>::nArray2DT(void): fMajorDim(0), fMinorDim(0) { }

template <class nTYPE>
inline nArray2DT<nTYPE>::nArray2DT(int majordim, int minordim)
{
	Dimension(majordim, minordim);
}

template <class nTYPE>
inline nArray2DT<nTYPE>::nArray2DT(int majordim, int minordim, const nTYPE* MATHTYPEPtr):
	nArrayT<nTYPE>(majordim*minordim, MATHTYPEPtr),
	fMajorDim(majordim),
	fMinorDim(minordim)
{

}

template <class nTYPE>
inline nArray2DT<nTYPE>::nArray2DT(const nArray2DT& source)
{
	operator=(source);	
}

/* destructors */
template <class nTYPE>
inline nArray2DT<nTYPE>::~nArray2DT(void)
{
	fMajorDim = 0;
	fMinorDim = 0;
}

/* set fields - convert to shallow object */
template <class nTYPE>
inline void nArray2DT<nTYPE>::Alias(int majordim, int minordim, const nTYPE* MATHTYPEPtr)
{
	/* inherited */
	nArrayT<nTYPE>::Alias(majordim*minordim, MATHTYPEPtr);

	/* set dimensions */
	fMajorDim = majordim;
	fMinorDim = minordim;
}

template <class nTYPE>
inline void nArray2DT<nTYPE>::Alias(const nArray2DT& RHS)
{
	/* inherited */
	nArrayT<nTYPE>::Alias(RHS);

	/* set dimensions */
	fMajorDim = RHS.fMajorDim;
	fMinorDim = RHS.fMinorDim;
}

template <class nTYPE>
inline void nArray2DT<nTYPE>::Set(int majordim, int minordim, nTYPE* MATHTYPEPtr)
{
	Alias(majordim, minordim, MATHTYPEPtr);
}

/*
* Allocate an array of the specified size - works only with
* nArray2DT's created by default construction.
*/
template <class nTYPE>
inline void nArray2DT<nTYPE>::Dimension(int majordim, int minordim)
{
	/* zero dimensions */
	fMajorDim = fMinorDim = 0;

	/* inherited */
	nArrayT<nTYPE>::Dimension(majordim*minordim);

	/* set dimensions */
	fMajorDim = majordim;
	fMinorDim = minordim;
}

/* free memory (if allocated) and set size to zero */
template <class nTYPE>
inline void nArray2DT<nTYPE>::Free(void)
{
	/* inherited */
	nArrayT<nTYPE>::Free();
	
	/* reset dimensions */
	fMajorDim = 0;
	fMinorDim = 0;
}

/* resize to new major dimension, copying in at most what fits.
* extra space is initialized by specifying the fill. */
template <class nTYPE>
void nArray2DT<nTYPE>::Resize(int new_majordim)
{
	/* inherited */
	nArrayT<nTYPE>::Resize(new_majordim*fMinorDim, true);

	/* new dimension */
	fMajorDim = new_majordim;
}

template <class nTYPE>
void nArray2DT<nTYPE>::Resize(int new_majordim, const nTYPE& fill)
{
	/* inherited */
	nArrayT<nTYPE>::Resize(new_majordim*fMinorDim, fill);

	/* new dimension */
	fMajorDim = new_majordim;
}

/* dimensions */
template <class nTYPE>
inline int nArray2DT<nTYPE>::MajorDim(void) const { return fMajorDim; }
template <class nTYPE>
inline int nArray2DT<nTYPE>::MinorDim(void) const { return fMinorDim; }

/* accessors */
template <class nTYPE>
inline nTYPE& nArray2DT<nTYPE>::operator()(int majordim, int minordim)
{
#if __option (extended_errorcheck)
	/* range checking */
	if (majordim < 0 || majordim >= fMajorDim ||
		minordim < 0 || minordim >= fMinorDim) {
			ExceptionT::OutOfRange("nArray2DT<TYPE>::operator(int,int)",
				"majordim ! (0 <= %d <= %d), minordim ! (0 <= %d <= %d)", 
				majordim, fMajorDim-1,
				minordim, fMinorDim-1);
		}
#endif
	return fArray[majordim*fMinorDim + minordim];
}
template <class nTYPE>
inline const nTYPE& nArray2DT<nTYPE>::operator()(int majordim, int minordim) const
{
#if __option (extended_errorcheck)
	/* range checking */
	if (majordim < 0 || majordim >= fMajorDim ||
		minordim < 0 || minordim >= fMinorDim) {
			ExceptionT::OutOfRange("nArray2DT<TYPE>::operator(int,int)",
				"majordim ! (0 <= %d <= %d), minordim ! (0 <= %d <= %d)", 
				majordim, fMajorDim-1,
				minordim, fMinorDim-1);
		}
#endif
	return fArray[majordim*fMinorDim + minordim];
}

template <class nTYPE>
inline nTYPE* nArray2DT<nTYPE>::operator()(int majordim)
{
#if __option (extended_errorcheck)
	/* range checking */
	if (majordim < 0 || majordim >= fMajorDim) {
		ExceptionT::OutOfRange("nArray2DT<TYPE>::operator(int)",
			"majordim ! (0 <= %d <= %d)", majordim, fMajorDim-1);
	}
#endif
	return fArray + majordim*fMinorDim;
}
template <class nTYPE>
inline const nTYPE* nArray2DT<nTYPE>::operator()(int majordim) const
{
#if __option (extended_errorcheck)
	/* range checking */
	if (majordim < 0 || majordim >= fMajorDim) {
		ExceptionT::OutOfRange("nArray2DT<TYPE>::operator(int)",
			"majordim ! (0 <= %d <= %d)", majordim, fMajorDim-1);
	}
#endif
	return fArray + majordim*fMinorDim;
}

/*
* Copy/assignment operators - by a scalar or element by element
*
* Note: VC++ requires template argument in return type.
*/
template <class nTYPE>
inline nArray2DT<nTYPE>& nArray2DT<nTYPE>::operator=(const nArray2DT& RHS)
{
	/* inherited */
	nArrayT<nTYPE>::operator=(RHS);

	/* set dimensions */
	fMajorDim = RHS.fMajorDim;
	fMinorDim = RHS.fMinorDim;

	return *this;
}

template <class nTYPE>
inline nArray2DT<nTYPE>& nArray2DT<nTYPE>::operator=(const nTYPE& value)
{
	/* inherited */
	nArrayT<nTYPE>::operator=(value);

	return *this;
}

/* exchange data */
template <class nTYPE>
inline void nArray2DT<nTYPE>::Swap(nArray2DT<nTYPE>& source)
{
	/* inherited */
	nArrayT<nTYPE>::Swap(source);

	/* dimensions */
	int tmp = fMajorDim;
	fMajorDim = source.fMajorDim;
	source.fMajorDim = tmp;

	tmp = fMinorDim;
	fMinorDim = source.fMinorDim;
	source.fMinorDim = tmp;
}

/* return transposed list (re-dimensions) */
template <class nTYPE>
void nArray2DT<nTYPE>::Transpose(const nArray2DT<nTYPE>& source)
{
	/* called with data belonging to this*/
	if (source.Pointer() == Pointer()) 
	{
		nArray2DT<nTYPE> temp = source;
		Transpose(temp);
	}
	else
	{
		/* allocate memory */
		Dimension(source.fMinorDim, source.fMajorDim);
	
		const nTYPE* psrc = source.Pointer();	
		for (int i = 0; i < fMinorDim; i++)
		{
			nTYPE* pthis = Pointer() + i;
			for (int j = 0; j < fMajorDim; j++)
			{
				*pthis = *psrc++;
				pthis += fMinorDim;
			}	
		}
	}
}

/* create sublist with transposed indexing.*/
template <class nTYPE>
void nArray2DT<nTYPE>::SetLocal(const ArrayT<int>& rows,
	nArrayT<nTYPE>& sublist) const
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if (rows.Length()*fMinorDim != sublist.Length())
		ExceptionT::SizeMismatch();
#endif

//NOTE: could unroll loops for speed with
//      fMinorDim = 1,2,3
	int sublength = rows.Length();
	const int* prows = rows.Pointer();
	nTYPE* pSub = sublist.Pointer();	
	for (int i = 0; i < sublength; i++)
	{
		const nTYPE* parray = (*this)(prows[i]);
		nTYPE* psub = pSub;
		for (int j = 0; j < fMinorDim; j++)
		{
			*psub = *parray++;
			psub += sublength;
		}
		pSub++;	
	}
}

template <class nTYPE>
inline void nArray2DT<nTYPE>::RowCopy(int row, nTYPE* array) const
{
	/* safe wrapper for memcpy */
	MemCopy(array, (*this)(row), MinorDim());
}

template <class nTYPE>
inline void nArray2DT<nTYPE>::RowCopy(int row, nArrayT<nTYPE>& array) const
{
	/* redimension if needed */
	if (array.Length() != MinorDim()) array.Dimension(MinorDim());

	/* call wrapper */
	RowCopy(row, array.Pointer());
}

template <class nTYPE>
inline void nArray2DT<nTYPE>::ColumnCopy(int col, nTYPE* array) const
{
	/* don't do anything if empty */
	if (fMajorDim > 0) {
		nTYPE* pout = array;
		const nTYPE* pcol = &(*this)(0,col);
		for (int i = 0; i < fMajorDim; i++) {
			*pout++ = *pcol;
			pcol   += fMinorDim;
		}
	}
}

template <class nTYPE>
inline void nArray2DT<nTYPE>::ColumnCopy(int col, nArrayT<nTYPE>& array) const
{
#if __option (extended_errorcheck)
	if (array.Length() != fMajorDim) ExceptionT::SizeMismatch();
#endif

	/* call wrapper */
	ColumnCopy(col, array.Pointer());
}

/* copy the specified column out of the source array */
template <class nTYPE>
void nArray2DT<nTYPE>::ColumnCopy(int col, const nArray2DT<nTYPE>& a, int col_a)
{
#if __option (extended_errorcheck)
	if (a.MajorDim() != fMajorDim) ExceptionT::SizeMismatch("nArray2DT<nTYPE>::ColumnCopy");
#endif

	if (fMajorDim > 0) {
		int a_dim = a.MinorDim();
		nTYPE* pthis = Pointer(col);
		const nTYPE* pthat = a.Pointer(col_a);
		for (int i = 0; i < fMajorDim; i++) {
			*pthis = *pthat;
			pthis += fMinorDim;
			pthat += a_dim;
		}
	}
}

template <class nTYPE>
inline void nArray2DT<nTYPE>::RowAlias(int row,
	nArrayT<nTYPE>& array) const
{
	array.Alias(fMinorDim, (*this)(row));
}	

/* set values in batch */
template <class nTYPE>
inline void nArray2DT<nTYPE>::SetRow(int row, const nTYPE& value)
{
	nTYPE* prow = (*this)(row);
	for (int i = 0; i < fMinorDim; i++)
		*prow++ = value;
}			

template <class nTYPE>
inline void nArray2DT<nTYPE>::SetRow(int row, const nTYPE* array)
{
	/* copy */	
	MemCopy((*this)(row), array, fMinorDim);	
}

template <class nTYPE>
inline void nArray2DT<nTYPE>::SetRow(int row, const nArrayT<nTYPE>& array)
{
/* range checking */
#if __option (extended_errorcheck)
	if (array.Length() != fMinorDim) ExceptionT::SizeMismatch();
#endif
	
	/* copy */	
	MemCopy((*this)(row), array.Pointer(), fMinorDim);	
}

template <class nTYPE>
void nArray2DT<nTYPE>::SetColumn(int col, const nTYPE& value)
{
	if (fMajorDim > 0) {
		nTYPE* pcol = &(*this)(0,col);
		for (int i = 0; i < fMajorDim; i++) {
			*pcol = value;
			pcol += fMinorDim;
		}
	}
}

template <class nTYPE>
void nArray2DT<nTYPE>::SetColumn(int col, const nTYPE* array)
{
	if (fMajorDim > 0) {
		nTYPE* pcol = &(*this)(0,col);
		for (int i = 0; i < fMajorDim; i++) {
			*pcol = *array++;
			pcol += fMinorDim;
		}
	}
}

template <class nTYPE>
inline void nArray2DT<nTYPE>::SetColumn(int col, const nArrayT<nTYPE>& array)
{
/* range checking */
#if __option (extended_errorcheck)
	if (array.Length() != fMajorDim) ExceptionT::SizeMismatch();
#endif

	/* wrapper */
	SetColumn(col, array.Pointer());
}


/* row and column sums */
template <class nTYPE>
nTYPE nArray2DT<nTYPE>::RowSum(int row) const
{
	nTYPE sum = 0.0;
	const nTYPE* p = (*this)(row);
	for (int i = 0; i < fMinorDim; i++)
		sum += *p++;		
	return sum;
}

template <class nTYPE>
nTYPE nArray2DT<nTYPE>::ColumnSum(int col) const
{
	if (fMajorDim == 0)
		return 0.0;
	else 
	{
		nTYPE sum = 0.0;
		const nTYPE* p = &(*this)(0,col);
		for (int i = 0; i < fMajorDim; i++) {
			sum += *p;
			p   += fMinorDim;
		}
		return sum;
	}
}

/* row scaling */	
template <class nTYPE>
inline void nArray2DT<nTYPE>::ScaleRow(int row, const nTYPE& scale)
{
	nTYPE* p = (*this)(row);
	for (int i = 0; i < fMinorDim; i++)
		*p++ *= scale;
}

template <class nTYPE>
inline void nArray2DT<nTYPE>::ScaleColumn(int col, const nTYPE& scale)
{
	if (fMajorDim > 0) {
		nTYPE* p = &(*this)(0,col);
		for (int i = 0; i < fMajorDim; i++)	{
			*p *= scale;
			p += fMinorDim;
		}
	}
}

/* dot the specified row number with the array.  No
* check on the array dimensions. */
template <class nTYPE>
nTYPE nArray2DT<nTYPE>::DotRow(int row, const nTYPE* array) const
{
	const nTYPE *p = (*this)(row);
	register nTYPE sum = 0.0;
	register nTYPE temp;
	for (int i = 0; i < fMinorDim; i++)
	{
		temp  = *p++;
		temp *= *array++;
	
		sum  += temp;
	}
	return sum;
}

template <class nTYPE>
inline nTYPE nArray2DT<nTYPE>::DotRow(int row,
	const nArrayT<nTYPE>& array) const
{
#if __option (extended_errorcheck)
	/* check */
	if (array.Length() != fMinorDim) ExceptionT::SizeMismatch();
#endif

	return DotRow(row, array.Pointer());
}

template <class nTYPE>
nTYPE nArray2DT<nTYPE>::DotColumn(int col, const nArrayT<nTYPE>& array) const
{
#if __option (extended_errorcheck)
	/* check */
	if (array.Length() != fMajorDim) ExceptionT::SizeMismatch();
#endif
	
	if (fMajorDim == 0)
		return 0.0;
	else
	{
		const nTYPE *p = Pointer(col);
		const nTYPE *parray = array.Pointer();
		register nTYPE sum = 0.0;
		register nTYPE temp;
		for (int i = 0; i < fMajorDim; i++) {
			temp  = *p;
			temp *= *parray++;
			sum += temp;
			p += fMinorDim;
		}
		return sum;
	}
}

/* row to row operations */
template <class nTYPE>
void nArray2DT<nTYPE>::CopyRowFromRow(int copyrow, int fromrow)
{
	/* quick exit */
	if (copyrow == fromrow)
		return;
	else /* copy */
		MemCopy((*this)(copyrow), (*this)(fromrow), fMinorDim);	
}

template <class nTYPE>
void nArray2DT<nTYPE>::SwapRows(int row1, int row2)
{
	nTYPE* p1 = (*this)(row1);
	nTYPE* p2 = (*this)(row2);

	nTYPE temp;
	for (int i = 0; i < fMinorDim; i++)
	{
		temp = *p1;
		*p1++ = *p2;
		*p2++ = temp;
	}
}

/* deep and shallow row copies */
template <class nTYPE>
void nArray2DT<nTYPE>::RowCollect(const ArrayT<int>& rows,
	const nArray2DT<nTYPE>& array2D)
{
/* must have same minor dimension */
#if __option (extended_errorcheck)
	if (fMinorDim != array2D.fMinorDim || rows.Length() != fMajorDim)
		ExceptionT::SizeMismatch();
#endif

	/* shallow wrapper */
	const int* p = rows.Pointer();
	for (int i = 0; i < fMajorDim; i++)
		SetRow(i,array2D(*p++));
}

template <class nTYPE>
void nArray2DT<nTYPE>::RowCollect(const int* rows,
	const nArray2DT<nTYPE>& array2D)
{
/* must have same minor dimension */
#if __option (extended_errorcheck)
	if (fMinorDim != array2D.fMinorDim) ExceptionT::SizeMismatch();
#endif

	for (int i = 0; i < fMajorDim; i++)
		SetRow(i, array2D(*rows++));
}

/* assemble/add the rows of vals to the existing rows of data */
template <class nTYPE>
void nArray2DT<nTYPE>::Assemble(const ArrayT<int>& rows,
	const nArray2DT& vals)
{
	for (int i = 0; i < rows.Length(); i++)
	{	
		nTYPE* p = (*this)(rows[i]);
		const nTYPE* pvals = vals(i);
		for (int j = 0; j < vals.MinorDim(); j++)
			*p++ = *pvals++;
	}
}

template <class nTYPE>
void nArray2DT<nTYPE>::Accumulate(const ArrayT<int>& rows, const nArray2DT& vals)
{
	for (int i = 0; i < rows.Length(); i++)
	{	
		nTYPE* p = (*this)(rows[i]);
		const nTYPE* pvals = vals(i);
		for (int j = 0; j < vals.MinorDim(); j++)
			*p++ += *pvals++;
	}
}

template <class nTYPE>
inline void nArray2DT<nTYPE>::Accumulate(int col, const ArrayT<int>& rows, const nArrayT<nTYPE>& vals)
{
#if __option(extended_errorcheck)
	if (rows.Length() != vals.Length()) ExceptionT::GeneralFail();
#endif

	Accumulate(col, rows, vals.Pointer());
}

template <class nTYPE>
void nArray2DT<nTYPE>::Accumulate(int col, const ArrayT<int>& rows, const nTYPE* vals)
{
	for (int i = 0; i < rows.Length(); i++)
		(*this)(rows[i], col) += *vals++;
}

/*
* Commonly used single row operations;
*
*      this_i =  scale*RHS_i;
*		this_i += scale*RHS_i;
*
*  Note: this and RHS only need to have same minor dimension, as
*        long as row is in range for both.
*/		
template <class nTYPE>
void nArray2DT<nTYPE>::SetToRowScaled(int row, const nTYPE& scale,
	const nArray2DT& RHS)
{
	/* dimension checks */
#if __option(extended_errorcheck)	
	if (fMinorDim != RHS.fMinorDim) ExceptionT::SizeMismatch();
#endif

	nTYPE* pthis = (*this)(row);
	nTYPE* pRHS  = RHS(row);
	
	for (int i = 0; i < fMinorDim; i++)
	{
		*pthis    = scale;
		*pthis++ *= *pRHS++;
	}
}

template <class nTYPE>
void nArray2DT<nTYPE>::AddRowScaled(int row, const nTYPE& scale,
	const nArray2DT& RHS)
{
	/* dimension checks */
#if __option(extended_errorcheck)	
	if (fMinorDim != RHS.fMinorDim) ExceptionT::SizeMismatch();
#endif

	nTYPE* pthis = (*this)(row);
	nTYPE* pRHS  = RHS(row);
	
	nTYPE temp;
	
	for (int i = 0; i < fMinorDim; i++)
	{
		temp  = scale;
		temp *= *pRHS++;

		*pthis++ += temp;
	}
}

/* map addition on rows */
template <class nTYPE>
inline void nArray2DT<nTYPE>::AddToRowScaled(int row, const nTYPE& scale,
	const nArrayT<nTYPE>& array)
{
/* range checking */
#if __option (extended_errorcheck)
	if (array.Length() != fMinorDim) ExceptionT::SizeMismatch();
#endif

	AddToRowScaled(row, scale, array.Pointer());
}

template <class nTYPE>
void nArray2DT<nTYPE>::AddToRowScaled(int row, const nTYPE& scale, const nTYPE* array)
{
	nTYPE temp;
	const nTYPE* parray = array;	
	nTYPE* prow = (*this)(row);
	for (int i = 0; i < fMinorDim; i++)
	{
		temp = *parray++;
		temp *= scale;
		*prow++ += temp;	
	}
}	

template <class nTYPE>
void nArray2DT<nTYPE>::AddToRowsScaled(const nTYPE& scale,
	const nArrayT<nTYPE>& array)
{
/* range checking */
#if __option (extended_errorcheck)
	if (array.Length() != fMinorDim) ExceptionT::SizeMismatch();
#endif

	if (fMajorDim > fMinorDim)
	{
		nTYPE temp;
		const nTYPE* parray = array.Pointer();
		for (int j = 0; j < fMinorDim; j++)
		{
			nTYPE* pcol = Pointer(j);
			for (int i = 0; i < fMajorDim; i++)
			{
				temp  = *parray;
				temp *= scale;
			
				*pcol += temp;
				pcol += fMinorDim;
			}	
			parray++;
		}
	}
	else
	{
		nTYPE  temp;
		nTYPE* p = Pointer();
		for (int i = 0; i < fMajorDim; i++)
		{
			const nTYPE* parray = array.Pointer();
			for (int j = 0; j < fMinorDim; j++)
			{
				temp  = *parray++;
				temp *= scale;
				*p++ += temp;
			}
		}
	}
}

/* map addition on columns */
template <class nTYPE>
inline void nArray2DT<nTYPE>::AddToColumnScaled(int col, const nTYPE& scale,
	const nArrayT<nTYPE>& array)
{
/* range checking */
#if __option (extended_errorcheck)
	if (array.Length() != fMajorDim) ExceptionT::SizeMismatch();
#endif

	AddToColumnScaled(col, scale, array.Pointer());
}

template <class nTYPE>
void nArray2DT<nTYPE>::AddToColumnScaled(int col, const nTYPE& scale, const nTYPE* array)
{
	if (fMajorDim > 0) {
		nTYPE temp;
		const nTYPE* parray = array;	
		nTYPE* pcol = Pointer(col);
		for (int i = 0; i < fMajorDim; i++) {
			temp = *parray++;
			temp *= scale;
			*pcol += temp;
			pcol += fMinorDim;
		}
	}
}	

/* copy all rows/columns from source at start */
template <class nTYPE>
void nArray2DT<nTYPE>::BlockRowCopyAt(const nArray2DT& source, int start)
{
	/* quick exit */
	if (source.Length() == 0) return;

#if __option(extended_errorcheck)
	/* dimensions check */
	if (fMinorDim != source.fMinorDim) ExceptionT::SizeMismatch();
	if (start + source.fMajorDim > fMajorDim) ExceptionT::OutOfRange();
#endif

	/* copy */
	MemCopy((*this)(start), source.Pointer(), source.Length());
}

/* copy all rows/columns from source at start */
template <class nTYPE>
void nArray2DT<nTYPE>::BlockColumnCopyAt(const nArray2DT& source, int start)
{
	/* quick exit */
	if (source.Length() == 0) return;

#if __option(extended_errorcheck)
	/* dimension checks */
	if (fMajorDim != source.fMajorDim) ExceptionT::SizeMismatch();
	if (start + source.fMinorDim > fMinorDim) ExceptionT::OutOfRange();
#endif

	for (int i = 0; i < fMajorDim; i++)
	  MemCopy((*this)(i) + start, source(i), source.fMinorDim);
}

/* read/write with line numbers (1,...) */
template <class nTYPE>
void nArray2DT<nTYPE>::ReadNumbered(istream& in)
{
	for (int i = 0; i < fMajorDim; i++)
	{
		int row;
		in >> row;
		nTYPE* p = (*this)(--row);
		for (int j = 0; j < fMinorDim; j++)
			in >> *p++;
	}
}

template <class nTYPE>
void nArray2DT<nTYPE>::WriteNumbered(ostream& out) const
{
	const nTYPE *p = Pointer();
	int width = OutputWidth(out, p);
	for (int i = 0; i < fMajorDim; i++)
	{
		out << setw(kIntWidth) << i + 1;
		for (int j = 0; j < fMinorDim; j++)
			out << setw(width) << *p++;		
		out << '\n';
	}
}

/* print row */	
template <class nTYPE>
inline void nArray2DT<nTYPE>::PrintRow(int row, ostream& out) const
{
	PrintRow(row, fMinorDim, out);
}

template <class nTYPE>
void nArray2DT<nTYPE>::PrintRow(int row, int rowlength, ostream& out) const
{
#if __option(extended_errorcheck)
	/* no more than the whole row */
	if (rowlength > fMinorDim) ExceptionT::OutOfRange();
#endif

	const nTYPE* p = (*this)(row);
	int width = OutputWidth(out, p);
	for (int i = 0; i < rowlength; i++)
		out << setw(width) << *p++;
	out << '\n';
}

} // namespace Tahoe 
#endif /* _NARRAY2D_T_H_ */
