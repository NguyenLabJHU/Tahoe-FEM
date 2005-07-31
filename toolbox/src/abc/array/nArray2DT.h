/* $Id: nArray2DT.h,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (07/09/1996)                                          */
/* nArrayT with subdimension - row major storage                          */

#ifndef _NARRAY2D_T_H_
#define _NARRAY2D_T_H_

/* base class */
#include "nArrayT.h"

template <class nTYPE>
class nArray2DT: public nArrayT<nTYPE>
{
public:

	/* constructors */
	nArray2DT(void);
	nArray2DT(int majordim, int minordim);
	nArray2DT(int majordim, int minordim, nTYPE* MATHTYPEPtr);
	nArray2DT(const nArray2DT& source);

	/* destructors */
	~nArray2DT(void);

	/* set fields - convert to shallow object */
	void Set(int majordim, int minordim, nTYPE* MATHTYPEPtr);

	/* allocate an array of the specified size */
	void Allocate(int majordim, int minordim);

	/* free memory (if allocated) and set size to zero */
	void Free(void);

	/* resize to new major dimension, copying in at most what fits.
	 * extra space is initialized by specifying the fill. */
	void Resize(int new_majordim);
	void Resize(int new_majordim, const nTYPE& fill);

	/* dimensions */
	int MajorDim(void) const;
	int MinorDim(void) const;

	/* accessors */
	nTYPE& operator()(int majordim, int minordim) const;
	nTYPE* operator()(int majordim) const;

	/* copy/assignment operators - by a scalar or element by element */
	nArray2DT<nTYPE>& operator=(const nArray2DT& RHS);
	nArray2DT<nTYPE>& operator=(const nTYPE& value);
	void Alias(const nArray2DT& RHS);
	
	/* exchange data */
	void Swap(nArray2DT<nTYPE>& source);

	/* return transposed list (re-dimensions) - don't call with self! */
	void Transpose(const nArray2DT<nTYPE>& source);
	  	
	/* deep/shallow row copies */  	 							
	void RowCopy(int row, nArrayT<nTYPE>& array) const;
	void RowAlias(int row, nArrayT<nTYPE>& array) const;
	void ColumnCopy(int col, nArrayT<nTYPE>& array) const;
	
	/* set values in batch */
	void SetRow(int row, const nTYPE& value);
	void SetRow(int row, const nTYPE* array);
	void SetRow(int row, const nArrayT<nTYPE>& array);
	
	void SetColumn(int col, const nTYPE& value);
	void SetColumn(int col, const nTYPE* array);
	void SetColumn(int col, const nArrayT<nTYPE>& array);

	/* row and column sums */
	nTYPE RowSum(int row) const;
	nTYPE ColumnSum(int col) const;

	/* row scaling */	
	void ScaleRow(int row, const nTYPE& scale);

	/* dot the specified row number with the array */
	nTYPE DotRow(int row, nTYPE* array) const;
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

	/* assemble/add the rows of vals to the existing rows of data */
	void Assemble(const ArrayT<int>& rows, const nArray2DT& vals);
	void Accumulate(const ArrayT<int>& rows, const nArray2DT& vals);

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

	/* map addition on rows */
	void AddToRowScaled(int row, const nTYPE& scale, const nArrayT<nTYPE>& array);
	void AddToRowsScaled(const nTYPE& scale, const nArrayT<nTYPE>& array);

	/* copy all rows/columns from source at start */
	void BlockRowCopyAt(const nArray2DT& source, int start);
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
	nTYPE* p = array.Pointer();
	int width = OutputWidth(out, p);
	for (int i = 0; i < array.MajorDim(); i++)
	{
		if (i > 0) out << '\n';
		for (int j = 0; j < array.MinorDim(); j++)
			out << setw(width) << *p++;			
	}
		return out;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructor */
template <class nTYPE>
inline nArray2DT<nTYPE>::nArray2DT(void): fMajorDim(0), fMinorDim(0) { }

template <class nTYPE>
inline nArray2DT<nTYPE>::nArray2DT(int majordim, int minordim)
{
	Allocate(majordim, minordim);
}

template <class nTYPE>
inline nArray2DT<nTYPE>::nArray2DT(int majordim, int minordim,
	nTYPE* MATHTYPEPtr):
	nArrayT<nTYPE>(majordim*minordim, MATHTYPEPtr),
	fMajorDim(majordim), fMinorDim(minordim)
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
inline void nArray2DT<nTYPE>::Set(int majordim, int minordim,
	nTYPE* MATHTYPEPtr)
{
	/* inherited */
	nArrayT<nTYPE>::Set(majordim*minordim,MATHTYPEPtr);

	/* set dimensions */
	fMajorDim = majordim;
	fMinorDim = minordim;
}

/*
* Allocate an array of the specified size - works only with
* nArray2DT's created by default construction.
*/
template <class nTYPE>
inline void nArray2DT<nTYPE>::Allocate(int majordim, int minordim)
{
	/* inherited */
	nArrayT<nTYPE>::Allocate(majordim*minordim);

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
inline nTYPE& nArray2DT<nTYPE>::operator()(int majordim, int minordim) const
{
/* range checking */
#if __option (extended_errorcheck)
if (majordim < 0 || majordim >= fMajorDim ||
	minordim < 0 || minordim >= fMinorDim) throw(eOutOfRange);
#endif

	return fArray[majordim*fMinorDim + minordim];
}

template <class nTYPE>
inline nTYPE* nArray2DT<nTYPE>::operator()(int majordim) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (majordim < 0 || majordim >= fMajorDim) throw(eOutOfRange);
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

template <class nTYPE>
inline void nArray2DT<nTYPE>::Alias(const nArray2DT& RHS)
{
	/* inherited */
	nArrayT<nTYPE>::Alias(RHS);

	/* set dimensions */
	fMajorDim = RHS.fMajorDim;
	fMinorDim = RHS.fMinorDim;
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

/* return transposed list (re-dimensions) - don't call with self! */
template <class nTYPE>
void nArray2DT<nTYPE>::Transpose(const nArray2DT<nTYPE>& source)
{
	/* allocate memory */
	Allocate(source.fMinorDim, source.fMajorDim);
	
	nTYPE* psrc   = source.Pointer();	
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

/* create sublist with transposed indexing.*/
template <class nTYPE>
void nArray2DT<nTYPE>::SetLocal(const ArrayT<int>& rows,
	nArrayT<nTYPE>& sublist) const
{
//NOTE: could unroll loops for speed with
//      fMinorDim = 1,2,3

	int	 sublength = rows.Length();
	int*	 prows = rows.Pointer();
	nTYPE* pSub = sublist.Pointer();	
		
	for (int i = 0; i < sublength; i++)
	{
		nTYPE* parray = Pointer() + prows[i]*fMinorDim;
		nTYPE* psub	 = pSub;

		for (int j = 0; j < fMinorDim; j++)
		{
			*psub = *parray++;
			psub += sublength;
		}
		
		pSub++;	
	}
}

template <class nTYPE>
void nArray2DT<nTYPE>::RowCopy(int row, nArrayT<nTYPE>& array) const
{
	/* temp wrapper */
	nArrayT<nTYPE> temp(fMinorDim, (*this)(row) );

	/* deep copy */
	array = temp;
}

template <class nTYPE>
inline void nArray2DT<nTYPE>::RowAlias(int row,
	nArrayT<nTYPE>& array) const
{
	array.Set(fMinorDim, (*this)(row));
}	

template <class nTYPE>
void nArray2DT<nTYPE>::ColumnCopy(int col,
	nArrayT<nTYPE>& array) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (array.Length() != fMajorDim) throw eSizeMismatch;
#endif

	nTYPE* pout = array.Pointer();
	nTYPE* pcol = &(*this)(0,col);

	for (int i = 0; i < fMajorDim; i++)
	{
		*pout++ = *pcol;
		pcol   += fMinorDim;
	}
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
	if (array.Length() != fMinorDim) throw eSizeMismatch;
#endif
	
	/* copy */	
	MemCopy((*this)(row), array.Pointer(), fMinorDim);	
}

template <class nTYPE>
void nArray2DT<nTYPE>::SetColumn(int col, const nTYPE& value)
{
	nTYPE* pcol = &(*this)(0,col);

	for (int i = 0; i < fMajorDim; i++)
	{
		*pcol = value;
		pcol += fMinorDim;
	}
}

template <class nTYPE>
void nArray2DT<nTYPE>::SetColumn(int col,
	const nArrayT<nTYPE>& array)
{
	/* dimension check */
	if (array.Length() != fMajorDim)
		throw eSizeMismatch;
	else if (fMajorDim == 0)
		return;
	
	nTYPE* pcol = &(*this)(0,col);
	nTYPE* prhs = array.Pointer();
	for (int i = 0; i < fMajorDim; i++)
	{
		*pcol = *prhs++;
		pcol += fMinorDim;
	}
}

template <class nTYPE>
void nArray2DT<nTYPE>::SetColumn(int col, const nTYPE* array)
{
	nTYPE* pcol = &(*this)(0,col);
	for (int i = 0; i < fMajorDim; i++)
	{
		*pcol = *array++;
		pcol += fMinorDim;
	}
}

/* row and column sums */
template <class nTYPE>
nTYPE nArray2DT<nTYPE>::RowSum(int row) const
{
	nTYPE sum = 0.0;
	nTYPE* p  = (*this)(row);
	for (int i = 0; i < fMinorDim; i++)
		sum += *p++;		
	return sum;
}

template <class nTYPE>
nTYPE nArray2DT<nTYPE>::ColumnSum(int col) const
{
	nTYPE sum = 0.0;
	nTYPE* p  = &(*this)(0,col);
	for (int i = 0; i < fMajorDim; i++)
	{
		sum += *p;
		p   += fMinorDim;
	}
	return sum;
}

/* row scaling */	
template <class nTYPE>
inline void nArray2DT<nTYPE>::ScaleRow(int row, const nTYPE& scale)
{
	nTYPE* p = (*this)(row);
	for (int i = 0; i < fMinorDim; i++)
		*p++ *= scale;
}

/* dot the specified row number with the array.  No
* check on the array dimensions. */
template <class nTYPE>
nTYPE nArray2DT<nTYPE>::DotRow(int row, nTYPE* array) const
{
	nTYPE *p = (*this)(row);
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
	if (array.Length() != fMinorDim) throw eSizeMismatch;
#endif

	return DotRow(row, array.Pointer());
}

template <class nTYPE>
nTYPE nArray2DT<nTYPE>::DotColumn(int col, const nArrayT<nTYPE>& array) const
{
#if __option (extended_errorcheck)
	/* check */
	if (array.Length() != fMajorDim) throw eSizeMismatch;
#endif
	
	nTYPE *p = Pointer(col);
	nTYPE *parray = array.Pointer();
	register nTYPE sum = 0.0;
	register nTYPE temp;
	for (int i = 0; i < fMajorDim; i++)
	{
		temp  = *p;
		temp *= *parray++;
		sum += temp;
		p += fMinorDim;
	}
	return sum;
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
		throw eSizeMismatch;
#endif

	/* shallow wrapper */
	int* p = rows.Pointer();
	for (int i = 0; i < fMajorDim; i++)
		SetRow(i,array2D(*p++));
}

template <class nTYPE>
void nArray2DT<nTYPE>::RowCollect(const int* rows,
	const nArray2DT<nTYPE>& array2D)
{
/* must have same minor dimension */
#if __option (extended_errorcheck)
	if (fMinorDim != array2D.fMinorDim) throw eSizeMismatch;
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
		nTYPE*     p = (*this)(rows[i]);
		nTYPE* pvals = vals(i);
	
		for (int j = 0; j < vals.MinorDim(); j++)
			*p++ = *pvals++;
	}
}

template <class nTYPE>
void nArray2DT<nTYPE>::Accumulate(const ArrayT<int>& rows,
	const nArray2DT& vals)
{
	for (int i = 0; i < rows.Length(); i++)
	{	
		nTYPE*     p = (*this)(rows[i]);
		nTYPE* pvals = vals(i);
	
		for (int j = 0; j < vals.MinorDim(); j++)
			*p++ += *pvals++;
	}
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
	if (fMinorDim != RHS.fMinorDim) throw eSizeMismatch;
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
	if (fMinorDim != RHS.fMinorDim) throw eSizeMismatch;
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
void nArray2DT<nTYPE>::AddToRowScaled(int row, const nTYPE& scale,
	const nArrayT<nTYPE>& array)
{
/* range checking */
#if __option (extended_errorcheck)
	if (array.Length() != fMinorDim) throw eSizeMismatch;
#endif

	nTYPE  temp;
	nTYPE* parray = array.Pointer();	
	nTYPE* prow   = (*this)(row);
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
	if (array.Length() != fMinorDim) throw eSizeMismatch;
#endif

	if (fMajorDim > fMinorDim)
	{
		nTYPE  temp;
		nTYPE* parray = array.Pointer();
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
			nTYPE* parray = array.Pointer();
			for (int j = 0; j < fMinorDim; j++)
			{
				temp  = *parray++;
				temp *= scale;
				*p++ += temp;
			}
		}
	}
}

/* copy all rows/columns from source at start */
template <class nTYPE>
void nArray2DT<nTYPE>::BlockRowCopyAt(const nArray2DT& source, int start)
{
#if __option(extended_errorcheck)
	/* dimensions check */
if (fMinorDim != source.fMinorDim) throw eSizeMismatch;
	if (start + source.fMajorDim > fMajorDim) throw eOutOfRange;
#endif

	/* copy */
	if (source.Length() > 0)
		MemCopy((*this)(start), source.Pointer(), source.Length());
}

/* copy all rows/columns from source at start */
template <class nTYPE>
void nArray2DT<nTYPE>::BlockColumnCopyAt(const nArray2DT& source, int start)
{
#if __option(extended_errorcheck)
/* dimension checks */
if (fMajorDim != source.fMajorDim) throw eSizeMismatch;
	if (start + source.fMinorDim > fMinorDim) throw eOutOfRange;
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
	nTYPE *p = Pointer();
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
	if (rowlength > fMinorDim) throw eOutOfRange;
#endif

	nTYPE* p = (*this)(row);
	int width = OutputWidth(out, p);
	for (int i = 0; i < rowlength; i++)
		out << setw(width) << *p++;
	out << '\n';
}

#endif /* _NARRAY2D_T_H_ */
