/* $Id: nMatrixT.h,v 1.32 2005-07-29 03:09:33 paklein Exp $ */
/* created: paklein (05/24/1996) */
#ifndef _NMATRIX_T_H_
#define _NMATRIX_T_H_

/* base class */
#include "nArrayT.h"

namespace Tahoe {

/** two-dimensional matrix mathematics template object */
template <class nTYPE>
class nMatrixT: public nArrayT<nTYPE>
{
public:

	/** \name control flags */
	/*@{*/
	enum SymmetryFlagT {kWhole = 0, kUpperOnly = 1};
	enum AssemblyModeT {kOverwrite = 0, kAccumulate = 1};
	/*@}*/

	/** \name constructors */
	/*@{*/
	nMatrixT(void);
	nMatrixT(int numrows, int numcols);
	explicit nMatrixT(int squaredim);
	nMatrixT(const nMatrixT& source);	

	/** construct alias */
	nMatrixT(int numrows, int numcols, const nTYPE* p);	
	/*@}*/

	/** destructor*/
	~nMatrixT(void);

	/** \name dimensioning methods 
	 * Also see copy/assignment operators */
	/*@{*/
	/** dimension matrix */
	void Dimension(int numrows, int numcols);

	/** dimension square matrix */
	void Dimension(int squaredim);

	/** dimensions this matrix to the same dimensions as the source, but no data is copied */
	void Dimension(const nMatrixT& source) { Dimension(source.Rows(), source.Cols()); };

	/** \name convert to a shallow object */
	/*@{*/
	void Alias(int numrows, int numcols, const nTYPE* p);
	void Alias(int dim, const nTYPE* p) { Alias(dim, dim, p); };
	void Alias(const nMatrixT& RHS);
	/*@}*/

	/** \deprecated replaced by nMatrixT::Alias on 09/04/2003 */
	void Set(int numrows, int numcols, nTYPE* p);

	/** free memory (if allocated) and set size to zero */
	void Free(void);

	/** number of rows */
	int Rows(void) const;

	/** number of columns */
	int Cols(void) const;
	/*@}*/

	/** \name deprecated dimensioning methods */
	/*@{*/
	/** \deprecated replaced by nMatrixT::Dimension on 02/13/2002 */
	void Allocate(int numrows, int numcols) { Dimension(numrows, numcols); } ;
	/** \deprecated replaced by nMatrixT::Dimension on 02/13/2002 */
	void Allocate(int squaredim) { Dimension(squaredim); } ;
	/*@}*/

	/** \name accessors */
	/*@{*/
	/** return a single element of the matrix */
	nTYPE& operator()(int nrow, int ncol);

	/** return a single element of the matrix */
	const nTYPE& operator()(int nrow, int ncol) const;

	/** return a pointer a column in the matrix */
	nTYPE* operator()(int ncol);

	/** return a pointer a column in the matrix */
	const nTYPE* operator()(int ncol) const;
	/*@}*/

	/** \name accessing blocks
	 * row and col in the upper left */
	/*@{*/
	void AddBlock(int row, int col, const nMatrixT<nTYPE>& block);
	void SetBlock(int row, int col, const nMatrixT<nTYPE>& block);

	void CopyBlock(int row, int col, nMatrixT<nTYPE>& block) const;
	void CopyBlock(const ArrayT<int>& rc, nMatrixT<nTYPE>& block) const; //block must be square
	void CopyBlock(const ArrayT<int>& r, const ArrayT<int>& c, nMatrixT<nTYPE>& block) const;
	/*@}*/

	/** \name copy/assignment operators */
	/*@{*/
	nMatrixT<nTYPE>& operator=(const nMatrixT& RHS);
	nMatrixT<nTYPE>& operator=(const nTYPE& value);

	/** exchange data */
	void Swap(nMatrixT<nTYPE>& source);
	/*@}*/
	
	/** \name accessing rows and columns */
	/*@{*/
	void CopyRow(int rownum, ArrayT<nTYPE>& row) const;
	void CopyFromRow(int rownum, int start_col, ArrayT<nTYPE>& row) const;
	void CopyRows(const ArrayT<int>& rows,
	              const nMatrixT<nTYPE>& source);
	void CopyColumn(int colnum, ArrayT<nTYPE>& col) const;
	void CopyColumns(const ArrayT<int>& cols,
	                 const nMatrixT<nTYPE>& source);
	void ColumnAlias(int colnum, ArrayT<nTYPE>& col) const;
	/*@}*/

	/** create a symmetric matrix. assumes the data is stored
	 * in the upper triangle of the matrix.  Setting IsUpper = 0,
	 * copies the data from the lower triangle */
	void CopySymmetric(int IsUpper = 1);

	/** set this to the matrix transpose.
	 * \param matrix source matrix to transpose
	 * \return reference to *this */
	nMatrixT<nTYPE>& Transpose(const nMatrixT<nTYPE>& matrix, int fillmode = kOverwrite);

	/** set this its matrix transpose.
	 * \return reference to *this */
	nMatrixT<nTYPE>& Transpose(int fillmode = kOverwrite);
	
	/** \name matrix-matrix multiplication
	 * operates on this using \e a and \e b. Operations allowed on entire matrices 
	 * only, all matrix dimensions must be consistent */
	/*@{*/
	void MultAB(const nMatrixT& a, const nMatrixT& b, int upper = 0);
	void MultATB(const nMatrixT& a, const nMatrixT& b, int upper = 0);
	void MultABT(const nMatrixT& a, const nMatrixT& b, int upper = 0);
	void MultATBT(const nMatrixT& a, const nMatrixT& b);
	/*@}*/

	/** \name matrix-matrix-matrix operations */
	/*@{*/
	void MultABC(const nMatrixT& a, const nMatrixT& b, const nMatrixT& c,
		int range = kWhole, int fillmode = kOverwrite);
	void MultABCT(const nMatrixT& a, const nMatrixT& b, const nMatrixT& c,
		int range = kWhole, int fillmode = kOverwrite);
	void MultATBC(const nMatrixT& a, const nMatrixT& b, const nMatrixT& c,
		int range = kWhole, int fillmode = kOverwrite);
	/*@}*/

	/** \name symmetric matrix-matrix-matrix operations
	 * Useful for tensor basis transformations */
	/*@{*/	 
	void MultQBQT(const nMatrixT& q, const nMatrixT& b,
		int range = kWhole, int fillmode = kOverwrite);
	void MultQTBQ(const nMatrixT& q, const nMatrixT& b,
		int range = kWhole, int fillmode = kOverwrite);
	/*@}*/	 
	
	/** \name matrix-vector multiplication
	 * Methods taking pointer arguments assume arrays have correct dimensions
	 * \param x vector contracted with this matrix
	 * \param b returns with results of product */
	/*@{*/
	void Multx(const nArrayT<nTYPE>& x, nArrayT<nTYPE>& b, const nTYPE& scale = nTYPE(1.0), int fillmode = kOverwrite) const;
	void Multx(const nTYPE* x, nTYPE* b, const nTYPE& scale = nTYPE(1.0), int fillmode = kOverwrite) const;
	void MultTx(const nArrayT<nTYPE>& x, nArrayT<nTYPE>& b, const nTYPE& scale = nTYPE(1.0), int fillmode = kOverwrite) const;
	void MultTx(const nTYPE* x, nTYPE* b, const nTYPE& scale = nTYPE(1.0), int fillmode = kOverwrite) const;
	/*@}*/

	/** vector-matrix-vector product */
	nTYPE MultmBn(const nArrayT<nTYPE>& m, const nArrayT<nTYPE>& n) const;
	   		
	/** \name dyadic product
	 * Set this to the the outer product of the two vectors, or
	 * in dyadic notation \f$ \mathbf{v}_1 \otimes \mathbf{v}_2 \f$. */
	/*@{*/
	void Outer(const nArrayT<nTYPE>& v1, const nArrayT<nTYPE>& v2, 
		const nTYPE& scale = nTYPE(1.0), int fillmode = kOverwrite);
	void Outer(const nTYPE* v1, const nTYPE* v2, 
		const nTYPE& scale = nTYPE(1.0), int fillmode = kOverwrite);
	/*@}*/

	/** \name identity operations
	 * For use with square matrices \e only */
	/*@{*/
	void PlusIdentity(const nTYPE& value = nTYPE(1.0));
	nMatrixT<nTYPE>& Identity(const nTYPE& value = nTYPE(1.0));
	/*@}*/

	/** \name writing to rows */
	/*@{*/	
	void SetRow(int row, const nArrayT<nTYPE>& vec);
	void SetRow(int row, const nTYPE* vec);
	void SetRow(int row, const nTYPE& value);
	/*@}*/	

	/** \name writing to cols */
	/*@{*/	
	void SetCol(int col, const nArrayT<nTYPE>& vec);
	void SetCol(int col, const nTYPE* vec);
	void SetCol(int col, const nTYPE& value);
	/*@}*/	

	/** \name row/column products
	 * dot the specified row/column number with the array */
	/*@{*/
	nTYPE DotRow(int rownum, const nArrayT<nTYPE>& array) const;
	nTYPE DotRow(int rownum, const nTYPE* array) const;
	nTYPE DotCol(int colnum, const nArrayT<nTYPE>& array) const;
	nTYPE DotCol(int colnum, const nTYPE* array) const;
	/*@}*/

protected:
	
	/** \name dimensions */
	/*@{*/
	int	fRows;
	int	fCols;
	/*@}*/
};

/* I/O operators */
template <class nTYPE>
istream& operator>>(istream& in, nMatrixT<nTYPE>& matrix)
{
	for (int j = 0; j < matrix.Rows(); j++)
		for (int i = 0; i < matrix.Cols(); i++)
				in >> matrix(j,i);
	return in;
}

template <class nTYPE>
ostream& operator<<(ostream& out, const nMatrixT<nTYPE>& matrix)
{
	int width = OutputWidth(out, matrix.Pointer());	
	for (int j = 0; j < matrix.Rows(); j++)	
	{
		if (j > 0) out << '\n';
		for (int i = 0; i < matrix.Cols(); i++)
			out << setw(width) << matrix(j,i);		
	}
	return out;
}

/*************************************************************************
* Implementation
*************************************************************************/

/* constructor */
template <class nTYPE>
inline nMatrixT<nTYPE>::nMatrixT(void): fRows(0), fCols(0) { }

template <class nTYPE>
inline nMatrixT<nTYPE>::nMatrixT(int numrows, int numcols)
{
	Dimension(numrows,numcols);
}

template <class nTYPE>
inline nMatrixT<nTYPE>::nMatrixT(int squaredim)
{
	Dimension(squaredim);
}

template <class nTYPE>
inline nMatrixT<nTYPE>::nMatrixT(int numrows, int numcols, const nTYPE* p)
{
	Alias(numrows, numcols, p);
}

template <class nTYPE>
inline nMatrixT<nTYPE>::nMatrixT(const nMatrixT& source)
{
	operator=(source);
}

/* destructor*/
template <class nTYPE>
inline nMatrixT<nTYPE>::~nMatrixT(void)
{
	fRows = 0;
	fCols = 0;
}

/* post construction dimensioning */
template <class nTYPE>
inline void nMatrixT<nTYPE>::Dimension(int numrows, int numcols)
{
	/* zero dimensions */
	fRows = fCols = 0;

	/* inherited */
	nArrayT<nTYPE>::Dimension(numrows*numcols);

	/* set dimensions */
	fRows = numrows;
	fCols = numcols;
}

template <class nTYPE>
inline void nMatrixT<nTYPE>::Dimension(int squaredim)
{
	Dimension(squaredim,squaredim);
}

template <class nTYPE>
inline void nMatrixT<nTYPE>::Alias(int numrows, int numcols, const nTYPE* p)
{
	/* inherited */
	nArrayT<nTYPE>::Alias(numrows*numcols, p);

	/* set dimensions */
	fRows = numrows;
	fCols = numcols;
}

template <class nTYPE>
inline void nMatrixT<nTYPE>::Alias(const nMatrixT& RHS)
{
	Alias(RHS.Rows(), RHS.Cols(), RHS.Pointer() );	
}

template <class nTYPE>
inline void nMatrixT<nTYPE>::Set(int numrows, int numcols, nTYPE* p)
{
	Alias(numrows, numcols, p);
}

/* free memory (if allocated) and set size to zero */
template <class nTYPE>
inline void nMatrixT<nTYPE>::Free(void)
{
	/* inherited */
	nArrayT<nTYPE>::Free();
	
	/* set size parameters */
	fRows = 0;
	fCols = 0;
}

/* element accessor */
template <class nTYPE>
inline nTYPE& nMatrixT<nTYPE>::operator()(int nrow, int ncol)
{
#if __option (extended_errorcheck)
	/* range checking */
	if (nrow < 0 || nrow >= fRows || ncol < 0 || ncol >= fCols)
		ExceptionT::OutOfRange("nMatrixT<TYPE>::operator(int,int)",
			"row ! (0 <= %d <= %d), col ! (0 <= %d <= %d)",
			nrow, fRows-1, ncol, fCols-1);
#endif
	
	return(this->fArray[ncol*fRows + nrow]);
}

template <class nTYPE>
inline const nTYPE& nMatrixT<nTYPE>::operator()(int nrow, int ncol) const
{
#if __option (extended_errorcheck)
	/* range checking */
	if (nrow < 0 || nrow >= fRows || ncol < 0 || ncol >= fCols)
		ExceptionT::OutOfRange("nMatrixT<TYPE>::operator(int,int)",
			"row ! (0 <= %d <= %d), col ! (0 <= %d <= %d)",
			nrow, fRows-1, ncol, fCols-1);
#endif
	
	return(this->fArray[ncol*fRows + nrow]);
}

/* returns a pointer to the top of the specified column */
template <class nTYPE>
inline nTYPE* nMatrixT<nTYPE>::operator()(int ncol)
{
#if __option (extended_errorcheck)
/* range checking */
	if (ncol < 0 || ncol >= fCols)
		ExceptionT::OutOfRange("nMatrixT<TYPE>::operator(int)",
			"col ! (0 <= %d <= %d)", ncol, fCols-1);
#endif
	
	return(this->fArray + ncol*fRows);
}

template <class nTYPE>
inline const nTYPE* nMatrixT<nTYPE>::operator()(int ncol) const
{
#if __option (extended_errorcheck)
/* range checking */
	if (ncol < 0 || ncol >= fCols)
		ExceptionT::OutOfRange("nMatrixT<TYPE>::operator(int)",
			"col ! (0 <= %d <= %d)", ncol, fCols-1);
#endif
	
	return(this->fArray + ncol*fRows);
}

/* assemble beginning with row and col in the upper left. */
template <class nTYPE>
void nMatrixT<nTYPE>::AddBlock(int row, int col,
const nMatrixT<nTYPE>& block)
{
/* range checking */
#if __option (extended_errorcheck)
	if (row + block.Rows() > fRows ||
	    col + block.Cols() > fCols) ExceptionT::SizeMismatch("nMatrixT<nTYPE>::AddBlock");
#endif

	double* pstart = &(*this)(row,col);
	const double* pblock = block.Pointer();
	for (int i = 0; i < block.Cols(); i++)
	{
		double* pcol = pstart;
		for (int j = 0; j < block.Rows(); j++)
			*pcol++ += *pblock++;
			
		pstart += fRows;
	}
}

template <class nTYPE>
void nMatrixT<nTYPE>::SetBlock(int row, int col,
const nMatrixT<nTYPE>& block)
{
/* range checking */
#if __option (extended_errorcheck)
	if (row + block.Rows() > fRows ||
	    col + block.Cols() > fCols) ExceptionT::SizeMismatch("nMatrixT<nTYPE>::SetBlock");
#endif

	double* pstart = &(*this)(row,col);
	const double* pblock = block.Pointer();
	for (int i = 0; i < block.Cols(); i++)
	{
		double* pcol = pstart;
		for (int j = 0; j < block.Rows(); j++)
			*pcol++ = *pblock++;
			
		pstart += fRows;
	}
}

template <class nTYPE>
void nMatrixT<nTYPE>::CopyBlock(int row, int col,
	nMatrixT<nTYPE>& block) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (row + block.Rows() > fRows ||
	    col + block.Cols() > fCols) ExceptionT::SizeMismatch("nMatrixT<nTYPE>::CopyBlock");
#endif

	const double* pstart = &(*this)(row,col);
	double* pblock = block.Pointer();
	for (int i = 0; i < block.Cols(); i++)
	{
		const double* pcol = pstart;
		for (int j = 0; j < block.Rows(); j++)
			*pblock++ = *pcol++;
		pstart += fRows;
	}
}

template <class nTYPE>
void nMatrixT<nTYPE>::CopyBlock(const ArrayT<int>& rc,
	nMatrixT<nTYPE>& block) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (block.Rows() != block.Cols() ||
	    block.Cols() != rc.Length()) ExceptionT::SizeMismatch("nMatrixT<nTYPE>::CopyBlock");
#endif
	
	nTYPE* pblock = block.Pointer();
	for (int j = 0; j < rc.Length(); j++)
	{
		const nTYPE* pcol = (*this)(rc[j]);
		const int* prc = rc.Pointer();
		for (int i = 0; i < rc.Length(); i++)
			*pblock++ = pcol[*prc++];
	}
}

template <class nTYPE>
void nMatrixT<nTYPE>::CopyBlock(const ArrayT<int>& r,
	const ArrayT<int>& c, nMatrixT<nTYPE>& block) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (block.Rows() != r.Length() ||
	    block.Cols() != c.Length()) ExceptionT::SizeMismatch("nMatrixT<nTYPE>::CopyBlock");
#endif
	
	nTYPE* pblock = block.Pointer();
	const int* pc = c.Pointer();
	for (int j = 0; j < c.Length(); j++)
	{
		const nTYPE* pcol = (*this)(*pc++);
		const int* pr = r.Pointer();
		for (int i = 0; i < r.Length(); i++)
			*pblock++ = pcol[*pr++];
	}
}

/* dimensions */
template <class nTYPE>
inline int nMatrixT<nTYPE>::Rows(void) const { return(fRows); }

template <class nTYPE>
inline int nMatrixT<nTYPE>::Cols(void) const { return(fCols); }

/* copy/assignment operators */
template <class nTYPE>
inline nMatrixT<nTYPE>& nMatrixT<nTYPE>::operator=(const nMatrixT& RHS)
{
	/* inherited */
	nArrayT<nTYPE>::operator=(RHS);

	/* set dimensions */
	fRows = RHS.fRows;
	fCols = RHS.fCols;

	return(*this);
}

template <class nTYPE>
inline nMatrixT<nTYPE>& nMatrixT<nTYPE>::operator=(const nTYPE& value)
{
	/* inherited */
	nArrayT<nTYPE>::operator=(value);
	return(*this);
}

/* exchange data */
template <class nTYPE>
inline void nMatrixT<nTYPE>::Swap(nMatrixT<nTYPE>& source)
{
	/* inherited */
	nArrayT<nTYPE>::Swap(source);

	/* dimensions */
	int tmp = fRows;
	fRows = source.fRows;
	source.fRows = tmp;

	tmp = fCols;
	fCols = source.fCols;
	source.fCols = tmp;
}

/* selected row(s) or column(s) */
template <class nTYPE>
void nMatrixT<nTYPE>::CopyRow(int rownum, ArrayT<nTYPE>& row) const
{
/* dimension check */
#if __option(extended_errorcheck)
	if (row.Length() != fCols) ExceptionT::SizeMismatch("nMatrixT");
	if (rownum < 0 || rownum >= fRows) ExceptionT::OutOfRange("nMatrixT");
#endif

	nTYPE* prow = row.Pointer();
	const nTYPE* pthis = this->Pointer() + rownum;
	for (int i = 0; i < fCols; i++)
	{
		*prow++ = *pthis;
		pthis += fRows;
	}
}

template <class nTYPE>
void nMatrixT<nTYPE>::CopyFromRow(int rownum, int start_col,
	ArrayT<nTYPE>& row) const
{
/* dimension check */
#if __option(extended_errorcheck)
	if (start_col < 0) ExceptionT::OutOfRange("nMatrixT");
	if (start_col + row.Length() > fCols) ExceptionT::SizeMismatch("nMatrixT");
	if (rownum < 0 || rownum >= fRows) ExceptionT::OutOfRange("nMatrixT");
#endif

	int num_vals = row.Length();
	nTYPE* prow  = row.Pointer();
	const nTYPE* pthis = this->Pointer(start_col) + rownum;
	for (int i = 0; i < fCols; i++)
	{
		*prow++ = *pthis;
		pthis += fRows;
	}
}

template <class nTYPE>
void nMatrixT<nTYPE>::CopyRows(const ArrayT<int>& rows,
	const nMatrixT<nTYPE>& source)
{
/* dimension check */
#if __option(extended_errorcheck)
	if (fCols != source.Cols() || fRows != rows.Length())
		ExceptionT::SizeMismatch("nMatrixT");
#endif

	int* prows = rows.Pointer();
	for (int i = 0; i < rows.Length(); i++)
	{
		double* psrc  = source.Pointer(*prows++);
	    double* pthis = this->Pointer(i);
		
		for (int j = 0; j < fCols; j++)
		{
			*pthis = *psrc;
		
			pthis += fRows;
			psrc  += source.Rows();
		}
	}
}

template <class nTYPE>
void nMatrixT<nTYPE>::CopyColumn(int colnum, ArrayT<nTYPE>& col) const
{
/* dimension check */
#if __option(extended_errorcheck)
	if (col.Length() != fRows) ExceptionT::SizeMismatch("nMatrixT");
#endif
	
	/* byte copy */
	MemCopy(col.Pointer(), (*this)(colnum), fRows);
}
	
template <class nTYPE>
void nMatrixT<nTYPE>::CopyColumns(const ArrayT<int>& cols,
	const nMatrixT<nTYPE>& source)
{
/* dimension check */
#if __option(extended_errorcheck)
	if (fRows != source.Rows || fCols != cols.Length())
		ExceptionT::SizeMismatch("nMatrixT");
#endif

	int* prows = cols.Pointer();
	for (int i = 0; i < cols.Length(); i++)
	{
		double* psrc  = source(*prows++);
	    double* pthis = (*this)(i);
		
		for (int j = 0; j < fRows; j++)
			*pthis++ = *psrc++;
	}
}

template <class nTYPE>
void nMatrixT<nTYPE>::ColumnAlias(int colnum, ArrayT<nTYPE>& col) const
{
	col.Alias(fRows, (*this)(colnum));
}	

/*
* Creates a symmetric matrix, assuming the data is stored
* in the upper triangle of the matrix.  Setting IsUpper = 0,
* copies the data from the lower triangle.
*/
template <class nTYPE>
void nMatrixT<nTYPE>::CopySymmetric(int IsUpper)
{
/* must be square */
#if __option (extended_errorcheck)
	if (fRows != fCols) ExceptionT::GeneralFail("nMatrixT");
#endif

	int col, dex, row;

	/* copy from the upper triangle */
	if (IsUpper)
	{
		for (col = 1; col < fCols; col++)
		{
			dex = col*fRows; /* top of the column */
			
			for (row = 0; row < col; row++)
				this->fArray[row*fCols + col] = this->fArray[dex++];
		}
	}
	/* copy from the lower triangle */
	else
	{
		for (col = 0; col < fCols; col++)
		{
			dex = col*fRows + col + 1; /* just below the diagonal */
			
			for (row = col+1; row < fRows; row++)
				this->fArray[row*fCols + col] = this->fArray[dex++];
		}
	}
}

/* tranposition */
template <class nTYPE>
nMatrixT<nTYPE>& nMatrixT<nTYPE>::Transpose(const nMatrixT<nTYPE>& matrix, int fillmode)
{
#if __option (extended_errorcheck)	
	if (fRows != matrix.fCols ||
	    fCols != matrix.fRows) ExceptionT::SizeMismatch("nMatrixT");
#endif

	/* selve transposition */
	if (this->fArray == matrix.fArray) return Transpose(fillmode);

	nTYPE *pthis = this->fArray;
	nTYPE *pm    = matrix.fArray;
	if (fillmode == kOverwrite)
	{
		for (int i = 0; i < matrix.fRows; i++)
		{
			nTYPE* pmj = pm++;
			for (int j = 0; j < matrix.fCols; j++)
			{
				*pthis++ = *pmj;
				    pmj += matrix.fRows;
			}
		}
	}
	else if (fillmode == kAccumulate)
	{
		for (int i = 0; i < matrix.fRows; i++)
		{
			nTYPE* pmj = pm++;
			for (int j = 0; j < matrix.fCols; j++)
			{
				*pthis++ += *pmj;
				    pmj += matrix.fRows;
			}
		}
	}
	else ExceptionT::GeneralFail("nMatrixT<nTYPE>::Transpose", "unrecognized fill mode: %d", fillmode);

	return *this;
}

template <class nTYPE>
nMatrixT<nTYPE>& nMatrixT<nTYPE>::Transpose(int fillmode)
{
	if (fillmode == kOverwrite)
	{
		register nTYPE temp;
		for (int i = 0; i < fRows - 1; i++)
		{
			nTYPE* prow = (*this)(i+1) + i;
			nTYPE* pcol = (*this)(i) + i + 1;
			for (int j = i + 1; j < fCols; j++)
			{
				temp  = *prow;
				*prow = *pcol;
				*pcol = temp;
	
				pcol++;
				prow += fRows;
			}
		}
	}
	else if (fillmode == kAccumulate)
	{
		register nTYPE temp;
		for (int i = 0; i < fRows - 1; i++)
		{
			nTYPE* prow = (*this)(i+1) + i;
			nTYPE* pcol = (*this)(i) + i + 1;
			for (int j = i + 1; j < fCols; j++)
			{
				temp  = *prow;
				*prow += *pcol;
				*pcol += temp;
	
				pcol++;
				prow += fRows;
			}
		}
	}
	else ExceptionT::GeneralFail("nMatrixT<nTYPE>::Transpose", "unrecognized fill mode: %d", fillmode);

	return *this;
}

/*
* Matrix multiplication - operates on this using a and b.
* Operations allowed on entire matrices only, ie no apparent
* dimensions.
*/

/* this(i,j) = A(i,k)*B(k,j) , (this = A * B) 	*/
template <class nTYPE>
void nMatrixT<nTYPE>::MultAB(const nMatrixT& A, const nMatrixT& B, int upper)
{
	/* dimension checks */
#if __option (extended_errorcheck)
	if (fRows != A.fRows ||
	    fCols != B.fCols ||
	  A.fCols != B.fRows) ExceptionT::SizeMismatch("nMatrixT");
#endif

	/* 2 x 2 specialization */
	if (fRows == 2 && fCols == 2 && A.fCols == 2) {
		nTYPE* c = this->Pointer();
		const nTYPE* a = A.Pointer();
		const nTYPE* b = B.Pointer();

		c[0] = a[0]*b[0] + a[2]*b[1];
		c[1] = a[1]*b[0] + a[3]*b[1];
		c[2] = a[0]*b[2] + a[2]*b[3]; 
		c[3] = a[1]*b[2] + a[3]*b[3];
	}	
	/* 3 x 3 specialization */
	else if (fRows == 3 && fCols == 3 && A.fCols == 3) {
		nTYPE* c = this->Pointer();
		const nTYPE* a = A.Pointer();
		const nTYPE* b = B.Pointer();

		c[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
		c[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
		c[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
		c[3] = a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
		c[4] = a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
		c[5] = a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
		c[6] = a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
		c[7] = a[1]*b[6] + a[4]*b[7] + a[7]*b[8]; 
		c[8] = a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
	}
	else if (!upper)	/* entire matrix */		
	{		
		int dotcount = A.fCols;
		nTYPE* c = this->Pointer();
		const nTYPE* BCol = B.Pointer();

		register nTYPE temp;
		register nTYPE sum;
		
		for (int Bcol = 0; Bcol < fCols; Bcol++)
		{
			const nTYPE* ARow = A.Pointer();
		 		 	
			for (int Arow = 0; Arow < fRows; Arow++)
			{
				sum = 0.0;
				
				const nTYPE* AR = ARow;
				const nTYPE* BC = BCol;
				
				for (int i = 0; i < dotcount; i++)
				{
					temp  = *AR;
					temp *= *BC++;
					sum  += temp;
					
					AR += fRows;
				}
				
				*c++ = sum;
				ARow++;
			}
				
			BCol += dotcount;
		}
	}
	else /* upper triangle only - must be square */	
	{
		/* dimension checks */
#if __option (extended_errorcheck)
		if (fRows != fCols) ExceptionT::GeneralFail("nMatrixT");
#endif
		
		int dotcount = A.fCols;
		const nTYPE* BCol = B.Pointer();
		
		register nTYPE temp;
		register nTYPE sum;
		
		for (int Bcol = 0; Bcol < fCols; Bcol++)
		{
			nTYPE*  c   = this->Pointer() + Bcol*fRows;
		 	const nTYPE* ARow = A.Pointer();
		 	 	
			for (int Arow = 0; Arow <= Bcol; Arow++)
			{
				sum = 0.0;
				const nTYPE* AR = ARow;
				const nTYPE* BC = BCol;
				
				for (int i = 0; i < dotcount; i++)
				{
					temp  = *AR;
					temp *= *BC++;
					sum  += temp;
					
					AR += fRows;
				}
				
				*c++ = sum;
				ARow++;
			}
			
			BCol += dotcount;
		}
	}	
}

/*                                    T		    */
/* this(i,j) = A(k,i)*B(k,j) (this = A  * B)	*/
template <class nTYPE>
void nMatrixT<nTYPE>::MultATB(const nMatrixT& A, const nMatrixT& B, int upper)
{
	/* dimension checks */
#if __option (extended_errorcheck)
	if (fRows != A.fCols ||
		fCols != B.fCols ||
	  A.fRows != B.fRows) ExceptionT::SizeMismatch("nMatrixT");
#endif

	if (!upper)	/* entire matrix */
	{
		int dotcount = A.fRows;
		nTYPE* c = this->Pointer();
		const nTYPE* BCol = B.Pointer();
		
		register nTYPE temp;
		register nTYPE sum;
		
		for (int Bcol = 0; Bcol < fCols; Bcol++)
		{
			const nTYPE* ACol = A.Pointer();
		 		 	
			for (int Acol = 0; Acol < fRows; Acol++)
			{
				sum = 0.0;
				const nTYPE* AC = ACol;
				const nTYPE* BC = BCol;
				
				for (int i = 0; i < dotcount; i++)
				{
					temp  = *AC++;
					temp *= *BC++;
				
					sum += temp;
				}
					
				*c++ = sum;
				ACol += dotcount;
			}
				
			BCol += dotcount;
		}
	}	
	else /* upper triangle only - must be square */
	{
		/* dimension checks */
#if __option (extended_errorcheck)
		if(fRows != fCols) ExceptionT::GeneralFail("nMatrixT");
#endif

	  	int dotcount = A.fRows;
		const nTYPE* BCol = B.Pointer();

		register nTYPE temp;
		register nTYPE sum;
		
		for (int Bcol = 0; Bcol < fCols; Bcol++)
		{
		 	nTYPE* c = this->Pointer() + Bcol*fRows;
		 	const nTYPE* ACol = A.Pointer();
		 	 	
			for (int Acol = 0; Acol <= Bcol; Acol++)
			{
				sum = 0.0;
				const nTYPE* AC = ACol;
				const nTYPE* BC = BCol;
				
				for (int i = 0; i < dotcount; i++)
				{
					temp  = *AC++;
					temp *= *BC++;
				
					sum  += temp;
				}
				
				*c++ = sum;
				ACol += dotcount;
			}
			
			BCol += dotcount;
		}
	}
}	

/*                               T		*/
/* (*this) = A(i,k)*B(j,k) (A * B )		*/
template <class nTYPE>
void nMatrixT<nTYPE>::MultABT(const nMatrixT& A, const nMatrixT& B, int upper)
{	
	/* dimension checks */
#if __option (extended_errorcheck)
	if (fRows != A.fRows ||
	    fCols != B.fRows ||
	  A.fCols != B.fCols) ExceptionT::SizeMismatch("nMatrixT");
#endif

	if (!upper) /* entire matrix */
	{	
		int dotcount = A.fCols;
		nTYPE* c = this->Pointer();
		const nTYPE* BRow = B.Pointer();
	
		register nTYPE temp;
		register nTYPE sum;

		for (int Brow = 0; Brow < fCols; Brow++)
		{
			const nTYPE* ARow = A.Pointer();
		 		 	
			for (int Arow = 0; Arow < fRows; Arow++)
			{
				sum = 0.0;
				const nTYPE* AR = ARow;
				const nTYPE* BR = BRow;
					
				for (int i = 0; i < dotcount; i++)
				{
					temp  = *AR;
					temp *= *BR;
					sum  += temp;
					
					AR  += fRows;
					BR  += fCols;
				}
					
				*c++ = sum;
				ARow++;
			}
				
			BRow++;
		}
	}
	else /* upper triangle only - must be square */
	{
		/* dimension checks */
#if __option (extended_errorcheck)
		if (fRows != fCols) ExceptionT::GeneralFail("nMatrixT");
#endif

		int dotcount = A.fCols;
		const nTYPE* BRow = B.Pointer();

		register nTYPE temp;
		register nTYPE sum;
		
		for (int Brow = 0; Brow < fCols; Brow++)
		{
		 	nTYPE* c = this->Pointer() + Brow*fRows;
		 	const nTYPE* ARow = A.Pointer();
		 		 	
			for (int Arow = 0; Arow <= Brow; Arow++)
			{
				sum  = 0.0;
				const nTYPE* AR = ARow;
				const nTYPE* BR = BRow;
				
				for (int i = 0; i < dotcount; i++)
				{
					temp  = *AR;
					temp *= *BR;				
					sum  += temp;
					
					AR  += fRows;
					BR  += fCols;
				}
				
				*c++ = sum;
				ARow++;
			}
				
			BRow++;
		}
	}
}
	
/*                             T    T	*/
/* (*this) = A(k,i)*B(j,k) = (A  * B )	*/
template <class nTYPE>
void nMatrixT<nTYPE>::MultATBT(const nMatrixT& A, const nMatrixT& B)
{
	/* dimension checks */
#if __option (extended_errorcheck)
	if (fRows != A.fCols &&
		fCols != B.fRows &&
A.fRows != B.fCols) ExceptionT::SizeMismatch("nMatrixT");
#endif

	int dotcount = A.fRows;
	nTYPE* cRow = this->Pointer();
	const nTYPE* ACol = A.Pointer();

	register nTYPE temp;
	register nTYPE sum;

	for (int Acol = 0; Acol < fRows; Acol++)
	{
		const nTYPE* BRow = B.Pointer();
		nTYPE* cR = cRow;
		 		 	
		for (int Brow = 0; Brow < fCols; Brow++)
		{
			sum = 0.0;
			const nTYPE* AC = ACol;
			const nTYPE* BR = BRow;
					
			for (int i = 0; i < dotcount; i++)
			{
				temp  = *AC++;
				temp *= *BR;			
				sum  += temp;
				
				BR  += fCols;
			}
			
			*cR = sum;		
			cR += fRows;
			BRow++;
		}
				
			ACol += dotcount;
			cRow++;
	}
}

/* matrix-matrix-matrix operations, ie. tensor basis transformations */

/* this_ij = p_iI b_IJ q_Ji */
template <class nTYPE>
void nMatrixT<nTYPE>::MultABC(const nMatrixT& p, const nMatrixT& b, const nMatrixT& q,
	int range, int fillmode)
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (fRows != p.fRows ||
	    fCols != q.fCols ||
	  b.fRows != p.fCols ||
	  b.fCols != q.fRows) ExceptionT::SizeMismatch("nMatrixT");
#endif

	/* initialize */
	if (fillmode == kOverwrite) *this = 0;

	register nTYPE temp;
	register nTYPE bq_Ji;
	
	nTYPE* pthisj = this->Pointer();
	const nTYPE* pqj = q.Pointer();
	for (int j = 0; j < fCols; j++)
	{
		int istop = (range == kUpperOnly) ? j + 1 : fRows;
	
		const nTYPE* pbI = b.Pointer();
		const nTYPE* ppI = p.Pointer();
		for (int I = 0; I < b.fRows; I++)
		{			
			const nTYPE* pbJ = pbI++;
			const nTYPE* pqJ = pqj;
			bq_Ji = 0.0;
			for (int J = 0; J < b.fCols; J++)
			{
				temp  = *pbJ;
				temp *= *pqJ;
				
				bq_Ji += temp;
				
				pbJ += b.fRows;
				pqJ++;
			}

			const nTYPE* ppi = ppI;
			nTYPE* pthisi = pthisj;
			for (int i = 0; i < istop; i++)
			{
				temp  = *ppi++;
				temp *= bq_Ji;
				
				*pthisi++ += temp;
			}
			ppI += p.fRows;
		}
		pthisj += fRows;
		pqj    += q.fRows;
	}
}

/* this_ij = p_iI b_IJ q_jJ */
template <class nTYPE>
void nMatrixT<nTYPE>::MultABCT(const nMatrixT& p, const nMatrixT& b, const nMatrixT& q,
	int range, int fillmode)
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (fRows != p.fRows ||
	    fCols != q.fRows ||
	  b.fRows != p.fCols ||
	  b.fCols != q.fCols) ExceptionT::SizeMismatch("nMatrixT");
#endif

	/* initialize */
	if (fillmode == kOverwrite) *this = 0;

	register nTYPE temp;
	register nTYPE bqT_Ij;
	
	nTYPE* pthisj = this->Pointer();
	const nTYPE* pqj = q.Pointer();
	for (int j = 0; j < fCols; j++)
	{
		int istop = (range == kUpperOnly) ? j + 1 : fRows;
	
		const nTYPE* pbI = b.Pointer();
		const nTYPE* ppI = p.Pointer();
		for (int I = 0; I < b.fRows; I++)
		{			
			const nTYPE* pbJ = pbI++;
			const nTYPE* pqJ = pqj;
			bqT_Ij = 0.0;
			for (int J = 0; J < b.fCols; J++)
			{
				temp  = *pbJ;
				temp *= *pqJ;
				
				bqT_Ij += temp;
				
				pbJ += b.fRows;
				pqJ += q.fRows;
			}

			const nTYPE* ppi = ppI;
			nTYPE* pthisi = pthisj;
			for (int i = 0; i < istop; i++)
			{
				temp  = *ppi++;
				temp *= bqT_Ij;
				
				*pthisi++ += temp;
			}
			ppI += p.fRows;
		}
		pthisj += fRows;
		pqj++;
	}
}

/* this_IJ = p_iI b_ij q_jJ */
template <class nTYPE>
void nMatrixT<nTYPE>::MultATBC(const nMatrixT& p, const nMatrixT& b, const nMatrixT& q,
	int range, int fillmode)
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (fRows != p.fCols ||
	    fCols != q.fCols ||
	  b.fRows != p.fRows ||
	  b.fCols != q.fRows) ExceptionT::SizeMismatch("nMatrixT");
#endif

	/* initialize */
	if (fillmode == kOverwrite) *this = 0;

	register nTYPE temp;
	register nTYPE bq_iJ;
	
	nTYPE* pthisJ = this->Pointer();
	const nTYPE* pqJ = q.Pointer();
	for (int J = 0; J < fCols; J++)
	{
		int Istop = (range == kUpperOnly) ? J + 1 : fRows;
	
		const nTYPE* ppi = p.Pointer();
		const nTYPE* pbi = b.Pointer();
		for (int i = 0; i < b.fRows; i++)
		{			
			const nTYPE* pbj = pbi++;
			const nTYPE* pqj = pqJ;
			bq_iJ = 0.0;
			for (int j = 0; j < b.fCols; j++)
			{
				temp  = *pbj;
				temp *= *pqj++;
				
				bq_iJ += temp;
				
				pbj += b.fRows;
			}

			const nTYPE* ppI = ppi++;
			nTYPE* pthisI = pthisJ;
			for (int I = 0; I < Istop; I++)
			{
				temp  = *ppI;
				temp *= bq_iJ;
				
				*pthisI++ += temp;
				
				ppI += p.fRows;
			}
		}
		pthisJ += fRows;
		pqJ    += q.fRows;
	}
}

//TEMP - debugging
#if 0
/* this_ij = q_iI b_IJ q_jJ */
template <class nTYPE>
inline void nMatrixT<nTYPE>::MultQBQT(const nMatrixT& q,
	const nMatrixT& b, int range, int fillmode) 
{
	MultABCT(q, b, q, range, fillmode);
}

/* this_IJ = q_iI b_ij q_jJ */
template <class nTYPE>
inline void nMatrixT<nTYPE>::MultQTBQ(const nMatrixT& q,
	const nMatrixT& b, int range, int fillmode)
{
	MultATBC(q, b, q, range, fillmode);
}
#else
/* this_ij = q_iI b_IJ q_jJ */
template <class nTYPE>
void nMatrixT<nTYPE>::MultQBQT(const nMatrixT& q,
	const nMatrixT& b, int range, int fillmode)
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (fRows != q.fRows ||
	    fCols != q.fRows ||
	  b.fRows != q.fCols ||
	  b.fCols != q.fCols) ExceptionT::SizeMismatch("nMatrixT");
#endif

	/* initialize */
	if (fillmode == kOverwrite) *this = 0;

	register nTYPE temp;
	register nTYPE bqT_Ij;
	
	nTYPE* pthisj = this->Pointer();
	const nTYPE* pqj = q.Pointer();

	for (int j = 0; j < fCols; j++)
	{
		int istop = (range == kUpperOnly) ? j + 1 : fRows;
	
		const nTYPE* pbI = b.Pointer();
		const nTYPE* pqI = q.Pointer();
	
		for (int I = 0; I < b.fRows; I++)
		{			
			const nTYPE* pbJ = pbI++;
			const nTYPE* pqJ = pqj;
			
			bqT_Ij = 0.0;
			for (int J = 0; J < b.fCols; J++)
			{
				temp  = *pbJ;
				temp *= *pqJ;
				
				bqT_Ij += temp;
				
				pbJ += b.fRows;
				pqJ += q.fRows;
			}

			const nTYPE* pqi = pqI;
			nTYPE* pthisi = pthisj;
			
			for (int i = 0; i < istop; i++)
			{
				temp  = *pqi++;
				temp *= bqT_Ij;
				
				*pthisi++ += temp;
			}
			
			pqI += q.fRows;
		}
		
		pthisj += fRows;
		pqj++;
	}
}

/* this_IJ = q_iI b_ij q_jJ */
template <class nTYPE>
void nMatrixT<nTYPE>::MultQTBQ(const nMatrixT& q,
	const nMatrixT& b, int range, int fillmode)
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (fRows != q.fCols ||
	    fCols != q.fCols ||
	  b.fRows != q.fRows ||
	  b.fCols != q.fRows) ExceptionT::SizeMismatch("nMatrixT");
#endif

	if (fRows == 3 && fCols == 3 && b.fCols == 3)
	{
		nTYPE* A = this->Pointer();
		const nTYPE* Q = q.Pointer();
		const nTYPE* B = b.Pointer();

		/* initialize */
		if (fillmode == kOverwrite) {
			*A++ = 0.0;
			*A++ = 0.0;
			*A++ = 0.0;
			*A++ = 0.0;
			*A++ = 0.0;
			*A++ = 0.0;
			*A++ = 0.0;
			*A++ = 0.0;
			*A   = 0.0;
			A = this->Pointer();
		}

		/* diagonal and upper triangle */
		A[0] += Q[0]*(B[0]*Q[0] + B[1]*Q[1] + B[2]*Q[2]) + Q[1]*(B[3]*Q[0] + B[4]*Q[1] + B[5]*Q[2]) + Q[2]*(B[6]*Q[0] + B[7]*Q[1] + B[8]*Q[2]);
		A[3] += (B[0]*Q[0] + B[1]*Q[1] + B[2]*Q[2])*Q[3] + (B[3]*Q[0] + B[4]*Q[1] + B[5]*Q[2])*Q[4] + (B[6]*Q[0] + B[7]*Q[1] + B[8]*Q[2])*Q[5];
		A[4] += Q[3]*(B[0]*Q[3] + B[1]*Q[4] + B[2]*Q[5]) + Q[4]*(B[3]*Q[3] + B[4]*Q[4] + B[5]*Q[5]) + Q[5]*(B[6]*Q[3] + B[7]*Q[4] + B[8]*Q[5]);
		A[6] += (B[0]*Q[0] + B[1]*Q[1] + B[2]*Q[2])*Q[6] + (B[3]*Q[0] + B[4]*Q[1] + B[5]*Q[2])*Q[7] + (B[6]*Q[0] + B[7]*Q[1] + B[8]*Q[2])*Q[8];
		A[7] += (B[0]*Q[3] + B[1]*Q[4] + B[2]*Q[5])*Q[6] + (B[3]*Q[3] + B[4]*Q[4] + B[5]*Q[5])*Q[7] + (B[6]*Q[3] + B[7]*Q[4] + B[8]*Q[5])*Q[8];
		A[8] += Q[6]*(B[0]*Q[6] + B[1]*Q[7] + B[2]*Q[8]) + Q[7]*(B[3]*Q[6] + B[4]*Q[7] + B[5]*Q[8]) + Q[8]*(B[6]*Q[6] + B[7]*Q[7] + B[8]*Q[8]);
	
		/* lower triangle */
		if (range == kWhole) {
			A[1] += Q[0]*(B[0]*Q[3] + B[1]*Q[4] + B[2]*Q[5]) + Q[1]*(B[3]*Q[3] + B[4]*Q[4] + B[5]*Q[5]) + Q[2]*(B[6]*Q[3] + B[7]*Q[4] + B[8]*Q[5]);
			A[2] += Q[0]*(B[0]*Q[6] + B[1]*Q[7] + B[2]*Q[8]) + Q[1]*(B[3]*Q[6] + B[4]*Q[7] + B[5]*Q[8]) + Q[2]*(B[6]*Q[6] + B[7]*Q[7] + B[8]*Q[8]);
			A[5] += Q[3]*(B[0]*Q[6] + B[1]*Q[7] + B[2]*Q[8]) + Q[4]*(B[3]*Q[6] + B[4]*Q[7] + B[5]*Q[8]) + Q[5]*(B[6]*Q[6] + B[7]*Q[7] + B[8]*Q[8]);
		}
	}
	else
	{
		/* initialize */
		if (fillmode == kOverwrite) *this = 0;

		register nTYPE temp;
		register nTYPE bq_iJ;
	
		nTYPE* pthisJ = this->Pointer();
		const nTYPE* pqJ = q.Pointer();

		for (int J = 0; J < fCols; J++)
		{
			int Istop = (range == kUpperOnly) ? J + 1 : fRows;
	
			const nTYPE* pqi = q.Pointer();
			const nTYPE* pbi = b.Pointer();
	
			for (int i = 0; i < b.fRows; i++)
			{			
				const nTYPE* pbj = pbi++;
				const nTYPE* pqj = pqJ;
			
				bq_iJ = 0.0;
				for (int j = 0; j < b.fCols; j++)
				{
					temp  = *pbj;
					temp *= *pqj++;
					bq_iJ += temp;
				
					pbj += b.fRows;
				}

				const nTYPE* pqI = pqi++;
				nTYPE* pthisI = pthisJ;
			
				for (int I = 0; I < Istop; I++)
				{
					temp  = *pqI;
					temp *= bq_iJ;
				
					*pthisI++ += temp;
				
					pqI += q.fRows;
				}
			}
		
			pthisJ += fRows;
			pqJ    += q.fRows;
		}
	}
}
#endif

/* matrix-vector multiplication */

/* b_i = A_ij*x_j */
template <class nTYPE>
inline void nMatrixT<nTYPE>::Multx(const nArrayT<nTYPE>& x, nArrayT<nTYPE>& b, 
	const nTYPE& scale, int fillmode) const
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (fRows != b.Length() || fCols != x.Length()) ExceptionT::SizeMismatch("nMatrixT::Multx");
#endif

	Multx(x.Pointer(), b.Pointer(), scale, fillmode);
}

template <class nTYPE>
void nMatrixT<nTYPE>::Multx(const nTYPE* x, nTYPE* b, const nTYPE& scale, int fillmode) const
{
	const nTYPE* ARow = this->Pointer();
	const nTYPE* px0 = x;
	nTYPE* pb = b;

	register nTYPE temp;
	register nTYPE sum;

	if (fillmode == kOverwrite)
	{
		for (int i = 0; i < fRows; i++)
		{
			sum = 0.0;
			const nTYPE *px = px0;
			const nTYPE *AR = ARow;
			for (int j = 0; j < fCols; j++)
			{
				temp  = *px++;
				temp *= *AR;		
				 sum += temp;
			
				AR += fRows;	
			}
			sum *= scale;
			*pb++ = sum;
			ARow++;
		}
	}
	else if (fillmode == kAccumulate)
	{
		for (int i = 0; i < fRows; i++)
		{
			sum = 0.0;
			const nTYPE *px = px0;
			const nTYPE *AR = ARow;
			for (int j = 0; j < fCols; j++)
			{
				temp  = *px++;
				temp *= *AR;		
				 sum += temp;
			
				AR += fRows;	
			}
			sum *= scale;
			*pb++ += sum;
			ARow++;
		}
	}
	else ExceptionT::GeneralFail("nMatrixT::Multx");
}

/* b_i = A_ji*x_j */
template <class nTYPE>
inline void nMatrixT<nTYPE>::MultTx(const nArrayT<nTYPE>& x, nArrayT<nTYPE>& b,
	const nTYPE& scale, int fillmode) const
{
	/* dimension checks */
#if __option (extended_errorcheck)
	if (fRows != x.Length() || fCols != b.Length()) ExceptionT::SizeMismatch("nMatrixT::MultTx");
#endif

	MultTx(x.Pointer(), b.Pointer(), scale, fillmode);
}

template <class nTYPE>
void nMatrixT<nTYPE>::MultTx(const nTYPE* x, nTYPE* b, const nTYPE& scale, int fillmode) const
{
	const nTYPE* ARow = this->Pointer();
	const nTYPE* px0 = x;
	nTYPE* pb = b;

	register nTYPE temp;
	register nTYPE sum;

	if (fillmode == kOverwrite)
	{
		for (int i = 0; i < fCols; i++)
		{
			sum = 0.0;
			const nTYPE *px = px0;
			const nTYPE *AR = ARow;
			for (int j = 0; j < fRows; j++)
			{
				temp  = *px++;
				temp *= *AR++;
				sum  += temp;
			}
			sum *= scale;
			*pb++ = sum;
			ARow += fRows;
		}
	}
	else if (fillmode == kAccumulate)
	{
		for (int i = 0; i < fCols; i++)
		{
			sum = 0.0;
			const nTYPE *px = px0;
			const nTYPE *AR = ARow;
			for (int j = 0; j < fRows; j++)
			{
				temp  = *px++;
				temp *= *AR++;
				sum  += temp;
			}
			sum *= scale;
			*pb++ += sum;
			ARow += fRows;
		}
	}
	else ExceptionT::GeneralFail("nMatrixT::Multx");
}

/* vector-matrix-vector product */
template <class nTYPE>
nTYPE nMatrixT<nTYPE>::MultmBn(const nArrayT<nTYPE>& m,
	const nArrayT<nTYPE>& n) const
{
	/* dimension checks */
#if __option (extended_errorcheck)
	if (fRows != m.Length() ||
	    fCols != n.Length()) ExceptionT::SizeMismatch("nMatrixT");
#endif

	register nTYPE product = 0.0;
	register nTYPE mB;
	register nTYPE temp;
	
	const nTYPE* pnj = n.Pointer();
	const nTYPE* pBj = this->Pointer();

	for (int j = 0; j < fCols; j++)
	{
		mB = 0.0;
		const nTYPE* pBi = pBj;
		const nTYPE* pmi = m.Pointer();
	
		for (int i = 0; i < fRows; i++)
		{
			temp  = *pBi++;
			temp *= *pmi++;
		
			mB   += temp;		
		}
		
		mB      *= *pnj++;
		product += mB;
		
		pBj += fRows;
	}

	return(product);
}	

/* dyadic product */
template <class nTYPE>
inline void nMatrixT<nTYPE>::Outer(const nArrayT<nTYPE>& v1, const nArrayT<nTYPE>& v2,
	const nTYPE& scale, int fillmode)
{
	/* dimension checks */
#if __option (extended_errorcheck)
	if (v1.Length() != fRows || v2.Length() != fCols) ExceptionT::SizeMismatch("nMatrixT");
#endif

	/* inherited */
	nMatrixT<nTYPE>::Outer(v1.Pointer(), v2.Pointer(), scale, fillmode);
}

template <class nTYPE>
void nMatrixT<nTYPE>::Outer(const nTYPE* v1, const nTYPE* v2, const nTYPE& scale, int fillmode)
{
	if (fillmode == kOverwrite)
	{
		if (fCols == 2 && fRows == 2)
		{
			nTYPE* v1v2 = this->Pointer();
			nTYPE v10 = scale*v1[0];
			nTYPE v11 = scale*v1[1];

			v1v2[0] = v10*v2[0];
			v1v2[1] = v11*v2[0];
			v1v2[2] = v10*v2[1];
			v1v2[3] = v11*v2[1];
		}
		else if (fCols == 3 && fRows == 3)
		{
			nTYPE* v1v2 = this->Pointer();
			nTYPE v10 = scale*v1[0];
			nTYPE v11 = scale*v1[1];
			nTYPE v12 = scale*v1[2];

			v1v2[0] = v10*v2[0];
			v1v2[1] = v11*v2[0];
			v1v2[2] = v12*v2[0];

			v1v2[3] = v10*v2[1];
			v1v2[4] = v11*v2[1];
			v1v2[5] = v12*v2[1];

			v1v2[6] = v10*v2[2];
			v1v2[7] = v11*v2[2];
			v1v2[8] = v12*v2[2];
		}
		else
		{
			nTYPE* pthis = this->Pointer();
			const nTYPE* pv1 = v1;
			const nTYPE* pv2 = v2;
			for (int j = 0; j < fCols; j++)
			{
				const nTYPE* pcol = pv1;
				for (int i = 0; i < fRows; i++)
				{
					*pthis    = scale;
					*pthis   *= *pcol++;
					*pthis++ *= *pv2;
				}
				pv2++;
			}
		}
	}
	else if (fillmode == kAccumulate)
	{
		if (fCols == 2 && fRows == 2)
		{
			nTYPE* v1v2 = this->Pointer();
			nTYPE v10 = scale*v1[0];
			nTYPE v11 = scale*v1[1];

			v1v2[0] += v10*v2[0];
			v1v2[1] += v11*v2[0];
			v1v2[2] += v10*v2[1];
			v1v2[3] += v11*v2[1];
		}
		else if (fCols == 3 && fRows == 3)
		{
			nTYPE* v1v2 = this->Pointer();
			nTYPE v10 = scale*v1[0];
			nTYPE v11 = scale*v1[1];
			nTYPE v12 = scale*v1[2];

			v1v2[0] += v10*v2[0];
			v1v2[1] += v11*v2[0];
			v1v2[2] += v12*v2[0];

			v1v2[3] += v10*v2[1];
			v1v2[4] += v11*v2[1];
			v1v2[5] += v12*v2[1];

			v1v2[6] += v10*v2[2];
			v1v2[7] += v11*v2[2];
			v1v2[8] += v12*v2[2];
		}
		else
		{
			nTYPE* pthis = this->Pointer();
			const nTYPE* pv1 = v1;
			const nTYPE* pv2 = v2;
			register nTYPE temp;
			for (int j = 0; j < fCols; j++)
			{
				const nTYPE* pcol = pv1;
				for (int i = 0; i < fRows; i++)
				{
					temp  = scale;
					temp *= *pcol++;
					temp *= *pv2;
				
					*pthis++ += temp;
				}
				pv2++;
			}
		}
	}
	else ExceptionT::GeneralFail("nMatrixT<nTYPE>::Outer");
}

/* identity operations - square matrices ONLY */
template <class nTYPE>
inline void nMatrixT<nTYPE>::PlusIdentity(const nTYPE& value)
{
/* must be square */
#if __option (extended_errorcheck)
	if (fRows != fCols) ExceptionT::GeneralFail("nMatrixT");
#endif

	if (fRows == 2)
	{
		this->fArray[0] += value;	
		this->fArray[3] += value;	
	}
	else if (fRows == 3)
	{
		this->fArray[0] += value;	
		this->fArray[4] += value;	
		this->fArray[8] += value;	
	}
	else
	{
		nTYPE* dex = this->Pointer();
		int inc = fRows + 1;
		for (int i = 0; i < fRows; i++)
		{
			*dex += value;
			 dex += inc;
		}	
	}
}

template <class nTYPE>
nMatrixT<nTYPE>& nMatrixT<nTYPE>::Identity(const nTYPE& value)
{
/* must be square */
#if __option (extended_errorcheck)
	if (fRows != fCols) ExceptionT::GeneralFail("nMatrixT");
#endif
	
	if (fRows == 2)
	{
		nTYPE* p = this->Pointer();
		*p++ = value;	
		*p++ = 0.0;	
		*p++ = 0.0;	
		*p   = value;
	}
	else if (fRows == 3)
	{
		nTYPE* p = this->Pointer();
		*p++ = value;	
		*p++ = 0.0;	
		*p++ = 0.0;	

		*p++ = 0.0;	
		*p++ = value;	
		*p++ = 0.0;	

		*p++ = 0.0;	
		*p++ = 0.0;	
		*p   = value;	
	}
	else
	{
		*this = 0.0;
		nTYPE* dex = this->Pointer();
		int inc = fRows + 1;
		for (int i = 0; i < fRows; i++)
		{
			*dex  = value;
			 dex += inc;
		}
	}
	
	return *this;	
}

/* writing to rows/columns */
template <class nTYPE>
inline void nMatrixT<nTYPE>::SetRow(int row, const nArrayT<nTYPE>& vec)
{
	/* dimension check */
#if __option (extended_errorcheck)
	if (vec.Length() != fCols) ExceptionT::SizeMismatch("nMatrixT");
	if (row < 0 || row >= fRows) ExceptionT::OutOfRange("nMatrixT");
#endif
	
	SetRow(row, vec.Pointer());
}

template <class nTYPE>
inline void nMatrixT<nTYPE>::SetRow(int row, const nTYPE* vec)
{
	nTYPE* pcol = this->Pointer() + row;	
	for (int i = 0; i < fCols; i++)
	{
		*pcol = *vec++;
		pcol += fRows;
	}
}

template <class nTYPE>
inline void nMatrixT<nTYPE>::SetRow(int row, const nTYPE& value)
{
	nTYPE* pcol = this->Pointer() + row;	
	for (int i = 0; i < fCols; i++)
	{
		*pcol = value;
		pcol += fRows;
	}
}

template <class nTYPE>
inline void nMatrixT<nTYPE>::SetCol(int col, const nArrayT<nTYPE>& vec)
{
/* dimension check */
#if __option (extended_errorcheck)
	if (vec.Length() != fRows) ExceptionT::OutOfRange("nMatrixT");
#endif

	SetCol(col, vec.Pointer());
}

template <class nTYPE>
inline void nMatrixT<nTYPE>::SetCol(int col, const nTYPE* vec)
{
	nTYPE* pcol = (*this)(col);	
	for (int i = 0; i < fRows; i++)
		*pcol++ = *vec++;
}

template <class nTYPE>
inline void nMatrixT<nTYPE>::SetCol(int col, const nTYPE& value)
{
	nTYPE* pcol = (*this)(col);	
	for (int i = 0; i < fRows; i++)
		*pcol++ = value;
}

/* dot the specified row/column number with the array */
template <class nTYPE>
inline nTYPE nMatrixT<nTYPE>::DotRow(int rownum,
	const nArrayT<nTYPE>& vec) const
{
#if __option (extended_errorcheck)
	/* dimension check */
	if (vec.Length() != fCols) ExceptionT::SizeMismatch("nMatrixT");
#endif
	return DotRow(rownum, vec.Pointer());
}

template <class nTYPE>
inline nTYPE nMatrixT<nTYPE>::DotRow(int rownum, const nTYPE* pvec) const
{
	const nTYPE *p = &(*this)(rownum,0);
	register nTYPE sum = 0.0;
	register nTYPE temp;
	for (int i = 0; i < fCols; i++)
	{
		temp  = *p;
		temp *= *pvec++;
		sum  += temp;
		p    += fRows;
	}
	return sum;
}

template <class nTYPE>
inline nTYPE nMatrixT<nTYPE>::DotCol(int colnum,
	const nArrayT<nTYPE>& vec) const
{
#if __option (extended_errorcheck)
	/* dimension check */
	if (vec.Length() != fRows) ExceptionT::SizeMismatch("nMatrixT");
#endif
	return DotCol(colnum, vec.Pointer());
}

template <class nTYPE>
inline nTYPE nMatrixT<nTYPE>::DotCol(int colnum, const nTYPE* pvec) const
{
	const nTYPE *p = (*this)(colnum);
	register nTYPE sum = 0.0;
	register nTYPE temp;
	for (int i = 0; i < fRows; i++)
	{
		temp  = *p++;
		temp *= *pvec++;
		sum  += temp;
	}
	return sum;
}

} // namespace Tahoe 
#endif /* _NMATRIX_T_H_ */
