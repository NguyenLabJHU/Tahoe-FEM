/* $Id: nMatrixT.h,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (05/24/1996)                                          */
/* 2 dimensional matrix mathematics template object.                      */

#ifndef _NMATRIX_T_H_
#define _NMATRIX_T_H_

/* base class */
#include "nArrayT.h"

/* control flags */

template <class nTYPE>
class nMatrixT: public nArrayT<nTYPE>
{
public:

	/* control flags */
	enum SymmetryFlagT {kWhole = 0, kUpperOnly = 1};

	enum AssemblyModeT {kOverwrite = 0, kAccumulate = 1};

	/* constructor*/
	nMatrixT(void);
	nMatrixT(int numrows, int numcols);
	nMatrixT(int squaredim);
	nMatrixT(int numrows, int numcols, nTYPE* p);	
	nMatrixT(const nMatrixT& source);	

	/* destructor*/
	~nMatrixT(void);

	/* post construction dimensioning */
	void Allocate(int numrows, int numcols);
	void Allocate(int squaredim);
	void Set(int numrows, int numcols, nTYPE* p);

	/* free memory (if allocated) and set size to zero */
	void Free(void);

	/* element and column accessor */
	nTYPE& operator()(int nrow, int ncol) const;
	nTYPE* operator()(int ncol) const;

	/* block accessing with row and col in the upper left */
	void AddBlock(int row, int col, const nMatrixT<nTYPE>& block);
	void SetBlock(int row, int col, const nMatrixT<nTYPE>& block);
	void CopyBlock(int row, int col, nMatrixT<nTYPE>& block) const;
	void CopyBlock(const ArrayT<int>& rc, nMatrixT<nTYPE>& block) const; //block must be square
	void CopyBlock(const ArrayT<int>& r, const ArrayT<int>& c, nMatrixT<nTYPE>& block) const;

	/* dimensions */
	int Rows(void) const;
	int Cols(void) const;

	/* copy/assignment operators */
	nMatrixT<nTYPE>& operator=(const nMatrixT& RHS);
	nMatrixT<nTYPE>& operator=(const nTYPE& value);
	void Alias(const nMatrixT& RHS);

	/* exchange data */
	void Swap(nMatrixT<nTYPE>& source);
	
	/* selected row(s) or column(s) */
	void CopyRow(int rownum, ArrayT<nTYPE>& row) const;
	void CopyFromRow(int rownum, int start_col, ArrayT<nTYPE>& row) const;
	void CopyRows(const ArrayT<int>& rows,
	              const nMatrixT<nTYPE>& source);
	void CopyColumn(int colnum, ArrayT<nTYPE>& col) const;
	void CopyColumns(const ArrayT<int>& cols,
	                 const nMatrixT<nTYPE>& source);
	void ColumnAlias(int colnum, ArrayT<nTYPE>& col) const;

	/* creates a symmetric matrix, assuming the data is stored
	 * in the upper triangle of the matrix.  Setting IsUpper = 0,
	 * copies the data from the lower triangle */
	void CopySymmetric(int IsUpper = 1);
	
	/* matrix-matrix multiplication - operates on this using a and b.
	 * Operations allowed on entire matrices only, ie no apparent
	 * dimensions */
	void MultAB(const nMatrixT& a, const nMatrixT& b, int upper = 0);
	void MultATB(const nMatrixT& a, const nMatrixT& b, int upper = 0);
	void MultABT(const nMatrixT& a, const nMatrixT& b, int upper = 0);
	void MultATBT(const nMatrixT& a, const nMatrixT& b);

	/* matrix-matrix-matrix operations, ie. tensor basis transformations */
	void MultQBQT(const nMatrixT& q, const nMatrixT& b,
		int range = kWhole, int fillmode = kOverwrite);
	void MultQTBQ(const nMatrixT& q, const nMatrixT& b,
		int range = kWhole, int fillmode = kOverwrite);
	
	/* matrix-vector multiplication - returns the result in b */
	void Multx(const nArrayT<nTYPE>& x, nArrayT<nTYPE>& b) const;
	void MultTx(const nArrayT<nTYPE>& x, nArrayT<nTYPE>& b) const;

	/* vector-matrix-vector product */
	nTYPE MultmBn(const nArrayT<nTYPE>& m,
	                 const nArrayT<nTYPE>& n) const;
	   		
	/* returns the outer product of the 2 vectors, or
	 * in dyadic notation:
	 *
	 *        v1 (x) v2 */
	void Outer(const nArrayT<nTYPE>& v1, const nArrayT<nTYPE>& v2);

	/* identity operations - square matrices ONLY */
	void PlusIdentity(const nTYPE& value = nTYPE(1.0));
	nMatrixT<nTYPE>& Identity(const nTYPE& value = nTYPE(1.0));

	/* writing to rows/columns */
	void SetRow(int row, const nArrayT<nTYPE>& vec);
	void SetRow(int row, const nTYPE& value);
	void SetCol(int col, const nArrayT<nTYPE>& vec);
	void SetCol(int col, const nTYPE& value);

	/* dot the specified row/column number with the array */
	nTYPE DotRow(int rownum, const nArrayT<nTYPE>& array) const;
	nTYPE DotCol(int colnum, const nArrayT<nTYPE>& array) const;

protected:
	
	int	fRows;
	int	fCols;
};

/* I/O operators */
template <class nTYPE>
istream& operator>>(istream& in, nMatrixT<nTYPE>& matrix)
{
	for (int j = 0; j < matrix.Rows(); j++)
		for (int i = 0; i < matrix.Cols(); i++)
				in >> matrix(j,i);
	return in;
};

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
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructor */
template <class nTYPE>
inline nMatrixT<nTYPE>::nMatrixT(void): fRows(0), fCols(0) { }

template <class nTYPE>
inline nMatrixT<nTYPE>::nMatrixT(int numrows, int numcols)
{
	Allocate(numrows,numcols);
}

template <class nTYPE>
inline nMatrixT<nTYPE>::nMatrixT(int squaredim)
{
	Allocate(squaredim);
}

template <class nTYPE>
inline nMatrixT<nTYPE>::nMatrixT(int numrows, int numcols, nTYPE* p)
{
	Set(numrows, numcols, p);
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
inline void nMatrixT<nTYPE>::Allocate(int numrows, int numcols)
{
	/* inherited */
	nArrayT<nTYPE>::Allocate(numrows*numcols);

	/* set dimensions */
	fRows = numrows;
	fCols = numcols;
}

template <class nTYPE>
inline void nMatrixT<nTYPE>::Allocate(int squaredim)
{
	Allocate(squaredim,squaredim);
}

template <class nTYPE>
inline void nMatrixT<nTYPE>::Set(int numrows, int numcols, nTYPE* p)
{
	/* inherited */
	nArrayT<nTYPE>::Set(numrows*numcols, p);

	/* set dimensions */
	fRows = numrows;
	fCols = numcols;
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
inline nTYPE& nMatrixT<nTYPE>::operator()(int nrow, int ncol) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (nrow < 0 ||
	    nrow >= fRows ||
	    ncol < 0 ||
	    ncol >= fCols) throw eOutOfRange;
#endif
	
	return(fArray[ncol*fRows + nrow]);
}

/* returns a pointer to the top of the specified column */
template <class nTYPE>
inline nTYPE* nMatrixT<nTYPE>::operator()(int ncol) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (ncol < 0 || ncol >= fCols) throw eOutOfRange;
#endif
	
	return(fArray + ncol*fRows);
}

/* assemble beginning with row and col in the upper left. */
template <class nTYPE>
void nMatrixT<nTYPE>::AddBlock(int row, int col,
const nMatrixT<nTYPE>& block)
{
/* range checking */
#if __option (extended_errorcheck)
	if (row + block.Rows() > fRows ||
	    col + block.Cols() > fCols) throw eSizeMismatch;
#endif

	double* pstart = &(*this)(row,col);
	double* pblock = block.Pointer();
	
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
	    col + block.Cols() > fCols) throw eSizeMismatch;
#endif

	double* pstart = &(*this)(row,col);
	double* pblock = block.Pointer();
	
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
	    col + block.Cols() > fCols) throw eSizeMismatch;
#endif

	double* pstart = &(*this)(row,col);
	double* pblock = block.Pointer();
	for (int i = 0; i < block.Cols(); i++)
	{
		double* pcol = pstart;
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
	    block.Cols() != rc.Length()) throw eSizeMismatch;
#endif
	
	nTYPE* pblock = block.Pointer();
	for (int j = 0; j < rc.Length(); j++)
	{
		nTYPE* pcol = (*this)(rc[j]);
		int* prc = rc.Pointer();
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
	    block.Cols() != c.Length()) throw eSizeMismatch;
#endif
	
	nTYPE* pblock = block.Pointer();
	int* pc = c.Pointer();
	for (int j = 0; j < c.Length(); j++)
	{
		nTYPE* pcol = (*this)(*pc++);
		int* pr = r.Pointer();
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

template <class nTYPE>
inline void nMatrixT<nTYPE>::Alias(const nMatrixT& RHS)
{
	Set(RHS.Rows(), RHS.Cols(), RHS.Pointer() );	
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
	if (row.Length() != fCols) throw eSizeMismatch;
	if (rownum < 0 || rownum >= fRows) throw eOutOfRange;
#endif

	nTYPE* prow  = row.Pointer();
	nTYPE* pthis = Pointer() + rownum;
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
	if (start_col < 0) throw eOutOfRange;
	if (start_col + row.Length() > fCols) throw eSizeMismatch;
	if (rownum < 0 || rownum >= fRows) throw eOutOfRange;
#endif

	int num_vals = row.Length();
	nTYPE* prow  = row.Pointer();
	nTYPE* pthis = Pointer(start_col) + rownum;
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
		throw eSizeMismatch;
#endif

	int* prows = rows.Pointer();
	for (int i = 0; i < rows.Length(); i++)
	{
		double* psrc  = source.Pointer(*prows++);
	    double* pthis = Pointer(i);
		
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
	if (col.Length() != fRows) throw eSizeMismatch;
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
		throw eSizeMismatch;
#endif

	int* prows = rows.Pointer();
	for (int i = 0; i < rows.Length(); i++)
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
	col.Set(fRows, (*this)(colnum));
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
	if (fRows != fCols) throw eGeneralFail;
#endif

	int col, dex, row;

	/* copy from the upper triangle */
	if (IsUpper)
	{
		for (col = 1; col < fCols; col++)
		{
			dex = col*fRows; /* top of the column */
			
			for (row = 0; row < col; row++)
				fArray[row*fCols + col] = fArray[dex++];
		}
	}
	/* copy from the lower triangle */
	else
	{
		for (col = 0; col < fCols; col++)
		{
			dex = col*fRows + col + 1; /* just below the diagonal */
			
			for (row = col+1; row < fRows; row++)
				fArray[row*fCols + col] = fArray[dex++];
		}
	}
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
	  A.fCols != B.fRows) throw eSizeMismatch;
#endif		

	if (!upper)	/* entire matrix */		
	{		
		int        dotcount = A.fCols;
		nTYPE*	c       = Pointer();
		nTYPE*	BCol    = B.Pointer();

		register nTYPE temp;
		register nTYPE sum;
		
		for (int Bcol = 0; Bcol < fCols; Bcol++)
		{
			nTYPE*ARow = A.Pointer();
		 		 	
			for (int Arow = 0; Arow < fRows; Arow++)
			{
				sum = 0.0;
				
				nTYPE* AR = ARow;
				nTYPE* BC = BCol;
				
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
		if (fRows != fCols) throw eGeneralFail;
#endif
		
		int 	  dotcount = A.fCols;
		nTYPE* BCol     = B.Pointer();
		
		register nTYPE temp;
		register nTYPE sum;
		
		for (int Bcol = 0; Bcol < fCols; Bcol++)
		{
			nTYPE*  c   = Pointer() + Bcol*fRows;
		 	nTYPE* ARow = A.Pointer();
		 	 	
			for (int Arow = 0; Arow <= Bcol; Arow++)
			{
				sum = 0.0;
				nTYPE* AR = ARow;
				nTYPE* BC = BCol;
				
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
	  A.fRows != B.fRows) throw eSizeMismatch;
#endif

	if (!upper)	/* entire matrix */
	{
		int		  dotcount = A.fRows;
		nTYPE*  c       = Pointer();
		nTYPE*  BCol    = B.Pointer();
		
		register nTYPE temp;
		register nTYPE sum;
		
		for (int Bcol = 0; Bcol < fCols; Bcol++)
		{
			nTYPE* ACol = A.Pointer();
		 		 	
			for (int Acol = 0; Acol < fRows; Acol++)
			{
				sum = 0.0;
				nTYPE* AC = ACol;
				nTYPE* BC = BCol;
				
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
		if(fRows != fCols) throw eGeneralFail;
#endif

	  	int 	  dotcount = A.fRows;
		nTYPE*  BCol    = B.Pointer();

		register nTYPE temp;
		register nTYPE sum;
		
		for (int Bcol = 0; Bcol < fCols; Bcol++)
		{
		 	nTYPE* c    = Pointer() + Bcol*fRows;
		 	nTYPE* ACol = A.Pointer();
		 	 	
			for (int Acol = 0; Acol <= Bcol; Acol++)
			{
				sum = 0.0;
				nTYPE* AC = ACol;
				nTYPE* BC = BCol;
				
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
	  A.fCols != B.fCols) throw eSizeMismatch;
#endif

	if (!upper) /* entire matrix */
	{	
		int		  dotcount = A.fCols;
		nTYPE*  c       = Pointer();
		nTYPE*  BRow    = B.Pointer();
	
		register nTYPE temp;
		register nTYPE sum;

		for (int Brow = 0; Brow < fCols; Brow++)
		{
			nTYPE* ARow = A.Pointer();
		 		 	
			for (int Arow = 0; Arow < fRows; Arow++)
			{
				sum = 0.0;
				nTYPE* AR = ARow;
				nTYPE* BR = BRow;
					
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
		if (fRows != fCols) throw eGeneralFail;
#endif

		int 	  dotcount = A.fCols;
		nTYPE*  BRow    = B.Pointer();

		register nTYPE temp;
		register nTYPE sum;
		
		for (int Brow = 0; Brow < fCols; Brow++)
		{
		 	nTYPE*  c   = Pointer() + Brow*fRows;
		 	nTYPE* ARow = A.Pointer();
		 		 	
			for (int Arow = 0; Arow <= Brow; Arow++)
			{
				sum  = 0.0;
				nTYPE* AR = ARow;
				nTYPE* BR = BRow;
				
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
A.fRows != B.fCols) throw eSizeMismatch;
#endif

	int 	  dotcount = A.fRows;
	nTYPE* cRow     = Pointer();
	nTYPE* ACol     = A.Pointer();

	register nTYPE temp;
	register nTYPE sum;

	for (int Acol = 0; Acol < fRows; Acol++)
	{
		nTYPE* BRow  = B.Pointer();
		nTYPE* cR    = cRow;
		 		 	
		for (int Brow = 0; Brow < fCols; Brow++)
		{
			sum = 0.0;
			nTYPE* AC = ACol;
			nTYPE* BR = BRow;
					
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
	  b.fCols != q.fCols) throw eSizeMismatch;
#endif

	/* initialize */
	if (fillmode == kOverwrite) *this = 0;

	register nTYPE temp;
	register nTYPE bqT_Ij;
	
	nTYPE* pthisj =   Pointer();
	nTYPE* pqj    = q.Pointer();

	for (int j = 0; j < fCols; j++)
	{
		int istop = (range == kUpperOnly) ? j + 1 : fRows;
	
		nTYPE* pbI = b.Pointer();
		nTYPE* pqI = q.Pointer();
	
		for (int I = 0; I < b.fRows; I++)
		{			
			nTYPE* pbJ = pbI++;
			nTYPE* pqJ = pqj;
			
			bqT_Ij = 0.0;
			for (int J = 0; J < b.fCols; J++)
			{
				temp  = *pbJ;
				temp *= *pqJ;
				
				bqT_Ij += temp;
				
				pbJ += b.fRows;
				pqJ += q.fRows;
			}

			nTYPE* pqi    = pqI;
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
	  b.fCols != q.fRows) throw eSizeMismatch;
#endif

	/* initialize */
	if (fillmode == kOverwrite) *this = 0;

	register nTYPE temp;
	register nTYPE bq_iJ;
	
	nTYPE* pthisJ =   Pointer();
	nTYPE* pqJ    = q.Pointer();

	for (int J = 0; J < fCols; J++)
	{
		int Istop = (range == kUpperOnly) ? J + 1 : fRows;
	
		nTYPE* pqi = q.Pointer();
		nTYPE* pbi = b.Pointer();
	
		for (int i = 0; i < b.fRows; i++)
		{			
			nTYPE* pbj = pbi++;
			nTYPE* pqj = pqJ;
			
			bq_iJ = 0.0;
			for (int j = 0; j < b.fCols; j++)
			{
				temp  = *pbj;
				temp *= *pqj++;
				
				bq_iJ += temp;
				
				pbj += b.fRows;
			}

			nTYPE* pqI    = pqi++;
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

/* matrix-vector multiplication */

/* b_i = A_ij*x_j */
template <class nTYPE>
void nMatrixT<nTYPE>::Multx(const nArrayT<nTYPE>& x,
	nArrayT<nTYPE>& b) const
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (fRows != b.Length() || fCols != x.Length()) throw eSizeMismatch;
#endif

	nTYPE* ARow = Pointer();
	nTYPE* px0  = x.Pointer();
	nTYPE* pb   = b.Pointer();

	register nTYPE temp;
	register nTYPE sum;

	for (int i = 0; i < fRows; i++)
	{
		sum = 0.0;
		nTYPE *px = px0;
		nTYPE *AR = ARow;
	
		for (int j = 0; j < fCols; j++)
		{
			temp  = *px++;
			temp *= *AR;		
			 sum += temp;
			
			AR += fRows;	
		}

		*pb++ = sum;
		ARow++;
	}
}

/* b_i = A_ji*x_j */
template <class nTYPE>
void nMatrixT<nTYPE>::MultTx(const nArrayT<nTYPE>& x,
	nArrayT<nTYPE>& b) const
{
	/* dimension checks */
#if __option (extended_errorcheck)
	if (fRows != x.Length() && fCols != b.Length()) throw eSizeMismatch;
#endif

	nTYPE* ARow = Pointer();
	nTYPE* px0  = x.Pointer();
	nTYPE* pb   = b.Pointer();

	register nTYPE temp;
	register nTYPE sum;

	for (int i = 0; i < fCols; i++)
	{
		sum = 0.0;
		nTYPE *px = px0;
		nTYPE *AR = ARow;
	
		for (int j = 0; j < fRows; j++)
		{
			temp  = *px++;
			temp *= *AR++;
			sum  += temp;
			
		}

		*pb++ = sum;
		ARow += fRows;
	}
}

/* vector-matrix-vector product */
template <class nTYPE>
nTYPE nMatrixT<nTYPE>::MultmBn(const nArrayT<nTYPE>& m,
	const nArrayT<nTYPE>& n) const
{
	/* dimension checks */
#if __option (extended_errorcheck)
	if (fRows != m.Length() ||
	    fCols != n.Length()) throw eSizeMismatch;
#endif

	register nTYPE product = 0.0;
	register nTYPE mB;
	register nTYPE temp;
	
	nTYPE* pnj = n.Pointer();
	nTYPE* pBj = (*this).Pointer();
	
	for (int j = 0; j < fCols; j++)
	{
		mB = 0.0;
		nTYPE* pBi = pBj;
		nTYPE* pmi = m.Pointer();
	
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

/*
* Returns the outer product of the 2 vectors, or
* in dyadic notation:
*
*        v1 (x) v2
*
*/
template <class nTYPE>
void nMatrixT<nTYPE>::Outer(const nArrayT<nTYPE>& v1,
	const nArrayT<nTYPE>& v2)
{
	/* dimension checks */
#if __option (extended_errorcheck)
	if (v1.Length() != fRows || v2.Length() != fCols) throw eSizeMismatch;
#endif

	nTYPE* pthis = Pointer();
	nTYPE* pv1 = v1.Pointer();
	nTYPE* pv2 = v2.Pointer();

	for (int j = 0; j < fCols; j++)
	{
		nTYPE* pcol = pv1;

		for (int i = 0; i < fRows; i++)
		{
			*pthis    = *pcol++;
			*pthis++ *= *pv2;
		}
			
		pv2++;
	}
}

/* identity operations - square matrices ONLY */
template <class nTYPE>
void nMatrixT<nTYPE>::PlusIdentity(const nTYPE& value)
{
/* must be square */
#if __option (extended_errorcheck)
	if (fRows != fCols) throw eGeneralFail;
#endif

	if (fRows == 2)
	{
		fArray[0] += value;	
		fArray[3] += value;	
	}
	else if (fRows == 3)
	{
		fArray[0] += value;	
		fArray[4] += value;	
		fArray[8] += value;	
	}
	else
	{
		nTYPE* dex = Pointer();
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
	if (fRows != fCols) throw eGeneralFail;
#endif
	
	if (fRows == 2)
	{
		nTYPE* p = Pointer();
		*p++ = value;	
		*p++ = 0.0;	
		*p++ = 0.0;	
		*p   = value;
	}
	else if (fRows == 3)
	{
		nTYPE* p = Pointer();
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
		nTYPE* dex = Pointer();
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
void nMatrixT<nTYPE>::SetRow(int row, const nArrayT<nTYPE>& vec)
{
	/* dimension check */
#if __option (extended_errorcheck)
	if (vec.Length() != fCols) throw eSizeMismatch;
	if (row < 0 || row >= fRows) throw(eOutOfRange);
#endif

	nTYPE* pvec = vec.Pointer();
	nTYPE* pcol = Pointer() + row;
	
	for (int i = 0; i < fCols; i++)
	{
		*pcol = *pvec++;
		pcol += fRows;
	}
}

template <class nTYPE>
void nMatrixT<nTYPE>::SetRow(int row, const nTYPE& value)
{
	nTYPE* pcol = Pointer() + row;	
	for (int i = 0; i < fCols; i++)
	{
		*pcol = value;
		pcol += fRows;
	}
}

template <class nTYPE>
void nMatrixT<nTYPE>::SetCol(int col, const nArrayT<nTYPE>& vec)
{
/* dimension check */
#if __option (extended_errorcheck)
	if (vec.Length() != fRows) throw(eOutOfRange);
#endif

	nTYPE* pvec = vec.Pointer();
	nTYPE* pcol = (*this)(col);
	
	for (int i = 0; i < fRows; i++)
		*pcol++ = *pvec++;
}

template <class nTYPE>
void nMatrixT<nTYPE>::SetCol(int col, const nTYPE& value)
{
	nTYPE* pcol = (*this)(col);	
	for (int i = 0; i < fRows; i++)
		*pcol++ = value;
}

/* dot the specified row/column number with the array */
template <class nTYPE>
nTYPE nMatrixT<nTYPE>::DotRow(int rownum,
	const nArrayT<nTYPE>& vec) const
{
/* dimension check */
#if __option (extended_errorcheck)
	if (vec.Length() != fCols) throw eSizeMismatch;
#endif

	nTYPE *p    = &(*this)(rownum,0);
	nTYPE *pvec = vec.Pointer();

	register nTYPE sum = 0.0;
	register nTYPE temp;

	for (int i = 0; i < fCols; i++)
	{
		temp  = *p;
		temp *= *pvec++;
	
		sum  += temp;
		p    += fRows;
	}
		
	return(sum);
}

template <class nTYPE>
nTYPE nMatrixT<nTYPE>::DotCol(int colnum,
	const nArrayT<nTYPE>& vec) const
{
/* dimension check */
#if __option (extended_errorcheck)
	if (vec.Length() != fRows) throw eSizeMismatch;
#endif

	nTYPE *p    = (*this)(colnum);
	nTYPE *pvec = vec.Pointer();

	register nTYPE sum = 0.0;
	register nTYPE temp;

	for (int i = 0; i < fRows; i++)
	{
		temp  = *p++;
		temp *= *pvec++;
	
		sum  += temp;
	}
		
	return(sum);
}

#endif /* _NMATRIX_T_H_ */
