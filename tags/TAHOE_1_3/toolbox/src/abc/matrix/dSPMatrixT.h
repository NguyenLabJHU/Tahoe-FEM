/* $Id: dSPMatrixT.h,v 1.4 2003-11-21 22:41:36 paklein Exp $ */
/* created MLK 10/3/00 */
#ifndef _DSPMATRIX_T_H_
#define _DSPMATRIX_T_H_

/* direct members */
#include "nAutoFill2DT.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;
class dMatrixT;
#include "ios_fwd_decl.h"
//class dLACOOMatrixT;

/** sparse matrix class */
class dSPMatrixT
{
  public:

	/** \name constructors */
	/*@{*/
	dSPMatrixT(void);
	dSPMatrixT(int numrows, int numcols, int max_cols);
	dSPMatrixT(int squaredim, int max_cols);
	dSPMatrixT(const dSPMatrixT& source);
	/*@}*/
	
	/** dimensioning */
	void Dimension(int numrows, int numcols, int max_cols);
	
	/** clears matrix */
	void Clear(void);
	
	/** matrix element input (*this)(row, col) = element */
	void SetElement(int row, int col, double element);
	
	/** matrix element input (*this)(row, col) += element */
	void AddElement(int row, int col, double element);

	/** element accessor */
	double operator()(int row, int col) const;
	
	/** \name dimensions */
	/*@{*/
	int Rows(void) const;
	int Cols(void) const;
	int Nnz(void) const;
	int RowNnz(int row) const;
	/*@}*/
	
	/** \name operators */
	/*@{*/
	/** assignment operators */
	dSPMatrixT& operator=(const dSPMatrixT& RHS);
	dSPMatrixT& operator=(double value);

	/** multiplication by a scalar */
	dSPMatrixT& operator*=(double value);
	
	/** matrix addition */  	 
  	dSPMatrixT& operator+=(const dSPMatrixT& RHS);

	/** matrix subtraction */  	 
  	dSPMatrixT& operator-=(dSPMatrixT& RHS);
	/*@}*/

	/** matrix-diagonal matrix multiplication, M = M.D or M = D.M */
	/*@{*/
	void MultMD(const dArrayT& D);
	void MultDM(const dArrayT& D);
	/*@}*/
	
	/** \name adding blocks
	 * Add a block matrix (dMatrixT) to dSPMatrixT structure with an offset
	 * of rowstart rows and colstart cols in this */
	/*@{*/
	void SetBlock(int row, int col, const dMatrixT& block);
	void SetBlock(int row, int col, const dSPMatrixT& block);
	void AddBlock(int row, int col, const dMatrixT& block);
	void AddBlock(int row, int col, const dSPMatrixT& matrix);
	void AddBlockT(int row, int col, const dSPMatrixT& matrix);
	/*@}*/

	/** I/O operators */
	friend ostream& operator<<(ostream& out, const dSPMatrixT& matrix);
			
	/** transpose */
	void Transpose(const dSPMatrixT& RHS);
	
	/* matrix-vector multiplication - returns the result in b */
	void Multx(const dArrayT& x, dArrayT& b) const;
	void MultTx(const dArrayT& x, dArrayT& b) const;

	/** \name matrix-matrix multiplication */
	/*@{*/
	void MultAB(const dSPMatrixT& A, const dSPMatrixT& B);
	void MultABT(const dSPMatrixT& A, const dSPMatrixT& B);
	void MultATB(const dSPMatrixT& A, const dSPMatrixT& B);
	/*@}*/
	
	/** conversion to dLACOOMatrixT */
//	void Convert2dLACOOMatrixT(dLACOOMatrixT& dLACOOmatrix) const;
	
	/** copy row, copy column, copy column of M from transposed M */
	/*@{*/
	void CopyRow(int i, dArrayT& row) const;
	void CopyColumn(int j, dArrayT& col) const;
	void CopyColumnT(int j, dArrayT& col) const;
	/*@}*/
	
	/** this function is accessed by dCOOTensor3DT in the contraction with a matrix */
	int Width(void) const;
	
  private:
  	
  	int ColNum(int row, int j) const;
	double Value(int row, int j) const;
	
	int* ColPointer(int row);
	const int* ColPointer(int row) const;
	double* ValPointer(int row);
	const double* ValPointer(int row) const;

	//iAutoFill2DT& ColMatrix(void) const;
	//dAutoFill2DT& ValMatrix(void) const;

  protected:

	/** \name sparse matrix */
	/*@{*/
	AutoFill2DT<int>     fCol_Matrix; /**< row and column data */
	nAutoFill2DT<double> fVal_Matrix; /**< row and value data */
	/*@}*/
		
	/** \name dimensions */
	/*@{*/
	int	fRows; /**< number of rows */
	int	fCols; /**< number of columns */
	int fNnz;  /**< number of nonzero values */
	/*@}*/	
};

/******************** 
 * inline functions *
 ********************/

/* dimensions */
inline int dSPMatrixT::Rows(void) const { return fRows; }
inline int dSPMatrixT::Cols(void) const { return fCols; }
inline int dSPMatrixT::Nnz(void) const { return fNnz; }
inline int dSPMatrixT::RowNnz(int row) const { return fCol_Matrix.MinorDim(row); }

/* returns pointer to first element of column numbers matrix */
inline int* dSPMatrixT::ColPointer(int row) {
	return fCol_Matrix(row);
}
inline const int* dSPMatrixT::ColPointer(int row) const {
	return fCol_Matrix(row);
}

/* returns pointer to first element of values matrix */
inline double* dSPMatrixT::ValPointer(int row) {
	return fVal_Matrix(row);
}
inline const double* dSPMatrixT::ValPointer(int row) const {
	return fVal_Matrix(row);
}

/* matrix element accessor */
inline double dSPMatrixT::operator()(int row, int col) const
{
#if __option (extended_errorcheck)
	/* range checking */
	if (row < 0 || row >= fRows || col < 0 || col >= fCols)
		ExceptionT::OutOfRange();
#endif
	
	int num_cols_in_row = fCol_Matrix.MinorDim(row);
	
	/* search through column numbers in row */
	const int* pC = fCol_Matrix(row);				
	for( int i = 0; i < num_cols_in_row; i++)
	{
		int colnum = *pC++;
			
		/* gone past column number---return 0 */
		if( colnum > col )
			return(0.0);
		
		/* if element exists */
		if( colnum == col )
			return (fVal_Matrix(row, i));
	}
			
	/* if element doesn't exist */
	return 0.0;
}

/* matrix element input */
inline void dSPMatrixT::SetElement(int row, int col, double element)
{
#if __option (extended_errorcheck)
	/* range checking */
	if (row < 0 || row >= fRows || col < 0 || col >= fCols)
		ExceptionT::OutOfRange();
#endif
	
	if( element != 0.0 )
	{
		int num_cols_in_row = fCol_Matrix.MinorDim(row);
				
		/* if row is empty or col is greater than last column number in row */
		/* then append column and element values to end of row */
		if( num_cols_in_row == 0 || col > fCol_Matrix(row, num_cols_in_row - 1) )
		{
			fCol_Matrix.Append(row, col);
			fVal_Matrix.Append(row, element);
			fNnz++;
			return;
		}
		
		/* search through column numbers in row to determine where to insert or overwrite */
		int* pC = fCol_Matrix(row);
		for( int i = 0; i < num_cols_in_row; i++)
		{
			int colnum = *pC++;
			
			/* gone past column number---insert element at position i */
			if( colnum > col )
			{
				fCol_Matrix.Insert(row, col, i);
				fVal_Matrix.Insert(row, element, i);
				fNnz++;
				break;
			}
			
			/* if element already exists */
			if( colnum == col )
			{
				fVal_Matrix(row, i) = element;
				break;
			}	
		}
	}
}

/* matrix element input */
inline void dSPMatrixT::AddElement(int row, int col, double element)
{
#if __option (extended_errorcheck)
	/* range checking */
	if (row < 0 || row >= fRows || col < 0 || col >= fCols)
		ExceptionT::OutOfRange();
#endif
	
	if( element != 0.0 )
	{
		int num_cols_in_row = fCol_Matrix.MinorDim(row);
				
		/* if row is empty or col is greater than last column number in row */
		/* then append column and element values to end of row */
		if( num_cols_in_row == 0 || col > fCol_Matrix(row, num_cols_in_row - 1) )
		{
			fCol_Matrix.Append(row, col);
			fVal_Matrix.Append(row, element);
			fNnz++;
			return;
		}
		
		/* search through column numbers in row to determine where to insert or add */
		int* pC = fCol_Matrix(row);
		for( int i = 0; i < num_cols_in_row; i++)
		{
			int colnum = *pC++;
			
			/* gone past column number---insert element at position i */
			if( colnum > col )
			{
				fCol_Matrix.Insert(row, col, i);
				fVal_Matrix.Insert(row, element, i);
				fNnz++;
				break;
			}
			
			/* if element already exists */
			if( colnum == col )
			{
				fVal_Matrix(row, i) += element;
				break;
			}	
		}
	}
}

} /* namespace Tahoe */

#endif /* _DSPMATRIX_T_H_ */
