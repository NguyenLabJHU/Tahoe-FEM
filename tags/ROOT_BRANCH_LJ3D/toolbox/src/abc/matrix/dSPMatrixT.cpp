/* $Id: dSPMatrixT.cpp,v 1.8 2004-03-16 05:37:16 paklein Exp $ */
/* created MLK 10/3/00 */
#include "dSPMatrixT.h"

#include <iostream.h>
#include <iomanip.h>
#include <fstream.h> 

#include "dMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"

const int kheadroom = 25;

using namespace Tahoe;

/* return the partition size for the given number */
static int NumPartitions(int total) {
	int np = total/250;
	return (np < 1) ? 1 : np;
};

/* constructors */ 
dSPMatrixT::dSPMatrixT(void): 
	fRows(0), 
	fCols(0), 
	fNnz(0)
{ 
	
}

dSPMatrixT::dSPMatrixT(int numrows, int numcols, int max_cols): 
	fCol_Matrix(numrows, NumPartitions(numrows), kheadroom, max_cols),
	fVal_Matrix(numrows, NumPartitions(numrows), kheadroom, max_cols),
	fRows(numrows),
	fCols(numcols),
	fNnz(0)
{ 
	
}

dSPMatrixT::dSPMatrixT(int squaredim, int max_cols):
	fCol_Matrix(squaredim, NumPartitions(squaredim), kheadroom, max_cols),
	fVal_Matrix(squaredim, NumPartitions(squaredim), kheadroom, max_cols),
	fRows(squaredim),
	fCols(squaredim),
	fNnz(0)
{ 
	
}

/*dSPMatrixT::dSPMatrixT(const dSPMatrixT& source):
	fRows(source.Rows()),
	fCols(source.Cols()),
	fNnz(source.Nnz()),
	
	fCol_Matrix(source.ColMatrix()),
	fVal_Matrix(source.ValMatrix())
{

}*/

/* allocate */
void dSPMatrixT::Dimension(int numrows, int numcols, int max_cols)
{ 
	fRows = numrows;
	fCols = numcols;
	fNnz = 0;
	
	fCol_Matrix.SetHeadRoom(kheadroom);
	fCol_Matrix.SetSize(fRows, NumPartitions(fRows));
	fCol_Matrix.SetMaxMinorDim(max_cols);

	fVal_Matrix.SetHeadRoom(kheadroom);
	fVal_Matrix.SetSize(fRows, NumPartitions(fRows));
	fVal_Matrix.SetMaxMinorDim(max_cols);
}

/* clear matrix */
void dSPMatrixT::Clear(void)
{
	fCol_Matrix.Reset();
	fVal_Matrix.Reset();
	fNnz = 0;
}

namespace Tahoe {

/* I/O operators */
ostream& operator<<(ostream& out, const dSPMatrixT& matrix)
{
	out << "row\t" << "col\t" << "val\n\n"; 
	for( int i = 0; i < matrix.Rows(); i++ )
	{
		for( int j = 0; j < matrix.RowNnz(i); j++ )
		{
			out << i << "\t";
			out << matrix.ColNum(i,j) << "\t";
			out << matrix.Value(i,j) << "\n";
		}
	}
	out << endl;
	return (out);
}
}

dSPMatrixT& dSPMatrixT::operator=(const dSPMatrixT& RHS)
{
#if __option (extended_errorcheck)	
	/* dimension checks */
	if ( fRows != RHS.Rows() || fCols != RHS.Cols() )
		ExceptionT::SizeMismatch();
#endif

	SetBlock(0,0,RHS);
	return(*this);
}

dSPMatrixT& dSPMatrixT::operator=(double value)
{
	fVal_Matrix = value;
	return(*this);	
}

/* multiplication by a scalar */
dSPMatrixT& dSPMatrixT::operator*=(double value)
{
	fVal_Matrix *= value;
	return(*this);	
}

/* matrix addition M += RHS */
dSPMatrixT& dSPMatrixT::operator+=(const dSPMatrixT& RHS)
{
#if __option (extended_errorcheck)	
		/* dimension checks */
		if ( fRows != RHS.Rows() || fCols != RHS.Cols() )
			ExceptionT::SizeMismatch();
#endif

	int RHSnumrows = RHS.Rows();
	for( int i = 0; i < RHSnumrows; i++ )
	{
		const double* pV = RHS.ValPointer(i);
		const int* pC = RHS.ColPointer(i);
		int RHSnnzrowi = RHS.RowNnz(i);
		for( int j = 0; j < RHSnnzrowi; j++ )
			AddElement(i, *pC++, *pV++);
	}
	
	return(*this);
}

/* matrix subtraction M -= RHS */
dSPMatrixT& dSPMatrixT::operator-=(dSPMatrixT& RHS)
{
#if __option (extended_errorcheck)
	/* dimension checks */	
	if ( fRows != RHS.Rows() || fCols != RHS.Cols() )
		ExceptionT::SizeMismatch();
#endif
	
	RHS *= -1.0;
	(*this) += RHS;
	return(*this);
}

/* copies row i of matrix into vector "row" */
void dSPMatrixT::CopyRow(int i, dArrayT& row) const
{
#if __option (extended_errorcheck)	
	/* dimension checks */
	if (fCols != row.Length() || i < 0 || i >= fRows )
		ExceptionT::OutOfRange();
#endif
	
	row = 0.0;
	
	/* pointers to beginning of rows */
	const double* pV = fVal_Matrix(i);
	const int* pC = fCol_Matrix(i);
	
	/* record each nonzero element of the row */
	int numrownnz = RowNnz(i);
	for(int k = 0; k < numrownnz; k++)
		row[*pC++] = *pV++;
}

/* copies col j of matrix into vector "col" */
void dSPMatrixT::CopyColumn(int j, dArrayT& col) const
{
#if __option (extended_errorcheck)	
	/* dimension checks */
	if (fRows != col.Length() || j < 0 || j >= fCols )
		ExceptionT::OutOfRange();
#endif
	
	dSPMatrixT MT( fCols, fRows, Width() );
	MT.Transpose((*this));
	
	col = 0.0;
	
	/* pointers to beginning of rows of MT */
	double* pV = MT.ValPointer(j);
	int* pC = MT.ColPointer(j);
	
	/* record each nonzero element of the row of MT */
	int numrownnz = MT.RowNnz(j);
	for(int k = 0; k < numrownnz; k++)
		col[*pC++] = *pV++;
}

/* copies col j of matrix into vector "col", with input M^T */
void dSPMatrixT::CopyColumnT(int j, dArrayT& col) const
{
	CopyRow(j, col);
}

/* add a block matrix to dSPMatrixT structure */
void dSPMatrixT::SetBlock(int row, int col, const dMatrixT& block)
{
#if __option (extended_errorcheck)
	/* range checking */
	if (row + block.Rows() > fRows || col + block.Cols() > fCols)
		ExceptionT::OutOfRange();
#endif
	
	for (int i = 0; i < block.Rows(); i++)
	{
		for (int j = 0; j < block.Cols(); j++)
		{
			SetElement( i+row, j+col, block(i,j) );
		}
	}
}

void dSPMatrixT::SetBlock(int row, int col, const dSPMatrixT& block)
{
#if __option (extended_errorcheck)
	/* range checking */
	if (row + block.Rows() > fRows || col + block.Cols() > fCols)
		ExceptionT::OutOfRange();
#endif

	int blocknumrows = block.Rows();
	for( int i = 0; i < blocknumrows; i++ )
	{
		const double* pV = block.ValPointer(i);
		const int* pC = block.ColPointer(i);
		int blocknnzrowi = block.RowNnz(i);
		for( int j = 0; j < blocknnzrowi; j++ )
		{
			SetElement(i + row, (*pC++) + col, *pV++);
		}
	}
}

void dSPMatrixT::AddBlock(int row, int col, const dMatrixT& block)
{
#if __option (extended_errorcheck)
	/* range checking */
	if (row + block.Rows() > fRows || col + block.Cols() > fCols)
		ExceptionT::OutOfRange();
#endif
	
	for (int i = 0; i < block.Rows(); i++)
	{
		for (int j = 0; j < block.Cols(); j++)
		{
			AddElement( i+row, j+col, block(i,j) );
		}
	}
}

void dSPMatrixT::AddBlock(int row, int col, const dSPMatrixT& block)
{
#if __option (extended_errorcheck)
	/* range checking */
	if (row + block.Rows() > fRows || col + block.Cols() > fCols)
		ExceptionT::OutOfRange();
#endif

	int blocknumrows = block.Rows();
	for( int i = 0; i < blocknumrows; i++ )
	{
		const double* pV = block.ValPointer(i);
		const int* pC = block.ColPointer(i);
		int blocknnzrowi = block.RowNnz(i);
		for( int j = 0; j < blocknnzrowi; j++ )
		{
			AddElement(i + row, (*pC++) + col, *pV++);
		}
	}
}

void dSPMatrixT::AddBlockT(int row, int col, const dSPMatrixT& block)
{
	/* range checking */
#if __option (extended_errorcheck)
	/* range checking */
	if (row + block.Cols() > fRows || col + block.Rows() > fCols)
		ExceptionT::OutOfRange();
#endif

	int blocknumrows = block.Rows();
	for( int i = 0; i < blocknumrows; i++ )
	{
		const double* pV = block.ValPointer(i);
		const int* pC = block.ColPointer(i);
		int blocknnzrowi = block.RowNnz(i);
		for( int j = 0; j < blocknnzrowi; j++ )
		{
			AddElement((*pC++) + row, i + col, *pV++);
		}
	}
}

/* Pre multiplication by a diagonal matrix, D.  M = D.M
 * (*this) = D.(*this)	*/
void dSPMatrixT::MultDM(const dArrayT& D)
{	
#if __option (extended_errorcheck)
	/* dimension checking */
	if ( fRows != D.Length() || fCols != D.Length() )
		ExceptionT::SizeMismatch();
#endif
	
	// Mij = Dii Mij
	for( int i = 0; i < fRows; i++ )
	{
		double* pV = fVal_Matrix(i);
		for( int j = 0; j < fCol_Matrix.MinorDim(i); j++ )
		{
			fVal_Matrix(i, j) = (*pV++)*D[i];
		}
	}
}

/* Post multiplication by a diagonal matrix, D.  M = M.D
 * (*this) = (*this).D	*/
void dSPMatrixT::MultMD(const dArrayT& D)
{	
#if __option (extended_errorcheck)
	/* dimension checking */
	if ( fRows != D.Length() || fCols != D.Length() )
		ExceptionT::SizeMismatch();
#endif
			
	// Mij = Dii Mij
	for( int i = 0; i < fRows; i++ )
	{
		double* pV = fVal_Matrix(i);
		int* pC = fCol_Matrix(i);
		for( int j = 0; j < fCol_Matrix.MinorDim(i); j++ )
		{
			fVal_Matrix(i, j) = (*pV++)*D[*pC++];
		}
	}
}

/* A_ij*x_j = b_i	 		 A.x = b */
void dSPMatrixT::Multx(const dArrayT& x, dArrayT& b) const
{
#if __option (extended_errorcheck)	
	/* dimension checks */
	if (fRows != b.Length() || fCols != x.Length())
		ExceptionT::SizeMismatch();
#endif
	
	b = 0.0;
	
	/* for each nonzero element of the matrix, perform a multiplication
		with the appropriate element of x and place in the appropriate
		element of b */
	for(int i = 0; i < fRows; i++)
	{
		/* pointers to beginning of rows */
		const double* pV = fVal_Matrix(i);
		const int* pC = fCol_Matrix(i);
		for(int j = 0; j < fCol_Matrix.MinorDim(i); j++)
		{
			b[i] += (*pV++)*x[*pC++];
		}
	}
}

/*								 T	 */
/* A_ji*x_j = b_i			b = A .x  */
void dSPMatrixT::MultTx(const dArrayT& x, dArrayT& b) const
{	
#if __option (extended_errorcheck)
	/* dimension checks */
	if (fRows != x.Length() || fCols != b.Length())
		ExceptionT::SizeMismatch();
#endif

	b = 0.0;
	
	/* for each nonzero element of the matrix, perform a multiplication
		with the appropriate element of x and place in the appropriate
		element of b */
	for(int i = 0; i < fRows; i++)
	{
		/* pointers to beginning of rows */
		const double* pV = fVal_Matrix(i);
		const int* pC = fCol_Matrix(i);
		for(int j = 0; j < fCol_Matrix.MinorDim(i); j++)
		{
			b[*pC++] += (*pV++)*x[i];
		}
	}
}

/* Transpose of RHS */   /* this = RHS^T */
void dSPMatrixT::Transpose(const dSPMatrixT& RHS)
{
	int RHSnumrows = RHS.Rows();
	for( int i = 0; i < RHSnumrows; i++ )
	{
		const double* pV = RHS.ValPointer(i);
		const int* pC = RHS.ColPointer(i);
		int RHSnnzrowi = RHS.RowNnz(i);
		for( int j = 0; j < RHSnnzrowi; j++ )
		{
			// efficient since always appending to end of rows in new matrix
			SetElement(*pC++, i, *pV++);
			////SetElement(RHS.ColNum(i,j), i, RHS.Value(i,j));
		}
	}
	
	fRows = RHS.Cols();
	fCols = RHS.Rows();
	fNnz = RHS.Nnz();
}

/*                              	  T		*/
/* (*this) = A(i,k)*B(j,k) 		(A * B )	*/
void dSPMatrixT::MultABT(const dSPMatrixT& A, const dSPMatrixT& B)
{
#if __option (extended_errorcheck)	
	/* dimension checks */
	if (A.Rows() != fRows || B.Rows() != fCols || A.Cols() != B.Cols())
		ExceptionT::SizeMismatch();
#endif
	
	int numrows_A = A.Rows();
	int numrows_B = B.Rows();
	dArrayT arowvec(A.Cols());
	iArrayT Bnumrowvals(B.Rows());
	/* record number of nonzero entries in each column of B */
	for( int jj = 0; jj < numrows_B; jj++ )
		Bnumrowvals[jj] = B.RowNnz(jj);
	
	/* iterate over rows of A */
	for(int i = 0; i < numrows_A; i++)
	{	
		int Anumrowvals = A.RowNnz(i);
		
		/* record elements of row i of A in vector, arowvec */
		arowvec = 0.0;
		const int* pC_A = A.ColPointer(i);    // pointer to column numbers of row i of A
		const double* pV_A = A.ValPointer(i); // pointer to values of entries in row i of A
		for(int kk = 0; kk < Anumrowvals; kk++ )
			arowvec[*pC_A++] = *pV_A++;
	
		/* iterate over rows of B */
		for(int j = 0; j < numrows_B; j++)
		{	
			double sum = 0.0;
			const int* pC_B = B.ColPointer(j);    // pointer to column numbers of row j of B
			const double* pV_B = B.ValPointer(j); // pointer to values of entries in row j of B
			
			/* iterate over nonzero values of B in row j */
			for(int k = 0; k < Bnumrowvals[j]; k++)
			{
				////double aval = arowvec[*pR_B++];
				////if( aval != 0.0 )
				////	sum += (*pV_B)*aval;
				////pV_B++;
				sum += (*pV_B++)*arowvec[*pC_B++];
			}
			
			////if( sum != 0.0 )
				SetElement(i, j, sum);
		}
	}
}

/*                               	     T		    */
/* this(i,j) = A(k,i)*B(k,j) 	(this = A  * B)		*/
void dSPMatrixT::MultATB(const dSPMatrixT& A, const dSPMatrixT& B)
{
#if __option (extended_errorcheck)
	/* dimension checks */	
	if (A.Cols() != fRows || B.Cols() != fCols || A.Rows() != B.Rows())
		ExceptionT::SizeMismatch();
#endif
	
	/* take transpose of A */
	dSPMatrixT AT( A.Cols(), A.Rows(), A.Width() );
	AT.Transpose( A );

	MultAB(AT, B);
}

/* this(i,j) = A(i,k)*B(k,j) , (this = A * B) 	*/
void dSPMatrixT::MultAB(const dSPMatrixT& A, const dSPMatrixT& B)
{
/* dimension checks */
#if __option (extended_errorcheck)	
	if (A.Rows() != fRows || B.Cols() != fCols || A.Cols() != B.Rows())
		ExceptionT::SizeMismatch();		
#endif
	
	/* take transpose of B */
	dSPMatrixT BT( B.Cols(), B.Rows(), B.Width() );
	BT.Transpose( B );
	
	int numrows_A = A.Rows();
	int numcols_B = B.Cols();
	dArrayT arowvec(A.Cols());
	iArrayT Bnumcolvals(B.Cols());
	/* record number of nonzero entries in each column of B */
	for( int jj = 0; jj < numcols_B; jj++ )
		Bnumcolvals[jj] = BT.RowNnz(jj);
	
	/* iterate over rows of A */
	for(int i = 0; i < numrows_A; i++)
	{	
		int Anumrowvals = A.RowNnz(i);
		
		/* record elements of row i of A in vector, arowvec */
		arowvec = 0.0;
		const int* pC_A = A.ColPointer(i);    // pointer to column numbers of row i of A
		const double* pV_A = A.ValPointer(i); // pointer to values of entries in row i of A
		for(int kk = 0; kk < Anumrowvals; kk++ )
			arowvec[*pC_A++] = *pV_A++;
	
		/* iterate over columns of B (rows of B^T) */
		for(int j = 0; j < numcols_B; j++)
		{	
			double sum = 0.0;
			int* pR_B = BT.ColPointer(j);    // pointer to row numbers of column j of B
			double* pV_B = BT.ValPointer(j); // pointer to values of entries in column j of B
			
			/* iterate over nonzero values of B in column j */
			for(int k = 0; k < Bnumcolvals[j]; k++)
			{
				////double aval = arowvec[*pR_B++];
				////if( aval != 0.0 )
				////	sum += (*pV_B)*aval;
				////pV_B++;
				sum += (*pV_B++)*arowvec[*pR_B++];
			}
			
			////if( sum != 0.0 )
				SetElement(i, j, sum);
		}
	}
}

#if 0
/* conversion to dLACOOMatrixT */
void dSPMatrixT::Convert2dLACOOMatrixT(dLACOOMatrixT& dLACOOmatrix) const
{
	/* extract rows and columns */
	iArrayT row_array(fNnz);
	iArrayT col_array(fNnz);
	dArrayT val_array(fNnz);

	int* p_row_array = row_array.Pointer();
	int* p_col_array = col_array.Pointer();
	double* p_val_array = val_array.Pointer();
	
	for(int i = 0; i < fRows; i++)
	{
		/* pointers to beginning of rows */
		double* pV = fVal_Matrix(i);
		int* pC = fCol_Matrix(i);
		for(int j = 0; j < fCol_Matrix.MinorDim(i); j++)
		{
			*p_row_array++ = i;
			*p_col_array++ = *pC++;
			*p_val_array++ = *pV++;
		}
	}
	
	/* construct dLACOOMatrixT */
	dLACOOmatrix.Form(fRows, fCols, row_array, col_array, val_array, fNnz);
}
#endif

#if 0
/* this(i,j) = A(i,k)*B(k,j) , (this = A * B) 	*/
void dSPMatrixT::MultAB(dSPMatrixT& A, dSPMatrixT& B)
{
	/* dimension checks */
#if __option (extended_errorcheck)	
		if (A.Rows() != fRows || B.Cols() != fCols || A.Cols() != B.Rows()) throw(eGeneralFail);
#endif
	
	dSPMatrixT BT( B.Cols(), B.Rows(), B.Width() );
	BT.Transpose( B );

	int numrows_A = A.Rows();
	int numcols_B = B.Cols();
	dArrayT bcolvec(B.Rows());
	iArrayT Bnumcolvals(B.Cols());
	/* record number of nonzero entries in each column of B */
	for( int jj = 0; jj < numcols_B; jj++ )
	{
		Bnumcolvals[jj] = BT.RowNnz(jj);
	}
	
	/* iterate over rows of A */
	for(int i = 0; i < numrows_A; i++)
	{	
		int Anumrowvals = A.RowNnz(i);
		
		/* iterate over columns of B (rows of B^T) */
		for(int j = 0; j < numcols_B; j++)
		{	
			double sum = 0.0;
			int* pC_A = A.ColPointer(i);    // pointer to col nums of nonzero entries of row i in A
			double* pV_A = A.ValPointer(i); // pointer to values of nonzero entries of row i in A
			
			/* record elements of column j of B in vector, bcolvec */
			bcolvec = 0;
			int* pR_B = BT.ColPointer(j);    // pointer to row numbers of column j of B
			double* pV_B = BT.ValPointer(j); // pointer to values of entries in column j of B
			//int Bnumcolvals = BT.RowNnz(j);
			//for(int kk = 0; kk < Bnumcolvals; kk++ )
			for(int kk = 0; kk < Bnumcolvals[j]; kk++ )
			{
				bcolvec[*pR_B++] = *pV_B++;
			}
			
			/* iterate over nonzero values of A in row i */
			for(int k = 0; k < Anumrowvals; k++)
			{
				double bval = bcolvec[*pC_A++];
				if( bval != 0.0 )
					sum += (*pV_A)*bval;
				pV_A++;
				
				///////////sum += (*pV_A++)*bcolvec[*pC_A++];
			}
			
			if( sum != 0.0 )
				SetElement(i, j, sum);
		}
	}
}
#endif

#if 0
/* this(i,j) = A(i,k)*B(k,j) , (this = A * B) 	*/
/*void dSPMatrixT::MultAB(dSPMatrixT& A, dSPMatrixT& B) */
{
	/* dimension checks */
#if __option (extended_errorcheck)	
		if (A.Rows() != fRows || B.Cols() != fCols || A.Cols() != B.Rows()) throw(eGeneralFail);
#endif
	
	dSPMatrixT BT( B.Cols(), B.Rows(), B.Width() );
	BT.Transpose( B );
	
	/* iterate over rows of A */
	for(int i = 0; i < A.Rows(); i++)
	{
		/* iterate over columns of B (rows of B^T) */
		for(int j = 0; j < B.Cols(); j++)
		{
			int maxrow_B;  // maximum row number of column j of B (max col num of row j of BT)
			double sum = 0.0;
			if( BT.RowNnz(j) == 0 )
				maxrow_B = -1;
			else
				maxrow_B = BT.ColNum( j, BT.RowNnz(j)-1 );
			
			int* pC_A = A.ColPointer(i);    // pointer to col nums of nonzero entries of row i in A
			double* pV_A = A.ValPointer(i); // pointer to values of nonzero entries of row i in A
			
			/* iterate over nonzero values of A in row i */
			int startrow = 0;
			for(int k = 0; k < A.RowNnz(i); k++)
			{
				
				/* if the current column number of A is greater than the maximum row number
					of the B column, then skip to next column of B */
				if( *pC_A > maxrow_B )
					break;
				
				/* iterate over values of rows for nonzero entries in column j of B 
					to see if any row number matches column number of i row of A */
				int* pR_B = BT.ColPointer(j);    // pointer to row numbers of column j of B
				double* pV_B = BT.ValPointer(j); // pointer to values of entries in column j of B
				pR_B += startrow;
				pV_B += startrow;
				for( int kk = startrow; kk < BT.RowNnz(j); kk++)
				{
					/* if the column number of A is less than the current row number of B,
						then skip to next column number in row i of A */
					if( *pC_A < *pR_B )
						break;
						
					/* if column number of A is greater than current row number of B,
						then set startrow to kk+1 */
					if( *pC_A > *pR_B )
						startrow = kk+1;
					
					/* if column number of A matches row number of B, augment sum */
					if( *pC_A == *pR_B )
					{
						sum += (*pV_A)*(*pV_B);
						startrow = kk+1;
						break;
					}
				
					pR_B++;
					pV_B++;
				}
	
				pC_A++;
				pV_A++;
			}
			
			if( sum != 0.0 )
				SetElement(i, j, sum);
		}
	}
}
#endif

/***********/
/* private */
/***********/

/* returns column number of jth element in row i */
int dSPMatrixT::ColNum(int row, int j) const
{
	return (fCol_Matrix(row,j));
}

/* returns matrix entry value of jth element in row i */
double dSPMatrixT::Value(int row, int j) const
{
	return (fVal_Matrix(row,j));
}

/* returns width of storage matrices */
int dSPMatrixT::Width(void) const
{
	return ( fCol_Matrix.MaxMinorDim() );
}

#if 0
/* returns fCol_Matrix */
iAutoFill2DT& dSPMatrixT::ColMatrix(void) const
{
	return( fCol_Matrix );
}

/* returns fVal_Matrix */
dAutoFill2DT& dSPMatrixT::ValMatrix(void) const
{
	return( fVal_Matrix );
}
#endif
