/* $Id: LAdMatrixT.cpp,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (12/05/1996)                                          */
/* Matrix2D with some linear algebra functions                            */

#include "LAdMatrixT.h"
#include <math.h>
#include "Constants.h"
#include "dArrayT.h"

/* constructor */
LAdMatrixT::LAdMatrixT(void) { }
LAdMatrixT::LAdMatrixT(int squaredim): dMatrixT(squaredim) { }
LAdMatrixT::LAdMatrixT(const LAdMatrixT& source): dMatrixT(source)
{ 	
	/* must be square */
	if (fRows != fCols) throw eGeneralFail;
}

/* pivoting functions */
void LAdMatrixT::RowPivot(int row1, int row2)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (row1 < 0 || row1 >= fRows ||
	    row2 < 0 || row2 >= fRows) throw eOutOfRange;
#endif

	double* p1 = (*this)(0) + row1;
	double* p2 = (*this)(0) + row2;
	
	for (int i = 0; i < fCols; i++)
	{
		double temp = *p1;
		*p1 = *p2;
		*p2 = temp;
		
		p1 += fRows;
		p2 += fRows;
	}
}

void LAdMatrixT::ColumnPivot(int col1, int col2)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (col1 < 0 || col1 >= fCols ||
	    col2 < 0 || col2 >= fCols) throw eOutOfRange;
#endif

	double* p1 = (*this)(col1);
	double* p2 = (*this)(col2);
	
	for (int i = 0; i < fRows; i++)
	{
		double temp = *p1;
		*p1++ = *p2;
		*p2++ = temp;
	}
}

/* square matrices only */
void LAdMatrixT::SymmetricPivot(int dex1, int dex2)
{
	RowPivot(dex1, dex2);
	ColumnPivot(dex1, dex2);
}

void LAdMatrixT::LinearSolve(dArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (RHS.Length() != fRows) throw eSizeMismatch;
#endif

	/* mean matrix value */
	double mean = Average();
			
	/* forward reduction */
	for (int col = 0; col < fCols-1; col++)
	{
		/* partial pivoting */
		int maxrow = col; //initialize to max on column	
		double currmax = fabs( (*this)(col,col) );
		double* pcol = &(*this)(col+1,col);
		for (int row = col + 1; row < fRows; row++)
			if ( fabs( *pcol++ ) > currmax )
			{
				currmax = fabs( *(pcol - 1) );
				maxrow = row;
			}
				
		if (maxrow != col)
		{
			/* swap row elements right of the diagonal */
			double* pfromrow = &(*this)(col,col);
			double* ptorow   = &(*this)(maxrow,col);
			for (int col2 = col; col2 < fCols; col2++)
			{
				double temp = *pfromrow;
				*pfromrow   = *ptorow;
				*ptorow     = temp;
			
				pfromrow += fRows;
				ptorow   += fRows;
			}
			
			/* swap RHS components */
			double temp = RHS[col];
			RHS[col]    = RHS[maxrow];
			RHS[maxrow] = temp;	
		}
			
		/* forward reduction */
		double diagvalue = (*this)(col,col);
		if (fabs( diagvalue/mean ) < kSmall) throw eGeneralFail;

		for (int row1 = col + 1; row1 < fRows; row1++)
		{
			double fact = (*this)(row1,col)/diagvalue;
			if (fabs(fact) > kSmall)
			{
				double* prow1 = &(*this)(row1,col+1);
				double* prow2 = &(*this)(col ,col+1);
				for (int col2 = col + 1; col2 < fCols; col2++)
				{
					*prow1 -= fact*(*prow2);
					
					prow1 += fRows;
					prow2 += fRows;
				}

				/* RHS */
				RHS[row1] -= fact*RHS[col];
			}
		}		
	}
	
	/* back substitution */
	if (fabs( (*this)(fRows-1,fCols-1)/mean ) < kSmall)
		throw eGeneralFail;

	RHS[fRows-1] /= (*this)(fRows-1,fCols-1); 		
	for (int row = fRows-2; row > -1; row--)
	{
		double sum = RHS[row];
		
		double* pRHS = &RHS[row+1];
		double* pcol = &(*this)(row,row+1);
		for (int col = row + 1; col < fCols; col++)
		{
			sum -= (*pcol)*(*pRHS++);
		
			pcol += fRows;
		}
			
		RHS[row] = sum/(*this)(row,row);
	}
}

void LAdMatrixT::LinearSolve2(dArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (RHS.Length() != fRows) throw eSizeMismatch;
#endif

	/* mean matrix value */
	double mean = Average();
			
	/* forward reduction */
	for (int col = 0; col < fCols-1; col++)
	{
		double diagvalue = (*this)(col,col);
		if (fabs( diagvalue/mean ) < kSmall) throw eGeneralFail;
		
		for (int row = col + 1; row < fRows; row++)
		{
			double fact = (*this)(row,col)/diagvalue;
			
			if (fabs(fact) > kSmall)
			{			
				double* prow1 = &(*this)(row,col+1);
				double* prow2 = &(*this)(col,col+1);
				
				for (int col2 = col + 1; col2 < fCols; col2++)
				{
					*prow1 -= fact*(*prow2);
					
					prow1 += fRows;
					prow2 += fRows;
				}

				/* RHS */
				RHS[row] -= fact*RHS[col];
			}
		}		
	}
	
	/* back substitution */
	if (fabs( (*this)(fRows-1,fCols-1)/mean ) > kSmall) throw eGeneralFail;

	RHS[fRows-1] /= (*this)(fRows-1,fCols-1); 		
	for (int row = fRows-2; row > -1; row--)
	{
		double sum = RHS[row];
		
		double* pRHS = &RHS[row+1];
		double* pcol = &(*this)(row,row+1);
		
		for (int col = row + 1; col < fCols; col++)
		{
			sum -= (*pcol)*(*pRHS++);
		
			pcol += fRows;
		}
			
		RHS[row] = sum/(*this)(row,row);
	}
}
