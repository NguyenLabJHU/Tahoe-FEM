/* $Id: BandedLAdMatrixT.cpp,v 1.6 2005-07-29 03:09:33 paklein Exp $ */
/* created: MLK (05/21/1997)                                              */
/* square banded matrix operations                                        */
/* banded matrix elements stored in columns                               */
#include "BandedLAdMatrixT.h"
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include "toolboxConstants.h"
#include "dArrayT.h"
#include "dMatrixT.h"

using namespace Tahoe;
const char caller[] = "BandedLAdMatrixT";

/* constructor */
BandedLAdMatrixT::BandedLAdMatrixT(int squaredim, int leftbandsize, int rightbandsize):
	nArrayT<double>((leftbandsize+rightbandsize+1)*squaredim),
	fRows(squaredim),
	fCols(squaredim),
	fLband(leftbandsize),
	fRband(rightbandsize),
	fColumnHeight(fLband + fRband + 1) //plus the diagonal element
{

}

/* element accessor */
double BandedLAdMatrixT::GetElement(int row, int col) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (row < 0 || row >= fRows || col < 0 || col >= fCols)
		ExceptionT::OutOfRange(caller);
#endif
	
	if(row-col > fLband || col-row > fRband)
		return 0.0;
	else
		return (fArray[col*(fLband + fRband + 1) + fRband + row - col]);
};

namespace Tahoe {

/* I/O operators */
ostream& operator<<(ostream& out, const BandedLAdMatrixT& matrix)
{
	int d_width = OutputWidth(out, matrix.Pointer());

	for (int j = 0; j < matrix.fRows; j++)
	{
		for (int i = 0; i < matrix.fCols; i++)
				out << setw(d_width) << matrix.GetElement(j,i);
		
		out << '\n';
	}
	
	return (out);
}

} // namespace Tahoe

using namespace Tahoe;

/* assemble beginning with row and col in the upper left. */
void BandedLAdMatrixT::AddBlock(int row, int col, const dMatrixT& block)
{
/* range check */
#if __option(extended_errorcheck)
	/* within bounds */
	if (row < 0 || row + block.Rows() > fRows ||
	    col < 0 || row + block.Cols() > fCols) ExceptionT::OutOfRange(caller);

	/* within the band */
	if ((row + block.Rows() - col) > fLband + 1 ||
	    (col + block.Cols() - row) > fRband + 1 ) ExceptionT::OutOfRange(caller);
#endif

	const double* pblock = block.Pointer();
	double* pstart = fArray + col*fColumnHeight + fRband + row - col;
	
	for (int i = 0; i < block.Cols(); i++)
	{
		double* pcol = pstart;
		for (int j = 0; j < block.Rows(); j++)
			*pcol++ += *pblock++;
			
		pstart += (fColumnHeight - 1);
	}
}

/* transpose copy */
void BandedLAdMatrixT::Transpose(const BandedLAdMatrixT& matrix)
{
#if __option(extended_errorcheck)
	/* dimension checks */
	if( fRband != matrix.Lband() || fLband != matrix.Rband() )
		ExceptionT::SizeMismatch(caller);
#endif

	int mincol, maxcol;
	for(int row = 0; row < fRows; row++) {
		mincol = row - fLband;
		if(mincol < 0)
			mincol = 0;
			
		maxcol = row + fRband + 1;
		if(maxcol > fCols)
			maxcol = fCols;
		
		for(int col = mincol; col < maxcol; col++) {
			(*this)(row,col) = matrix(col,row);
		}
	}
}

/* Gaussian elimination with the given RHS vector or RHS matrix */
void BandedLAdMatrixT::LinearSolve(dArrayT& RHS)
{
	/* compute mean value of elements contained in bands */
	double* pA = (*this)(0);
	double  mean   = 0.0;
	double  length = (fLband + fRband + 1)*fCols;
	
	for (int i = 0; i < length; i++)
		mean += *pA++;
	mean /= length;	
			
	/* forward reduction */
	for (int col = 0; col < fCols-1; col++)
	{
		double diagvalue = (*this)(col,col);
		if(fabs( diagvalue/mean ) < 1.0e-12)
			ExceptionT::GeneralFail(caller);
		
		int maxrow = fRows-1;
		if(col + fLband < maxrow)
			maxrow = col + fLband;
		
		for (int row = col + 1; row <= maxrow; row++)
		{
			double fact = (*this)(row,col)/diagvalue;
			
			if (fabs(fact) > 1.0e-12)
			{
				double* prow1 = &(*this)(row,col+1);
				double* prow2 = &(*this)(col,col+1);
				
				int maxcol = fCols-1;
				if(col + fRband < maxcol)
					maxcol = col + fRband;
				
				for (int col2 = col + 1; col2 <= maxcol; col2++)
				{
					*prow1 -= fact*(*prow2);
					
					prow1 += fLband + fRband;
					prow2 += fLband + fRband;
				}

				/* RHS */
				RHS[row] -= fact*RHS[col];
			}
		}		
	}
	
	/* back substitution */
	if(fabs( (*this)(fRows-1,fCols-1)/mean ) < 1.0e-12)
		ExceptionT::GeneralFail(caller);

	RHS[fRows-1] /= (*this)(fRows-1,fCols-1); 	
			
	for (int row = fRows-2; row >= 0; row--)
	{
		double sum = RHS[row];
		
		double* pRHS = &RHS[row+1];
		double* pcol = &(*this)(row,row+1);
		
		int maxcol = fCols-1;
		if(row + fRband < maxcol)
			maxcol = row + fRband;
			
		for (int col = row + 1; col <= maxcol; col++)
		{
			sum -= (*pcol)*(*pRHS++);
		
			pcol += fLband + fRband;
		}
			
		RHS[row] = sum/(*this)(row,row);
	}
}

/*
*	LinearSolve
*
*	for Matrix RHS
*/
void BandedLAdMatrixT::LinearSolve(dMatrixT& RHS)
{
	/* dimension checks */
	if(fRows != fCols || RHS.Rows() != fRows)
		ExceptionT::SizeMismatch(caller);

	/* compute mean value of elements contained in bands */
	double* pA = (*this)(0);
	double  mean   = 0.0;
	double  length = (fLband + fRband + 1)*fCols;
	
	for (int i = 0; i < length; i++)
		mean += *pA++;
	mean /= length;	
			
	/* forward reduction */
	for (int col = 0; col < fCols-1; col++)
	{
		double diagvalue = (*this)(col,col);
		if(fabs( diagvalue/mean ) < 1.0e-12)
			ExceptionT::GeneralFail(caller);
		
		int maxrow = fRows-1;
		if(col + fLband < maxrow)
			maxrow = col + fLband;
		
		for (int row = col + 1; row <= maxrow; row++)
		{
			double fact = (*this)(row,col)/diagvalue;
			
			if (fabs(fact) > 1.0e-12)
			{
				double* prow1 = &(*this)(row,col+1);
				double* prow2 = &(*this)(col,col+1);
				
				int maxcol = fCols-1;
				if(col + fRband < maxcol)
					maxcol = col + fRband;
				
				for (int col2 = col + 1; col2 <= maxcol; col2++)
				{
					*prow1 -= fact*(*prow2);
					
					prow1 += fLband + fRband;
					prow2 += fLband + fRband;
				}

				/* RHS */
				for (int colRHS = 0; colRHS < RHS.Cols(); colRHS++)
				{
					RHS(row,colRHS) -= fact*RHS(col,colRHS);
				}
				
			}
		}		
	}
	
	/* back substitution */
	if(fabs( (*this)(fRows-1,fCols-1)/mean ) < 1.0e-12)
		ExceptionT::GeneralFail(caller);

	for (int k = 0; k < RHS.Cols(); k++)
	{
		RHS(fRows-1,k) /= (*this)(fRows-1,fCols-1);
	}	
	
	for (int colRHS = 0; colRHS < RHS.Cols(); colRHS++)
	{
			
		for (int row = fRows-2; row >= 0; row--)
		{
			double sum = RHS(row,colRHS);
			
			double* pRHS = &RHS(row+1,colRHS);
			double* pcol = &(*this)(row,row+1);
			
			int maxcol = fCols-1;
			if(row + fRband < maxcol)
				maxcol = row + fRband;
				
			for (int col = row + 1; col <= maxcol; col++)
			{
				sum -= (*pcol)*(*pRHS++);
			
				pcol += fLband + fRband;
			}
				
			RHS(row,colRHS) = sum/(*this)(row,row);
		}
	}
}


/*
*	BandedInverse
*
*/
void BandedLAdMatrixT::BandedInverse(dMatrixT& RHS)
{
	/* dimension checks */
	if(fRows != fCols || RHS.Rows() != fRows)
		ExceptionT::SizeMismatch(caller);

	/* compute mean value of elements contained in bands */
	double* pA = (*this)(0);
	double  mean   = 0.0;
	double  length = (fLband + fRband + 1)*fCols;
	
	for (int i = 0; i < length; i++)
		mean += *pA++;
	mean /= length;	
	
	RHS = 0.0;
	RHS.PlusIdentity();
			
	/* forward reduction */
	for (int col = 0; col < fCols-1; col++)
	{
		double diagvalue = (*this)(col,col);
		if(fabs( diagvalue/mean ) < 1.0e-12)
			ExceptionT::GeneralFail(caller);
		
		int maxrow = fRows-1;
		if(col + fLband < maxrow)
			maxrow = col + fLband;
		
		for (int row = col + 1; row <= maxrow; row++)
		{
			double fact = (*this)(row,col)/diagvalue;
			
			if (fabs(fact) > 1.0e-12)
			{
				double* prow1 = &(*this)(row,col+1);
				double* prow2 = &(*this)(col,col+1);
				
				int maxcol = fCols-1;
				if(col + fRband < maxcol)
					maxcol = col + fRband;
				
				for (int col2 = col + 1; col2 <= maxcol; col2++)
				{
					*prow1 -= fact*(*prow2);
					
					prow1 += fLband + fRband;
					prow2 += fLband + fRband;
				}

				/* RHS */
				for (int colRHS = 0; colRHS < RHS.Cols(); colRHS++)
				{
					RHS(row,colRHS) -= fact*RHS(col,colRHS);
				}
				
			}
		}		
	}
	
	/* back substitution */
	if(fabs( (*this)(fRows-1,fCols-1)/mean ) < 1.0e-12)
		ExceptionT::GeneralFail(caller);

	for (int colRHS1 = 0; colRHS1 < RHS.Cols(); colRHS1++)
	{
		RHS(fRows-1,colRHS1) /= (*this)(fRows-1,fCols-1);
	}	
	
	for (int colRHS = 0; colRHS < RHS.Cols(); colRHS++)
	{
			
		for (int row = fRows-2; row >= 0; row--)
		{
			double sum = RHS(row,colRHS);
			
			double* pRHS = &RHS(row+1,colRHS);
			double* pcol = &(*this)(row,row+1);
			
			int maxcol = fCols-1;
			if(row + fRband < maxcol)
				maxcol = row + fRband;
				
			for (int col = row + 1; col <= maxcol; col++)
			{
				sum -= (*pcol)*(*pRHS++);
			
				pcol += fLband + fRband;
			}
				
			RHS(row,colRHS) = sum/(*this)(row,row);
		}
	}
}



/* b_i = A_ij*x_j */
void BandedLAdMatrixT::Multx(const dArrayT& x,
	dArrayT& b) const
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (fRows != b.Length() || fCols != x.Length()) ExceptionT::SizeMismatch(caller);
#endif

	//double* ARow = Pointer() + fRband;
	//double* ARow;
	const double* px0  = x.Pointer();
	double* pb   = b.Pointer();

	register double temp;
	register double sum;

	for (int i = 0; i < fRows; i++)
	{
		sum = 0.0;
		
		// start column and stop column
		int startcol, stopcol;
		
		if(i >= fLband)
			startcol = i - fLband;
		else
			startcol = 0;
		
		if(i+fRband < fCols)
			stopcol = i + fRband;
		else
			stopcol = fCols-1;
		
		// set pointer ARow to first element of i row of A
		//ARow = &(*this)(i,startcol);
		const double *AR = &(*this)(i,startcol);

		// set pointer to x vector to appropriate row
		const double *px = px0 + startcol;
		
		//double *AR = ARow;
		
		for (int j = startcol; j <= stopcol; j++)
		{
			temp  = *px++;
			temp *= *AR;		
			 sum += temp;
			
			AR += fRband + fLband;	
		}

		*pb++ = sum;
		
	}
}



/* B_ij = this_ik*M_kj */
void BandedLAdMatrixT::MultM(const dMatrixT& M,
	dMatrixT& B) const
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if ( ( fRows != B.Rows() || B.Cols() != M.Cols() ) || fCols != M.Rows() )
		ExceptionT::SizeMismatch(caller);
#endif

	for (int col = 0; col < B.Cols(); col++) {

		const double* pM0  = M.Pointer() + col*(M.Rows());
		double* pB   = B.Pointer() + col*(B.Rows());

		register double temp;
		register double sum;

		for (int i = 0; i < fRows; i++)
		{
			sum = 0.0;
			
			// start column and stop column
			int startcol, stopcol;
			
			if(i >= fLband)
				startcol = i - fLband;
			else
				startcol = 0;
			
			if(i+fRband < fCols)
				stopcol = i + fRband;
			else
				stopcol = fCols-1;
			
			// set pointer ARow to first element of i row of A
			//ARow = &(*this)(i,startcol);
			const double *AR = &(*this)(i,startcol);	

			// set pointer to M matrix to appropriate row
			const double *pM = pM0 + startcol;
		
			//double *AR = ARow;
		
			for (int j = startcol; j <= stopcol; j++)
			{
				temp  = *pM++;
				temp *= *AR;		
				 sum += temp;
				
				AR += fRband + fLband;	
			}

			*pB++ = sum;
		}
	}
}



/* b_i = A_ji*x_j */
void BandedLAdMatrixT::MultTx(const dArrayT& x,
	dArrayT& b) const
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (fRows != x.Length() || fCols != b.Length()) ExceptionT::SizeMismatch(caller);
#endif

	const double* px0  = x.Pointer();
	double* pb   = b.Pointer();

	register double temp;
	register double sum;

	for (int i = 0; i < fCols; i++)
	{
		sum = 0.0;
		
		// start row and stop row
		int startrow, stoprow;
		
		if(i >= fRband)
			startrow = i - fRband;
		else
			startrow = 0;
		
		if(i+fLband < fRows)
			stoprow = i + fLband;
		else
			stoprow = fRows-1;
		
		// set pointer ARow to first element of i row of A
		//ARow = &(*this)(i,startcol);
		const double *AR = &(*this)(startrow,i);

		// set pointer to x vector to appropriate row
		const double *px = px0 + startrow;
		
		//double *AR = ARow;
		
		for (int j = startrow; j <= stoprow; j++)
		{
			temp  = *px++;
			temp *= *AR;		
			 sum += temp;
			
			AR ++;	
		}

		*pb++ = sum;
		
	}
}




/*                             T    T	*/
/* B_ij = this_ki*M_jk  = (this  * M )  */
void BandedLAdMatrixT::MultTM(const dMatrixT& M,
	dMatrixT& B) const
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if ( ( fRows != B.Cols() || B.Rows() != M.Cols() ) || fCols != M.Rows() )
		ExceptionT::SizeMismatch(caller);
#endif

	for (int row = 0; row < B.Rows(); row++) {

		const double* pM0  = M.Pointer() + row*(M.Rows());
		double* pB   = B.Pointer() + row;

		register double temp;
		register double sum;

		for (int i = 0; i < fRows; i++)
		{
			sum = 0.0;
			
			// start column and stop column
			int startcol, stopcol;
			
			if(i >= fLband)
				startcol = i - fLband;
			else
				startcol = 0;
			
			if(i+fRband < fCols)
				stopcol = i + fRband;
			else
				stopcol = fCols-1;
			
			// set pointer ARow to first element of i row of A
			//ARow = &(*this)(i,startcol);
			const double *AR = &(*this)(i,startcol);	

			// set pointer to M matrix to appropriate row
			const double *pM = pM0 + startcol;
		
			//double *AR = ARow;
		
			for (int j = startcol; j <= stopcol; j++)
			{
				temp  = *pM++;
				temp *= *AR;		
				 sum += temp;
				
				AR += fRband + fLband;	
			}
			*pB = sum;
			pB += (B.Rows());
		}
	}
}
