/* $Id: LAdMatrixT.cpp,v 1.8 2003-11-10 22:14:05 cjkimme Exp $ */
/* created: paklein (12/05/1996)                                          */
/* Matrix2D with some linear algebra functions                            */

#include "LAdMatrixT.h"
#include <math.h>
#include "toolboxConstants.h"
#include "dArrayT.h"

/* constructor */

using namespace Tahoe;

LAdMatrixT::LAdMatrixT(void) { }
LAdMatrixT::LAdMatrixT(int squaredim): dMatrixT(squaredim) { }
LAdMatrixT::LAdMatrixT(const LAdMatrixT& source): dMatrixT(source)
{ 	
	/* must be square */
	if (fRows != fCols) throw ExceptionT::kGeneralFail;
}

/* pivoting functions */
void LAdMatrixT::RowPivot(int row1, int row2)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (row1 < 0 || row1 >= fRows ||
	    row2 < 0 || row2 >= fRows) throw ExceptionT::kOutOfRange;
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
	    col2 < 0 || col2 >= fCols) throw ExceptionT::kOutOfRange;
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
	const char caller[] = "LAdMatrixT::LinearSolve";
/* dimension checks */
#if __option (extended_errorcheck)
	if (RHS.Length() != fRows) ExceptionT::SizeMismatch(caller);
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
		if (fabs(diagvalue/mean) < kSmall) 
			ExceptionT::BadJacobianDet(caller, "small pivot %g", diagvalue/mean);

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
	if (fabs((*this)(fRows-1,fCols-1)/mean) < kSmall)
		ExceptionT::BadJacobianDet(caller, "small pivot %g", (*this)(fRows-1,fCols-1)/mean);

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
	if (RHS.Length() != fRows) throw ExceptionT::kSizeMismatch;
#endif

	/* mean matrix value */
	double mean = Average();
			
	/* forward reduction */
	for (int col = 0; col < fCols-1; col++)
	{
		double diagvalue = (*this)(col,col);
		if (fabs( diagvalue/mean ) < kSmall) throw ExceptionT::kGeneralFail;
		
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
	if (fabs( (*this)(fRows-1,fCols-1)/mean ) > kSmall) throw ExceptionT::kGeneralFail;

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




int LAdMatrixT::BiCGStab(dArrayT& x ,const dArrayT& RHS ,const double M , int max_iter , double tol)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (RHS.Length() != fRows) throw ExceptionT::kSizeMismatch;
#endif


 	double resid;
 	double rho_1, rho_2 = 1., alpha, beta, omega;
    dArrayT p(fRows) , phat(fRows), s(fRows), shat(fRows), t(fRows), v(fRows);
	dArrayT r(fRows) , rtilde(fRows), t1(fRows), t2(fRows);

    double normb = RHS.Magnitude();
    
    //r = b - A * x;
    (*this).Multx(x,r);
    r *= -1.0;
    r += RHS;
    
    rtilde = r;
    
    if (normb == 0.0)
    normb = 1;
    
    if ((resid = r.Magnitude()/ normb) <= tol)
    {
    	tol = resid;
    	max_iter = 0;
    	return 0;
    }

  	for (int i = 1; i <= max_iter; i++) 
  	{
    	rho_1 = Dot(rtilde, r);
    	if (rho_1 == 0) 
    	{
    		tol = r.Magnitude() / normb;
      		return 2;
    	}

        if (i == 1)
      		p = r;
    	else 
    	{
      		beta = (rho_1/rho_2) * (alpha/omega);      		
 	     	//p = r + beta * (p - omega * v);
 	     	t1 = v;
	      	t1 *= -1.0*omega;
	      	t1 += p;
	      	t1 *= beta;
 	     	t1 += r;
	      	p = t1;
    	}
    	
    	//phat = M.solve(p);
    	//assume M := M*Identity
    	phat = p;
    	phat /= M;
   		//v = A * phat;
   		(*this).Multx(phat,v);
   		
    	alpha = rho_1 / Dot(rtilde, v);
    	
    	//s = r - alpha * v;
    	s.SetToCombination(1, r,-1.0*alpha,v);

    	
   		if ((resid = s.Magnitude()/normb) < tol) 
   		{
      		//x += alpha*phat;
      		t1 = phat;
      		t1 *= alpha;
      		x += t1;
      		
      		tol = resid;
      		return 0;
    	}

    	//shat = M.solve(s);
    	//assume M := M*identity
    	shat = s;
    	shat /=M;
     	(*this).Multx(shat,t);
   	 	omega = Dot(t,s) / Dot(t,t);
 
 
    	//x += alpha * phat + omega * shat;
    	x.AddCombination(alpha , phat ,omega , shat);
    	
    	//r = s - omega * t;
		r.SetToCombination(1.0, s,-1.0*omega,t);


    	rho_2 = rho_1;
    	if ((resid = r.Magnitude() / normb) < tol) 
    	{
      		tol = resid;
     		max_iter = i;
      		return 0;
    	}
   		if (omega == 0) 
   		{
      		tol = r.Magnitude() / normb;
      		return 3;
    	}
	}// end iteration loop
	
  tol = resid;
  return 1;

	
}


int LAdMatrixT::BiCGStab(dArrayT& x ,const dArrayT& RHS ,const dArrayT& M , int max_iter , double tol)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (RHS.Length() != fRows) throw ExceptionT::kSizeMismatch;
	if (  M.Length() != fRows) throw ExceptionT::kSizeMismatch;
#endif

	int j;
 	
 	double resid;
 	double rho_1, rho_2 = 1., alpha, beta, omega;
    dArrayT p(fRows) , phat(fRows), s(fRows), shat(fRows), t(fRows), v(fRows);
	dArrayT r(fRows) , rtilde(fRows), t1(fRows), t2(fRows);

    double normb = RHS.Magnitude();
    
    //r = b - A * x;
    (*this).Multx(x,r);
    r *= -1.0;
    r += RHS;
    
    rtilde = r;
    
    if (normb == 0.0)
    normb = 1;
    
    if ((resid = r.Magnitude()/ normb) <= tol)
    {
    	tol = resid;
    	max_iter = 0;
    	return 0;
    }

  	for (int i = 1; i <= max_iter; i++) 
  	{
    	rho_1 = Dot(rtilde, r);
    	if (rho_1 == 0) 
    	{
    		tol = r.Magnitude() / normb;
      		return 2;
    	}

        if (i == 1)
      		p = r;
    	else 
    	{
      		beta = (rho_1/rho_2) * (alpha/omega);      		
 	     	//p = r + beta * (p - omega * v);
 	     	t1 = v;
	      	t1 *= -1.0*omega;
	      	t1 += p;
	      	t1 *= beta;
 	     	t1 += r;
	      	p = t1;
    	}
    	
    	//phat = M.solve(p);
    	//assume M = diagonal
    	for (  j = 0; j <  M.Length(); j++) 
    		phat[j] = p[j]/M[j];
    		
   		//v = A * phat;
   		(*this).Multx(phat,v);
   		
    	alpha = rho_1 / Dot(rtilde, v);
    	
    	//s = r - alpha * v;
    	s.SetToCombination(1, r,-1.0*alpha,v);

    	
   		if ((resid = s.Magnitude()/normb) < tol) 
   		{
      		//x += alpha*phat;
      		t1 = phat;
      		t1 *= alpha;
      		x += t1;
      		
      		tol = resid;
      		return 0;
    	}

    	//shat = M.solve(s);
    	//assume M = diagonal
    	for ( j = 0; j <  M.Length(); j++) 
    		shat[j] = s[j]/M[j];
    	
    	
     	(*this).Multx(shat,t);
   	 	omega = Dot(t,s) / Dot(t,t);
 
 
    	//x += alpha * phat + omega * shat;
    	x.AddCombination(alpha , phat ,omega , shat);
    	
    	//r = s - omega * t;
		r.SetToCombination(1.0, s,-1.0*omega,t);


    	rho_2 = rho_1;
    	if ((resid = r.Magnitude() / normb) < tol) 
    	{
      		tol = resid;
     		max_iter = i;
      		return 0;
    	}
   		if (omega == 0) 
   		{
      		tol = r.Magnitude() / normb;
      		return 3;
    	}
	}// end iteration loop
	
  tol = resid;
  return 1;

	
}








