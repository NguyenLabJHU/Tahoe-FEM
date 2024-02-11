// $Id: FEA_dScalarT.cpp,v 1.4 2003-04-23 23:34:22 creigh Exp $
#include "FEA.h"

using namespace Tahoe; 

//----------------------------------------------------

void FEA_dScalarT::Dot ( FEA_dVectorT &a, FEA_dVectorT &b ) 
{
  int i,l, n_rows = a.Rows(), n_ip = a.IPs(); 
  register double dot = 0.0;
	double *pa = a[0].Pointer ();
	double *pb = b[0].Pointer ();

	for (l=0; l<n_ip; l++) {
 	  for (i=0; i<n_rows; i++) 
     	dot += (*pa++)*(*pb++); 
		(*this)[l] = dot;
		dot = 0.0;
	}

}

//----------------------------------------------------

void FEA_dScalarT::Double_Dot ( FEA_dMatrixT &A, FEA_dMatrixT &B ) 
{
  int i,l, n_rows_x_n_cols = A.Rows()*A.Cols(); 
  //register double dot = 0.0;
  double dot = 0.0;
	double *pA = A[0].Pointer ();
	double *pB = B[0].Pointer ();

	for (l=0; l<fLength; l++) {
 	  for (i=0; i<n_rows_x_n_cols; i++) 
     	dot += (*pA++)*(*pB++); 
		(*this)[l] = dot;
		dot = 0.0;
	}

}

//----------------------------------------------------

void FEA_dScalarT::L1_Norm (FEA_dMatrixT &A) 
{
	for (int i=0; i<fLength; i++) (*this)[i] = A[i].AbsSum (); 
}

//----------------------------------------------------

void FEA_dScalarT::L1_Diagonal ( FEA_dMatrixT &A )  //--- Get sum of abs ( diagonal terms )
{
  int i,l, n_rows = A.Rows(); 
  (*this) = 0.0;

	for (l=0; l<fLength; l++) {
		double *pA = A[l].Pointer ();
 	  for (i=0; i<n_rows; i++) {
     	(*this)[l] += fabs((*pA)); 
			pA += (n_rows+1);
		}
	}
}

//----------------------------------------------------

void FEA_dScalarT::Diagonal_Bias ( FEA_dMatrixT &A ) 
{
	FEA_dScalarT L1_A(A.n_ip);

	(*this).L1_Diagonal ( A );			
	L1_A.L1_Norm ( A );
	(*this) /= 	L1_A;
}

//----------------------------------------------------

void FEA_dScalarT::Max (FEA_dVectorT &a) 
{
	for (int i=0; i<fLength; i++) (*this)[i] = a[i].Max();
}

//----------------------------------------------------

void FEA_dScalarT::Max (FEA_dMatrixT &A) 
{
	for (int i=0; i<fLength; i++) (*this)[i] = A[i].Max();
}

//----------------------------------------------------

void FEA_dScalarT::AbsMax (FEA_dVectorT &a) 
{
	for (int i=0; i<fLength; i++) (*this)[i] = a[i].AbsMax();
}

//----------------------------------------------------

void FEA_dScalarT::AbsMax (FEA_dMatrixT &A) 
{
	for (int i=0; i<fLength; i++) (*this)[i] = A[i].AbsMax();
}


//----------------------------------------------------

void FEA_dScalarT::Determinant (FEA_dMatrixT &A) 
{
	for (int i=0; i<fLength; i++) (*this)[i] = A[i].Det();
}

//----------------------------------------------------

void FEA_dScalarT::Trace (FEA_dMatrixT &A) 
{
	for (int i=0; i<fLength; i++) (*this)[i] = A[i].Trace();
}

//----------------------------------------------------

