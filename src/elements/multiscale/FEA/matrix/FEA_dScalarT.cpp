// $Id: FEA_dScalarT.cpp,v 1.3 2003-02-03 04:40:24 paklein Exp $
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
  int i,l, n_rows_x_n_cols = A.Rows()*A.Cols(), n_ip = A.IPs(); 
  //register double dot = 0.0;
  double dot = 0.0;
	double *pA = A[0].Pointer ();
	double *pB = B[0].Pointer ();

	for (l=0; l<n_ip; l++) {
 	  for (i=0; i<n_rows_x_n_cols; i++) 
     	dot += (*pA++)*(*pB++); 
		(*this)[l] = dot;
		dot = 0.0;
	}

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

