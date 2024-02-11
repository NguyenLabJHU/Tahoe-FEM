// $Id: FEA_dVectorT.cpp,v 1.6 2003-11-21 22:54:44 paklein Exp $
#include "FEA.h"

using namespace Tahoe; 

//------------- Constructors -------------------------

FEA_dVectorT::FEA_dVectorT(void) : 
ArrayT<dArrayT>(), 
ip_components(0) 
{
	n_ip = 0;
	n_sd = 0;
}

//----------------------------------------------------

FEA_dVectorT::FEA_dVectorT(int ips) : 
ArrayT<dArrayT>(ips),
ip_components(ips)
{
	n_ip = ips;
	n_sd = 0;
}

//----------------------------------------------------

FEA_dVectorT::FEA_dVectorT(int ips,int length) :
ip_components(ips) // Allocates, sets length
{
	FEA_Dimension(ips,length);
}

//------------ Utilities -----------------------------

void FEA_dVectorT::FEA_Dimension (int ips,int length) 
{

  if (ip_components.length==0)  ip_components.Allocate(ips);
  n_ip   = ips;  // this also is/equals fLength
	n_sd   = length;

	Block_Memory.Allocate(n_ip*length); // Dont need to do this way but just keep
  Allocate(n_ip); // Allocate fArray
	for (int i=0; i<n_ip; i++)
		fArray[i].Set(length, Block_Memory.Pointer(i*length));

}		

//----------------------------------------------------

void FEA_dVectorT::FEA_Dimension (const FEA_dVectorT &a) // Copy Dimension of another
{
	if (a.fLength == 0)
    fLength = 0; 
	else
		FEA_Dimension(a.Length(), a[0].Length());
}

//----------------------------------------------------

void FEA_dVectorT::Print(char *c) { // overload << later

  if (fLength==0)
		cout << "...ERROR >> FEA_dVectorT::Print() : "<<c<<" Unallocated \n\n";

  cout << "\n Vector "<<c<<" evaluated at "<<fLength<<" inegration points (ip): \n"; 
  for (int i=0; i<fLength; i++) 
    cout <<"\n "<< c <<" @ ip "<<i<<": \n\n"<< (*this)[i] << "\n";

}

//------------ Vector Operations ----------------------

//----------------------------------------------------

void FEA_dVectorT::Magnitude (FEA_dScalarT &s)
{
	for (int l=0; l<n_ip; l++) 
		s[l] = (*this)[l].Magnitude(); 
}

 
void FEA_dVectorT::Dot (const FEA_dVectorT& b, FEA_dScalarT& c) 
{
  if (fLength==0 || b.Length()==0)
		cout << "...ERROR >> FEA_dVectorT::Mult_aTb : Unallocated a \n\n";

	for (int i=0; i<fLength; i++)
    c[i] = (*this)[i].Dot((*this)[i],b[i]); // <--- Will this work ?

}

//----------------------------------------------------

void FEA_dVectorT::Dot (const FEA_dMatrixT &B, FEA_dVectorT &c)
{
  if (fLength==0 || B.Length()==0)
		cout << "...ERROR >> FEA_dVectorT::Mult_aTb : Unallocated a \n\n";

	int j,l;

	for (l=0; l<n_ip; l++)
	  for (j=0; j<B[l].Cols(); j++)
    	c[l][j] = B[l].DotCol(j,(*this)[l]); 
}	

//----------------------------------------------------
 
void FEA_dVectorT::Dot (const FEA_dMatrixT &B, const FEA_dVectorT &c, FEA_dScalarT& d) 
{
  if (n_ip==0 || B.Length()==0 || c.n_ip==0)
		cout << "...ERROR >> FEA_dVectorT::Mult_aTb : Unallocated a \n\n";

	FEA_dVectorT temp(n_ip,n_sd);  // Faster way exists for aBc w/o temp
	Dot(B,temp); 
	temp.Dot(c,d); 

}	

//----------------------------------------------------
		
void FEA_dVectorT::Outer (const FEA_dVectorT& b, FEA_dMatrixT& C) 
{
  if (fLength==0 || b.Length()==0)
		cout << "...ERROR >> FEA_dVectorT::Mult_aTb : Unallocated a \n\n";

	for (int i=0; i<fLength; i++)
    C[i].Outer((*this)[i],b[i]); 

}

 
//----------------------------------------------------

void FEA_dVectorT::SumOf(const FEA_dVectorT &a, const FEA_dVectorT &b) {

  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++) 
    (*this)[i].SumOf(a[i],b[i]);
		
}

//----------------------------------------------------

void FEA_dVectorT::DiffOf(const FEA_dVectorT &a, const FEA_dVectorT &b) {

  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++) 
    (*this)[i].DiffOf(a[i],b[i]); 
	
}

//----------------------------------------------------

void FEA_dVectorT::MultAb  (const FEA_dMatrixT &A, const FEA_dVectorT &b) // Untested 
{
	for (int i=0; i<n_ip; i++)
		A[i].Multx ( b[i], (*this)[i] );
};

//----------------------------------------------------

void FEA_dVectorT::MultAb  (const FEA_dMatrixT &A, const dArrayT &b) // Untested 
{
	for (int i=0; i<n_ip; i++)
		A[i].Multx ( b, (*this)[i] );
};

//----------------------------------------------------
		
void FEA_dVectorT::MultATb (const FEA_dMatrixT &A, const FEA_dVectorT &b) // Untested 
{
	for (int i=0; i<n_ip; i++)
		A[i].MultTx ( b[i], (*this)[i] );
};

//----------------------------------------------------

void FEA_dVectorT::operator = (const FEA_dVectorT &a) {

  if (a.n_ip==0) 
		cout << "...ERROR >> FEA_dVectorT::operator = : rhs 0 length \n";
	if (n_ip==0)
  	FEA_Dimension (a);
	else if ( n_ip!=a.n_ip || n_sd!=a.n_sd) 
  	FEA_Dimension (a);

	for (int i=0; i<n_ip; i++) 
		(*this)[i] = a[i]; 
	
}

//----------------------------------------------------

void FEA_dVectorT::operator +=  (const FEA_dVectorT &a) { 
  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++) 
		(*this)[i] += a[i]; 
}

//----------------------------------------------------

 void FEA_dVectorT::operator -=  (const FEA_dVectorT &a) { 
  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++) 
		(*this)[i] -= a[i]; 
}

//----------------------------------------------------

void FEA_dVectorT::operator =  (const double &a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dVectorT: Vector unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] = a; 
}

//----------------------------------------------------

void FEA_dVectorT::operator *=  (const double &a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dVectorT: Vector unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] *= a; 
}

//----------------------------------------------------

void FEA_dVectorT::operator /=  (const double &a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dVectorT: Vector unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] /= a; 
}

/*

//----------------------------------------------------

void FEA_dVectorT::operator =  (const double *a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dVectorT: Vector unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] = a; 
}

//----------------------------------------------------

void FEA_dVectorT::operator *=  (const double *a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dVectorT: Vector unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] *= a; 
}

//----------------------------------------------------

void FEA_dVectorT::operator /=  (const double *a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dVectorT: Vector unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] /= a; 
}

*/

//----------------------------------------------------
//Un-tested
void FEA_dVectorT::operator *=  (const FEA_dScalarT &s) { 
	for (int i=0; i<fLength; i++) 
		(*this)[i] *= s[i]; 
}

//----------------------------------------------------
//Un-tested
void FEA_dVectorT::operator /=  (const FEA_dScalarT &s) { 
	for (int i=0; i<fLength; i++) 
		(*this)[i] /= s[i]; 
}

//----------------------------------------------------

FEA_EquateT& FEA_dVectorT::operator()(const int i) 
{
	double *p = FEA_Pointer (i);
	extern FEA_StackT* fStack;
	int n = fStack->Next_Shallow_Stack();
	fStack->Shallow_Stack[n].length = n_ip; 

	for (int l=0; l<n_ip; l++) { 
    fStack->Shallow_Stack[n].vec_ptrs[l] = p;  // Shallow copy allows data of (*this) LHS matrix to change
		p += n_sd;
	}

  return	fStack->Shallow_Stack[n];
}

#if 0 
//-- Antiquated method : Changed to new one 28FEB03
FEA_EquateT& FEA_dVectorT::operator()(const int i) {
  if (fLength==0) cout <<"..ERROR>> FEA_dVectorT: Vector unallocated"; 

	for (int l=0; l<n_ip; l++) 
		ip_components.vec_ptrs[l] = (*this)[l].Pointer(i);

  return(ip_components);
}
#endif


