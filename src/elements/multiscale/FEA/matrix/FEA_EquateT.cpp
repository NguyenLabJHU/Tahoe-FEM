// $Id: FEA_EquateT.cpp,v 1.5 2004-07-30 17:21:20 paklein Exp $
#include "FEA.h"

using namespace Tahoe;

//--------------------

FEA_EquateT::FEA_EquateT(void):
	length(0),
	vec_ptrs(NULL)
{

}

//--------------------

FEA_EquateT::FEA_EquateT(const int len):
	length(0),
	vec_ptrs(NULL)
{
	Allocate(len);
}

//--------------------

/* destructor */
FEA_EquateT::~FEA_EquateT(void) 
{
	delete [] vec_ptrs; 
}

//--------------------

void FEA_EquateT::Allocate(const int len)
{
	if (len != length)
	{
		/* free existing */
		delete[] vec_ptrs;
		
		length = len;
		if (length > 0)
			vec_ptrs = new double*[length];
		else
			vec_ptrs = NULL;
	} 
}

//##################################################################################
// NOTE: For these methods, don't error check length vs. length of stack (100)
//##################################################################################
//##################################################################################

void FEA_EquateT::operator = (const FEA_EquateT& a)
{

	for (int i=0; i<length; i++)  { 
  	double* p = vec_ptrs[i]; 
  	double* q = a.vec_ptrs[i];
    *p = *q;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//-------------------- 

void FEA_EquateT::operator += (const FEA_EquateT& a)  
{

	for (int i=0; i<length; i++)  { 
  	double* p = vec_ptrs[i]; 
  	double* q = a.vec_ptrs[i];
    *p += *q;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//-------------------- 

void FEA_EquateT::operator -= (const FEA_EquateT& a)  
{

	for (int i=0; i<length; i++)  { 
  	double* p = vec_ptrs[i]; 
  	double* q = a.vec_ptrs[i];
    *p -= *q;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//-------------------- 

void FEA_EquateT::operator *= (const FEA_EquateT& a)  
{

	for (int i=0; i<length; i++)  { 
  	double* p = vec_ptrs[i]; 
  	double* q = a.vec_ptrs[i];
    *p *= *q;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//-------------------- 

void FEA_EquateT::operator /= (const FEA_EquateT& a)  
{

	for (int i=0; i<length; i++)  { 
  	double* p = vec_ptrs[i]; 
  	double* q = a.vec_ptrs[i];
		cout << " p = "<<*p<<"\n\n";
		cout << " q = "<<*q<<"\n\n";
    *p /= *q;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//##################################################################################
//##################################################################################

void FEA_EquateT::operator = (const FEA_dScalarT& a)  
{
	const double *q = a.Pointer();
	for (int i=0; i<length; i++)  { 
  	double* p = vec_ptrs[i]; 
    *p = *q++;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//-------------------- 

void FEA_EquateT::operator += (const FEA_dScalarT& a)  
{
	const double *q = a.Pointer();
	for (int i=0; i<length; i++)  { 
  	double* p = vec_ptrs[i]; 
    *p += *q++;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//-------------------- 

void FEA_EquateT::operator -= (const FEA_dScalarT& a)  
{
	const double *q = a.Pointer();
	for (int i=0; i<length; i++)  { 
  	double* p = vec_ptrs[i]; 
    *p -= *q++;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//-------------------- 

void FEA_EquateT::operator *= (const FEA_dScalarT& a)  
{
	const double *q = a.Pointer();
	for (int i=0; i<length; i++)  { 
  	double* p = vec_ptrs[i]; 
    *p *= *q++;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//-------------------- 

void FEA_EquateT::operator /= (const FEA_dScalarT& a)  
{
	const double *q = a.Pointer();
	for (int i=0; i<length; i++)  { 
  	double* p = vec_ptrs[i]; 
    *p /= *q++;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//##################################################################################
//##################################################################################

void FEA_EquateT::operator = (const double *vector)
{

	for (int i=0; i<length; i++)  {
  	double* p = vec_ptrs[i]; 
    *p = vector[i]; 
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//##################################################################################
//##################################################################################
// CAUTION: This sets the value at ALL ips to value 

void FEA_EquateT::operator = (const double& value)
{

	for (int i=0; i<length; i++)  {
  	double* p = vec_ptrs[i]; 
    *p = value;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}
//-------------------- 

void FEA_EquateT::operator += (const double& a)
{

	for (int i=0; i<length; i++)  {
  	double* p = vec_ptrs[i]; 
    *p += a;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}
//-------------------- 

void FEA_EquateT::operator -= (const double& a)
{

	for (int i=0; i<length; i++)  {
  	double* p = vec_ptrs[i]; 
    *p -= a;
	}

extern FEA_StackT* fStack;
fStack->Reset();
}

//-------------------- 

void FEA_EquateT::operator *= (const double& a)
{

	for (int i=0; i<length; i++)  {
  	double* p = vec_ptrs[i]; 
    *p *= a;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//-------------------- 

void FEA_EquateT::operator /= (const double& a)
{

	for (int i=0; i<length; i++)  {
  	double* p = vec_ptrs[i]; 
    *p /= a;
	}	

extern FEA_StackT* fStack;
fStack->Reset();
}

//##################################################################################
//##################################################################################

FEA_EquateT& FEA_EquateT::operator + (const FEA_EquateT& a) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = length; // Set so != 100

	for (int i=0; i<length; i++)  
    (*fStack->Stack[n].vec_ptrs[i]) = (*vec_ptrs[i]) + (*a.vec_ptrs[i]); 

return fStack->Stack[n]; 
}

//-------------------- 

FEA_EquateT& FEA_EquateT::operator - (const FEA_EquateT& a) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = length; // Set so != 100

	for (int i=0; i<length; i++)  
    (*fStack->Stack[n].vec_ptrs[i]) = (*vec_ptrs[i]) - (*a.vec_ptrs[i]); 

return fStack->Stack[n]; 
}

//-------------------- 

FEA_EquateT& FEA_EquateT::operator * (const FEA_EquateT& a) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = length; // Set so != 100

	for (int i=0; i<length; i++)  
    (*fStack->Stack[n].vec_ptrs[i]) = (*vec_ptrs[i]) * (*a.vec_ptrs[i]); 

return fStack->Stack[n]; 
}

//-------------------- 

FEA_EquateT& FEA_EquateT::operator / (const FEA_EquateT& a) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = length; // Set so != 100

	for (int i=0; i<length; i++)  
    (*fStack->Stack[n].vec_ptrs[i]) = (*vec_ptrs[i]) / (*a.vec_ptrs[i]); 

return fStack->Stack[n]; 
}

//##################################################################################
//##################################################################################

FEA_EquateT& FEA_EquateT::operator + (const FEA_dScalarT& a) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = length; // Set so != 100

	for (int i=0; i<length; i++)  
    (*fStack->Stack[n].vec_ptrs[i]) = (*vec_ptrs[i]) + a[i]; 

return fStack->Stack[n]; 

}

//-------------------- 

FEA_EquateT& FEA_EquateT::operator - (const FEA_dScalarT& a) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = length; // Set so != 100

	for (int i=0; i<length; i++)  
    (*fStack->Stack[n].vec_ptrs[i]) = (*vec_ptrs[i]) - a[i]; 

return fStack->Stack[n]; 

}

//-------------------- 

FEA_EquateT& FEA_EquateT::operator * (const FEA_dScalarT& a) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = length; // Set so != 100

	for (int i=0; i<length; i++)  
    (*fStack->Stack[n].vec_ptrs[i]) = (*vec_ptrs[i]) * a[i]; 

return fStack->Stack[n]; 

}

//-------------------- 

FEA_EquateT& FEA_EquateT::operator / (const FEA_dScalarT& a) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = length; // Set so != 100

	for (int i=0; i<length; i++)  
    (*fStack->Stack[n].vec_ptrs[i]) = (*vec_ptrs[i]) / a[i]; 

return fStack->Stack[n]; 

}

//##################################################################################
//##################################################################################

FEA_EquateT& FEA_EquateT::operator + (const double& a) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = length; // Set so != 100

	for (int i=0; i<length; i++)  
    (*fStack->Stack[n].vec_ptrs[i]) = (*vec_ptrs[i]) + a; 

return fStack->Stack[n]; 

}

//-------------------- 

FEA_EquateT& FEA_EquateT::operator - (const double& a) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = length; // Set so != 100

	for (int i=0; i<length; i++)  
    (*fStack->Stack[n].vec_ptrs[i]) = (*vec_ptrs[i]) - a; 

return fStack->Stack[n]; 

}

//-------------------- 

FEA_EquateT& FEA_EquateT::operator * (const double& a) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = length; // Set so != 100

	for (int i=0; i<length; i++)  
    (*fStack->Stack[n].vec_ptrs[i]) = (*vec_ptrs[i]) * a; 

return fStack->Stack[n]; 

}

//-------------------- 

FEA_EquateT& FEA_EquateT::operator / (const double& a) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = length; // Set so != 100

	for (int i=0; i<length; i++)  
    (*fStack->Stack[n].vec_ptrs[i]) = (*vec_ptrs[i]) / a; 

return fStack->Stack[n]; 

}

//##################################################################################
//##################################################################################

void FEA_EquateT::Print() const { Print(" "); } 

void FEA_EquateT::Print(char *c) const // overload << later
{

	if (length==0)
		cout << "...ERROR >> FEA_EquateT::Print() : "<<c<<" Unallocated \n\n";

  cout << "\n "<<c<<" EquateT vec_ptrs[] evaluated at "<<length<<" inegration points (ip): \n"; 
  for (int i=0; i<length; i++) 
    cout <<"\n "<< c <<" @ ip "<<i<<": \n\n"<< (*vec_ptrs[i]) << "\n";
	cout << "\n";
}


