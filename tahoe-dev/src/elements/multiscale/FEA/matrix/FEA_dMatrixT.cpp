// $Id: FEA_dMatrixT.cpp,v 1.13 2003-11-21 22:54:43 paklein Exp $
#include "FEA.h"

using namespace Tahoe; 

//------------- Constructors -------------------------

 FEA_dMatrixT::FEA_dMatrixT(void) : 
ArrayT<dMatrixT>() 
{
	n_ip = n_rows = n_cols = n_rows_x_n_cols = 0;
}

//----------------------------------------------------

FEA_dMatrixT::FEA_dMatrixT(int ips) : 
ArrayT<dMatrixT>(ips)
{
  n_ip = ips;
	n_rows = n_cols = n_rows_x_n_cols = 0;
}

//----------------------------------------------------

FEA_dMatrixT::FEA_dMatrixT(int ips,int n) 
{
	FEA_Dimension(ips,n);
}

//----------------------------------------------------

FEA_dMatrixT::FEA_dMatrixT(int ips,int rows,int cols) 
{
	FEA_Dimension(ips,rows,cols);
}

//------------ Utilities -----------------------------

void FEA_dMatrixT::FEA_Dimension (int ips,int rows,int cols) 
{
  n_ip   = ips;  // this also is/equals fLength set w/ Allocate
	n_rows = rows;
	n_cols = cols;
	n_rows_x_n_cols = n_rows * n_cols;
	n_ip_x_n_rows_x_n_cols = n_ip * n_rows * n_cols;

	Block_Memory.Allocate(n_ip*rows*cols); 
  Allocate(n_ip); // Allocate fArray
	for (int i=0; i<n_ip; i++)
		fArray[i].Set(rows, cols, Block_Memory.Pointer(i*n_rows_x_n_cols));

}		

//----------------------------------------------------

void FEA_dMatrixT::FEA_Dimension (int nip,int n) { FEA_Dimension(nip,n,n); }

//----------------------------------------------------

void FEA_dMatrixT::FEA_Dimension (const FEA_dMatrixT &a) // Copy Dimension of another
{
	if (a.fLength == 0)
    fLength = 0; 
	else
		FEA_Dimension(a.Length(), a[0].Rows(), a[0].Cols());
}

//----------------------------------------------------

void FEA_dMatrixT::FEA_Set ( int rows,int cols,const FEA_dMatrixT &A, int i,int j) 
{
  n_ip   = A.n_ip;  // this also is/equals fLength set w/ Allocate
	n_rows = rows;
	n_cols = cols;
	n_rows_x_n_cols = n_rows * n_cols;
	n_ip_x_n_rows_x_n_cols = n_ip * n_rows * n_cols;
	int A_rows = A.n_rows;	

  Allocate(n_ip); // Allocate fArray to put dMatricies in 
	for (int l=0; l<n_ip; l++)
		fArray[l].Alias(n_rows, n_cols, A[l].Pointer(A_rows*j + i));

}		
//----------------------------------------------------

void FEA_dMatrixT::FEA_UnSet ( void ) { for (int l=0; l<n_ip; l++) fArray[l].Set(n_rows, n_cols, 0); }		

//----------------------------------------------------

void FEA_dMatrixT::FEA_Delete ( void ) 
{
  //Block_Memory.Free();
	for (int i=0; i<n_ip; i++)
		(*this)[i].Free();

	Free();
	n_ip = n_rows = n_cols = n_rows_x_n_cols = n_ip_x_n_rows_x_n_cols = 0;
}	

//----------------------------------------------------

void FEA_dMatrixT::Print() const { Print(" "); } 

void FEA_dMatrixT::Print(char *c) const
{ 

	if (fLength==0)
		cout << "...ERROR >> FEA_dMatrixT::Print() : "<<c<<" Unallocated \n\n";

  cout << "\n "<<c<<" Matrix evaluated at "<<fLength<<" inegration points (ip): \n"; 
  for (int i=0; i<fLength; i++) 
    cout <<"\n "<< c <<" @ ip "<<i<<": \n\n"<< (*this)[i] << "\n";
	cout << "\n";
}

//------------------------------------------------------------

void FEA_dMatrixT::print(char *c) const
{
	int _ncl=0;		 int _nch=n_cols-1;
	int _nrl=0;		 int _nrh=n_rows-1;

	if (fLength==0)
		cout << "...ERROR >> FEA_dMatrixT::Print() : "<<c<<" Unallocated \n\n";

	for (int l=0; l<n_ip; l++) {

  	cout << "\n "<<c<<" Matrix evaluated at "<<fLength<<" inegration points (ip): \n"; 

		printf("  "); 
		for (int j=_ncl; j<=_nch; j++) { if (j!=_nch) printf("%9d",j);  if(j==_nch) printf("%9d\n",j); }

		for (int i=_nrl; i<=_nrh; i++) {
    	printf("%3d |",i);
    	for (int j=_ncl; j<=_nch; j++) {
        if ( j<=_nch )  printf(" %8.1e", (*this)[l](i,j) );
        if ( j==_nch )  printf("\n"); 
    	}
		}
	}

	printf("\n");
}

//----------------------------------------------------

void FEA_dMatrixT::Random (int seed) 
{
	for (int i=0; i<fLength; i++) (*this)[i].Random (seed++); 
}

//----------------------------------------------------

void FEA_dMatrixT::Random	(double high_val,double low_val, int seed) 
{
	Random(seed);
	double *p = (*this)[0].Pointer();
	for (int i=0; i<n_ip_x_n_rows_x_n_cols; i++) {
		while (fabs(*p) < .1) // Get value:  .1 < x < 1.0
			(*p) *= 10.0;
		if (*p > 0)
			(*p) = low_val + (*p)*(high_val - low_val);
		else
			(*p) = low_val - (*p)*(high_val - low_val);

		p++;
	}

}

//------------ Matrix Operations ----------------------

FEA_dMatrixT& FEA_dMatrixT::T(void) { Transpose(); return (*this); }

void FEA_dMatrixT::Transpose(void) {
				
	for (int i=0; i<fLength; i++)
    (*this)[i].Transpose(); 
     
}

//----------------------------------------------------

void FEA_dMatrixT::Transpose(const FEA_dMatrixT& a) {

  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++)
    (*this)[i].Transpose(a[i]); 

}

//----------------------------------------------------

void FEA_dMatrixT::Inverse(void) {

	for (int i=0; i<fLength; i++)
    (*this)[i].Inverse(); 
}

//----------------------------------------------------

void FEA_dMatrixT::Inverse(const FEA_dMatrixT& a) {

  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++)
    (*this)[i].Inverse(a[i]); 
}

//----------------------------------------------------

void FEA_dMatrixT::SumOf(const FEA_dMatrixT &a, const FEA_dMatrixT &b) {

  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++) 
    (*this)[i].SumOf(a[i],b[i]);
		
}

//----------------------------------------------------

void FEA_dMatrixT::DiffOf(const FEA_dMatrixT &a, const FEA_dMatrixT &b) {

  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++) 
    (*this)[i].DiffOf(a[i],b[i]); 
}

//----------------------------------------------------
#if 0
void FEA_dMatrixT::Symmetrize (const FEA_dMatrixT &a) // ######### DOESNT WORK PROPERLY 
{
	for (int i=0; i<fLength; i++)
    (*this)[i].Symmetrize(a[i]); 
}
#endif

//----------------------------------------------------

void FEA_dMatrixT::Sym ( const FEA_dMatrixT &A ) 
{
	FEA_dMatrixT AT; 	AT.FEA_Dimension (A);
	AT.Transpose ( A );

	(*this)  = A;
	(*this) += AT;
	(*this) *= 0.5;

}

//----------------------------------------------------

void FEA_dMatrixT::Skew ( const FEA_dMatrixT &A ) 
{
	FEA_dMatrixT AT; 	AT.FEA_Dimension (A);
	AT.Transpose ( A );

	(*this)  = A;
	(*this) -= AT;
	(*this) *= 0.5;

}

//----------------------------------------------------

void FEA_dMatrixT::Identity (double value) {

	for (int i=0; i<fLength; i++)
    (*this)[i].Identity(value);
}

//----------------------------------------------------

void FEA_dMatrixT::PlusIdentity (double value) {

	for (int i=0; i<fLength; i++)
    (*this)[i].PlusIdentity(value);
}

//----------------------------------------------------

void FEA_dMatrixT::Determinant (FEA_dScalarT &det) {

	for (int i=0; i<fLength; i++)
    det[i] = (*this)[i].Det();
}

//----------------------------------------------------

void FEA_dMatrixT::Trace (FEA_dScalarT &trace) {

	for (int i=0; i<fLength; i++)
    trace[i] = (*this)[i].Trace();
}

//----------------------------------------------------

void FEA_dMatrixT::operator =  (const FEA_dMatrixT &a) {

  if ( fLength==0 && a.Length()==0 ) 
		return;
  if ( fLength!=0 && a.Length()==0 ) 
		cout << "...ERROR >> FEA_dMatrixT::operator = : rhs 0 length \n";
  if ( fLength==0 && a.Length()!=0 ) 
  	FEA_Dimension (a);
	else if ( fLength!=a.Length() || (*this)[0].Rows()!=a[0].Rows() || (*this)[0].Cols()!=a[0].Cols() ) 
  	FEA_Dimension (a);

	for (int i=0; i<fLength; i++) 
		(*this)[i] = a[i]; 
	
}
//----------------------------------------------------

void FEA_dMatrixT::operator +=  (const FEA_dMatrixT &a) { 
  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++) 
		(*this)[i] += a[i]; 
}

//----------------------------------------------------

void FEA_dMatrixT::operator -=  (const FEA_dMatrixT &a) { 
  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++) 
		(*this)[i] -= a[i]; 
}
/*
//----------------------------------------------------

void FEA_dMatrixT::operator =  (const double &a) {  

  if (Length()==0) cout << "...ERROR >> FEA_dMatrixT::operator = : Matrix Unallocated \n";
	for (int i=0; i<fLength; i++) 
		(*this)[i] = a;	
}
*/
//----------------------------------------------------

void FEA_dMatrixT::operator =  (const double &a) { // Un-tested 

  double *p = (*this)[0].Pointer();
  if (Length()==0) cout << "...ERROR >> FEA_dMatrixT::operator = : Matrix Unallocated \n";
	for (int i=0; i< n_ip_x_n_rows_x_n_cols; i++) 
		(*p++) = a;	
}

//----------------------------------------------------

void FEA_dMatrixT::operator +=  (const double &a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dMAtrixT: Matrix unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] += a; 
}

//----------------------------------------------------

void FEA_dMatrixT::operator -=  (const double &a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dMAtrixT: Matrix unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] -= a; 
}

//----------------------------------------------------

void FEA_dMatrixT::operator *=  (const double &a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dMAtrixT: Matrix unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] *= a; 
}

//----------------------------------------------------

void FEA_dMatrixT::operator /=  (const double &a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dMAtrixT: Matrix unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] /= a; 
}




//----------------------------------------------------
/*
void FEA_dMatrixT::operator =  (const double *a) { // Un-tested 

  double *p = (*this)[0].Pointer();
  if (Length()==0) cout << "...ERROR >> FEA_dMatrixT::operator = : Matrix Unallocated \n";
	for (int i=0; i< n_ip_x_n_rows_x_n_cols; i++) 
		(*p++) = a;	
}

//----------------------------------------------------

void FEA_dMatrixT::operator +=  (const double *a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dMAtrixT: Matrix unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] += a; 
}

//----------------------------------------------------

void FEA_dMatrixT::operator -=  (const double *a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dMAtrixT: Matrix unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] -= a; 
}

//----------------------------------------------------

void FEA_dMatrixT::operator *=  (const double *a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dMAtrixT: Matrix unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] *= a; 
}

//----------------------------------------------------

void FEA_dMatrixT::operator /=  (const double *a) { 
  if (fLength==0) cout <<"..ERROR>> FEA_dMAtrixT: Matrix unallocated"; 
	for (int i=0; i<fLength; i++) 
		(*this)[i] /= a; 
}

*/


//----------------------------------------------------
//Un-tested
void FEA_dMatrixT::operator *=  (const FEA_dScalarT &s) { 
	for (int i=0; i<fLength; i++) 
		(*this)[i] *= s[i]; 
}

//----------------------------------------------------
//Un-tested
void FEA_dMatrixT::operator /=  (const FEA_dScalarT &s) { 
	for (int i=0; i<fLength; i++) 
		(*this)[i] /= s[i]; 
}

//----------------------------------------------------

FEA_EquateT& FEA_dMatrixT::operator()(const int i, const int j) 
{
	double *p = FEA_Pointer (n_rows*j + i);
	extern FEA_StackT* fStack;
	int n = fStack->Next_Shallow_Stack();
	fStack->Shallow_Stack[n].length = n_ip; 

	for (int l=0; l<n_ip; l++) { 
    fStack->Shallow_Stack[n].vec_ptrs[l] = p;  // Shallow copy allows data of (*this) LHS matrix to change
		p += n_rows_x_n_cols;
	}

  return	fStack->Shallow_Stack[n];
}

//------------ Matrix-Matrix Operations ---------------

void FEA_dMatrixT::MultAB (const FEA_dMatrixT &a, const FEA_dMatrixT &b, int upper) {
				
  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++)
    (*this)[i].MultAB (a[i], b[i], upper); 
}

//----------------------------------------------------

void FEA_dMatrixT::MultAB (const dMatrixT &a, const FEA_dMatrixT &b, int upper) {
				
	for (int i=0; i<fLength; i++)
    (*this)[i].MultAB (a, b[i], upper); 
}

//----------------------------------------------------
void FEA_dMatrixT::MultAB (const FEA_dMatrixT &a, const dMatrixT &b, int upper) {
				
  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++)
    (*this)[i].MultAB (a[i], b, upper); 
}

//----------------------------------------------------
void FEA_dMatrixT::MultATB (const FEA_dMatrixT &a, const FEA_dMatrixT &b, int upper) {
				
  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++)
    (*this)[i].MultATB (a[i], b[i], upper); 
     
}

//----------------------------------------------------

void FEA_dMatrixT::MultABT (const FEA_dMatrixT &a, const FEA_dMatrixT &b, int upper) {
				
  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++)
    (*this)[i].MultABT (a[i], b[i], upper); 
     
}

//----------------------------------------------------

void FEA_dMatrixT::MultATBT (const FEA_dMatrixT &a, const FEA_dMatrixT &b) {
				
  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++)
    (*this)[i].MultATBT (a[i], b[i]); 
     
}

//----------------------------------------------------
void FEA_dMatrixT::Outer (const FEA_dVectorT &a, const FEA_dVectorT &b) {
				
	if (fLength==0) {
		int nip = a.IPs();
		int dim = a.Rows();
		FEA_Dimension (nip, dim);
	}
	for (int i = 0; i < fLength; i++)
    	(*this)[i].Outer (a[i], b[i]); 
}

//----------------------------------------------------

//------------- Matrix-Matrix-Matrix Operations ------

void FEA_dMatrixT::MultABC  (const FEA_dMatrixT &a, const FEA_dMatrixT &b, 
								const FEA_dMatrixT &c, int range, int fillmode) {
				
  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++)
    (*this)[i].MultABC (a[i], b[i], c[i], range, fillmode); 
     
}

//----------------------------------------------------

void FEA_dMatrixT::MultABCT (const FEA_dMatrixT &a, const FEA_dMatrixT &b, 
								const FEA_dMatrixT &c, int range, int fillmode) {
				
  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++)
    (*this)[i].MultABCT (a[i], b[i], c[i], range, fillmode); 
     
}

//----------------------------------------------------

void FEA_dMatrixT::MultATBC (const FEA_dMatrixT &a, const FEA_dMatrixT &b, 
								const FEA_dMatrixT &c, int range, int fillmode) {
				
  if (fLength==0) FEA_Dimension (a);
	for (int i=0; i<fLength; i++)
    (*this)[i].MultATBC (a[i], b[i], c[i], range, fillmode); 
     
}

//------------- Misc Dot Products ---------------------

void FEA_dMatrixT::Double_Dot (const FEA_dMatrixT &a, FEA_dScalarT &s) 
{
  int i,l; 
  register double dot = 0.0;
	double *p  = (*this)[0].Pointer (); 
	const double *pa = a[0].Pointer ();

	for (l=0; l<n_ip; l++) {
 	  for (i=0; i<n_rows_x_n_cols; i++) 
     	dot += (*p++)*(*pa++); 
		s[l] = dot;
		dot = 0.0;
	}

}

//----------------------------------------------------

void FEA_dMatrixT::Match_Signs (const FEA_dMatrixT &A) {
				
  if (fLength==0) FEA_Dimension (A);
	double *p  = (*this)[0].Pointer (); 
	const double *q  = A[0].Pointer (); 

	for (int i=0; i<n_ip_x_n_rows_x_n_cols-1; i++) {
		if ( (*p)*(*q) < 0.0 ) // If signs are opposite
			  (*p) *= -1.0;				// Change sign of (*this) to that of A
		*p++; *q++;
	}
	if ( (*p)*(*q) < 0.0 ) // Check last item of array (no ++ing)
			 (*p) *= -1.0;				
}

//----------------------------------------------------

void FEA_dMatrixT::AddBlock (const int &i,const int &j,const FEA_dMatrixT &block) 
{
	for (int l=0; l<n_ip; l++) (*this)[l].AddBlock(i,j,block[l]); 
}

//----------------------------------------------------

void FEA_dMatrixT::CopyBlock (const int &i,const int &j, FEA_dMatrixT &block) 
{
	for (int l=0; l<n_ip; l++) (*this)[l].CopyBlock(i,j,block[l]); 
}

//----------------------------------------------------

void FEA_dMatrixT::SetBlock (const int &i,const int &j,const FEA_dMatrixT &block) 
{
	for (int l=0; l<n_ip; l++) (*this)[l].SetBlock(i,j,block[l]); 
}

//----------------------------------------------------

void FEA_dMatrixT::Swap_Rows (const int &row1, const int &row2) 
{
	FEA_dMatrixT row1_hold (n_ip,1,n_cols);  // Fake as an Array
	FEA_dMatrixT row2_hold (n_ip,1,n_cols);

	CopyBlock	(row1,0, row1_hold); 
	CopyBlock	(row2,0, row2_hold); 
	SetBlock 	(row1,0, row2_hold); 
	SetBlock 	(row2,0, row1_hold); 
}

//----------------------------------------------------

void FEA_dMatrixT::Swap_Cols (const int &col1, const int &col2) 
{
	FEA_dMatrixT col1_hold (n_ip,n_rows,1);  // Fake as an Array
	FEA_dMatrixT col2_hold (n_ip,n_rows,1);

	CopyBlock	(0,col1, col1_hold); 
	CopyBlock	(0,col2, col2_hold); 
	SetBlock 	(0,col1, col2_hold); 
	SetBlock 	(0,col2, col1_hold); 
}

//----------------------------------------------------

void FEA_dMatrixT::Mag_and_Dir (FEA_dScalarT &mag, FEA_dMatrixT &N) 
{ 
	Magnitude (mag);
	N = (*this);
	if (mag!=0.0)
		N /= mag;
}

//----------------------------------------------------

void FEA_dMatrixT::Direction (FEA_dMatrixT &N) 
{ 
	FEA_dScalarT s(n_ip);
	Magnitude (s); 
	N = (*this);
	if (s!=0.0)
		N /= s;
}

//----------------------------------------------------

void FEA_dMatrixT::Magnitude (FEA_dScalarT &s) { Magnitude_Squared (s); s.Sqrt(); }

//----------------------------------------------------

void FEA_dMatrixT::Magnitude_Squared (FEA_dScalarT &s) 
{
	for (int l=0; l<n_ip; l++) 
		s[l] = (*this)[l].ScalarProduct(); 
}

//-------------- Special Operations -------------------

FEA_EquateT& FEA_dMatrixT::Dot (int rc,int ij,const FEA_dMatrixT &a,int a_rc,int a_ij)
{
	extern FEA_StackT* fStack;
	int n = fStack->Next_Stack();
	fStack->Stack[n].length = n_ip;
	int l,i,j,dim_check=1;

	//------------------ error checking 
  if (dim_check) {
  	if (fLength==0 || a.fLength==0) 
			cout <<"..ERROR>> FEA_dMAtrixT: Matrix unallocated"; 
		if (rc==FEA::kCol && a_rc==FEA::kCol)  
			if (n_rows != a.n_rows)
   			cout <<"...ERROR >> FEA_dMatrixT::Dot incompatible dimensions. ";
	  else if (rc==FEA::kRow && a_rc==FEA::kRow) 
			if (n_cols != a.n_cols)
   			cout <<"...ERROR >> FEA_dMatrixT::Dot incompatible dimensions. ";
		else if (rc==FEA::kCol && a_rc==FEA::kRow) 
			if (n_rows != a.n_cols)
   			cout <<"...ERROR >> FEA_dMatrixT::Dot incompatible dimensions. ";
		else if (rc==FEA::kRow && a_rc==FEA::kCol) 
			if (n_cols != a.n_rows)
   			cout <<"...ERROR >> FEA_dMatrixT::Dot incompatible dimensions. ";
	}
	//------------------

	if (rc==FEA::kCol && a_rc==FEA::kCol)  
		for (l=0; l<n_ip; l++) {
			double *p = (*this)[l].Pointer(n_rows*ij); // ij = { 0 < ij < n-1 }
			const double *q = a[l].Pointer(a.n_rows*a_ij);
      (*fStack->Stack[n].vec_ptrs[l]) = 0.0;
   	  for (i=0; i<n_rows; i++) 
      	(*fStack->Stack[n].vec_ptrs[l]) += (*p++)*(*q++);   // think of as two ops *p; then p++;
		}
	else if (rc==FEA::kRow && a_rc==FEA::kRow) 
		for (l=0; l<n_ip; l++) {
			double *p = (*this)[l].Pointer(ij);
			const double *q = a[l].Pointer(a_ij);
      (*fStack->Stack[n].vec_ptrs[l]) = 0.0;
   	  for (j=0; j<n_cols; j++) {
      	(*fStack->Stack[n].vec_ptrs[l]) += (*p)*(*q);
				p+=n_rows;
				q+=a.n_rows;
			}
		}
	else if (rc==FEA::kCol && a_rc==FEA::kRow) 
		for (l=0; l<n_ip; l++) {
			double *p = (*this)[l].Pointer(n_rows*ij);
			const double *q = a[l].Pointer(a_ij);
      (*fStack->Stack[n].vec_ptrs[l]) = 0.0;
   	  for (j=0; j<n_rows; j++) {
      	(*fStack->Stack[n].vec_ptrs[l]) += (*p)*(*q);
				p++; 
				q+=a.n_rows;
			}
		}
	else if (rc==FEA::kRow && a_rc==FEA::kCol) 
		for (l=0; l<n_ip; l++) {
			double *p = (*this)[l].Pointer(ij);
			const double *q = a[l].Pointer(a.n_rows*a_ij);
      (*fStack->Stack[n].vec_ptrs[l]) = 0.0;
   	  for (j=0; j<n_cols; j++) {
      	(*fStack->Stack[n].vec_ptrs[l]) += (*p)*(*q);
				p+=n_rows;
				q++; 
			}
		}
	else 	
		cout <<"...ERROR >> FEA_dMatrixT::Dot() : either FEA::kRow or FEA::kCol must be sent. ";

return fStack->Stack[n]; 
} 

//----------------------------------------------------

FEA_EquateT& FEA_dMatrixT::Dot_Aij (int rc,int ij,const FEA_dMatrixT &a,int a_rc,int a_ij, FEA_dMatrixT &c,int i,int j)
{
	extern FEA_StackT* fStack;
	int n = fStack->Next_Stack();
	fStack->Stack[n].length = n_ip;
	int k,l,dim_check=1;

	//------------------ ERROR CHECKING
  if (dim_check) {
  	if (fLength==0 || a.fLength==0) 
			cout <<"..ERROR>> FEA_dMAtrixT: Matrix unallocated"; 
		if (rc==FEA::kCol && a_rc==FEA::kCol)  
			if (n_rows != a.n_rows)
   			cout <<"...ERROR >> FEA_dMatrixT::Dot incompatible dimensions. ";
	  else if (rc==FEA::kRow && a_rc==FEA::kRow) 
			if (n_cols != a.n_cols)
   			cout <<"...ERROR >> FEA_dMatrixT::Dot incompatible dimensions. ";
		else if (rc==FEA::kCol && a_rc==FEA::kRow) 
			if (n_rows != a.n_cols)
   			cout <<"...ERROR >> FEA_dMatrixT::Dot incompatible dimensions. ";
		else if (rc==FEA::kRow && a_rc==FEA::kCol) 
			if (n_cols != a.n_rows)
   			cout <<"...ERROR >> FEA_dMatrixT::Dot incompatible dimensions. ";
	}
	//------------------

	double *r = c[0].Pointer(n_rows*j + i);

	if (rc==FEA::kCol && a_rc==FEA::kCol)  
		for (l=0; l<n_ip; l++) {
			double *p = (*this)[l].Pointer(n_rows*ij);
			const double *q = a[l].Pointer(a.n_rows*a_ij);
      (*fStack->Stack[n].vec_ptrs[l]) = 0.0;
   	  for (k=0; k<n_rows; k++) 
      	(*fStack->Stack[n].vec_ptrs[l]) += (*p++)*(*q++);  
			(*fStack->Stack[n].vec_ptrs[l]) *= (*r);
			r += c.n_rows_x_n_cols;
		}
	else if (rc==FEA::kRow && a_rc==FEA::kRow) 
		for (l=0; l<n_ip; l++) {
			double *p = (*this)[l].Pointer(ij);
			const double *q = a[l].Pointer(a_ij);
      (*fStack->Stack[n].vec_ptrs[l]) = 0.0;
   	  for (k=0; k<n_cols; k++) {
      	(*fStack->Stack[n].vec_ptrs[l]) += (*p)*(*q);
				p+=n_rows;
				q+=a.n_rows;
			}
			(*fStack->Stack[n].vec_ptrs[l]) *= (*r);
			r += c.n_rows_x_n_cols;
		}
	else if (rc==FEA::kCol && a_rc==FEA::kRow) 
		for (l=0; l<n_ip; l++) {
			double *p = (*this)[l].Pointer(n_rows*ij);
			const double *q = a[l].Pointer(a_ij);
      (*fStack->Stack[n].vec_ptrs[l]) = 0.0;
   	  for (k=0; k<n_rows; k++) {
      	(*fStack->Stack[n].vec_ptrs[l]) += (*p)*(*q);
				p++; 
				q+=a.n_rows;
			}
			(*fStack->Stack[n].vec_ptrs[l]) *= (*r);
			r += c.n_rows_x_n_cols;
		}
	else if (rc==FEA::kRow && a_rc==FEA::kCol) 
		for (l=0; l<n_ip; l++) {
			double *p = (*this)[l].Pointer(ij);
			const double *q = a[l].Pointer(a.n_rows*a_ij);
      (*fStack->Stack[n].vec_ptrs[l]) = 0.0;
   	  for (k=0; k<n_cols; k++) {
      	(*fStack->Stack[n].vec_ptrs[l]) += (*p)*(*q);
				p+=n_rows;
				q++; 
			}
			(*fStack->Stack[n].vec_ptrs[l]) *= (*r);
			r += c.n_rows_x_n_cols;
		}
	else 	
		cout <<"...ERROR >> FEA_dMatrixT::Dot() : either FEA::kRow or FEA::kCol must be sent. ";

return fStack->Stack[n]; 
} 

//----------------------------------------------------

FEA_EquateT& FEA_dMatrixT::ij_x_Aij (int i,int j,const FEA_dMatrixT &a,int ii,int jj) 
{
extern FEA_StackT* fStack;
int n = fStack->Next_Stack();
fStack->Stack[n].length = n_ip;

	int l,dim_check=1;

  if (dim_check) {
    if (fLength==0) 
			cout <<"..ERROR>> FEA_dMAtrixT: Matrix unallocated"; 
		if (n_rows != a.n_rows)
   		cout <<"...ERROR >> FEA_dMatrixT::Dot incompatible dimensions. ";
	}

	double *p = (*this)[0].Pointer (n_rows*j +i);  // ij Component of 1st (0th) dMatrix
	const double *q = a[0].Pointer (n_rows*jj +ii);      // ij Component of 1st (0th) dMatrix

	for (l=0; l<n_ip; l++) {
    (*fStack->Stack[n].vec_ptrs[l]) = (*p)*(*q);   
		p += n_rows_x_n_cols; 
		q += a.n_rows_x_n_cols; 
	}

return fStack->Stack[n]; 
}



