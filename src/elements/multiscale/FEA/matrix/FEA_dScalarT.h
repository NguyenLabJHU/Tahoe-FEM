// $Id: FEA_dScalarT.h,v 1.10 2004-03-21 21:18:03 raregue Exp $
#ifndef _FEA_DSCALART_H_
#define _FEA_DSCALART_H_

namespace Tahoe {

/** In FEA analysis, if gauss quadrature integration is used, a scalar, 
 *  (say c for example) must be evaluated at the integration points.  This
 *  class is a way of masking the fact that we must do operations on all the 
 *  c's (i.e. one per int point).  For example, to calculate a.b = c, simply write:
 *  FEA_dVectorT a(n_ip) dotted with FEA_dVectorT b(n_ip,n_sd) equals FEA_ScalarT c(n_ip) 
 *  No need to loop through ip's, it's done automatically. An instance of 
 *  FEA_dScalarT is an Array of scalar values at different ips. **/

//#########################################################

class FEA_dVectorT; // Forward Declarations
class FEA_dMatrixT;

class FEA_dScalarT: public dArrayT // For the name sake only
{
  public:

		FEA_dScalarT(void) : dArrayT() { }  // Call Base Class Constructor 
		FEA_dScalarT(int length): dArrayT(length) { } // Call Base Class Constructor 
		void Print(void) { Print(" "); }
		void Print(char*);
    	int IPs(void) { return fLength; }
		void FEA_Dimension(int n_ip) { Dimension(n_ip); }

		void Sin			( FEA_dScalarT &s ) { for (int i=0; i<fLength; i++) (*this)[i] = sin  ( s[i] ); }
		void Cos			( FEA_dScalarT &s ) { for (int i=0; i<fLength; i++) (*this)[i] = cos  ( s[i] ); }
		void Tan			( FEA_dScalarT &s ) { for (int i=0; i<fLength; i++) (*this)[i] = tan  ( s[i] ); }
		void Exp  		( FEA_dScalarT &s ) { for (int i=0; i<fLength; i++) (*this)[i] = exp  ( s[i] ); }
		void Sinh 		( FEA_dScalarT &s ) { for (int i=0; i<fLength; i++) (*this)[i] = sinh ( s[i] ); }
		void ArcSinh 	( FEA_dScalarT &s ) { for (int i=0; i<fLength; i++) (*this)[i] = asinh ( s[i] ); }
		void Cosh 		( FEA_dScalarT &s ) { for (int i=0; i<fLength; i++) (*this)[i] = cosh ( s[i] ); }
		void Tanh 		( FEA_dScalarT &s ) { for (int i=0; i<fLength; i++) (*this)[i] = tanh ( s[i] ); }
		void Sech 		( FEA_dScalarT &s ) { for (int i=0; i<fLength; i++) (*this)[i] = 1.0/cosh ( s[i] ); }
		void Sqrt		( FEA_dScalarT &s ) { for (int i=0; i<fLength; i++) (*this)[i] = sqrt ( s[i] ); }
		void Macaulay 	( FEA_dScalarT &s ) { for (int i=0; i<fLength; i++) if (s[i]<0.0) (*this)[i] = s[i]; else (*this)[i] = 0.0; }
		//void Macaulay 	( FEA_dScalarT &s ) { for (int i=0; i<fLength; i++) if (s[i]<0.0) (*this)[i] = 0.0; else (*this)[i] = s[i]; }
		void Pow		( FEA_dScalarT &s, FEA_dScalarT &t ) { for (int i=0; i<fLength; i++) (*this)[i] = pow ( s[i],t[i] ); }
		void Pow		( FEA_dScalarT &s, double &t ) { for (int i=0; i<fLength; i++) (*this)[i] = pow ( s[i],t ); }
		void Pow		( double &s, FEA_dScalarT &t ) { for (int i=0; i<fLength; i++) (*this)[i] = pow ( s,t[i] ); }
		//void Abs		( FEA_dScalarT &t ) { for (int i=0; i<fLength; i++) (*this)[i] = abs ( t[i] ); }
		void Abs		( FEA_dScalarT &t ) { for (int i=0; i<fLength; i++) (*this)[i] = fabs ( t[i] ); }

		void Sin			( void ) 	{ for (int i=0; i<fLength; i++) (*this)[i] = sin  	( (*this)[i] ); }
		void Cos			( void ) 	{ for (int i=0; i<fLength; i++) (*this)[i] = cos  	( (*this)[i] ); }
		void Tan			( void ) 	{ for (int i=0; i<fLength; i++) (*this)[i] = tan  	( (*this)[i] ); }
		void Exp			( void ) 	{ for (int i=0; i<fLength; i++) (*this)[i] = exp  	( (*this)[i] ); }
		void Sinh			( void ) 	{ for (int i=0; i<fLength; i++) (*this)[i] = sinh 	( (*this)[i] ); }
		void ArcSinh		( void ) 	{ for (int i=0; i<fLength; i++) (*this)[i] = asinh 	( (*this)[i] ); }
		void Cosh			( void ) 	{ for (int i=0; i<fLength; i++) (*this)[i] = cosh 	( (*this)[i] ); }
		void Tanh			( void ) 	{ for (int i=0; i<fLength; i++) (*this)[i] = tanh ( (*this)[i] ); }
		void Sech			( void ) 	{ for (int i=0; i<fLength; i++) (*this)[i] = 1.0/cosh ( (*this)[i] ); }
		void Sqrt			( void ) 	{ for (int i=0; i<fLength; i++) (*this)[i] = sqrt ( (*this)[i] ); }
		void Macaulay 		( void )	{ for (int i=0; i<fLength; i++) if ( (*this)[i]<0.0 ) (*this)[i]=0.0; }

		void Squared  		( void )	{ for (int i=0; i<fLength; i++) (*this)[i] *= (*this)[i]; }

		void operator  = (const double *a ) 				{ for (int i=0; i<fLength; i++) (*this)[i]  = a[i]; }

		void operator  = (const FEA_dScalarT &a ) 	{ for (int i=0; i<fLength; i++) (*this)[i]  = a[i]; }
		void operator += (const FEA_dScalarT &a ) 	{ for (int i=0; i<fLength; i++) (*this)[i] += a[i]; }
		void operator -= (const FEA_dScalarT &a ) 	{ for (int i=0; i<fLength; i++) (*this)[i] -= a[i]; }
		void operator *= (const FEA_dScalarT &a ) 	{ for (int i=0; i<fLength; i++) (*this)[i] *= a[i]; }
		void operator /= (const FEA_dScalarT &a ) 	{ for (int i=0; i<fLength; i++) (*this)[i] /= a[i]; }  

		void operator  = (const FEA_EquateT &a )	  { for (int i=0; i<fLength; i++) (*this)[i]  = (*a.vec_ptrs[i]); }
		void operator += (const FEA_EquateT &a )	  { for (int i=0; i<fLength; i++) (*this)[i] += (*a.vec_ptrs[i]); }
		void operator -= (const FEA_EquateT &a )	  { for (int i=0; i<fLength; i++) (*this)[i] -= (*a.vec_ptrs[i]); }
		void operator *= (const FEA_EquateT &a )	  { for (int i=0; i<fLength; i++) (*this)[i] *= (*a.vec_ptrs[i]); }
		void operator /= (const FEA_EquateT &a )	  { for (int i=0; i<fLength; i++) (*this)[i] /= (*a.vec_ptrs[i]); }  

		void operator  = (const double &a ) 				{ for (int i=0; i<fLength; i++) (*this)[i]  = a; }
		void operator += (const double &a ) 				{ for (int i=0; i<fLength; i++) (*this)[i] += a; }
		void operator -= (const double &a ) 				{ for (int i=0; i<fLength; i++) (*this)[i] -= a; }
		void operator *= (const double &a ) 				{ for (int i=0; i<fLength; i++) (*this)[i] *= a; }
		void operator /= (const double &a ) 				{ for (int i=0; i<fLength; i++) (*this)[i] /= a; }

		bool operator == (const double  a) {   for(int i=0; i<fLength; i++) { if((*this)[i] == a){return 1;} }  return 0;  }
		bool operator != (const double  a) {   for(int i=0; i<fLength; i++) { if((*this)[i] == a){return 0;} }  return 1;  }
		bool operator >  (const double  a) {   for(int i=0; i<fLength; i++) { if((*this)[i] > a){return 1;} }  return 0;  }
		bool operator <  (const double  a) {   for(int i=0; i<fLength; i++) { if((*this)[i] < a){return 1;} }  return 0;  }

		void Dot 				 		( FEA_dVectorT 	&a, FEA_dVectorT &b ); 
		void Double_Dot  		( FEA_dMatrixT 	&A, FEA_dMatrixT &B ); 
		void L1_Norm 				( FEA_dMatrixT 	&A );
		void L1_Diagonal  	( FEA_dMatrixT 	&A );  //--- Get sum of abs ( diagonal terms )
		void Diagonal_Bias 	( FEA_dMatrixT 	&A );
		void Max			  		( FEA_dVectorT 	&a ); 
		void Max			  		( FEA_dMatrixT 	&A ); 
		void AbsMax			  	( FEA_dVectorT 	&a ); 
		void AbsMax			  	( FEA_dMatrixT 	&A ); 

		// Tensor Invariants
		void Trace 					( FEA_dMatrixT 	&A ); 
		void Determinant 		( FEA_dMatrixT 	&A ); 

};
    // RECALL: operator = already overloaded by base class 

//----------------------------------------------------

inline void FEA_dScalarT::Print(char *c) { // overload << later

  cout << "\n Scalar "<<c<<" evaluated at "<<fLength<<" inegration points (ip): \n"; 
  for (int i=0; i<fLength; i++) 
    cout <<"\n "<< c <<" @ ip "<<i<<" = "<< (*this)[i] << "\n";

}

} /* namespace Tahoe */


#endif
