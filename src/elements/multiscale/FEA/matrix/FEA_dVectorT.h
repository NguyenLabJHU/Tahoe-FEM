// $Id: FEA_dVectorT.h,v 1.8 2003-11-21 22:54:44 paklein Exp $
#ifndef _FEA_DVECTORT_H_
#define _FEA_DVECTORT_H_

namespace Tahoe {

/** In FEA analysis, if gauss quadrature integration is used, a vector, 
 *  (say v for example) must be evaluated a the integration points.  This
 *  class is a way of masking the fact that we must do operations on all the 
 *  v's (i.e. one per int point).  For example, to calculate v.u, simply write:
 *  FEA_dMatrixT v(n_ip) dotted with FEA_dVectorT u(n_ip,n_sd) equals
 *  FEA_dScalarT c(n_ip) by:  v.Dot(u,c); 
 *  No need to loop through ip's, it's done automatically. An instance ou
 *  FEA_dVectorT is an Array of vectors, each evaluated at different ips. **/

class FEA_dVectorT: public ArrayT <dArrayT>
{
	public:

		// constructors
					
		FEA_dVectorT (void);
		FEA_dVectorT (int n_ip);
		FEA_dVectorT (int n_ip,int length); 
		
    // utilities -- name mangling problem when tried to just call it "Dimension()"
		
		void FEA_Dimension (const FEA_dVectorT &a);
		void FEA_Dimension (int ips,int length);
		void Print(void) { Print(" "); }
		void Print(char*);

		double* FEA_Pointer   (int offset) { return (*this)[0].Pointer(offset); }

		int IPs  (void) const { return n_ip; } 
		int Rows (void) const { return n_sd; } 

   	// vector operations
/*		
		void One_Norm      	(void);
		void Two_Norm      	(void);
		void Inf_Norm      	(void);
		void Projector  		(void);
		void Orthog_Projector (void);
*/
   		// vector-vector / vector-Tensor operations
   	
   		void Magnitude  (FEA_dScalarT &s);
		
    	/** a.b = c <==> aTb = c  (1x3)(3x1) = (1x1) */  
		void Dot (const FEA_dVectorT& b, FEA_dScalarT& c);  
		
		/** a.B = c <==>   aTB = cT  (1x3).(3x3) = (1x3) */
		void Dot (const FEA_dMatrixT &B, FEA_dVectorT &c); 

		/** a.Bc = c <==>  aTBc (1x3).(3x3).(3x1) = (1x1) */
    	void Dot (const FEA_dMatrixT &B, const FEA_dVectorT &c, FEA_dScalarT& d); 

		/** a o b = C <==> abT (3x1)(1x3) = (3x3) */
		void Outer 	 (const FEA_dVectorT& b, FEA_dMatrixT& C);  

		void SumOf   (const FEA_dVectorT& a, const FEA_dVectorT& b);
		void DiffOf  (const FEA_dVectorT& a, const FEA_dVectorT& b);

		void MultAb  (const FEA_dMatrixT &A, const FEA_dVectorT &b); 
		void MultAb  (const FEA_dMatrixT &A, const dArrayT &b); 			// Use for E(l)= B(l)*d
		void MultATb (const FEA_dMatrixT &A, const FEA_dVectorT &b); 

		// overloaded operators   NOTE: no way to do C=A*B or C=A+B w/o an extra deep copy (slower)
	
	  	void operator  =  (const FEA_dVectorT &a); 
		void operator +=  (const FEA_dVectorT &a); 
		void operator -=  (const FEA_dVectorT &a); 
	  	void operator  =  (const double &a);   
		void operator *=  (const double &a);   
		void operator /=  (const double &a);
		/* void operator  =  (const double *a);   
		void operator *=  (const double *a);   
		void operator /=  (const double *a); */

		void operator *=  (const FEA_dScalarT &s);   
		void operator /=  (const FEA_dScalarT &s);   

    FEA_EquateT& operator () (const int i); 
		
	protected:

		dArrayT  Block_Memory; // Don't really need this w/ new way of FEA_EquateT setup	
    FEA_EquateT  ip_components; 

   	int n_ip; // this is a redundant value of fLength 
	int n_sd; // this is a redundant value of (*this)[0].fLength
};

}

#endif /* _FEA_DVECTORT_H_ */
