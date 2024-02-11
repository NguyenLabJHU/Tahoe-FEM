// $Id: FEA_EquateT.h,v 1.3 2003-02-03 04:40:24 paklein Exp $
#ifndef _FEA_EQUATET_H_
#define _FEA_EQUATET_H_

namespace Tahoe {

class FEA_dScalarT; // Forward Declaration

/** This class is used to equate components of FEA_dMatrixT, FEA_dVectorT, 
 *  and FEA_dScalarT */
class FEA_EquateT {

	public:

		double **vec_ptrs; // A vector of pointers to scalar components
		int length;

    FEA_EquateT(void);
    FEA_EquateT(const int len);
    
	/** destructor */
	~FEA_EquateT(void);

		void Allocate(const int len);

		void operator  = (const FEA_EquateT& a);
		void operator += (const FEA_EquateT& a); 	
		void operator -= (const FEA_EquateT& a); 	
		void operator *= (const FEA_EquateT& a); 	
		void operator /= (const FEA_EquateT& a); 	

		void operator  = (const FEA_dScalarT& a); 
		void operator += (const FEA_dScalarT& a); 
		void operator -= (const FEA_dScalarT& a); 
		void operator *= (const FEA_dScalarT& a); 
		void operator /= (const FEA_dScalarT& a); 


		void operator  = (const double& a);
		void operator += (const double& a); 
		void operator -= (const double& a); 
		void operator *= (const double& a); 
		void operator /= (const double& a); 

		void operator  = (const double *vector);

		FEA_EquateT& operator + (const FEA_EquateT& a); 
		FEA_EquateT& operator - (const FEA_EquateT& a); 
		FEA_EquateT& operator * (const FEA_EquateT& a); 
		FEA_EquateT& operator / (const FEA_EquateT& a); 

		FEA_EquateT& operator + (const FEA_dScalarT& a); 
		FEA_EquateT& operator - (const FEA_dScalarT& a); 
		FEA_EquateT& operator * (const FEA_dScalarT& a); 
		FEA_EquateT& operator / (const FEA_dScalarT& a); 

		FEA_EquateT& operator + (const double& a); 
		FEA_EquateT& operator - (const double& a); 
		FEA_EquateT& operator * (const double& a); 
		FEA_EquateT& operator / (const double& a); 

		void Print(void) const;
		void Print(char*) const;

};

}

#endif // _FEA_EQUATET_H_
