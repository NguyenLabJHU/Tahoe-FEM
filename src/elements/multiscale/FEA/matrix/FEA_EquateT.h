//DEVELOPMENT

#ifndef _FEA_EQUATET_H_
#define _FEA_EQUATET_H_

/** This class is used to equate components of FEA_dMatrixT, FEA_dVectorT, 
 *  and FEA_dScalarT */

namespace Tahoe {

class FEA_dScalarT; // Forward Declaration

class FEA_EquateT {

	public:

		double **vec_ptrs; // A vector of pointers to scalar components
		int length;

    FEA_EquateT(void);
    FEA_EquateT(const int len);

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
