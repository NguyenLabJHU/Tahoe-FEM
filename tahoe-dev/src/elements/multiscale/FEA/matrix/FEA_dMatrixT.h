// $Id: FEA_dMatrixT.h,v 1.11 2003-11-21 22:54:44 paklein Exp $
#ifndef _FEA_DMATRIXT_H_
#define _FEA_DMATRIXT_H_

namespace Tahoe {

/** In FEA analysis, if gauss quadrature integration is used, a tensor, 
 *  (say F for example) must be evaluated a the integration points.  This
 *  class is a way of masking the fact that we must do operations on all the 
 *  F's (i.e. one per int point).  For example, to calculate C, simply write:
 *  FEA_dMatrixT F(n_ip,n_sd,n_sd), C(n_ip,n_sd,n_sd);  C.MultATB(F,F); 
 *  No need to loop through ip's, it's done automatically. An instance of
 *  FEA_dMatrixT is an Array of tensors evaluated at different ips. **/

class FEA_dMatrixT: public ArrayT <dMatrixT>
{
	public:

		// constructors
					
		FEA_dMatrixT (void);
		FEA_dMatrixT (int n_ip);
		FEA_dMatrixT (int n_ip,int n); 
		FEA_dMatrixT (int n_ip,int rows,int cols);
	
	  // data access routines
	 
	  int IPs(void)  const { return n_ip;   } 
	  int Rows(void) const { return n_rows; }	
	  int Cols(void) const { return n_cols; } 
		
		double* FEA_Pointer   (int offset) { return (*this)[0].Pointer(offset); }

    // utilities -- name mangling problem when tried to just call it "Dimension()"
		
		void FEA_Dimension (const FEA_dMatrixT &a);
		void FEA_Dimension (int ips,int n);
		void FEA_Dimension (int ips,int rows,int cols);
		void FEA_Set			 (int rows,int cols,const FEA_dMatrixT &A,int i,int j);
		void FEA_UnSet		 (void); 
		void FEA_Delete		 (void);  
		void Print  (void)  const;
		void Print  (char*) const;
		void print  (char*) const;
		void Random (int seed=1);
		void Random (double high_val, double low_val=0, int seed=1);

   	// matrix operations
		
		FEA_dMatrixT&    T(void); // convert (*this) into its transpose, then return it 
		void Transpose 		(void); 
		void Transpose 		(const FEA_dMatrixT& fea_matrix);
		void Inverse   		(void);
		void Inverse   		(const FEA_dMatrixT& fea_matrix);
		void SumOf     		(const FEA_dMatrixT &a, const FEA_dMatrixT &b);
		void DiffOf    		(const FEA_dMatrixT &a, const FEA_dMatrixT &b);
		//void Symmetrize   (const FEA_dMatrixT &a );
		//void Symmetrize   (void);
		void Sym				  (const FEA_dMatrixT &A );	
		void Skew				  (const FEA_dMatrixT &A );	
		void Identity 		(double value=1.0);
		void PlusIdentity (double value=1.0);
		void Determinant	(FEA_dScalarT &det);
		void Det					(FEA_dScalarT &det) { Determinant(det); }
		void Trace 				(FEA_dScalarT &trace);

   	// matrix-matrix operations
		
		void Magnitude_Squared (FEA_dScalarT &s);
		void Magnitude  		(FEA_dScalarT &s);
		void Direction  		(FEA_dMatrixT &N);
		void Mag_and_Dir		(FEA_dScalarT &mag, FEA_dMatrixT &N);
		void Double_Dot 		(const FEA_dMatrixT &a, FEA_dScalarT &s);
		void Match_Signs 		(const FEA_dMatrixT &A);
		void Insert_SubMatrix  (const int &i,const int &j,const FEA_dMatrixT &block) 	{ SetBlock(i,j,block); 	}
		void Extract_SubMatrix (const int &i,const int &j,FEA_dMatrixT &block) 	{ CopyBlock(i,j,block); }
		void Add_SubMatrix 	(const int &i,const int &j,const FEA_dMatrixT &a) 				{ AddBlock(i,j,a); 			}
		void AddBlock				(const int &i,const int &j,const FEA_dMatrixT &block);
		void CopyBlock			(const int &i,const int &j, FEA_dMatrixT &block);
		void SetBlock				(const int &i,const int &j,const FEA_dMatrixT &block);
		void Swap_Rows			(const int &row1,const int &row2); 
		void Swap_Cols			(const int &col1,const int &col2); 


		void MultAB   		(const FEA_dMatrixT &a, const FEA_dMatrixT &b, 	int upper=0);
		void MultAB   		(const dMatrixT &a, 		const FEA_dMatrixT &b, 	int upper=0);
		void MultAB   		(const FEA_dMatrixT &a, const dMatrixT &b, 			int upper=0);
    	void MultATB  		(const FEA_dMatrixT &a, const FEA_dMatrixT &b, 	int upper=0);
		void MultABT  		(const FEA_dMatrixT &a, const FEA_dMatrixT &b, 	int upper=0);
		void MultATBT 		(const FEA_dMatrixT &a, const FEA_dMatrixT &b);
		
		void Outer 			(const FEA_dVectorT &a, const FEA_dVectorT &b);

   	// matrix-matrix-matrix operations 
		
		void MultABC  (const FEA_dMatrixT &a, const FEA_dMatrixT &b, const FEA_dMatrixT &c,
								   int range=dMatrixT::kWhole,int fillmode=dMatrixT::kOverwrite); 
	 	void MultABCT (const FEA_dMatrixT &a, const FEA_dMatrixT &b, const FEA_dMatrixT &c,
								   int range=dMatrixT::kWhole,int fillmode=dMatrixT::kOverwrite); 
	 	void MultATBC (const FEA_dMatrixT &a, const FEA_dMatrixT &b, const FEA_dMatrixT &c,
								   int range=dMatrixT::kWhole,int fillmode=dMatrixT::kOverwrite); 

		// overloaded operators   NOTE: no way to do C=A*B or C=A+B w/o an extra deep copy (slower)
	
		void operator  =  (const FEA_dMatrixT &a); 
		void operator +=  (const FEA_dMatrixT &a); 
		void operator -=  (const FEA_dMatrixT &a); 
		//void operator *=  (const FEA_dMatrixT &a); 
		// Not implemented because Tahoe's A *= B is a component by component multiplication 
	  	void operator  =  (const double &a); 
		void operator +=  (const double &a);   
		void operator -=  (const double &a); 
		void operator *=  (const double &a);   
		void operator /=  (const double &a); 
		/*void operator  =  (const double *a); 
		void operator +=  (const double *a);   
		void operator -=  (const double *a); 
		void operator *=  (const double *a);   
		void operator /=  (const double *a);  */

		void operator *=  (const FEA_dScalarT &s);   
		void operator /=  (const FEA_dScalarT &s);  

    	FEA_EquateT& operator () (const int p, const int q);

		// special operations
		
	  	/** rc:  FEA::kRow or FEA::kCol, it indicates whether a row or column is selected for dotting. */
		/** ij:  i or j. parameter is the "index" specifying which row i or column j is desired for dotting. */

	  	FEA_EquateT&       Dot (int rc,int ij,const FEA_dMatrixT &a,int a_rc,int a_ij); 
		FEA_EquateT&   Dot_Aij (int rc,int ij,const FEA_dMatrixT &a,int a_rc,int a_ij, FEA_dMatrixT &c,int i,int j); 
		FEA_EquateT&  ij_x_Aij (int i, int  j,const FEA_dMatrixT &a,int   ii,int  jj); 

		void ij_EQ_Akl_x_Bmn (int i,int j,FEA_dMatrixT &A,int k,int l,FEA_dMatrixT &B,int m,int n)
				{ for (int ip=0; ip<n_ip; ip++) (*this)[ip](i,j) = A[ip](k,l)*B[ip](m,n); }

		//protected:

    	int n_ip, n_rows, n_cols, n_rows_x_n_cols, n_ip_x_n_rows_x_n_cols;
		
		dArrayT  Block_Memory; // Can use in-lieu-of FEA_Pointer() 

};

}

#endif /* _FEA_DMATRIXT_H_ */
