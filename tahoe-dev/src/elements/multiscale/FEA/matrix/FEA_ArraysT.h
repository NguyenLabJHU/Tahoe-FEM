// $Id: FEA_ArraysT.h,v 1.2 2003-02-03 04:40:24 paklein Exp $
#ifndef _FEA_ARRAY_T__ 
#define _FEA_ARRAY_T__ 

namespace Tahoe {

//-------------------------------------------------
				
class FEA_dMatrix_ArrayT : public ArrayT <FEA_dMatrixT>  
{

	public:

    FEA_dMatrix_ArrayT 	(int n_mat,int n_ip,int n_rows,int n_cols); 
 		void Construct			(int n_mat,int n_ip,int n_rows,int n_cols);

    FEA_dMatrix_ArrayT 	(void) : ArrayT <FEA_dMatrixT> () { }
 		void Print  (void);
		void Print  (char*);
  
};

//-------------------------------------------------
				
class FEA_dVector_ArrayT : public ArrayT <FEA_dVectorT>  
{

	public:

    FEA_dVector_ArrayT 	(int n_mat,int n_ip,int n_rows); 
 		void Construct			(int n_mat,int n_ip,int n_rows);
		 
    FEA_dVector_ArrayT (void) : ArrayT <FEA_dVectorT> () { }
 		void Print  (void);
		void Print  (char*);

};

//-------------------------------------------------
				
class FEA_dScalar_ArrayT : public ArrayT <FEA_dScalarT>  
{

	public:

    FEA_dScalar_ArrayT 	(int n_mat,int n_ip); 
 		void Construct			(int n_mat,int n_ip);

    FEA_dScalar_ArrayT (void) : ArrayT <FEA_dScalarT> () { }
 		void Print  (void);
		void Print  (char*);
  
};

//-------------------------------------------------

} // namespace Tahoe 
#endif 

