// $Id: FEA_Data_ProcessorT.h,v 1.3 2003-02-03 04:40:22 paklein Exp $
#ifndef _FEA_DATAPROCESSRT_H_
#define _FEA_DATAPROCESSRT_H_

namespace Tahoe {

class FEA_Data_ProcessorT  
{
	public:

		enum Spatial_DirectionT { dx1, dx2, dx3 };
		enum NumberT { k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 }; // <-- numbers are one less then face value
					
		FEA_Data_ProcessorT 			(); 
		FEA_Data_ProcessorT 			( FEA_dMatrixT &fdNdx );
		void Construct 						( FEA_dMatrixT &fdNdx ); 

		void grad_u     			(FEA_dMatrixT &B,int T_flag=FEA::kNonSymmetric);

		void grad_u_A   			(FEA_dMatrixT &A, 	FEA_dMatrixT &B);
		void A_grad_u_T  		  (FEA_dMatrixT &A, 	FEA_dMatrixT &B);

		void A_grad_u 				(FEA_dMatrixT &A, 	FEA_dMatrixT &B,int T_flag=FEA::kNonSymmetric);
		void grad_u_T_A 			(FEA_dMatrixT &A, 	FEA_dMatrixT &B) { A_grad_u (A,B,FEA::kNonSymTranspose); }


		void A_grad_u_B 			(FEA_dMatrixT &A, 	FEA_dMatrixT &B, FEA_dMatrixT &B_3hat) { AikBlj_Ukl(A,B,B_3hat); }
		void A_grad_u_T_B 		(FEA_dMatrixT &A, 	FEA_dMatrixT &B, FEA_dMatrixT &B_3hat) { AikBlj_Ukl(A,B,B_3hat, FEA::kNonSymTranspose); }
		void AikBlj_Ukl 			(FEA_dMatrixT &A, 	FEA_dMatrixT &B, FEA_dMatrixT &B_3hat, int T_flag=FEA::kNonSymmetric); 


		void A_o_B_grad_u 		(FEA_dMatrixT &A, 	FEA_dMatrixT &B, FEA_dMatrixT &B_sharp); 
		void A_o_B 						(FEA_dMatrixT &A, 	FEA_dMatrixT &B, FEA_dMatrixT &C); 
		void A_o_1 						(FEA_dMatrixT &A, 	FEA_dMatrixT &C) { A_o_B (A,I,C); }  // Slow 0's of I

		void II_minus_A_o_B		(FEA_dMatrixT &A, 	FEA_dMatrixT &A1, 	FEA_dMatrixT &P); 
		void IIsym						(FEA_dMatrixT &II);
		void Identity					(FEA_dMatrixT &I_mat) { I_mat=I; }
		void Identity					(FEA_dVectorT &I_vec, double scale=1.0); 
	 	
	  void C_IJKL 					(const double &lamda,const double &mu,FEA_dMatrixT &D,int kine=FEA::kPlaneStrain);
		void C_IJKL_E_KL			(double &lamda,double &mu, FEA_dMatrixT &E, FEA_dMatrixT &S); // Hooke's Law
		void c_ijkl						(double &lamda,double &mu, FEA_dScalarT &J, FEA_dMatrixT &F, FEA_dMatrixT &D);
		void c_ijkl_Alt				(double &lamda,double &mu, FEA_dScalarT &J, FEA_dMatrixT &F, FEA_dMatrixT &D);

		void Curl			    		(const ArrayT<dMatrixT> &T,	FEA_dMatrixT &curl) const;
		//void Curl			    	(FEA_dVectorT &a2,	FEA_dVectorT &c);
		//void A2_o_A1_grad_u (FEA_dMatrixT &A2, 	FEA_dMatrixT &A1, 	FEA_dMatrixT &B);
	
		void Reduce_Order			(	FEA_dMatrixT &A, 	FEA_dVectorT &a );
		//void Reduce_Order			(	FEA_dTensorO4T &AA, 	FEA_dMatrixT &A );


		void Form_Order_Reduction_Map(void);
    nMatrixT<int> Map;		
	  FEA_dMatrixT dN;	
	  FEA_dMatrixT 	I;	
		int n_ip, n_en, n_sd, n_sd_x_n_sd, n_sd_x_n_en;

		//--------------------------------------------------------------------------------------------------------
		//-- Alternative methods, often slower Methods Due to use of I and multiplication by 0's ... good for diagnostics
		//   These methods were original formulations, replaced by faster/ more streamlined ones/
		
		void A_grad_U_Alt 		(FEA_dMatrixT &A, 	FEA_dMatrixT &B_1bar) { AikBlj_Ukl ( A,I, B_1bar); } 
		void A_grad_U_T_Alt 	(FEA_dMatrixT &A, 	FEA_dMatrixT &B_1bar) { AikBlj_Ukl ( A,I, B_1bar, FEA::kNonSymTranspose); } 
		void grad_U_A_Alt 		(FEA_dMatrixT &A, 	FEA_dMatrixT &B_2hat) { AikBlj_Ukl ( I,A,	B_2hat); } 
		void grad_U_T_A_Alt 	(FEA_dMatrixT &A, 	FEA_dMatrixT &B_2hat) { AikBlj_Ukl ( I,A,	B_2hat, FEA::kNonSymTranspose); } 

		void sym_grad_u 				(FEA_dMatrixT &B);
		void grad_u_T    				(FEA_dMatrixT &B);
		void A2_grad_u_A1 			(FEA_dMatrixT &A2, 	FEA_dMatrixT &A1, 	FEA_dMatrixT &B);
		void A2_grad_u_T_A1			(FEA_dMatrixT &A2, 	FEA_dMatrixT &A1, 	FEA_dMatrixT &B);
		void A2_grad_u_T_A1_Alt	(FEA_dMatrixT &A2, 	FEA_dMatrixT &A1, 	FEA_dMatrixT &B);
		void A2_grad_u 					(FEA_dMatrixT &A, 	FEA_dMatrixT &B);
};

}

// inline routines go here ...

#endif
