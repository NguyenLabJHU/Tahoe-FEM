// $Id: FEA_Data_Pro_Diagnostics.cpp,v 1.5 2003-04-23 23:34:21 creigh Exp $
#include "FEA.h"  

using namespace Tahoe;

//---------------------------------------------------------------------

void FEA_Data_ProcessorT::sym_grad_u(FEA_dMatrixT &B) // B, dim  3x8 || 6x24
{
   // Insert Ba into B  [See (3.10.5), (3.10.6)]
   // This is Ba at the lth ip 
				
int j0,j1,j2; 
			
//if (!B.IPs()) // Doesn't work !
//  B.FEA_Dimension(n_ip, (n_sd==2) ? 3 : 6	, n_sd_x_n_en);

if (n_sd==2) {

   for (int a=0; a<n_en; a++) {  
					 
		 j0 = n_sd*a;     j1 = n_sd*a +1;  // Columns of B for node "a"

     B(0,j0) = dN(dx1,a);   B(0,j1) = 0.0;           // u1,1 
		 B(1,j0) = 0.0;     		B(1,j1) = dN(dx2,a);     // u2,2
		 B(2,j0) = dN(dx2,a);   B(2,j1) = dN(dx1,a);     // u1,2 + u2,1

      // u1                    u 2         
   }
}

else if (n_sd==3) {

   for (int a=0; a<n_en; a++) {   // B = [[B0][B1]...[j0j1j2]...[Bnen-1]]

		 j0 = n_sd*a;     j1 = n_sd*a +1;     j2 = n_sd*a +2; // Columns of B for node "a"

   }
}

else
	cout << "...ERROR >> FEA_Data_ProcessorT::sym_grad_u : Bad n_sd \n";

}  

//---------------------------------------------------------------------

void FEA_Data_ProcessorT::grad_u_T (FEA_dMatrixT &B) // B_1hat, dim 4x8 || 9x24, grad_u_T (non-sym)
{
					
int j0,j1,j2;				

if (n_sd==2) {

   for (int a=0; a<n_en; a++) {

		 j0 = n_sd*a;     j1 = n_sd*a +1;  // Columns of B for node "a"

     B(0,j0) = dN(dx1,a);   B(0,j1) = 0.0;           // u1,1  
		 B(1,j0) = 0.0;     		B(1,j1) = dN(dx2,a);     // u2,2 
		 B(2,j0) = 0.0;     		B(2,j1) = dN(dx1,a);     // u2,1
		 B(3,j0) = dN(dx2,a);   B(3,j1) = 0.0;           // u1,2

     // u 1                    u 2             
   }
}

else if (n_sd==3) {

   for (int a=0; a<n_en; a++) {

		 j0 = n_sd*a;     j1 = n_sd*a +1;     j2 = n_sd*a +2; // Columns of B for node a

     B(0,j0) = dN(dx1,a);   B(0,j1) = 0.0;            B(0,j2) = 0.0;           // u1,1  -> u0,0
		 B(1,j0) = 0.0;     		B(1,j1) = dN(dx2,a);			B(1,j2) = 0.0;           // u2,2  -> u1,1
		 B(2,j0) = 0.0;     		B(2,j1) = 0.0;            B(2,j2) = dN(dx3,a);		 // u3,3  -> u2,2

		 B(3,j0) = 0.0;     		B(3,j1) = 0.0;				    B(3,j2) = dN(dx2,a);		 // u3,2  -> u2,1
		 B(4,j0) = 0.0;   			B(4,j1) = 0.0;      			B(4,j2) = dN(dx1,a);		 // u3,1  -> u2,0
		 B(5,j0) = 0.0;				  B(5,j1) = dN(dx1,a);			B(5,j2) = 0.0;           // u2,1  -> u1,0

		 B(6,j0) = 0.0;     		B(6,j1) = dN(dx3,a);    	B(6,j2) = 0.0;      		 // u2,3  -> u1,2
		 B(7,j0) = dN(dx3,a);   B(7,j1) = 0.0;      			B(7,j2) = 0.0;					 // u1,3  -> u0,2
		 B(8,j0) = dN(dx2,a);   B(8,j1) = 0.0;						B(8,j2) = 0.0;           // u1,2  -> u0,1

     // u 1                    u 2                      u 3
   }
}

else
	cout << "...ERROR >> FEA_Data_ProcessorT::grad_u : Bad n_sd \n";
}  

//---------------------------------------------------------------------

void FEA_Data_ProcessorT::A2_grad_u (FEA_dMatrixT &A2, FEA_dMatrixT &B) 
{
					
int j0,j1,j2;				

if (n_sd==2) {

   for (int a=0; a<n_en; a++) {

		 j0 = n_sd*a;     j1 = n_sd*a +1;  // Note: Can also use B(0,j0).Aij_EQ_Bkl_x_Cmn(dN,dx1,a, A2,k1,k1) : Slightly faster 

     B(0,j0) = dN.ij_x_Aij(dx1,a, A2,k1,k1);   		B(0,j1) = dN.ij_x_Aij(dx1,a, A2,k1,k2);		// u1,1 
     B(1,j0) = dN.ij_x_Aij(dx2,a, A2,k2,k1);   		B(1,j1) = dN.ij_x_Aij(dx2,a, A2,k2,k2);		// u2,2

     B(2,j0) = dN.ij_x_Aij(dx2,a, A2,k1,k1);   		B(2,j1) = dN.ij_x_Aij(dx2,a, A2,k1,k2);		// u1,2
     B(3,j0) = dN.ij_x_Aij(dx1,a, A2,k2,k1);   		B(3,j1) = dN.ij_x_Aij(dx1,a, A2,k2,k2);		// u2,1

   }
}

else if (n_sd==3) {

   for (int a=0; a<n_en; a++) {

		 j0 = n_sd*a;     j1 = n_sd*a +1;     j2 = n_sd*a +2;

     B(0,j0) = dN.ij_x_Aij(dx1,a, A2,k1,k1); 	B(0,j1) = dN.ij_x_Aij(dx1,a, A2,k1,k2);   B(0,j2) = dN.ij_x_Aij(dx1,a, A2,k1,k3);   // u1,1      
     B(1,j0) = dN.ij_x_Aij(dx2,a, A2,k2,k1); 	B(1,j1) = dN.ij_x_Aij(dx2,a, A2,k2,k2);   B(1,j2) = dN.ij_x_Aij(dx2,a, A2,k2,k3);   // u2,2      
     B(2,j0) = dN.ij_x_Aij(dx3,a, A2,k3,k1); 	B(2,j1) = dN.ij_x_Aij(dx3,a, A2,k3,k2);   B(2,j2) = dN.ij_x_Aij(dx3,a, A2,k3,k3);   // u3,3      

     B(3,j0) = dN.ij_x_Aij(dx3,a, A2,k2,k1);  B(3,j1) = dN.ij_x_Aij(dx3,a, A2,k2,k2);   B(3,j2) = dN.ij_x_Aij(dx3,a, A2,k2,k3);   // u2,3      
     B(4,j0) = dN.ij_x_Aij(dx3,a, A2,k1,k1);  B(4,j1) = dN.ij_x_Aij(dx3,a, A2,k1,k2);   B(4,j2) = dN.ij_x_Aij(dx3,a, A2,k1,k3);   // u1,3      
     B(5,j0) = dN.ij_x_Aij(dx2,a, A2,k1,k1);  B(5,j1) = dN.ij_x_Aij(dx2,a, A2,k1,k2);   B(5,j2) = dN.ij_x_Aij(dx2,a, A2,k1,k3);   // u1,2      

     B(6,j0) = dN.ij_x_Aij(dx2,a, A2,k3,k1);  B(6,j1) = dN.ij_x_Aij(dx2,a, A2,k3,k2);   B(6,j2) = dN.ij_x_Aij(dx2,a, A2,k3,k3);   // u3,2      
     B(7,j0) = dN.ij_x_Aij(dx1,a, A2,k3,k1);  B(7,j1) = dN.ij_x_Aij(dx1,a, A2,k3,k2);   B(7,j2) = dN.ij_x_Aij(dx1,a, A2,k3,k3);   // u3,1      
     B(8,j0) = dN.ij_x_Aij(dx1,a, A2,k2,k1);  B(8,j1) = dN.ij_x_Aij(dx1,a, A2,k2,k2);   B(8,j2) = dN.ij_x_Aij(dx1,a, A2,k2,k3);   // u2,1      

   }
}

else
	cout << "...ERROR >> FEA_Data_ProcessorT::grad_u : Bad n_sd \n";
}  

//---------------------------------------------------------------------

void FEA_Data_ProcessorT::A2_grad_u_A1(FEA_dMatrixT &A2, FEA_dMatrixT &A1, FEA_dMatrixT &B) 
{
					
int j0,j1,j2;				

if (n_sd==2) {

 	for (int a=0; a<n_en; a++) {

  	j0 = n_sd*a;     j1 = n_sd*a +1; 

  	B(0,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k1,k1);       B(0,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k1,k2);    
  	B(1,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k2,k1);       B(1,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k2,k2);    

  	B(2,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k1,k1);       B(2,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k1,k2);    
  	B(3,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k2,k1);       B(3,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k2,k2);    

   }
}

else if (n_sd==3) {

 for (int a=0; a<n_en; a++) {

	j0 = n_sd*a;     j1 = n_sd*a +1;     j2 = n_sd*a +2;
		 
  B(0,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k1,k1);  B(0,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k1,k2);  B(0,j2)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k1,k3);
  B(1,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k2,k1);  B(1,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k2,k2);  B(1,j2)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k2,k3);
  B(2,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k3, A2,k3,k1);  B(2,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k3, A2,k3,k2);  B(2,j2)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k3, A2,k3,k3);
 
  B(3,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k3, A2,k2,k1);  B(3,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k3, A2,k2,k2);  B(3,j2)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k3, A2,k2,k3);
  B(4,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k3, A2,k1,k1);  B(4,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k3, A2,k1,k2);  B(4,j2)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k3, A2,k1,k3);
  B(5,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k1,k1);  B(5,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k1,k2);  B(5,j2)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k1,k3);

	B(6,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k3,k1);  B(6,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k3,k2);  B(6,j2)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k2, A2,k3,k3);
  B(7,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k3,k1);  B(7,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k3,k2);  B(7,j2)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k3,k3);
  B(8,j0)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k2,k1);  B(8,j1)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k2,k2);  B(8,j2)= dN.Dot_Aij(FEA::kCol,a,A1,FEA::kCol,k1, A2,k2,k3);

 }

}

else
	cout << "...ERROR >> FEA_Data_ProcessorT::grad_u : Bad n_sd \n";
} 

//---------------------------------------------------------------------

void FEA_Data_ProcessorT::A2_grad_u_T_A1_Alt (FEA_dMatrixT &A2, FEA_dMatrixT &A1, FEA_dMatrixT &B) 
{
					
	int a,i,j; 
	FEA_dMatrixT Ba_tau_2hat	(n_ip, n_sd_x_n_sd, n_sd);
	FEA_dMatrixT Ba_tau_3hat	(n_ip, n_sd_x_n_sd, n_sd);
	FEA_dMatrixT A1T 					(n_ip, n_sd, n_sd);

	A1T.Transpose ( A1 );
	A_grad_u_T ( A2, B ); // Is actually now B_tau_2hat

 	for (a=0; a<n_en; a++) {

  	i = k1;  j = n_sd*a;     
 
		B.Extract_SubMatrix ( i,j, Ba_tau_2hat ); // Time Saver, write: B.Mult_SubMatrix_by ( i,j, Ba_tau_2hat )
		Ba_tau_3hat.MultAB 	( Ba_tau_2hat, A1T );
		B.Insert_SubMatrix 	( i,j, Ba_tau_3hat ); 

	}
// returns B_tau_3hat
}  

//---------------------------------------------------------------------

void FEA_Data_ProcessorT::c_ijkl_Alt (double &lamda,double &mu, FEA_dScalarT &J, FEA_dMatrixT &F, FEA_dMatrixT &D) 
{
  int m = (n_sd-1)*3;	
	FEA_dMatrixT CC 		(n_ip, n_sd_x_n_sd, n_sd_x_n_sd);
	FEA_dMatrixT CC_sym (n_ip, (n_sd-1)*3, (n_sd-1)*3 );  // <-- symetric reduction map here (n_sd-1)*3

  C_IJKL	(	lamda,mu, CC );
	CC.Extract_SubMatrix ( 0,0, CC_sym ); // <-- Must reduce to smaller size common-place for sym matricies

	TensorTransformT TTrans	(n_sd);

	for (int l=0; l<n_ip; l++) 
		D[l] = TTrans.PushForward ( F[l], CC_sym[l] ); 

	D /= J;

}

/*
 
//---------------------------------------------------------------------
// 4th Order Elasticity Tensor CC 

void FEA_Data_ProcessorT::C_IJKL_O4	(const double &lamda,const double &mu,FEA_dTensorO4T &CC) 
{
	FEA_dTensorO4T II(n_ip,n_sd);
	FEA_dTensorO4T I_o_I(n_ip,n_sd);
  II.II();	
	I_o_I.I_o_I(); 

	I_o_I *= lamda;
	II    *= 2.0*mu; 

	CC.SumOf(I_o_I,II);

}

//---------------------------------------------------------------------
// FUTURE DEVELOPMENT: Uses notion of FEA_dTensorO4
// 4th Order modlus cc or cijkl (Reduced to 2nd Order)
// F_o_F_o_F_o_F::CC |---> [cc]^matrix 

void FEA_Data_ProcessorT::c_ijkl (FEA_dMatrixT &F, FEA_dTensorO4T &CC, FEA_dMatrixT &cc) 
{
	for (i=0; i<n_sd; i++)
		for (j=0; j<n_sd; j++)
			for (k=0; k<n_sd; k++)
				for (l=0; l<n_sd; l++) {
          cc = 0.0; // Summations
					for (I=0; I<n_sd; I++)
						for (J=0; J<n_sd; J++)
							for (K=0; K<n_sd; K++)
								for (L=0; L<n_sd; L++)
									cc ( Map(i,j), Map(k,l) ) += F(i,I)*F(j,J)*F(k,K)*F(l,L)*CC(I,J,K,L);
				}

};

*/


