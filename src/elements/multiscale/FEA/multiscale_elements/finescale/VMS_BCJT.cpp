//DEVELOPMENT

#include "FEA.h" 
#include "VMS.h" 


using namespace Tahoe;

VMS_BCJT::VMS_BCJT	(FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
						double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes,BCJ_Matl,np1,n,fdelta_t,Integration_Scheme);
}

/* destructor */
//VMS_BCJT::~VMS_BCJT(void)
//{
//}

//---------------------------------------------------------------------

void VMS_BCJT::Construct (FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
							double fdelta_t, int Integration_Scheme) 
{
	int NUM_TEMP_TERMS = 6;
	delta_t = fdelta_t;
	n_ip 		= np1.fVars[VMS::kGRAD_ua].IPs(); 
	n_rows	= np1.fVars[VMS::kGRAD_ua].Rows(); 
	n_cols	= np1.fVars[VMS::kGRAD_ua].Cols();
	n_en    = Shapes.dNdx.Cols();
  n_sd 		= n_rows;	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en = n_sd * n_en;
 
	Data_Pro.Construct ( Shapes.dNdx );
  
	Back_Stress = kNoBackStress;
	Time_Integration_Scheme = Integration_Scheme;

  Form_C_List			(BCJ_Matl);
  Form_A_S_Lists	(np1,n);
  Form_T4_List		();
  Form_B_List			();

	//cout <<"C Constants: "<<C<<"\n\n"; 
	//B.Print("B Matricies");
	//A.Print("A Matricies");
	//S.Print("S FEA Scalars");
	//T4.Print("T4 Matricies");

	Integral.Construct (Shapes.j, Shapes.W); 
/*
*/
}

//---------------------------------------------------------------------
/**  Form stiffness matricies k^Alpha and k^Beta 
 *   The ON swithes in this function are for Alpha and Beta 
 *   contributions respectively */

void VMS_BCJT::Form_LHS_Ka_Kb ( dMatrixT &Ka, dMatrixT &Kb )
{
 /* del(grad_wa) 		*/ 	Ka  = Integral.of( B[kB_1hat], B[kB00_2hat] );  	
												Kb  = Integral.of( B[kB_1hat], B[kB00_2hat] ); 

 /* del(Dm+1) 			*/ 	if (Time_Integration_Scheme != FEA::kBackward_Euler)
													Ka += Integral.of( B[kB_1hat], B[kB_IIA] );   	

 /* del(Dm) 				*/ 	if (Time_Integration_Scheme != FEA::kForward_Euler)
													Ka += Integral.of( B[kB_1hat], B[kB_IIB] ); 		

 /* del(DEV(S-Zeta) */	Ka += Integral.of( B[kB_1hat], delta_t, S[kBeta4b], T4[kN_1hat0], B[kB_H_bar_prime_a] ); 	
 												Kb += Integral.of( B[kB_1hat], delta_t, S[kBeta4b], T4[kN_1hat0], B[kB_H_bar_prime_b] ); 
 /* del(Kappa)  		*/ 	Ka += Integral.of( B[kB_1hat], delta_t*C[kBeta3]*0.5*C[kc_zeta]*C[kMu]*C[kh], T4[kN_1hat0], B[kB21_2hat] ); 
 /* del(N) 				  */	Ka += Integral.of( B[kB_1hat], delta_t, S[kAlpha], T4[kP_O_N1hat], B[kB_H_bar_prime_a] ); 						
												Kb += Integral.of( B[kB_1hat], delta_t, S[kAlpha], T4[kP_O_N1hat], B[kB_H_bar_prime_b] ); 
}

//---------------------------------------------------------------------
// F internal (F_int) dimensions here won't actually be in terms of Force

void VMS_BCJT::Form_RHS_F_int ( dArrayT &F_int ) // Untested
{
	FEA_dVectorT Da_mp1_vec	( n_ip, n_sd_x_n_sd ); 
	FEA_dVectorT Da_m_vec		( n_ip, n_sd_x_n_sd ); 
	FEA_dVectorT N_vec			( n_ip, n_sd_x_n_sd ); 

	if (Time_Integration_Scheme != FEA::kBackward_Euler) 
		Data_Pro.Reduce_Order		(	A[kDa_mp1], Da_mp1_vec	); 
	else
  	Data_Pro.Identity ( Da_mp1_vec );	

	if (Time_Integration_Scheme != FEA::kForward_Euler) 
		Data_Pro.Reduce_Order		(	A[kDa_m], 	Da_m_vec		); 
	else
  	Data_Pro.Identity ( Da_m_vec );	

	Data_Pro.Reduce_Order		(	A[kN], 			N_vec				); 

	F_int  = Integral.of		( B[kB_1hat], Da_mp1_vec	);  	
	F_int -= Integral.of		( B[kB_1hat], Da_m_vec		);  	
	F_int += Integral.of		( B[kB_1hat], delta_t, S[kAlpha], N_vec	);  	
}

//=== Private =========================================================

//##################################################################################
//################ B_TERMS #########################################################
//##################################################################################

void VMS_BCJT::Form_B_List (void)
{

		B.Construct (kNUM_B_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_en);  // B = B(9,24)	

    // Generic B 
		
	 	Data_Pro.grad_u   			( B[kB], FEA::kSymmetric 					); // <-- Not used in VMS_BCJ
	 	Data_Pro.grad_u       	( B[kB_1hat], FEA::kNonSymmetric 	); 

    // del (grad_wa)
	
	 	Data_Pro.grad_u_A    		( A[kXI], B[kB00_2hat] 	  ); 
	
		if (Time_Integration_Scheme != FEA::kBackward_Euler) {
							
			// del(Dm+1)
		
	 		Data_Pro.A_grad_u_T_B 		( A[kFbT], 		A[kDa_mp1],  	B[kB01_tau_3hat] 	);
	 		Data_Pro.A_grad_u_B 			( A[kA02], 		A[kA03],  		B[kB02_tau_3hat] 	);

	 		Data_Pro.A_grad_u_B 			( A[kA04], 		A[kA05],  		B[kB03_3hat] 			);
	 		Data_Pro.A_grad_u_B 			( A[kDa_mp1], A[kFb],   		B[kB04_3hat] 			);

			B[kB_IIA]  = B[kB01_tau_3hat]; 
			B[kB_IIA] += B[kB02_tau_3hat]; 
			B[kB_IIA] += B[kB03_3hat]; 
			B[kB_IIA] += B[kB04_3hat]; 

		}

		if (Time_Integration_Scheme != FEA::kForward_Euler) {

			// del(Dm) 
		
	 		Data_Pro.A_grad_u_T_B 		( A[kFbT], 		A[kDa_m],  		B[kB05_tau_3hat] 	);
	 		Data_Pro.A_grad_u_B 			( A[kDa_m], 	A[kFb],   		B[kB06_3hat] 			);

			B[kB_IIB]  = B[kB05_tau_3hat]; 
			B[kB_IIB] += B[kB06_3hat]; 

		}

   	//---------------- del( DEV(H) )  &  del( N )
		
		//---- Preliminary Computation 		
		
		// del ( Cb )     	Notes:  14. - 17. 
	
	 	Data_Pro.A_grad_u_B 	( A[kA1],   A[kFb],   B[kBa_cb_3hat] 			);
	 	Data_Pro.A_grad_u_B 	( A[kFb],  	A[kFb],	  B[kBb_cb_3hat] 			);
	 	Data_Pro.A_grad_u_B 	( A[kFbT],  A[kA1T],  B[kBa_cb_tau_3hat] 	);
	 	Data_Pro.A_grad_u_B 	( A[kFbT],  A[kFb],  	B[kBb_cb_tau_3hat] 	);
	  	
		// del ( Cb^-1 )		Notes:  18. - 21. 
	
	 	Data_Pro.A_grad_u_B 	( A[kA2],   A[kA3],   B[kBa_cbi_3hat] 			);
	 	Data_Pro.A_grad_u_B 	( A[kA3T],  A[kA3],	  B[kBb_cbi_3hat] 			);
	 	Data_Pro.A_grad_u_B 	( A[kA3T],  A[kA2T],  B[kBa_cbi_tau_3hat] 	);
	 	Data_Pro.A_grad_u_B 	( A[kA3T],  A[kA3],  	B[kBb_cbi_tau_3hat] 	);

	  // del ( S )				Notes:	12. - 13.
		
		B[kB_H_bar_a].MultAB( T4[kD],  B[kBa_cb_3hat] );
		B[kB_H_bar_b].MultAB( T4[kD],  B[kBb_cb_3hat] );

		// del ( Zeta )			

		if (Back_Stress==kSteinmann) {
		
			//-- Preliminary Computation 			Notes:  22. - 26. 
		
	 		Data_Pro.A_grad_u_B 	( A[kA10],  A[kR2],   B[kBa31_3hat] );
		 	Data_Pro.A_grad_u_B 	( A[kFbi],  A[kR2],	  B[kBb31_3hat] );

			B[kBa_sharp].MultAB   	( T4[kAa12T_2bar], 		B[kB_1hat] );
			B[kBa_sharp] 	*= C[kGamma5];    
			B[kBa_sharp] 	*= S[kJb];       

			B[kBb_sharp].MultAB   	( T4[kAb12T_2bar], 		B[kB_1hat] );
			B[kBb_sharp] 	*= C[kGamma5];
			B[kBb_sharp] 	*= S[kJb];

			B[kBrE].MultAB       		( T4[krET_2bar], 			B[kB_1hat] );
			B[kBrE] 			*= S[kJb];
			B[kBrE] 			*= C[kGamma5];

			//-- Main Summation 		Notes:	12. - 13.
			
			B[kB_Temp0].MultAB	(	T4[kZ_sharp],	B[kB_1hat]);
			B[kB_H_bar_a] += B[kB_Temp0]; 
			B[kB_H_bar_a] += B[kBa_sharp];
			B[kB_H_bar_a] += B[kBrE];
			B[kB_H_bar_a] *= -1.0; 
			B[kB_H_bar_a] += B[kBa31_3hat];; 

			B[kB_Temp0].MultAB	(	T4[kZI],	B[kB_1hat]);
			B[kB_H_bar_b] += B[kB_Temp0]; 
			B[kB_H_bar_b] += B[kBb_sharp];
			B[kB_H_bar_b] += B[kBrE];
			B[kB_H_bar_b] -= B[kBb31_3hat];; 

		}	

    //---  del ( DEV( H ) ) -- Main Routine : Modify later for del ( DEV (*) )  
		
		//	Notes:	line 9. term a.
		T4[kT4_Temp0]  = T4[kC_bib];		
		T4[kT4_Temp0] *= -C[kOneThird]; 
		T4[kT4_Temp0].PlusIdentity(); 
		B[kB_Temp1].MultAB( T4[kT4_Temp0], B[kB_H_bar_a] );

		// Notes: line 9. term b.	
	  B[kB_Temp2].SumOf  ( 	B[kBa_cb_3hat], B[kBa_cb_tau_3hat]	) ;		
		B[kB_Temp2] *= 			 -C[kOneThird]; 
		B[kB_Temp3].MultAB ( 	T4[kC_biH], B[kB_Temp2]  );

		// Notes: line 9. term c.	
	  B[kB_Temp4].SumOf( 		B[kBa_cbi_3hat], B[kBa_cbi_tau_3hat]	) ;	
		B[kB_Temp4] *=  			S[kBeta_HC]; 
		B[kB_Temp4] *= 		 	 -C[kOneThird]; 

		// 	Notes: 	summation of line 9. terms
		B[kB_H_bar_prime_a]  = B[kB_Temp1];  
		B[kB_H_bar_prime_a] += B[kB_Temp3];
		B[kB_H_bar_prime_a] += B[kB_Temp4];

		//	Notes:	line 10. term a.
		T4[kT4_Temp0]  = 	T4[kC_bib];
		T4[kT4_Temp0] *=  -C[kOneThird]; 
		T4[kT4_Temp0].PlusIdentity(); 
		B[kB_Temp1].MultAB	( T4[kT4_Temp0], B[kB_H_bar_b] );

	 
		//	Notes:	line 10. term b.
	  B[kB_Temp2].SumOf  ( B[kBb_cb_3hat], B[kBb_cb_tau_3hat]	) ;	
		B[kB_Temp2] *= 			-C[kOneThird]; 
		B[kB_Temp3].MultAB ( T4[kC_biH], B[kB_Temp2] );

		//	Notes:	line 10. term c.
	  B[kB_Temp4].SumOf( 	B[kBb_cbi_3hat], B[kBb_cbi_tau_3hat]	) ;	
		B[kB_Temp4] *=  		S[kBeta_HC]; 
		B[kB_Temp4] *= 		 -C[kOneThird]; 

		// 	Notes: 	summation of line 10. terms
		B[kB_H_bar_prime_b]  = B[kB_Temp1];
		B[kB_H_bar_prime_b] += B[kB_Temp3];
		B[kB_H_bar_prime_b] += B[kB_Temp4];

		// del ( Kappa )

	 	Data_Pro.grad_u_A  ( A[kF], B[kB21_2hat] 	  ); 


}
						
//##################################################################################
//################ A_TERMS #########################################################
//##################################################################################

//NOTE: np1 := "n+1" time step; n := "n" time step; npt := "n+theta" time step
//      *** No subscript implies n+theta time step ***

void VMS_BCJT::Form_A_S_Lists (VMS_VariableT &np1,VMS_VariableT &n)
{
	//---- Developer cheat: put npt in function door in-lieu-of np1 or np for speed

	// cout << "FLAG 00 \n";

	//n.Print("n");

VMS_VariableT npt(	n.Get(VMS::kGRAD_ua), n.Get(VMS::kGRAD_ub)	);  

	// cout << "FLAG 0 \n";

	if 			(	Time_Integration_Scheme == FEA::kForward_Euler		)		npt = n;
	else if (	Time_Integration_Scheme == FEA::kBackward_Euler		)		npt = np1;
	else if (	Time_Integration_Scheme == FEA::kCrank_Nicholson	) { npt.SumOf(np1,n); npt *= 0.5; }
	else 	cout << " ...ERROR >> VMS_BCJT::Form_A_List() : Bad theta value for time stepping \n";
		
	//-----

  // npt.Print("npt");

	// cout << "FLAG 1 \n";

	A.Construct ( kNUM_A_TERMS, n_ip, n_sd, n_sd );	
	S.Construct ( kNUM_S_TERMS, n_ip );

	//-----
	
  A[kF]   			= npt.Get (	VMS::kF					);	 // NOTE: kF != VMS::kF 	
  A[kFa]  			= npt.Get (	VMS::kFa				);	 	
	A[kFai] 			= npt.Get (	VMS::kFai				);
  A[kFb]  			= npt.Get (	VMS::kFb				);	 	
  A[kFbi] 			= npt.Get (	VMS::kFbi				);	 	
	A[kgrad_ub] 	= npt.Get (	VMS::kgrad_ub		);

	A[kFbT].Transpose  			( A[kFb] );
	A[kCb].MultATB     			( A[kFb],  			A[kFb] 	);					
	A[kCbi].MultABT    			( A[kFbi], 			A[kFbi] );					
	A[kEb]  = A[kCb]; 
	A[kEb].PlusIdentity(-1.0); 
	A[kEb] *= 0.5; 
	A[kF_sharp].MultAB  		( A[kgrad_ub], 	A[kFb] 	);
	A[kF_sharpT].Transpose	( A[kF_sharp] );

	// B-Euler doesn't need following data since Da_mp1 = .5I thus del(Da_mp1) = 0
  if ( Time_Integration_Scheme != FEA::kBackward_Euler ) { 

		A[kFa_np1]  = np1.Get 	( VMS::kFa );
		A[kCa_np1].MultATB  		( A[kFa_np1], A[kFa_np1] 	);
		A[kDa_mp1].MultATBC			(	A[kFai], 		A[kCa_np1], 	 A[kFai]	);
		A[kDa_mp1] *= 0.5;
	
  	A[kA02].MultATBT 				(	A[kFai],		A[kFa_np1]	);
  	A[kA02] 	*= 0.5; 
  	A[kA03].MultAB   				(	A[kFa_np1], A[kFai]			);
  	A[kA04].MultATBT 				(	A[kFai], 		A[kFa_np1]	);
  	A[kA04] 	*= 0.5;
		A[kF_np1]  = np1.Get 		( VMS::kF );
  	A[kA05].MultAB   				(	A[kF_np1], 	A[kFai]			);
    if ( Time_Integration_Scheme == FEA::kForward_Euler ) { // Trapezoidal not yet available
			A[kDa_m] = 0.0;
			A[kDa_m].PlusIdentity(0.5); 
		}

	}

	// F-Euler doesn't need following data since Da_m = .5I thus del(Da_m) = 0
  if ( Time_Integration_Scheme != FEA::kForward_Euler ) { 

		A[kFa_n]    = n.Get			(	VMS::kFa );
		A[kCa_n].MultATB				( A[kFa_n], A[kFa_n]    );
		A[kDa_m].MultATBC 			(	A[kFai], 		A[kCa_n], 		A[kFai]	);
		A[kDa_m] 	*= 0.5;
    if ( Time_Integration_Scheme == FEA::kBackward_Euler ) { // Trapezoidal not yet available
			A[kDa_mp1] = 0.0;
			A[kDa_mp1].PlusIdentity(0.5); 
		}
	}	

  A[kA1].MultAB						(	A[kFbT], 		A[kF_sharp] 	);
  A[kA2].MultAB						( A[kCbi], 		A[kA1]	 			);
  A[kA3].MultAB						( A[kFb], 		A[kCbi] 			);
  A[kA1T].Transpose 			( A[kA1] ); 
  A[kA2T].Transpose 			( A[kA2] ); 
  A[kA3T].Transpose 			( A[kA3] ); 

  A[kA10].MultAB					( A[kFbi], 			A[kF_sharp] );

	//---- sE, Zeta, S, H, N

  Data_Pro.C_IJKL_E_KL 	( C[kMu], 	C[kLamda], A[kEb], A[kS] );	// Doesn't multiply excessive zeros 

	if (Back_Stress==kNoBackStress) 	
  	A[kH] = A[kS]; 

	if (Back_Stress==kSteinmann) {	
		// Calculate_sE    	  	( A[k_grad_ub], A[ksE] 			); // ##### Needs to be implemented 
  	A[kFbi_1hat]  = A[ksE]; 
  	A[kFbi_1hat] *= -1.0; 
  	A[kFbi_1hat].PlusIdentity(); 
		Data_Pro.Curl					( A[ksE],       A[kCurl_sE] );

		A[kR2].Transpose			( A[kCurl_sE] 	); 
		A[kR2] *= C[kGamma5];
		A[kR2] *= S[kJb]; 

		A[kZeta].MultAB				( A[kFbi],	A[kR2] );
  	A[kH].DiffOf     			( A[kS], 		A[kZeta] );
	}	
	
	S[kBeta_SCb].Double_Dot ( A[kS], A[kCb] ); 
	
	A[kDEV_H]  = 	A[kCbi];
	A[kDEV_H] *= 	S[kBeta_SCb]; 
	A[kDEV_H] *= -C[kOneThird]; 
	A[kDEV_H] += 	A[kH]; 

	A[kDEV_H].Mag_and_Dir 	( S[kMag_DEV_H], A[kN] );
  A[kN_1hat] = A[kN]; 

	//---- Betas and Alpha (i.e. S[i] )

	S[kKappa] = 1.0;

  S[kBeta]  = S[kMag_DEV_H]; 
  S[kBeta] *= C[kRoot3by2]; 
  S[kBeta] -= S[kKappa]; 
  S[kBeta] -= C[kY]; 
  S[kBeta] /= C[kV]; 

	S[kBeta4].Cosh ( S[kBeta] );	 
	//S[kBeta4].Macaulay_Bracket();  	// Leave unimplemented during testing
	S[kBeta4]  *= C[kGamma1];

	S[kAlpha].Sinh ( S[kBeta] );		
	//S[kAlpha].Macaulay_Bracket();		// Leave unimplemented during testing
	S[kAlpha]  *= C[kGamma1];

	S[kBeta3]   = S[kBeta4];
	S[kBeta3]  *= C[kGamma3p1];

	S[kBeta4b]  = S[kBeta4];
	S[kBeta4b] *= C[kGamma2];

	S[kBeta_HC].Double_Dot ( A[kH], A[kCb] );

  S[kJ].Determinant	  ( A[kF]		);  
	S[kJa].Determinant	( A[kFa]	); 
	S[kJb].Determinant	( A[kFb] 	); 

  //--- XI
	
  A[kXI]  =  A[kN_1hat];
  A[kXI] *=  S[kAlpha];
  A[kXI] *=  delta_t; 
  A[kXI] +=  A[kDa_mp1];
 	A[kXI] -=  A[kDa_m];

}

//##################################################################################
//################ 4th Order Tensor Terms (i.e. T4[i] ) ########################
//##################################################################################

void VMS_BCJT::Form_T4_List (void)  // These matricies are all 9x9 (not 3x3 like A)
{
	
	T4.Construct ( kNUM_T4_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_sd );

  Data_Pro.C_IJKL (	C[kLamda],  C[kMu],					T4[kD]				); 
  Data_Pro.A_o_B 	( A[kCbi], 		A[kCb], 				T4[kC_bib] 		);
  Data_Pro.A_o_B 	( A[kCbi], 		A[kH], 					T4[kC_biH] 		);
  Data_Pro.A_o_B 	( A[kZeta], 	A[kF_sharpT], 	T4[kZ_sharp] 	);
  Data_Pro.A_o_1 	( A[kZeta], 									T4[kZI] 			);
  Data_Pro.A_o_B 	( A[kN_1hat],	A[kN], 					T4[kN_1hat0] 	);

  Data_Pro.II_minus_A_o_B ( A[kN_1hat],	A[kN_1hat],	T4[kP_O_N1hat] ); // 4th Order Orthogonal Projector 

	//--------- Problem Specific Data Processing Routines
	
  /** For the following methods:
   * bb is Black Bold Font (4th Order Tensor) bbA_2bar_ijsr 
	 * where E is permutation symbol. See notes, terms with
	 * 4 free indicies are repacked into a 4th order tensor, then
	 * an order reduction (via mapping) is performed. rE is Roman Epsilon. */
	
 	/** A2ix*Ajs,m*Exrm*Bsr = bbA_2bar_ijsr*Bsr  |--> A_2bar_pq*Bq */
  //Data_Pro.bbA_2bar  ( A[kFbi], A[ksE],     	A[kAa12T_2bar] );     // ##### Needs to be implemented 
  //Data_Pro.bbA_2bar  ( A[kFbi], A[kFbi_1hat]  A[kAb12T_2bar] ); 		// ##### Needs to be implemented 
	
  // LATEST:: Data_Pro.Aikm_Elmj_to_Aijkl_to_AIJ  ( A[kFbi], A[kFbi_1hat]  A[kAb12T_2bar] ); 		// ##### Needs to be implemented 
	// Best if implement 3rd and 4th order tensors 
		
  /** A2ix*Ajm,s*Exmr*Bsr |--> bbAe_2bar_ijsr*Bsr |--> Ae_2bar_pq*Bq */
  //Data_Pro.bbAe_2bar ( A[kFbi], A[ksE]   			A[krET_2bar] );				// ##### Needs to be implemented 

}


//##################################################################################
//################ Material Constant Terms C[i] #####################################
//##################################################################################

void VMS_BCJT::Form_C_List (VMF_MaterialT *BCJ_Matl)
{

	C.Dimension (kNUM_C_TERMS);
	
	C[kLamda]    	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kLamda 	);
	C[kMu]    	 	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kMu 		);
	C[kl]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kl 			);
	C[kc_zeta]    = BCJ_Matl -> Retrieve ( BCJ_MatlT::kc_zeta );
	C[kh]    			= BCJ_Matl -> Retrieve ( BCJ_MatlT::kh 			);
  C[kf]        	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kf 			);
  C[kY]        	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kY 			);
  C[kV]        	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kV 			);

	C[kOneThird]  = 1.0/3.0;
	C[kRoot3by2]  = sqrt(1.5);
	C[kGamma1]   	= C[kRoot3by2]*C[kf];
	C[kGamma2]   	= C[kRoot3by2]/C[kV];
	C[kGamma3p1] 	= 1.0/C[kV];
	C[kGamma3p5] 	= C[kY]/C[kV];

	C[kGamma5] 	= C[kl]*C[kc_zeta]*C[kMu];

}


