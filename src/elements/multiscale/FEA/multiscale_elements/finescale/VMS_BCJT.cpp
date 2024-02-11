// $Id: VMS_BCJT.cpp,v 1.17 2003-09-16 16:41:26 raregue Exp $
#include "FEA.h" 
#include "VMS.h" 

using namespace Tahoe;

VMS_BCJT::VMS_BCJT	(FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
						int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes,BCJ_Matl,np1,n,fTime_Step,fdelta_t,Integration_Scheme);
}

/* destructor */
//VMS_BCJT::~VMS_BCJT(void)
//{
//}

//---------------------------------------------------------------------

void VMS_BCJT::Initialize (int &in_ip,int &in_sd,int &in_en, int Initial_Time_Step)
{
  n_ip = in_ip;  			// Note: Need to call Initialize() for each elmt set
  n_sd = in_sd;
	n_en = in_en;
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en = n_sd * n_en;

	C.Dimension 	(	kNUM_C_TERMS );
	S.Construct 	( kNUM_S_TERMS, 	n_ip 	);
	A.Construct 	( kNUM_A_TERMS, 	n_ip, n_sd, n_sd );
	T4.Construct 	( kNUM_T4_TERMS, 	n_ip, n_sd_x_n_sd, n_sd_x_n_sd );

	NP.Construct	( kNUM_NP_TERMS, n_en, n_sd, n_sd );  // Special one (not FEA at ip) rather values at Node Point (NP)
	NP[kIdentity].Identity( );

	S[kIV_Alpha_n] 	= 0.0; 		//-- Initialize for flow rule 
	S[kIV_Alpha] 	= 0.0;   	

	time_step = Initial_Time_Step;
}

//---------------------------------------------------------------------

void VMS_BCJT::Construct (FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
							int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Time_Integration_Scheme = Integration_Scheme;

	if ( fTime_Step != time_step) { 	// New time step
		S[kIV_Alpha_n] = S[kIV_Alpha];  // Eventually do this with all the kVar_n and make method Next_Step()
		time_step = fTime_Step;
	}

	delta_t = fdelta_t;


	Data_Pro.Construct ( Shapes.dNdx 	);
	Data_Pro.Insert_N  ( Shapes.N 		);
	Integral.Construct (Shapes.j, Shapes.W); 

  //-- The following don't allocate, that's done in Initialize()
	
  Form_C_List			( BCJ_Matl );
  Form_A_S_Lists	( np1,n );
  Form_T4_List		( );
  Form_B_List			( );

	//cout <<"C Constants: "<<C<<"\n\n"; 
	//B.Print("B Matricies");
	//A.Print("A Matricies");
	//S.Print("S Components");
	//T4.Print("T4 Matricies");

}

//---------------------------------------------------------------------
/**  Form stiffness matricies k^Alpha and k^Beta 
 *   The ON swithes in this function are for Alpha and Beta 
 *   contributions respectively */

void VMS_BCJT::Form_LHS_Ka_Kb ( dMatrixT &Ka, dMatrixT &Kb )
{
 	Ka  = Integral.of( B[kB_1hat], B[kB05_tau_3hat] ); 
 	Ka += Integral.of( B[kB_1hat], B[kB06_3hat]  		); 

	Ka += Integral.of( B[kB_1hat], C[kNeg_dt_Root3by2_f], T4[kMM], B[kBa_DEV_H] );
	Kb  = Integral.of( B[kB_1hat], C[kNeg_dt_Root3by2_f], T4[kMM], B[kBb_DEV_H] );

 	if 	(Iso_Hard_Type != kNo_Iso_Hard) 
		Ka += Integral.of( B[kB_1hat], S[kBetaK], T4[kN_o_Nea], B[kBa_Kappa] );
}

//---------------------------------------------------------------------
// F internal (F_int) dimensions here won't actually be in terms of Force

void VMS_BCJT::Form_RHS_F_int ( dArrayT &F_int ) // Untested
{
	FEA_dVectorT G2_vec		( n_ip, n_sd_x_n_sd ); 
	Data_Pro.Reduce_Order	(	A[kG2], G2_vec 		); 

	F_int  = Integral.of	( B[kB_1hat], G2_vec	);  	
}

//=== Private =========================================================

//##################################################################################
//################ B_TERMS #########################################################
//##################################################################################

void VMS_BCJT::Form_B_List (void)
{

		B.Construct (kNUM_B_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_en);  // B = B(9,24)	

	 	Data_Pro.grad_u     	( B[kB_1hat], FEA::kNonSymmetric 	); 
	 	Data_Pro.A_grad_u_T_B	( A[kFbT], 		A[kDa_m], B[kB05_tau_3hat] 		);
	 	Data_Pro.A_grad_u_B 	( A[kDa_m], 	A[kFb],   B[kB06_3hat] 				);

    //===================  del ( DEV( H ) ) 
		
		// NOTE: the term "del" is omitted before each variable for notational ease
		
		//--- Preparatory data
		
	 	Data_Pro.A_grad_u_B 	( A[kA1], 		A[kFb],   B[kBa_Cb_3hat] 			); 		B[kBa_Cb_3hat] 			*= -1.0;
	 	Data_Pro.A_grad_u_B 	( A[kFbT], 		A[kFb],   B[kBb_Cb_3hat] 			);
	 	Data_Pro.A_grad_u_T_B ( A[kFbT],  	A[kA1T],  B[kBa_Cb_tau_3hat] 	);		B[kBa_Cb_tau_3hat]	*= -1.0; 	
	 	Data_Pro.A_grad_u_T_B ( A[kFbT],  	A[kFb],  	B[kBb_Cb_tau_3hat] 	);
	  Data_Pro.A_grad_u_B 	( A[kA2],   	A[kA3],   B[kBa_Cbi_3hat] 		);
	 	Data_Pro.A_grad_u_B 	( A[kA3T],  	A[kA3],	  B[kBb_Cbi_3hat] 		);		B[kBb_Cbi_3hat]			*= -1.0; 	
	 	Data_Pro.A_grad_u_T_B ( A[kA3T],  	A[kA2T],  B[kBa_Cbi_tau_3hat] );
	 	Data_Pro.A_grad_u_T_B ( A[kA3T],  	A[kA3],  	B[kBb_Cbi_tau_3hat] );		B[kBb_Cbi_tau_3hat]	*= -1.0; 	

		//-- Calculation of del ( S )
		B[kBa_S].MultAB( T4[kCC],  B[kBa_Cb_3hat] );  
		B[kBb_S].MultAB( T4[kCC],  B[kBb_Cb_3hat] );

		//-- Calculation of del ( Sym_Zeta )

 		if 	( Back_Stress_Type == kNo_Back_Stress ) {
			B[kBa_H] = B[kBa_S];
			B[kBb_H] = B[kBb_S];
		}
		else {
			Form_del_Zeta_B ( ); // Calculates B[kBa_Z],  B[kBb_Z]	
			//B[kBa_Z].Match_Signs ( B[kBa_S] );
			//B[kBb_Z].Match_Signs ( B[kBb_S] );
			B[kBa_H].DiffOf (	B[kBa_S], B[kBa_Z] );
			B[kBb_H].DiffOf (	B[kBb_S], B[kBb_Z] );
		}	

		//--- Calculation of Ba_DEV_H : del (DEV(H)) = (1-Cbi_o_Cb):del(H) - 1/3(CbiT:del(Cb)) - 1/3.BetaTC.del(Cbi)
		
		T4[kT4_Temp0]  = T4[kCbi_o_Cb];		
		T4[kT4_Temp0] *= -C[k1by3]; 
		T4[kT4_Temp0].PlusIdentity(); 
		B[kBa_DEV_H].MultAB( 	T4[kT4_Temp0], B[kBa_H] );

	  B[kB_Temp0].SumOf  ( 	B[kBa_Cb_3hat], B[kBa_Cb_tau_3hat]	) ;		
		B[kB_Temp1].MultAB ( 	T4[kCbi_o_H], B[kB_Temp0]  );
		B[kB_Temp1]  	*= 		  C[k1by3]; 
		B[kBa_DEV_H] 	-=			B[kB_Temp1]; 

	  B[kB_Temp0].SumOf( 		B[kBa_Cbi_3hat], B[kBa_Cbi_tau_3hat]	) ;	
		B[kB_Temp0] 	*=  		S[kCb_i_H]; 
		B[kB_Temp0] 	*= 		 	C[k1by3]; 
		B[kBa_DEV_H] 	-=			B[kB_Temp0]; 

		//--- Calculation of Bb_DEV_H
		
		B[kBb_DEV_H].MultAB	( T4[kT4_Temp0], B[kBb_H] );

	  B[kB_Temp0].SumOf  ( 	B[kBb_Cb_3hat], B[kBb_Cb_tau_3hat]	) ;		
		B[kB_Temp1].MultAB ( 	T4[kCbi_o_H], B[kB_Temp0]  );
		B[kB_Temp1]  	*= 		  C[k1by3]; 
		B[kBb_DEV_H] 	-=			B[kB_Temp1]; 

	  B[kB_Temp0].SumOf( 		B[kBb_Cbi_3hat], B[kBb_Cbi_tau_3hat]	) ;	
		B[kB_Temp0] 	*=  		S[kCb_i_H]; 
		B[kB_Temp0] 	*= 		 	C[k1by3]; 
		B[kBb_DEV_H] 	-=			B[kB_Temp0]; 

    //===================  del ( DEV( H ) ) [END]
	
		//--- Calculation of Ba_Kappa_bar ( Iso Hardening )

 		if 	(Iso_Hard_Type != kNo_Iso_Hard) {

	 		Data_Pro.A_grad_u_B 	( A[kFaT], 	A[kF],   B[kBa_Kappa_3hat] 			);
	 		Data_Pro.A_grad_u_T_B ( A[kFT],  	A[kFa],  B[kBa_Kappa_tau_3hat] 	);	
			B[kBa_Kappa].SumOf 		( B[kBa_Kappa_tau_3hat], B[kBa_Kappa_3hat] 	);

		}
		
}
						
//##################################################################################
//################ A_TERMS #########################################################
//##################################################################################

//NOTE: np1 := "n+1" time step; n := "n" time step; npt := "n+theta" time step
//      *** No subscript implies n+theta time step ***

void VMS_BCJT::Form_A_S_Lists (VMS_VariableT &npt,VMS_VariableT &n)
{
	//---- Developer cheat: put npt in function door in-lieu-of np1 for speed

	A[kF] 				= npt.Get (	VMS::kF					);
	A[kFa] 				= npt.Get (	VMS::kFa				);
	A[kgrad_ub] 	= npt.Get (	VMS::kgrad_ub		);
	A[kFai] 			= npt.Get (	VMS::kFai				);
  A[kFb]  			= npt.Get (	VMS::kFb				);	 	
  A[kFbi]  			= npt.Get (	VMS::kFbi				);	 	
	A[kFbT].Transpose  			( A[kFb] 					);
	A[kFaT].Transpose  			( A[kFa] 					);
	A[kCa].MultATB     			( A[kFa], A[kFa] 	);	
	A[kCb].MultATB     			( A[kFb], A[kFb] 	);	
	A[kCbi].Inverse    			( A[kCb] ); 
	S[kJb].Determinant 			( A[kFb] );

	//--- LHS: del(Da) terms
	
	A[kFa_n]  = n.Get				(	VMS::kFa );
	A[kFb_n]  = n.Get				(	VMS::kFb );
	A[kCa_n].MultATB				( A[kFa_n], 	A[kFa_n] );
	A[kCb_n].MultATB				( A[kFb_n], 	A[kFb_n] );
	A[kDa_m].MultATBC 			(	A[kFai], 		A[kCa_n], 		A[kFai]	);
	A[kDa_m] *= 0.5;

	//--- LHS: del(sym(Eb)) terms

	A[kF_sharp].MultAB  		( A[kgrad_ub], 	A[kFb] 				);
  A[kA1].MultAB						(	A[kFbT], 			A[kF_sharp] 	);
  A[kA2].MultAB						(	A[kCbi], 			A[kA1] 				);
  A[kA3].MultAB						(	A[kFb],				A[kCbi] 			);
	A[kA1T].Transpose 			( A[kA1] ); 
	A[kA2T].Transpose 			( A[kA2] ); 
	A[kA3T].Transpose 			( A[kA3] ); 

	//--- LHS: del( DEV(Zeta) ) terms
	
	A[kF_sharp_T].Transpose ( A[kF_sharp] );

	//=============================== RHS: Fint terms
	
	A[kDa_mp1] = 0.0;
	A[kDa_mp1].PlusIdentity(0.5); 
	
	//--- Get Eb 

	A[kEb]  = A[kCb]; 
	A[kEb].PlusIdentity (-1.0); 
	A[kEb] *= 0.5; 

	//--- Get Zeta	

 	if 	( Back_Stress_Type != kNo_Back_Stress ) {		//-- Apply Back Stress 
		Get_Back_Stress ( ); //-- Builds  A[kZeta]
		A[kZetaT].Transpose ( A[kZeta]  );
		A[kSym_Zeta]  = A[kZeta];
		A[kSym_Zeta] += A[kZetaT];
		A[kSym_Zeta] /= 2.0; 		// Put all this in a Symmetrize() method
	}
	else
		A[kSym_Zeta] = 0.0;

	//--- Get S	

 	if 	(bControl_Eb)  {
		Form_CC ( ); // if (bControl_Eb) { CC = CCvmst, Lamda_ba, Mu_ba also formed	}
		Data_Pro.C_IJKL_E_KL 	( S[kLamda_ba], S[kMu_ba],	A[kEb], 	A[kS] 			);	
	}
	else {
  	Data_Pro.C_IJKL 			(	C[kLamda],  C[kMu],		T4[kCC]	  			);  
		Data_Pro.C_IJKL_E_KL 	( C[kLamda], 	C[kMu], 	A[kEb],  A[kS] 	);	// Doesn't multiply excessive zeros 
	}

	//-- Get Sigma (For plot and diagnostics) 
	//S[kJ].Determinant( A[kF] );
	A[kSigma].MultABCT 		(  A[kFb], 		A[kS], 	A[kFb] 	);
	A[kSigma] /= S[kJb]; // kJb or kJ (If Consv. of Plastic Volume holds) !!
	if ( S[kJb] == 0.0 ) 	cout << "VMS_BCJ::Form_A_S_Lists(): Diagnostic Varable J=0 \n";

	//A[kSym_Zeta].Match_Signs ( A[kS] );
	A[kH_bar].DiffOf ( A[kS], A[kSym_Zeta] );

	//--- Get DEV(S) : Put this as a method in DataPro ::DEV(Fb,S) ( use MultABCT(Fb,dev(S),Fb) )
	S[kCb_i_H].Double_Dot ( A[kCb], A[kH_bar] ); 
	A[kDEV_H]  = 	A[kCbi];
	A[kDEV_H] *= 	S[kCb_i_H]; 
	A[kDEV_H] *= -C[k1by3]; 
	A[kDEV_H] += 	A[kH_bar]; 

	//--- Get ||DEV(H)||, N
  A[kDEV_H].Mag_and_Dir ( S[kMag_DEV_H], A[kN] );

	//--- Get Beta
	S[kBeta]  = S[kMag_DEV_H];
	S[kBeta] *= C[kRoot3by2]; 

 	if 	( Iso_Hard_Type != kNo_Iso_Hard ) //-- Apply Iso Hard to RHS
		Get_Iso_Hard_Kappa_bar ( ); 
	else
		S[kKappa_bar] = 0.0;

	S[kS_Temp0]  = 	S[kKappa_bar];
	S[kS_Temp0] += 	C[kY];  //-- Y is the yield stress
	S[kBeta] 		-= 	S[kS_Temp0]; 
	S[kBeta] 		/= 	C[kV];

	//--- Get Sinh(Beta) 
	S[kMacaulay_Sinh_Beta].Sinh( S[kBeta] ); 
	S[kMacaulay_Sinh_Beta].Macaulay ( );

	S[kBeta2].Cosh( S[kBeta] ); 
 	S[kBeta2] *= C[kRoot3by2byV]; 

	//--- Include Hardening Coef. Kappa_Bar
 	if 	( Iso_Hard_Type != kNo_Iso_Hard ) {
		S[kBetaK] 	= S[kBeta2]; 
		S[kBetaK]  *= ( delta_t * C[kf] * C[k1by3] * C[kKappa] ); //-- For LHS
	}	

	//-- RHS Start ******
	
  A[kG2]  = A[kN];
  A[kG2] *= S[kMacaulay_Sinh_Beta]; 
  A[kG2] *= C[kNeg_dt_Root3by2_f]; 
  A[kG2] -= A[kDa_m];
  A[kG2] += A[kDa_mp1];

	//-- RHS Finished  ******

	//---- Stuff for purely Diagnostics NOT NEEDED IN CALCULATIONS ------ 

	if (bDiagnosis_Variables) {	

 		if 	(Iso_Hard_Type == kNo_Iso_Hard) {
			//-- Get Ea
			A[kEa]  = A[kCa]; 
			A[kEa].PlusIdentity (-1.0); 
			A[kEa] *= 0.5; 
		}

		//-- Get S_hat and S_tilde
		Data_Pro.C_IJKL_E_KL 	( C[kLamda1], 	C[kMu1], 		A[kEb], 	A[kS_tilde] );	
		Data_Pro.C_IJKL_E_KL 	( C[kLamda2], 	C[kMu2], 		A[kEb], 	A[kS_hat] 	);	
		A[kS_tilde] *= S[kPsi1];	// Diagnostics

		//-- Form La
		A[kFa_dot].DiffOf( A[kFa], A[kFa_n] ); 	// Back Euler	
		A[kFa_dot] /= delta_t; 
		A[kLa].MultAB ( A[kFa_dot], A[kFai] );
		//A[kDa].Sym ( A[kLa] );  // <-- gives bad results

		//-- Form Da	
	  A[kDa] = A[kLa];
	  A[kA_Temp0].Transpose( A[kLa] );
	  A[kDa] += A[kA_Temp0];
	  A[kDa] *= C[k1by2];

		//-- Form la
		A[kla].MultABC (	A[kFb],	 A[kLa], 	A[kFbi]	 );
	
		//-- Form da	
	  A[kda] = A[kla];
	  A[kA_Temp0].Transpose( A[kla] );
	  A[kda] += A[kA_Temp0];
	  A[kda] *= C[k1by2];
	
		//-- Form ArcSinh()	
  	A[kDa].Mag_and_Dir ( S[kMag_Da], A[kNda] ); 
		S[kArcBeta]   = S[kMag_Da];
		S[kArcBeta]  /= C[kf];
		S[kArcSinh].ArcSinh ( S[kArcBeta] );

		S[kpLamda]  =  S[kArcSinh];
		S[kpLamda] *=  C[kV];

 		if 	( Iso_Hard_Type != kNo_Iso_Hard ) 
			S[kpLamda] +=  S[kKappa_bar];

		if ( S[kMag_Da] != 0.0 )
			S[kpLamda] /= S[kMag_Da]; 
		else
			S[kpLamda] = 0.0; 

		/*
		S[kpLamda]  = S[kMag_DEV_H];
		S[kpLamda] /= S[kMacaulay_Sinh_Beta]; 
		S[kpLamda] /= ( C[kRoot3by2] * C[kf] );
		*/
			
		A[kC].MultATB  ( A[kF], A[kF] );	
		A[kE]  = A[kC]; 
		A[kE].PlusIdentity (-1.0); 
		A[kE] *= 0.5; 
	}
	
}

//##################################################################################
//################ 4th Order Tensor Terms (i.e. T4[i] ) ########################
//##################################################################################

void VMS_BCJT::Form_T4_List (void)  // These matricies are all 9x9 (not 3x3 like A)
{

	A[kA_Temp0] = A[kN];
	Data_Pro.A_o_B 					( A[kCbi], 			A[kH_bar], 			T4[kCbi_o_H] 	);
	Data_Pro.A_o_B 					( A[kCbi], 			A[kCb], 				T4[kCbi_o_Cb] );
	Data_Pro.A_o_B 					( A[kA_Temp0], 	A[kN], 					T4[kN_o_N] 		);  

 	Data_Pro.II_minus_A_o_B ( A[kA_Temp0],	A[kN],					T4[kPP] 			); // 4th Order Projector 
	if ( S[kMag_DEV_H] != 0.0 )
  	T4[kPP] /= S[kMag_DEV_H]; // See notes on del( N ) why we do this
	else
  	T4[kPP] = 0.0; 

 	if 	(Iso_Hard_Type != kNo_Iso_Hard) 
		Data_Pro.A_o_B 				( A[kA_Temp0], 	A[kNea], 				T4[kN_o_Nea] 	);  

 	if 	( Back_Stress_Type != kNo_Back_Stress ) { 
		Data_Pro.A_o_B 				( A[kSym_Zeta], 		A[kF_sharp_T], 	T4[kZ_o_F_sharp_T] 		);  
		Data_Pro.A_o_1 				( A[kSym_Zeta], 		T4[kZ_o_1] 		);  
	}

	//-- Build MM
	T4[kT4_Temp0]  = T4[kN_o_N];
	T4[kT4_Temp0] *= S[kBeta2]; 
  T4[kMM]  = T4[kPP];
  T4[kMM] *= S[kMacaulay_Sinh_Beta];
  T4[kMM] += T4[kT4_Temp0];

}

//---------------------------------------------------------------------------------------
//NOTE: the only place CC is used is in forming the LHS tangent, do it makes sense
//to pack CCvmst into CC.  Recall that on the RHS, S is formed by using Lame Constants.

void VMS_BCJT::Form_CC (void) //-- Puts either regular CC or CCvmst in CC 
{

		CCba_Scalars ( );  // Also creates Eb_dot and N_Eb_dot, Mag(Eb_dot)

		//cout << "Lamda1, Mu1 = "<<C[kLamda1]<<", "<<C[kMu1]<<endl;

		//-- Build CC1 & CC2 
  	Data_Pro.C_IJKL (	C[kLamda1],  C[kMu1],		T4[kCC1]	); 
  	Data_Pro.C_IJKL (	C[kLamda2],  C[kMu2],		T4[kCC2]  );  

		//-- Build CC_ba
		T4[kCC_ba]  = T4[kCC1];
		T4[kCC_ba] *= S[kPsi1];
		T4[kCC_ba] += T4[kCC2]; 	
	
		//-- Build CC_gamma	
		Data_Pro.A_o_B 					( A[kEb], 			A[kN_Eb_dot],	T4[kEEn] 			);  
		T4[kCC_gamma].MultAB		( T4[kCC1],			T4[kEEn]										);
		T4[kCC_gamma] *= S[kPsi]; 
		T4[kCC_gamma] /= ( delta_t * C[kGamma_b] ); 
		T4[kCC_gamma] *= C[kAlphaY_plus_1]; 

		//-- Build CC_psi	
		A[kA_Temp0] = A[kEb];
		Data_Pro.A_o_B 					( A[kA_Temp0], 	A[kEb],				T4[kEb_o_Eb] 	);  
		T4[kCC_psi].MultAB			( T4[kCC1],			T4[kEb_o_Eb]								);
		T4[kCC_psi] *= S[kGammaY_bar_Eta1_plus_Eta2];

		//-- Build CCvmst (pack it in CC, go on as usual with del(S) term)
		T4[kCC]  = T4[kCC_psi];
		T4[kCC] += T4[kCC_ba];
		T4[kCC] += T4[kCC_gamma]; // <-- Comment/Uncomment (Turn term off/on) for no rate dependency


}

//##################################################################################
//################ Material Constant Terms C[i] #####################################
//##################################################################################

void VMS_BCJT::Form_C_List (VMF_MaterialT *BCJ_Matl)
{

	//-- Elasticity 	
	C[kLamda]    	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kLamda 		);
	C[kMu]    	 	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kMu 			);

	//-- BCJ
	C[kf]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kf	 			);
	C[kV]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kV	 			);
	C[kY]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kY	 			);

	//-- CCba 
	C[kLamda1]    = BCJ_Matl -> Retrieve ( BCJ_MatlT::kLamda1 	);
	C[kMu1]    	 	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kMu1 			);
	C[kLamda2]    = BCJ_Matl -> Retrieve ( BCJ_MatlT::kLamda2		);
	C[kMu2]    	 	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kMu2 			);
	C[kPi]    		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kPi 			);
	C[kRho]    	 	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kRho 			);
	C[kGamma_b]   = BCJ_Matl -> Retrieve ( BCJ_MatlT::kGamma_b 	);
	C[kAlphaY]   	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kAlphaY 	);

	//-- Iso Hard
	C[kKappa]    	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kPlastic_Modulus_K );
	C[kH]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kH	 			);

	//-- Kine Hard (Back Stress)
	C[kc]    			= BCJ_Matl -> Retrieve ( BCJ_MatlT::kc_zeta 	);
	C[kl]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kl 				);

	//-- Combinations
	C[k2Gamma_b]						= 2.0 * C[kGamma_b]; 
	C[kAlphaY_plus_1]				= C[kAlphaY] + 1.0; 
	C[kRoot3by2]						= sqrt(1.5);
	C[kNeg_dt_Root3by2_f]  	= -1.0 * delta_t * C[kRoot3by2] * C[kf]; 
	C[kRoot3by2byV]					= C[kRoot3by2] / C[kV];
	C[k1by2] 								= 1.0/2.0;
	C[k1by3] 								= 1.0/3.0;
	C[k2by3] 								= 2.0/3.0;
	C[kRoot2by3]						= sqrt( C[k2by3] ); 
	C[kMu_c_l]							= C[kMu] * C[kc] * C[kl];

}

//##################################################################################
//################ Isotropic Hardening #############################################
//##################################################################################

void VMS_BCJT::Get_Iso_Hard_Kappa_bar ( void )
{
		A[kEa]  = A[kCa]; 
		A[kEa].PlusIdentity (-1.0); 
		A[kEa] *= 0.5; 

		A[kEa_n]  = A[kCa_n]; 
		A[kEa_n].PlusIdentity (-1.0); 
		A[kEa_n] *= 0.5; 

		A[kEa_dot].DiffOf ( A[kEa], A[kEa_n] );
		A[kEa_dot] /= delta_t; 
  	A[kEa_dot].Mag_and_Dir 		( S[kMag_Ea_dot], 	A[kNea] );

		if ( Iso_Hard_Type == kMethod1 || Iso_Hard_Type == kMethod2 ) {

			//-- Derive Kappa Bar Term 1
			S[kKappa_bar]		 = S[kIV_Alpha_n];    
			S[kKappa_bar]		*= ( C[kKappa]*C[kRoot2by3] );    

			//-- Derive Kappa Bar Term 2 
  		S[kS_Temp0]  = S[kMag_Ea_dot]; 
  		S[kS_Temp0] *= ( C[kRoot2by3] * C[kKappa] * C[kH] * delta_t ); 

			//-- Derive Kappa Bar: Sum Term 1 and Term 2 
  		S[kKappa_bar]  += S[kS_Temp0];  

			//-- Derive Alpha (This is needed to be used as Alpha_n next step)
			S[kIV_Alpha]		 = S[kMag_Ea_dot]; 
			S[kIV_Alpha]		*= ( C[kRoot2by3] * C[kH] * delta_t ); 
			S[kIV_Alpha]		+= S[kIV_Alpha_n];    

		}
		else
			cout << " ...ERROR >> VMS_BCJT::Get_Iso_Hard_Kappa() : Bad Iso_Hard_Type \n";

}

//##################################################################################
//################ Kinematic Hardening ( Back Stress ) #############################
//##################################################################################

//##################################################### Back Stress (RHS) 

void VMS_BCJT::Get_Back_Stress ( void )
{
	if ( Back_Stress_Type == kSteinmann ) {

		int n_ed = n_sd_x_n_sd;
		int n_ed_x_n_en = n_ed*n_en; 
																																				// -- 2D --
		FEA_dMatrixT 	me 					( n_ip, n_ed_x_n_en, n_ed_x_n_en );  			// (16x16)
		FEA_dMatrixT 	B_tilde 		( n_ip, n_ed, n_ed_x_n_en ); 							// (4x16)
		FEA_dVectorT 	grad_ub_vec ( n_ip, n_sd_x_n_sd ); 										// (4x1)

		dMatrixT 			M						( n_ed_x_n_en, n_ed_x_n_en ); 						// (16x16)
		dArrayT				fe					( n_ed_x_n_en ); 													// (16x1)
		dArrayT				sE_vec			( n_sd_x_n_sd * n_en ); 									// (16x1)

		/*sE_mat.Dimension ( n_en ); // sE_mat is a member variable
		for (int a=0; a<n_en; a++)
			sE_mat[a].Dimension ( n_sd ); */

		Data_Pro.Reduce_Order	(	A[kgrad_ub], grad_ub_vec ); 
		Data_Pro.Mass 				( n_ed, me );
		Data_Pro.Mass_B 			( n_ed, B_tilde );

 		M  = Integral.of  ( me );
		fe = Integral.of	( B_tilde, grad_ub_vec	);  	// (16x4)(4x1) = (16x1)

		M.Inverse ( );

		M.Multx ( fe,sE_vec ); // (16x16)(16x1) = (16x1)

		Data_Pro.Element_Nodal_Values_Expand_Order ( sE_vec , NP[ksE_mat] );	// (16x1) --> (4x2x2)

		//-- Static Condensation follows:
	
		Data_Pro.Curl ( NP[ksE_mat], A[kcurl_sE] );  // Recall Curl ( 1 - sE ) = Curl ( 1 ) - Curl (sE) = -Curl(sE);

		// Note: two negatives cancel each other: Zeta = - c*mu*l*Jb*Fbi*(Curl( -sE )^T) 

		A[kZeta].MultABT ( A[kFbi], A[kcurl_sE] );
		A[kZeta] *= S[kJb]; 
		A[kZeta] *= C[kMu_c_l]; 

		//A[kZeta].Print( "Zeta" );

	}
	else
		cout << " VMS_BCJT::Get_Back_Stress() >> Unknown Back_Stress_Type \n";
	
}

//##################################################### Back Stress (LHS) [ Linearized ]
// Note: all C[],S[],A[],and T4[] terms already formed

void VMS_BCJT::Form_del_Zeta_B ( void ) // Calculates B[kBa_Z],  B[kBb_Z]	
{

	bool bDEL_sE=1, bDEL_grad_ub=1;

	//-- Mu*c*l*del(Jb)*Fb^-1*curl(sE)^T  (Symmetry of Zeta accounted for in T4)
			
	B[kBa_Zeta_Jb].MultAB ( T4[kZ_o_F_sharp_T], 	B[kB_1hat] );  	// Needs neg sign
	B[kBb_Zeta_Jb].MultAB ( T4[kZ_o_1], 					B[kB_1hat] );		// Note: Zeta = Mu_c_l_Fbi_(curl(sE))^T

	//-- Mu*c*l*Jb*del(Fb^-1)*curl(sE)^T

	A[kA4].MultAB ( A[kFbi], A[kF_sharp] );	
	A[kA4T].Transpose 	( A[kA4] 	);	
	A[kFbiT].Transpose 	( A[kFbi] );	

	A[kMu_c_l_Jb_curl_sE_T].Transpose(  A[kcurl_sE] );		// Mu_c_l_Jb_curl_sE_T is R2 in the notes
	A[kMu_c_l_Jb_curl_sE_T] *= S[kJb]; 
	A[kMu_c_l_Jb_curl_sE_T] *= C[kMu_c_l]; 
	A[kR2T].Transpose ( 	A[kMu_c_l_Jb_curl_sE_T] 	);

  Data_Pro.A_grad_u_B ( A[kA4],         A[kMu_c_l_Jb_curl_sE_T],   	B[kBa_Zeta_Fbi_reg] );
  Data_Pro.A_grad_u_B ( A[kFbi],        A[kMu_c_l_Jb_curl_sE_T],   	B[kBb_Zeta_Fbi_reg] ); // Needs neg sign

	Data_Pro.A_grad_u_T_B ( A[kR2T], 			A[kA4T],   									B[kBa_Zeta_Fbi_tau] ); 
	Data_Pro.A_grad_u_T_B ( A[kR2T], 			A[kFbiT],   								B[kBb_Zeta_Fbi_tau] ); // Needs neg sign

	B[kBa_Zeta_Fbi].SumOf ( B[kBa_Zeta_Fbi_reg], 	B[kBa_Zeta_Fbi_tau] );
	B[kBa_Zeta_Fbi] *= C[k1by2]; 

	B[kBb_Zeta_Fbi].SumOf ( B[kBb_Zeta_Fbi_reg], 	B[kBb_Zeta_Fbi_tau] );
	B[kBb_Zeta_Fbi] *= C[k1by2]; // Needs neg sign

	if ( bDEL_sE == bDEL_grad_ub ) { // sE behaves as grad_ub, likewise del(sE) will behave as del(grad_ub)

		//-- Mu*c*l*Jb*Fb^-1*del(curl(sE)^T)
	
		NP[kFbi_hat_mat]  = NP[kIdentity];	// Fbi_hat = 1 - sE
		NP[kFbi_hat_mat] -= NP[ksE_mat];

		int i,j,s,r,m,p;
		double E_rmp;

		FEA_dScalarT Fip_sEjsm_Ermp ( n_ip );  // F,sE,E -->  Fbi, script Epsilon, permutation sybol
		FEA_dScalarT sE_jsm 				( n_ip );  // d(sE(j,s))/dm 

		FEA_dScalarT Fip_hFjsm_Ermp ( n_ip );  // F,hF,E -->  Fbi, Fbi_hat, permutation sybol
		FEA_dScalarT Fbi_hat_jsm		( n_ip );  // d(Fbi_hat(j,s))/dm 

		FEA_dScalarT Fip_sEjms_Emrp ( n_ip );  // F,sE,E -->  Fbi, script Epsilon, permutation sybol
		FEA_dScalarT sE_jms 				( n_ip );  // d(sE(j,m))/ds 

		for (i=0; i<n_sd; i++)		// Form 4th Order Tensors
			for (j=0; j<n_sd; j++)
				for (s=0; s<n_sd; s++)
					for (r=0; r<n_sd; r++) {

						T4[kAAe_2bar]( Data_Pro.Map(i,j) , Data_Pro.Map(s,r) ) = 0.0; 
						T4[kAAf_2bar]( Data_Pro.Map(i,j) , Data_Pro.Map(s,r) ) = 0.0; 
						T4[kEE_2bar] ( Data_Pro.Map(i,j) , Data_Pro.Map(s,r) ) = 0.0; 

						//-- Loop over dummy indicies m,p --> summation forms AA_2bar 4th order tensors
					
						for (m=0; m<n_sd; m++) {
							for (p=0; p<n_sd; p++) {

								E_rmp = Data_Pro.e_ijk[r][m][p];  // Permutation symbol

								Data_Pro.Grad_ij ( NP[ksE_mat],     	j,s,m, sE_jsm 			); 
								Data_Pro.Grad_ij ( NP[kFbi_hat_mat],	j,s,m, Fbi_hat_jsm 	);
								Data_Pro.Grad_ij ( NP[ksE_mat],     	j,m,s, sE_jms 			); 

								//-- AAe_2bar 
								Fip_sEjsm_Ermp  = A[kFbi](i,p); 
								Fip_sEjsm_Ermp *= sE_jsm;
								Fip_sEjsm_Ermp *= E_rmp; 
								T4[kAAe_2bar]( Data_Pro.Map(i,j) , Data_Pro.Map(s,r) ) +=  Fip_sEjsm_Ermp;

								//-- AAf_2bar 
								Fip_hFjsm_Ermp  = A[kFbi](i,p); 
								Fip_hFjsm_Ermp *= Fbi_hat_jsm;
								Fip_hFjsm_Ermp *=  E_rmp; 
								T4[kAAf_2bar]( Data_Pro.Map(i,j) , Data_Pro.Map(s,r) ) +=  Fip_hFjsm_Ermp;

								//-- EE_2bar 
								Fip_sEjms_Emrp  = A[kFbi](i,p); 
								Fip_sEjms_Emrp *= sE_jms;
								Fip_sEjms_Emrp *=  -E_rmp; // Recall: E_mrp = -E_rmp
								T4[kEE_2bar]( Data_Pro.Map(i,j) , Data_Pro.Map(s,r) ) +=  Fip_sEjms_Emrp;

					  	}	
						}
					}
	
		S[kS_Temp0]  =   S[kJb]; 
		S[kS_Temp0] *=   C[kMu_c_l]; 
 
		T4[kT4_Temp0].SumOf ( T4[kAAe_2bar], T4[kEE_2bar] );
		T4[kT4_Temp0] *=  S[kS_Temp0]; 
		B[kBa_Zeta_curl_sE_T].MultAB ( T4[kT4_Temp0], 	B[kB_1hat] );  // Needs neg sign

		T4[kT4_Temp0].DiffOf ( T4[kAAf_2bar], T4[kEE_2bar] );
		T4[kT4_Temp0] *=  S[kS_Temp0]; 
		B[kBb_Zeta_curl_sE_T].MultAB ( T4[kT4_Temp0], 	B[kB_1hat] ); 

	}	
	else {
		B[kBa_Zeta_curl_sE_T] = 0.0; 
		B[kBb_Zeta_curl_sE_T] = 0.0; 
	}

	//-- Summation of all contributing B terms

	B[kBa_Z]  = B[kBa_Zeta_Fbi];
	B[kBa_Z] -= B[kBa_Zeta_Jb];
	B[kBa_Z] -= B[kBa_Zeta_curl_sE_T];

	B[kBb_Z]  = B[kBb_Zeta_Jb];
	B[kBb_Z] -= B[kBb_Zeta_Fbi];				// Neg sign applied here
	B[kBb_Z] += B[kBb_Zeta_curl_sE_T];

}

//################################## CCba ################################################

void VMS_BCJT::CCba_Scalars( )  
{

#if 0
	int Mag_Eb_NOT_Zero=1;
	for (int l=0; l<n_ip; l++)
		if ( S[kMag_Eb][l] == 0.0 ) 
			Mag_Eb_NOT_Zero==0;
	if ( Mag_Eb_NOT_Zero) {
#endif

  A[kEb].Mag_and_Dir ( S[kMag_Eb], A[kNeb] ); 

	if ( S[kMag_Eb] != 0.0 ) {

  	S[kRho_Mag_Eb]  = S[kMag_Eb]; 
  	S[kRho_Mag_Eb] *= C[kRho]; 

		S[kPi_Tanh].Tanh( S[kRho_Mag_Eb] ); 
		S[kPi_Tanh] *= C[kPi];
		S[kPsi]  	 	 = S[kPi_Tanh]; 
		S[kPsi] 		/= S[kMag_Eb]; 

		S[kEta1]  = S[kPsi];			// Eta1 = Psi/(Mag_Eb^2)
		S[kEta1] /= S[kMag_Eb]; 
		S[kEta1] /= S[kMag_Eb]; 

		S[kEta2].Sech( S[kRho_Mag_Eb] ); 
		S[kEta2].Squared ( );
		S[kEta2] *= C[kPi]; 
		S[kEta2] /= S[kMag_Eb]; 
		S[kEta2] /= S[kMag_Eb]; 

		S[kEta1_plus_Eta2]  = S[kEta1];
		S[kEta1_plus_Eta2] += S[kEta2];

		S[kPsi1] 		= S[kPsi];
		S[kGammaY_bar_Eta1_plus_Eta2]  = S[kEta1_plus_Eta2];

		//------------- Add rate dependency to S_tilde (Modify Psi1 with rate term GammaY)
#if 1 
		//--- Get Eb_dot
		
		A[kEb_n]  = A[kCb_n]; 
		A[kEb_n].PlusIdentity (-1.0); 
		A[kEb_n] *= 0.5; 

		A[kEb_dot].DiffOf ( A[kEb], A[kEb_n] );
		A[kEb_dot] /= delta_t; 
  	A[kEb_dot].Mag_and_Dir 		( S[kMag_Eb_dot], 	A[kN_Eb_dot] );

		//--- Finish Psi1
	
		S[kGammaY]  	  = S[kMag_Eb_dot];	
		S[kGammaY] 		 /= C[kGamma_b];	
		S[kGammaY_bar]  = S[kGammaY];	
		S[kGammaY_bar] *= C[kAlphaY_plus_1];	
		S[kGammaY_bar] -= C[kAlphaY];	

		//-- Comment/Uncomment these for rate dependency
		S[kPsi1]   *= S[kGammaY_bar];
		S[kGammaY_bar_Eta1_plus_Eta2]  *= S[kGammaY_bar];
#endif 
		//-------------- End Modification
		
		//-- Build Lame_ba Constants
	  S[kLamda_ba]  = 	S[kPsi1];
	  S[kLamda_ba] *= 	C[kLamda1];
	  S[kLamda_ba] += 	C[kLamda2];
	  S[kMu_ba]  		= 	S[kPsi1];
	  S[kMu_ba] 	 *= 	C[kMu1];
	  S[kMu_ba] 	 += 	C[kMu2];

	}
	else { // So it doesn't blow up
		S[kPsi1] = 0.0;
		S[kEta1_plus_Eta2]  = 0.0;
	}

}

//##############################################################################################

void VMS_BCJT::Get ( StringT &Name, FEA_dMatrixT &tensor )
{
	if ( Name == "F" )
		tensor = A[kF];
	else if ( Name == "Fa" )
		tensor = A[kFa];
	else if ( Name == "Fb" )
		tensor = A[kFb];
	else if ( Name == "C" )
		tensor = A[kC];
	else if ( Name == "Ca" )
		tensor = A[kCa];
	else if ( Name == "Cb" )
		tensor = A[kCb];
	else if ( Name == "Cbi" )
		tensor = A[kCbi];
	else if ( Name == "S" )
		tensor = A[kS];
	else if ( Name == "S_hat" )
		tensor = A[kS_hat];
	else if ( Name == "S_tilde" )
		tensor = A[kS_tilde];
	else if ( Name == "Sigma" )
		tensor = A[kSigma];
	else if ( Name == "Da" )
		tensor = A[kDa];
	else if ( Name == "La" )
		tensor = A[kLa];
	else if ( Name == "da" )
		tensor = A[kda];
	else if ( Name == "la" )
		tensor = A[kla];
	else if ( Name == "Ea_dot" )
		tensor = A[kEa_dot];
	else if ( Name == "Eb_dot" )
		tensor = A[kEb_dot];
	else if ( Name == "E" )
		tensor = A[kE];
	else if ( Name == "Ea" )
		tensor = A[kEa];
	else if ( Name == "Eb" )
		tensor = A[kEb];
	else if ( Name == "grad_ub" )
		tensor = A[kgrad_ub];
	else if ( Name == "curl_sE" )
		tensor = A[kcurl_sE];
	else if ( Name == "Zeta" )
		tensor = A[kZeta];
	else if ( Name == "Sym_Zeta" )
		tensor = A[kSym_Zeta];
	else if ( Name == "H" )
		tensor = A[kH_bar];
	else
		cout << " ...ERROR: VMS_BCJT::Get() >> Unknown tensor '"<<Name<<"' requested. \n";
}

//##################################################################################

void VMS_BCJT::Get ( StringT &Name, FEA_dScalarT &scalar )
{
	if ( Name == "J" )
		scalar = S[kJ];
	else if ( Name == "Jb" )
		scalar = S[kJb];
	else if ( Name == "Rho_Mag_Eb" )
		scalar = S[kRho_Mag_Eb];
	else if ( Name == "Psi" )
		scalar = S[kPsi];
	else if ( Name == "Psi1" )
		scalar = S[kPsi1];
	else if ( Name == "Pi_Tanh" )
		scalar = S[kPi_Tanh];
	else if ( Name == "Beta" )
		scalar = S[kBeta];
	else if ( Name == "Sinh_Beta" )
		scalar = S[kMacaulay_Sinh_Beta];
	else if ( Name == "Mag_DEV_H" )
		scalar = S[kMag_DEV_H];
	else if ( Name == "Mag_Eb_dot" )
		scalar = S[kMag_Eb_dot];
	else if ( Name == "GammaY" )
		scalar = S[kGammaY];
	else if ( Name == "GammaY_bar" )
		scalar = S[kGammaY_bar];
	else if ( Name == "pLamda" )
		scalar = S[kpLamda];
	else if ( Name == "IV_Alpha" )
		scalar = S[kIV_Alpha];
	else if ( Name == "Kappa_bar" )
		scalar = S[kKappa_bar];
	else if ( Name == "ArcBeta" )
		scalar = S[kArcBeta];
	else if ( Name == "ArcSinh" )
		scalar = S[kArcSinh];
	else if ( Name == "Mag_Da" )
		scalar = S[kMag_Da];
	else
		cout << " ...ERROR: VMS_BCJT::Get() >> Unknown scalar '"<<Name<<"' requested. \n";
}

