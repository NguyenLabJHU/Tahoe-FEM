//DEVELOPMENT


#include "FEA.h" 
#include "VMS.h" // <-- Switch name to VMF later 

using namespace Tahoe;

VMF_Virtual_Work_EqT::VMF_Virtual_Work_EqT	( FEA_ShapeFunctionT &Shapes,VMF_MaterialT *Iso_Matl,VMS_VariableT &np1,VMS_VariableT &n, 
											int Integration_Scheme) 
{
	Construct (Shapes,Iso_Matl,np1,n,Integration_Scheme);
}

/* destructor */
//VMF_Virtual_Work_EqT::~VMF_Virtual_Work_EqT(void)
//{
//}

//---------------------------------------------------------------------

void VMF_Virtual_Work_EqT::Construct ( 	FEA_ShapeFunctionT &Shapes,VMF_MaterialT *Iso_Matl,VMS_VariableT &np1,VMS_VariableT &n, 
														int Integration_Scheme) 
{
	Time_Integration_Scheme = Integration_Scheme;

	n_ip 		= np1.fVars[0].IPs(); 
	n_rows	= np1.fVars[0].Rows(); 
	n_cols	= np1.fVars[0].Cols();
	n_en    = Shapes.dNdx.Cols();
  n_sd 		= n_rows;	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en = n_sd * n_en;

	lamda = Iso_Matl -> Retrieve ( Iso_MatlT::kLamda 	);
	mu 		= Iso_Matl -> Retrieve ( Iso_MatlT::kMu 		);
	Pi 		= Iso_Matl -> Retrieve ( Iso_MatlT::kPi 		);
	Rho 	= Iso_Matl -> Retrieve ( Iso_MatlT::kRho 		);
 
	Data_Pro.Construct ( Shapes.dNdx );

  Form_A_S_Lists 	(np1,n);
  Form_B_List			();

	//A.Print("A Matricies");
	//T4.Print("T4 Matrix");
	//S.Print("S Scalar");
	//B.Print("B Matricies");

	Integral.Construct ( Shapes.j, Shapes.W ); 

}

//---------------------------------------------------------------------

/** RECALL CLASSIC NEWTON-RAPHSON :
 
    K.delta_d = -[k.d - f]  : Where both K and k are element matricies.
    And where K is the tangent, drive it to zero for roots (RHS and delta_d will also vanish) 
		Note, k.d are the internal forces (F_int) and f are the external forces (F_ext).

	  MULTI-FIELD NEWTON-RAPHSON :	

    Ka.delta_da + Kb.delat_db = -[F_int - f]  : Where Ka,Kb, and k are element matricies.
    Ka and Kb are tangents, drive them to zero in a staggered scheme (RHS, delta_da, and delta_db
	 	will also vanish).  Note, dislocation glide (da) contributes to stress in a round-a-bout way.
		While it is true that the intermediate configuration is stress free, second PK S = S(Fb) 
		and Fb := dx/dX_bar. Recall X_bar = X + ua.  So changes in ua will ultimately affect S and
		sigma sinse sigma = 1/j FbSFb^T.  */	

//---------------------------------------------------------------------

void VMF_Virtual_Work_EqT::Form_LHS_Ka_Kb	( dMatrixT &Ka, dMatrixT &Kb, double delta_t)  // Untested
{
#pragma unused(delta_t)

	/* Term I. 		*/		Ka 	= Integral.of( 	B[kB_1hat], B[kBI_tau_3hat] 									);  	
	/* Term IIb. 	*/	 	Kb  = Integral.of( 	B[kB_1hat], B[kBbII_2hat] 										);  	
	/* Term IIa. 	*/	 	Ka -= Integral.of( 	B[kB_1hat], B[kBaII_3hat] 										);  	
	/* Term IIIb.	*/	 	Kb += Integral.of( 	B[kB_1hat], T4[kcc_b_tilde], B[kB_1hat] 			);  	
	/* Term IIIa.	*/	 	Ka -= Integral.of( 	B[kB_1hat], T4[kcc_b_tilde], B[kBaIII_2bar] 	);  	
}

//---------------------------------------------------------------------

void VMF_Virtual_Work_EqT::Form_RHS_F_int ( dArrayT &F_int, double delta_t) // Untested
{
#pragma unused(delta_t)

	FEA_dVectorT sigma_vec	( n_ip, n_sd_x_n_sd 		); // <-- Dimensionality problem
	Data_Pro.Reduce_Order		(	A[kSigma], sigma_vec 	); 

	F_int = Integral.of			( B[kB_1hat], sigma_vec	);  // <-- sigma_vec must be dim n_sd_x_n_en	
}

//=== Private =========================================================

//##################################################################################
//################ B_TERMS #########################################################
//##################################################################################

void VMF_Virtual_Work_EqT::Form_B_List (void)
{
		B.Construct ( kNUM_B_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_en );	

	 	Data_Pro.grad_u       		( 														B[kB_1hat]  			); 
	 	Data_Pro.A_grad_u_T_B 		( 	A[kSigma], 		A[kFbT],  	B[kBI_tau_3hat] 	);
	 	Data_Pro.grad_u_A			 		( 	A[kSigma], 								B[kBbII_2hat] 		);
	 	Data_Pro.A_grad_u_B			 	( 	A[kF_sharp],	A[kSigma], 	B[kBaII_3hat] 		);
	 	Data_Pro.A_grad_u				 	( 	A[kF_sharp],						 	B[kBaIII_2bar] 		);
}
						
//##################################################################################
//################ A_TERMS #########################################################
//##################################################################################

//NOTE: np1 := "n+1" time step; n := "n" time step; npt := "n+theta" time step
//      *** No subscript implies n+theta time step ***

void VMF_Virtual_Work_EqT::Form_A_S_Lists (VMS_VariableT &npt,VMS_VariableT &n,int Integration_Scheme)
{

#if 0
	//---- Developer cheat: put npt in function door in-lieu-of np1 or n for speed 
	VMS_VariableT npt(	n.Get(VMS::kGRAD_ua), n.Get(VMS::kGRAD_ub)	);  

	if 			(	Integration_Scheme == FEA::kForward_Euler		)		npt = n;
	else if (	Integration_Scheme == FEA::kBackward_Euler	)		npt = np1;
	else if (	Integration_Scheme == FEA::kCrank_Nicholson	) { npt.SumOf(np1,n); npt *= 0.5; }
	else 	cout << " ...ERROR >> VMF_Virtual_Work_EqT::Form_A_List() : Bad theta value for time stepping \n";
#endif

	A.Construct ( kNUM_A_TERMS, n_ip, n_sd, n_sd);
	S.Construct ( kNUM_S_TERMS, n_ip);

  A[kF]   		= npt.Get (	VMS::kF					);	 // NOTE: kF != VMS::kF 	
  A[kFb]  		= npt.Get (	VMS::kFb				);	 	
	A[kgrad_ub]	= npt.Get (	VMS::kgrad_ub		);

	A[kF].Determinant	 	  ( S[kJ] 								);  
	A[kFbT].Transpose  		( A[kFb] 								);
	A[kF_sharp].MultAB  	( A[kgrad_ub], 	A[kFb] 	);

	//----- Calculate Eb 
	
	A[kCb].MultATB   			( A[kFb],  			A[kFb] 	);					
	A[kEb]  = A[kCb]; 
	A[kEb].PlusIdentity(-1.0); 
	A[kEb] *= 0.5; 

	if 	( bControl_Eb ) {
		Control_Eb ( ); // Eb_tilde calculated
		//A[kEb_tilde].Print ("Eb_tilde Coarse");
	}
	else
		A[kEb_tilde] = A[kEb];

	//----- Calculate stresses S and Sigma

  Data_Pro.C_IJKL_E_KL	(  lamda, mu, A[kEb_tilde], A[kS] 	); // UT

	A[kSigma].MultABCT 		(  A[kFb], 		A[kS], 	A[kFb] 	);
	A[kSigma] /= S[kJ]; // not kJb !!

	//----- 4th Order VMF Finite Strain Elasticity Tensor

	T4.Construct ( kNUM_T4_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_sd );

  Data_Pro.c_ijkl				(	lamda,mu, S[kJ], A[kFb], T4[kcc_b]	);
  //Data_Pro.c_ijkl_Alt		(	lamda,mu, S[kJ], A[kFb], T4[kT4_Temp0]	);

 	if 	(0) {  //------- Build CC_tilde for Control of Eb
 	//if 	(bControl_Eb) {  //------- Build CC_tilde for Control of Eb

		A[kA_Temp0] = A[kNeb];
		Data_Pro.A_o_B 					( A[kA_Temp0], 	A[kNeb], 				T4[kNeb_o_Neb] 	);  
 		Data_Pro.II_minus_A_o_B ( A[kA_Temp0],	A[kNeb],				T4[kPPeb] 			); // 4th Order Projector 

		S[kPi_Sech_Rho_Mag_Eb_2].Sech ( S[kRho_Mag_Eb] );
		S[kPi_Sech_Rho_Mag_Eb_2].Squared ( );
		S[kPi_Sech_Rho_Mag_Eb_2] *= Pi; 

		T4[kT4_Temp0]  = T4[kNeb_o_Neb];
		T4[kT4_Temp0] *=  S[kPi_Sech_Rho_Mag_Eb_2];
		T4[kRR]  = T4[kPPeb]; 
		T4[kRR] *= S[kPi_Tanh_Rho_Mag_Eb]; 
		T4[kRR] += T4[kT4_Temp0];  

		Data_Pro.C_IJKL 	(	lamda, mu, T4[kCC]	); 
		Data_Pro.c_ijkl  	(	S[kJ], A[kFb], T4[kCC], T4[kRR], T4[kcc_b_tilde] );

	}
	else
		T4[kcc_b_tilde] =  T4[kcc_b];
	
}


//################################## Enhance Eb ################################################

void VMF_Virtual_Work_EqT::Control_Eb ( )  // Produces Eb_tilde
{
  A[kEb].Mag_and_Dir ( S[kRho_Mag_Eb], A[kNeb] ); 
  S[kRho_Mag_Eb] *= Rho; 

	S[kPi_Tanh_Rho_Mag_Eb].Tanh( S[kRho_Mag_Eb] ); 
	S[kPi_Tanh_Rho_Mag_Eb] *= Pi; 

	 A[kEb_tilde]  = A[kNeb]; 
	 A[kEb_tilde] *= S[kPi_Tanh_Rho_Mag_Eb];	
}

//##################################################################################

void VMF_Virtual_Work_EqT::Get ( StringT &Name, FEA_dMatrixT &tensor )
{
	if ( Name == "Sigma" )
		tensor = A[kSigma];
	else if ( Name == "S" )
		tensor = A[kS];
	else if ( Name == "Eb" )
		tensor = A[kEb];
	else if ( Name == "Eb_tilde" )
		tensor = A[kEb_tilde];
	else if ( Name == "F" )
		tensor = A[kF];
	else if ( Name == "Fb" )
		tensor = A[kFb];
	else if ( Name == "Cb" )
		tensor = A[kCb];
	else if ( Name == "grad_ub" )
		tensor = A[kgrad_ub];
	else
		cout << " ...ERROR: VMF_Virtual_Work_EqT::Get() >> Unknown tensor '"<<Name<<"' requested. \n";
}

//##################################################################################

void VMF_Virtual_Work_EqT::Get ( StringT &Name, FEA_dScalarT &scalar )
{
	if ( Name == "J" )
		scalar = S[kJ];
	else if ( Name == "Rho_Mag_Eb" )
		scalar = S[kRho_Mag_Eb];
	else if ( Name == "Pi_Tanh_Rho_Mag_Eb" )
		scalar = S[kPi_Tanh_Rho_Mag_Eb];
	else
		cout << " ...ERROR: VMF_Virtual_Work_EqT::Get() >> Unknown scalar '"<<Name<<"' requested. \n";
}

