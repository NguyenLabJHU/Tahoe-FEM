// $Id: APS_kappa_alphaT.cpp,v 1.2 2005-05-03 15:54:47 raregue Exp $
#include "APS_kappa_alphaT.h"

using namespace Tahoe;

APS_kappa_alphaT::APS_kappa_alphaT	(FEA_ShapeFunctionT &Shapes_displ, FEA_ShapeFunctionT &Shapes_plast, APS_MaterialT *APS_Matl, 
						APS_VariableT &np1, 
						APS_VariableT &n, int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes_displ, Shapes_plast, APS_Matl,np1,n,fTime_Step,fdelta_t,Integration_Scheme);
}

/* destructor */
//APS_kappa_alphaT::~APS_kappa_alphaT(void);


//---------------------------------------------------------------------

void APS_kappa_alphaT::Initialize (int &in_ip, int &in_sd, int &in_en_displ, int &in_en_plast, int &in_state, int &in_str, 
							int Initial_Time_Step )
{
  	n_ip = in_ip;  			// Note: Need to call Initialize() for each elmt set
  	n_sd = in_sd;
	n_en_displ = in_en_displ;
	n_en_plast = in_en_plast;
	n_state = in_state;
	n_str = in_str;
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en_plast = n_sd * n_en_plast;

	C.Dimension 	( kNUM_C_TERMS );
	
	S.Construct 	( kNUM_S_TERMS, 	n_ip 		);
	V.Construct 	( kNUM_V_TERMS, 	n_ip, n_sd 	);
	V_state.Construct ( kNUM_V_state_TERMS, n_ip, n_state );
	int dum=n_state+n_str;
	V_out.Construct ( kNUM_V_out_TERMS, n_ip, dum );
	VB_d.Construct 	( kNUM_VB_d_TERMS, n_ip, n_en_displ );
	VB_eps.Construct ( kNUM_VB_eps_TERMS, n_ip, n_sd_x_n_en_plast );
	B_d.Construct 	( kNUM_B_d_TERMS, n_ip, n_sd, n_en_displ);  
	B_eps.Construct ( kNUM_B_eps_TERMS, n_ip, n_sd, n_sd_x_n_en_plast);
	dum=1;
	B_gradu.Construct ( kNUM_B_gradu_TERMS, n_ip, dum, n_sd );  
	B_gradgammap.Construct ( kNUM_B_gradgammap_TERMS, n_ip, n_sd, n_sd);  	

	time_step = Initial_Time_Step;
}

//---------------------------------------------------------------------

void APS_kappa_alphaT::Construct (FEA_ShapeFunctionT &Shapes_displ, FEA_ShapeFunctionT &Shapes_plast, APS_MaterialT *APS_Matl, 
							APS_VariableT &np1, APS_VariableT &n, 
							int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Time_Integration_Scheme = Integration_Scheme;

	delta_t = fdelta_t;

	Data_Pro_Displ.Construct ( Shapes_displ.dNdx 	);
	Data_Pro_Displ.Insert_N  ( Shapes_displ.N 		);
	Integral.Construct ( Shapes_displ.j, Shapes_displ.W ); 
	Data_Pro_Plast.Construct ( Shapes_plast.dNdx 	);
	Data_Pro_Plast.Insert_N  ( Shapes_plast.N 		);

  	//-- The following don't allocate, that's done in Initialize()
	
  	Form_C_List			( APS_Matl );
  	Form_V_S_Lists		( np1,n );
  	Form_VB_List		( );
	Form_B_List			( );

}

//---------------------------------------------------------------------

void APS_kappa_alphaT::Form_LHS_Keps_Kd ( dMatrixT &Keps, dMatrixT &Kd )
{
 	Kd		= Integral.of( 	B_eps[kBgamma], B_d[kBmvgam1d] ); 
 	Kd		+= Integral.of( B_eps[kBgamma], B_d[kBmvgam2d] ); 
 	Kd		+= Integral.of( B_eps[kBgamma], B_d[kBmvgam3d] ); 
 	Kd  *= -1.0;
 	Keps	= Integral.of( 	B_eps[kBgamma], B_eps[kBgamma] 	);
 	Keps	-= Integral.of( B_eps[kBgamma], B_eps[kBmvgam1eps] );
 	Keps	-= Integral.of( B_eps[kBgamma], B_eps[kBmvgam2eps] );
 	Keps	-= Integral.of( B_eps[kBgamma], B_eps[kBmvgam3eps] );
}

//---------------------------------------------------------------------

void APS_kappa_alphaT::Form_RHS_F_int ( dArrayT &F_int ) 
{
	F_int   = Integral.of	( B_eps[kBgamma], V[kdel_gammap] );  	
	F_int  -= Integral.of	( B_eps[kBgamma], S[kdel_gamma1], V[km1_bar] ); 
	F_int  -= Integral.of	( B_eps[kBgamma], S[kdel_gamma2], V[km2_bar] ); 
	F_int  -= Integral.of	( B_eps[kBgamma], S[kdel_gamma3], V[km3_bar] ); 
}

//=== Private =========================================================


void APS_kappa_alphaT::Form_C_List (APS_MaterialT *APS_Matl)
{
	C[kMu]    	 		= APS_Matl -> Retrieve ( APS_MatlT::kMu 			);
	C[kl]    	 		= APS_Matl -> Retrieve ( APS_MatlT::kl 				);
	C[kgamma0_dot_1]    = APS_Matl -> Retrieve ( APS_MatlT::kgamma0_dot_1	 	);
	C[kgamma0_dot_2]    = APS_Matl -> Retrieve ( APS_MatlT::kgamma0_dot_2	 	);
	C[kgamma0_dot_3]    = APS_Matl -> Retrieve ( APS_MatlT::kgamma0_dot_3	 	);
	C[km_rate]    	 	= APS_Matl -> Retrieve ( APS_MatlT::km_rate	 		);
	C[km1_x]    	 	= APS_Matl -> Retrieve ( APS_MatlT::km1_x	 			);
	C[km1_y]    	 	= APS_Matl -> Retrieve ( APS_MatlT::km1_y	 			);
	C[km2_x]    	 	= APS_Matl -> Retrieve ( APS_MatlT::km2_x	 			);
	C[km2_y]    	 	= APS_Matl -> Retrieve ( APS_MatlT::km2_y	 			);
	C[km3_x]    	 	= APS_Matl -> Retrieve ( APS_MatlT::km3_x	 			);
	C[km3_y]    	 	= APS_Matl -> Retrieve ( APS_MatlT::km3_y	 			);
	C[kH]    	 		= APS_Matl -> Retrieve ( APS_MatlT::kH 				);
	C[kkappa0_1]    	= APS_Matl -> Retrieve ( APS_MatlT::kkappa0_1 		);
	C[kkappa0_2]    	= APS_Matl -> Retrieve ( APS_MatlT::kkappa0_2 		);
	C[kkappa0_3]    	= APS_Matl -> Retrieve ( APS_MatlT::kkappa0_3 		);
	C[k1by3] 			= 1.0/3.0;
	C[k1byRoot3] 		= 1.0/sqrt(3.0);
	C[kRoot3by3]		= sqrt( 3.0 )/3.0; 
	C[ksmall]			= 1.0e-8; 
}


void APS_kappa_alphaT::Form_V_S_Lists (  APS_VariableT &npt, APS_VariableT &n )
{
	// get state variables at times np1 and n
	//V_state[kstate] = npt.Get ( APS::kstate );
	V_state[kstate_n] = n.Get ( APS::kstate );
	
	// get and set grad_u
	B_gradu[kgrad_u] = npt.Get ( APS::kgrad_u );
	V[kV_Temp4](0)=B_gradu[kgrad_u](0,0);
	V[kV_Temp4](1)=B_gradu[kgrad_u](0,1);
	
	// get gamma_p and grad_gamma_p
	V[kgammap] = npt.Get ( APS::kgammap );
	V[kgammap_n] = n.Get ( APS::kgammap );
	B_gradgammap[kgrad_gammap] = npt.Get ( APS::kgrad_gammap );
	
	// calculate anti-plane curl component of gamma_p
	#pragma message("APS_kappa_alphaT::Form_V_S_Lists: are these grads what I think they are? ")
	S[kgammap_curl] = B_gradgammap[kgrad_gammap](1,0)-B_gradgammap[kgrad_gammap](0,1);
	// output mag of curl in out of plane direction
	V_out[kstressstate](2) = S[kgammap_curl];
		
	// calculate elastic gradient vector
	V[kV_Temp1].DiffOf( V[kV_Temp4], V[kgammap] ); 	
	
	// output displacement gradient vector
	V_out[kstressstate](0) = V[kV_Temp4](0);
	V_out[kstressstate](1) = V[kV_Temp4](1);	
	
	// calculate s13
	V_out[kstressstate](4) = V[kV_Temp1](0);
	V_out[kstressstate](4) *= C[kMu];
	// calculate s23
	V_out[kstressstate](5) = V[kV_Temp1](1);
	V_out[kstressstate](5) *= C[kMu];
	// calculate J2
	V[kV_Temp2](0) = V_out[kstressstate](4);
	V[kV_Temp2](1) = V_out[kstressstate](5);
	V[kV_Temp2].Magnitude(S[kS_Temp11]);
	V_out[kstressstate](6) = S[kS_Temp11];
	
	// calculate effective strain
	V_out[kV_out_Temp1](0) = V[kV_Temp4](0);
	V_out[kV_out_Temp1](0) *= V[kV_Temp4](0);
	V_out[kV_out_Temp1](1) = V[kV_Temp4](1);
	V_out[kV_out_Temp1](1) *= V[kV_Temp4](1);
	double dum2=sqrt(2.0);
	V_out[kV_out_Temp1](2) = V[kV_Temp4](0);
	V_out[kV_out_Temp1](2) *= V[kV_Temp4](1);
	V_out[kV_out_Temp1](2) *= dum2;
	V_out[kV_out_Temp1](3) = V[kV_Temp4](0);
	V_out[kV_out_Temp1](3) *= dum2;
	V_out[kV_out_Temp1](4) = V[kV_Temp4](1);
	V_out[kV_out_Temp1](4) *= dum2;
	V_out[kV_out_Temp1](5) = 0.0;
	V_out[kV_out_Temp1](6) = 0.0;
	V_out[kV_out_Temp1](7) = 0.0;
	V_out[kV_out_Temp1](8) = 0.0;
	V_out[kV_out_Temp1](9) = 0.0;
	V_out[kV_out_Temp1](10) = 0.0;
	V_out[kV_out_Temp1](11) = 0.0;
	V_out[kV_out_Temp1](12) = 0.0;
	V_out[kV_out_Temp1](13) = 0.0;
	V_out[kV_out_Temp1](14) = 0.0;
	V_out[kV_out_Temp1](15) = 0.0;
	V_out[kV_out_Temp1].Magnitude(S[kS_Temp8]);
	double dum3=sqrt(3.0);
	V_out[kstressstate](3) = S[kS_Temp8];
	V_out[kstressstate](3) /= dum3;

	// calculate \Delta \bgamma^p
	V[kdel_gammap].DiffOf(  V[kgammap], V[kgammap_n] ); 
	// calculate mag of \Delta \bgamma^p
	V[kdel_gammap].Magnitude(S[kmag_del_gammap]);
	
	//calculate constant to be used in solution of del_gamma^alpha of kappa^alpha
	S[kS_Temp13] = C[kRoot3by3];
	S[kS_Temp13] *= C[kH];
	S[kS_Temp13] *= C[kMu];
	
	
	/* slip system 1 */
	
	S[kIV_kappa1_n] = V_state[kstate_n](1);
	S[kdel_gamma1_n] = V_state[kstate_n](2);
	
	V[km1_bar](0) = C[km1_x];
	V[km1_bar](1) = C[km1_y];
	V[km1_bar].Magnitude(S[kmag_m1]);

	// calculate resolved stress on slip system
	V[kV_Temp1].Dot( V[km1_bar], S[kS_1] );
	S[kS_1]   *= C[kMu];
	
	// calculate relative stress
	S[kxi_1] = S[kS_1];
	S[kxi_1] *= S[kxi_curl_term];
	//S[kxi].Macaulay( S[kxi] );
	//if (S[kxi_1]<0.0) S[kxi_1] = 0.0;
	S[ksignxi_1] = S[kxi_1];
	S[kabs_xi_1].Abs( S[kxi_1] );
	if ( S[kabs_xi_1] > C[ksmall]) S[ksignxi_1] /= S[kabs_xi_1];
	
	// save relative stress
	V_state[kstate](0) = S[kxi_1];
	
	// output relative stress
	V_out[kstressstate](7) = S[kxi_1];
	
	// calculate sign of del_gamma1
	S[ksign_del_gamma1] = S[kdel_gamma1_n];
	S[kS_Temp10].Abs( S[kdel_gamma1_n] );
	if ( S[kS_Temp10] > C[ksmall]) S[ksign_del_gamma1] /= S[kS_Temp10];
	
	//calculate initial residual for solving for del_gamma1
	S[kS_Temp10] = S[kdel_gamma1_n];
	S[kS_Temp9] = C[kkappa0_1];
	S[kS_Temp9] += S[kIV_kappa1_n];
	S[kS_Temp1] = S[kabs_xi_1];
	S[kS_Temp1] /= S[kS_Temp9];
	S[kS_Temp12].Pow(S[kS_Temp1],C[km_rate]);
	S[kS_Temp12] *= C[kgamma0_dot_1];
	S[kS_Temp12] *= delta_t;
	S[kS_Temp12] *= S[ksignxi_1];
	S[kS_Temp10] -= S[kS_Temp12];
	S[kS_Temp11].Abs( S[kS_Temp10] );
	
	//initialize deldel_gamma1 and del_gamma1 before iteration
	S[kdeldel_gamma1] = 0.0;
	S[kdel_gamma1] = S[kdel_gamma1_n];
	
	/* iterate to solve for del_gamma1 */
	while (S[kS_Temp11] > C[ksmall]) {

	//calculate jacobian of Newton-Raphson iteration
	S[kS_Temp12] *= S[kS_Temp13];
	S[kS_Temp12] *= C[km_rate];
	S[kS_Temp12] *= S[kmag_m1];
	S[kS_Temp12] *= S[ksign_del_gamma1];
	S[kS_Temp12] /= S[kS_Temp9];
	S[kS_Temp14] = 1.0;
	S[kS_Temp14] += S[kS_Temp12];
	
	// calculate deldel_gamma1
	S[kdeldel_gamma1] = S[kS_Temp10];
	S[kdeldel_gamma1] *= -1.0;
	S[kdeldel_gamma1] /= S[kS_Temp14];
	
	// update del_gamma1
	S[kdel_gamma1] += S[kdeldel_gamma1];
	
	// calculate sign of del_gamma1
	S[ksign_del_gamma1] = S[kdel_gamma1];
	S[kS_Temp10].Abs( S[kdel_gamma1] );
	if ( S[kS_Temp10] > C[ksmall]) S[ksign_del_gamma1] /= S[kS_Temp10];
	
	// update kappa1
	S[kIV_kappa1] = S[kIV_kappa1_n];
	S[kS_Temp12] = S[kS_Temp13];
	S[kS_Temp12] *= S[kmag_m1];
	S[kS_Temp5].Abs( S[kdel_gamma1] );
	S[kS_Temp12] *= S[kS_Temp5];
	S[kIV_kappa1] += S[kS_Temp12];
	
	// update residual
	S[kS_Temp10] = S[kdel_gamma1];
	S[kS_Temp9] = C[kkappa0_1];
	S[kS_Temp9] += S[kIV_kappa1];
	S[kS_Temp1] = S[kabs_xi_1];
	S[kS_Temp1] /= S[kS_Temp9];
	S[kS_Temp12].Pow(S[kS_Temp1],C[km_rate]);
	S[kS_Temp12] *= C[kgamma0_dot_1];
	S[kS_Temp12] *= delta_t;
	S[kS_Temp12] *= S[ksignxi_1];
	S[kS_Temp10] -= S[kS_Temp12];
	S[kS_Temp11].Abs( S[kS_Temp10] );
	
	}
	/* end iteration */
	
	//output gamma1_dot
	V_out[kstressstate](9) = S[kdel_gamma1];
	V_out[kstressstate](9) /= delta_t;
	
	// save increment of slip along slip system
	V_state[kstate](2) = S[kdel_gamma1];
	
	// save kappa1
	V_state[kstate](1) = S[kIV_kappa1];
	
	// output kappa1
	V_out[kstressstate](8) = C[kkappa0_1];
	V_out[kstressstate](8) += S[kIV_kappa1];
	


	/* slip system 2 */
		
	S[kIV_kappa2_n] = V_state[kstate_n](4);
	S[kdel_gamma2_n] = V_state[kstate_n](5);
	
	V[km2_bar](0) = C[km2_x];
	V[km2_bar](1) = C[km2_y];
	V[km2_bar].Magnitude(S[kmag_m2]);

	// calculate resolved stress on slip system
	V[kV_Temp1].Dot( V[km2_bar], S[kS_2] );
	S[kS_2]   *= C[kMu];
	
	// calculate relative stress
	S[kxi_2] = S[kS_2];
	S[kxi_2] *= S[kxi_curl_term];
	//S[kxi].Macaulay( S[kxi] );
	//if (S[kxi_2]<0.0) S[kxi_2] = 0.0;
	S[ksignxi_2] = S[kxi_2];
	S[kabs_xi_2].Abs( S[kxi_2] );
	if ( S[kabs_xi_2] > C[ksmall]) S[ksignxi_2] /= S[kabs_xi_2];
	
	// save relative stress
	V_state[kstate](3) = S[kxi_2];
	
	// output relative stress
	V_out[kstressstate](10) = S[kxi_2];
	
	// calculate sign of del_gamma2
	S[ksign_del_gamma2] = S[kdel_gamma2_n];
	S[kS_Temp10].Abs( S[kdel_gamma2_n] );
	if ( S[kS_Temp10] > C[ksmall]) S[ksign_del_gamma2] /= S[kS_Temp10];
	
	//calculate initial residual for solving for del_gamma2
	S[kS_Temp10] = S[kdel_gamma2_n];
	S[kS_Temp9] = C[kkappa0_2];
	S[kS_Temp9] += S[kIV_kappa2_n];
	S[kS_Temp1] = S[kabs_xi_2];
	S[kS_Temp1] /= S[kS_Temp9];
	S[kS_Temp12].Pow(S[kS_Temp1],C[km_rate]);
	S[kS_Temp12] *= C[kgamma0_dot_2];
	S[kS_Temp12] *= delta_t;
	S[kS_Temp12] *= S[ksignxi_2];
	S[kS_Temp10] -= S[kS_Temp12];
	S[kS_Temp11].Abs( S[kS_Temp10] );
	
	//initialize deldel_gamma2 and del_gamma2 before iteration
	S[kdeldel_gamma2] = 0.0;
	S[kdel_gamma2] = S[kdel_gamma2_n];
	
	/* iterate to solve for del_gamma2 */
	while (S[kS_Temp11] > C[ksmall]) {

	//calculate jacobian of Newton-Raphson iteration
	S[kS_Temp12] *= S[kS_Temp13];
	S[kS_Temp12] *= C[km_rate];
	S[kS_Temp12] *= S[kmag_m2];
	S[kS_Temp12] *= S[ksign_del_gamma2];
	S[kS_Temp12] /= S[kS_Temp9];
	S[kS_Temp14] = 1.0;
	S[kS_Temp14] += S[kS_Temp12];
	
	// calculate deldel_gamma2
	S[kdeldel_gamma2] = S[kS_Temp10];
	S[kdeldel_gamma2] *= -1.0;
	S[kdeldel_gamma2] /= S[kS_Temp14];
	
	// update del_gamma2
	S[kdel_gamma2] += S[kdeldel_gamma2];
	
	// calculate sign of del_gamma2
	S[ksign_del_gamma2] = S[kdel_gamma2];
	S[kS_Temp10].Abs( S[kdel_gamma2] );
	if ( S[kS_Temp10] > C[ksmall]) S[ksign_del_gamma2] /= S[kS_Temp10];
	
	// update kappa2
	S[kIV_kappa2] = S[kIV_kappa2_n];
	S[kS_Temp12] = S[kS_Temp13];
	S[kS_Temp12] *= S[kmag_m2];
	S[kS_Temp5].Abs( S[kdel_gamma2] );
	S[kS_Temp12] *= S[kS_Temp5];
	S[kIV_kappa2] += S[kS_Temp12];
	
	// update residual
	S[kS_Temp10] = S[kdel_gamma2];
	S[kS_Temp9] = C[kkappa0_2];
	S[kS_Temp9] += S[kIV_kappa2];
	S[kS_Temp1] = S[kabs_xi_2];
	S[kS_Temp1] /= S[kS_Temp9];
	S[kS_Temp12].Pow(S[kS_Temp1],C[km_rate]);
	S[kS_Temp12] *= C[kgamma0_dot_2];
	S[kS_Temp12] *= delta_t;
	S[kS_Temp12] *= S[ksignxi_2];
	S[kS_Temp10] -= S[kS_Temp12];
	S[kS_Temp11].Abs( S[kS_Temp10] );
	
	}
	/* end iteration */
	
	//output gamma2_dot
	V_out[kstressstate](12) = S[kdel_gamma2];
	V_out[kstressstate](12) /= delta_t;
	
	// save increment of slip along slip system
	V_state[kstate](5) = S[kdel_gamma2];
	
	// save kappa2
	V_state[kstate](4) = S[kIV_kappa2];
	
	// output kappa2
	V_out[kstressstate](11) = C[kkappa0_2];
	V_out[kstressstate](11) += S[kIV_kappa2];
	

	
	
		
	/* slip system 3 */
	
	S[kIV_kappa3_n] = V_state[kstate_n](7);
	S[kdel_gamma3_n] = V_state[kstate_n](8);
	
	V[km3_bar](0) = C[km3_x];
	V[km3_bar](1) = C[km3_y];
	V[km3_bar].Magnitude(S[kmag_m3]);

	// calculate resolved stress on slip system
	V[kV_Temp1].Dot( V[km3_bar], S[kS_3] );
	S[kS_3]   *= C[kMu];
	
	// calculate relative stress
	S[kxi_3] = S[kS_3];
	S[kxi_3] *= S[kxi_curl_term];
	//S[kxi].Macaulay( S[kxi] );
	//if (S[kxi_3]<0.0) S[kxi_3] = 0.0;
	S[ksignxi_3] = S[kxi_3];
	S[kabs_xi_3].Abs( S[kxi_3] );
	if ( S[kabs_xi_3] > C[ksmall]) S[ksignxi_3] /= S[kabs_xi_3];
	
	// save relative stress
	V_state[kstate](6) = S[kxi_3];
	
	// output relative stress
	V_out[kstressstate](13) = S[kxi_3];
	
	// calculate sign of del_gamma3
	S[ksign_del_gamma3] = S[kdel_gamma3_n];
	S[kS_Temp10].Abs( S[kdel_gamma3_n] );
	if ( S[kS_Temp10] > C[ksmall]) S[ksign_del_gamma3] /= S[kS_Temp10];
	
	//calculate initial residual for solving for del_gamma3
	S[kS_Temp10] = S[kdel_gamma3_n];
	S[kS_Temp9] = C[kkappa0_3];
	S[kS_Temp9] += S[kIV_kappa3_n];
	S[kS_Temp1] = S[kabs_xi_3];
	S[kS_Temp1] /= S[kS_Temp9];
	S[kS_Temp12].Pow(S[kS_Temp1],C[km_rate]);
	S[kS_Temp12] *= C[kgamma0_dot_3];
	S[kS_Temp12] *= delta_t;
	S[kS_Temp12] *= S[ksignxi_3];
	S[kS_Temp10] -= S[kS_Temp12];
	S[kS_Temp11].Abs( S[kS_Temp10] );
	
	//initialize deldel_gamma3 and del_gamma3 before iteration
	S[kdeldel_gamma3] = 0.0;
	S[kdel_gamma3] = S[kdel_gamma3_n];
	
	/* iterate to solve for del_gamma1 */
	while (S[kS_Temp11] > C[ksmall]) {

	//calculate jacobian of Newton-Raphson iteration
	S[kS_Temp12] *= S[kS_Temp13];
	S[kS_Temp12] *= C[km_rate];
	S[kS_Temp12] *= S[kmag_m3];
	S[kS_Temp12] *= S[ksign_del_gamma3];
	S[kS_Temp12] /= S[kS_Temp9];
	S[kS_Temp14] = 1.0;
	S[kS_Temp14] += S[kS_Temp12];
	
	// calculate deldel_gamma3
	S[kdeldel_gamma3] = S[kS_Temp10];
	S[kdeldel_gamma3] *= -1.0;
	S[kdeldel_gamma3] /= S[kS_Temp14];
	
	// update del_gamma3
	S[kdel_gamma3] += S[kdeldel_gamma3];
	
	// calculate sign of del_gamma3
	S[ksign_del_gamma3] = S[kdel_gamma3];
	S[kS_Temp10].Abs( S[kdel_gamma3] );
	if ( S[kS_Temp10] > C[ksmall]) S[ksign_del_gamma3] /= S[kS_Temp10];
	
	// update kappa3
	S[kIV_kappa3] = S[kIV_kappa3_n];
	S[kS_Temp12] = S[kS_Temp13];
	S[kS_Temp12] *= S[kmag_m3];
	S[kS_Temp5].Abs( S[kdel_gamma3] );
	S[kS_Temp12] *= S[kS_Temp5];
	S[kIV_kappa3] += S[kS_Temp12];
	
	// update residual
	S[kS_Temp10] = S[kdel_gamma3];
	S[kS_Temp9] = C[kkappa0_3];
	S[kS_Temp9] += S[kIV_kappa3];
	S[kS_Temp1] = S[kabs_xi_3];
	S[kS_Temp1] /= S[kS_Temp9];
	S[kS_Temp12].Pow(S[kS_Temp1],C[km_rate]);
	S[kS_Temp12] *= C[kgamma0_dot_3];
	S[kS_Temp12] *= delta_t;
	S[kS_Temp12] *= S[ksignxi_3];
	S[kS_Temp10] -= S[kS_Temp12];
	S[kS_Temp11].Abs( S[kS_Temp10] );
	
	}
	/* end iteration */
	
	//output gamma3_dot
	V_out[kstressstate](15) = S[kdel_gamma3];
	V_out[kstressstate](15) /= delta_t;
	
	// save increment of slip along slip system
	V_state[kstate](8) = S[kdel_gamma3];
	
	// save kappa3
	V_state[kstate](7) = S[kIV_kappa3];
	
	// output kappa3
	V_out[kstressstate](14) = C[kkappa0_3];
	V_out[kstressstate](14) += S[kIV_kappa3];
	

	
	/* put current state variables */
	npt.Put ( APS::kstate, V_state[kstate] );
	
}


		
void APS_kappa_alphaT::Form_VB_List (void)
{					
		Data_Pro_Displ.APS_B(B_d[kB]);
		
 		Data_Pro_Plast.APS_Ngamma(B_eps[kBgamma]);
		Data_Pro_Plast.APS_Ngam1d2(VB_eps[kNgam_xdy]);
		Data_Pro_Plast.APS_Ngam2d1(VB_eps[kNgam_ydx]);
		
		/* calculate contribution of curl to relative stress */
		S[kxi_curl_term] = C[kl];
		S[kxi_curl_term] *= S[kgammap_curl];
		S[kxi_curl_term] /= 3.0;
		S[kxi_curl_term] += 1.0;
	
		/* slip system 1 */
		// calculate contributions to stiffness matricies
		S[kS_Temp8].Abs( S[kxi_1] );
		S[kS_Temp11] = C[kkappa0_1];
		S[kS_Temp11] += S[kIV_kappa1];
		S[kS_Temp8] /= S[kS_Temp11];
		S[kS_Temp9] = C[km_rate];
		S[kS_Temp9] -= 1.0;
		S[kS_Temp5].Pow(S[kS_Temp8], S[kS_Temp9]);
		S[kS_Temp6] = delta_t;	
		S[kS_Temp6] *= C[kgamma0_dot_1];
		S[kS_Temp6] *= C[km_rate];
		S[kS_Temp6] *= C[kMu];
		S[kS_Temp6] /= S[kS_Temp11];
		S[kA1] = S[kS_Temp5];
		S[kA1] *= S[kS_Temp6];
		
		S[kS_Temp8] = S[kxi_1];
		S[kS_Temp8] /= S[kS_Temp11];
		S[kB1] = C[kH];
		S[kB1] *= C[kRoot3by3];
		S[kB1] *= S[kmag_m1];
		S[kB1] *= S[ksign_del_gamma1];
		S[kB1] *= S[kS_Temp8];
		
		S[kS_Temp3] = S[kA1];
		S[kS_Temp3] *= S[kB1];
		S[kS_Temp3] += 1.0;
		S[kS_Temp2] = S[kA1];
		S[kS_Temp2] /= S[kS_Temp3];
		
 		V[km1_bar].Dot( B_d[kB], VB_d[kVdelgam1d] );
 		VB_d[kVdelgam1d] *= S[kS_Temp2];
 		VB_d[kVdelgam1d] *= S[kxi_curl_term];
 		
 		V[km1_bar].Dot( B_eps[kBgamma], VB_eps[kVB_eps_Temp4] );
		VB_eps[kVB_eps_Temp4] *= -1.0;
		VB_eps[kVB_eps_Temp4] *= S[kxi_curl_term];
 		S[kS_Temp3] = S[kS_1];
 		S[kS_Temp3] *= C[kl];
 		S[kS_Temp3] /= C[kMu];
 		S[kS_Temp3] /= 3.0; 
 		VB_eps[kVB_eps_Temp1].DiffOf( VB_eps[kNgam_ydx], VB_eps[kNgam_xdy] ); 
 		VB_eps[kVB_eps_Temp1] *= S[kS_Temp3];
 		VB_eps[kVdelgam1eps].SumOf(  VB_eps[kVB_eps_Temp4], VB_eps[kVB_eps_Temp1] ); 
 		VB_eps[kVdelgam1eps] *= S[kS_Temp2];
 		
 		
 		/* slip system 2 */
 		// calculate contributions to stiffness matricies
		S[kS_Temp8].Abs( S[kxi_2] );
		S[kS_Temp11] = C[kkappa0_2];
		S[kS_Temp11] += S[kIV_kappa2];
		S[kS_Temp8] /= S[kS_Temp11];
		S[kS_Temp9] = C[km_rate];
		S[kS_Temp9] -= 1.0;
		S[kS_Temp5].Pow(S[kS_Temp8], S[kS_Temp9]);
		S[kS_Temp6] = delta_t;	
		S[kS_Temp6] *= C[kgamma0_dot_2];
		S[kS_Temp6] *= C[km_rate];
		S[kS_Temp6] *= C[kMu];
		S[kS_Temp6] /= S[kS_Temp11];
		S[kA2] = S[kS_Temp5];
		S[kA2] *= S[kS_Temp6];
		
		S[kS_Temp8] = S[kxi_2];
		S[kS_Temp8] /= S[kS_Temp11];
		S[kB2] = C[kH];
		S[kB2] *= C[kRoot3by3];
		S[kB2] *= S[kmag_m2];
		S[kB2] *= S[ksign_del_gamma2];
		S[kB2] *= S[kS_Temp8];
		
		S[kS_Temp3] = S[kA2];
		S[kS_Temp3] *= S[kB2];
		S[kS_Temp3] += 1.0;
		S[kS_Temp2] = S[kA2];
		S[kS_Temp2] /= S[kS_Temp3];
		
 		V[km2_bar].Dot( B_d[kB], VB_d[kVdelgam2d] );
 		VB_d[kVdelgam2d] *= S[kS_Temp2];
 		VB_d[kVdelgam2d] *= S[kxi_curl_term];
 		
 		V[km2_bar].Dot( B_eps[kBgamma], VB_eps[kVB_eps_Temp4] );
		VB_eps[kVB_eps_Temp4] *= -1.0;
		VB_eps[kVB_eps_Temp4] *= S[kxi_curl_term];
 		S[kS_Temp3] = S[kS_2];
 		S[kS_Temp3] *= C[kl];
 		S[kS_Temp3] /= C[kMu];
 		S[kS_Temp3] /= 3.0; 
 		VB_eps[kVB_eps_Temp1].DiffOf( VB_eps[kNgam_ydx], VB_eps[kNgam_xdy] ); 
 		VB_eps[kVB_eps_Temp1] *= S[kS_Temp3];
 		VB_eps[kVdelgam2eps].SumOf(  VB_eps[kVB_eps_Temp4], VB_eps[kVB_eps_Temp1] ); 
 		VB_eps[kVdelgam2eps] *= S[kS_Temp2];
 			
 		 		
 		/* slip system 3 */
 		// calculate contributions to stiffness matricies
		S[kS_Temp8].Abs( S[kxi_3] );
		S[kS_Temp11] = C[kkappa0_3];
		S[kS_Temp11] += S[kIV_kappa3];
		S[kS_Temp8] /= S[kS_Temp11];
		S[kS_Temp9] = C[km_rate];
		S[kS_Temp9] -= 1.0;
		S[kS_Temp5].Pow(S[kS_Temp8], S[kS_Temp9]);
		S[kS_Temp6] = delta_t;	
		S[kS_Temp6] *= C[kgamma0_dot_3];
		S[kS_Temp6] *= C[km_rate];
		S[kS_Temp6] *= C[kMu];
		S[kS_Temp6] /= S[kS_Temp11];
		S[kA3] = S[kS_Temp5];
		S[kA3] *= S[kS_Temp6];
		
		S[kS_Temp8] = S[kxi_3];
		S[kS_Temp8] /= S[kS_Temp11];
		S[kB3] = C[kH];
		S[kB3] *= C[kRoot3by3];
		S[kB3] *= S[kmag_m3];
		S[kB3] *= S[ksign_del_gamma3];
		S[kB3] *= S[kS_Temp8];
		
		S[kS_Temp3] = S[kA3];
		S[kS_Temp3] *= S[kB3];
		S[kS_Temp3] += 1.0;
		S[kS_Temp2] = S[kA3];
		S[kS_Temp2] /= S[kS_Temp3];
		
 		V[km3_bar].Dot( B_d[kB], VB_d[kVdelgam3d] );
 		VB_d[kVdelgam3d] *= S[kS_Temp2];
 		VB_d[kVdelgam3d] *= S[kxi_curl_term];
 		
 		V[km3_bar].Dot( B_eps[kBgamma], VB_eps[kVB_eps_Temp4] );
		VB_eps[kVB_eps_Temp4] *= -1.0;
		VB_eps[kVB_eps_Temp4] *= S[kxi_curl_term];
 		S[kS_Temp3] = S[kS_3];
 		S[kS_Temp3] *= C[kl];
 		S[kS_Temp3] /= C[kMu];
 		S[kS_Temp3] /= 3.0; 
 		VB_eps[kVB_eps_Temp1].DiffOf( VB_eps[kNgam_ydx], VB_eps[kNgam_xdy] ); 
 		VB_eps[kVB_eps_Temp1] *= S[kS_Temp3];
 		VB_eps[kVdelgam3eps].SumOf(  VB_eps[kVB_eps_Temp4], VB_eps[kVB_eps_Temp1] ); 
 		VB_eps[kVdelgam3eps] *= S[kS_Temp2];
 		
}


void APS_kappa_alphaT::Form_B_List (void)
{
 		B_d[kBmvgam1d].Outer( V[km1_bar], VB_d[kVdelgam1d] );
 		B_d[kBmvgam2d].Outer( V[km2_bar], VB_d[kVdelgam2d] );
 		B_d[kBmvgam3d].Outer( V[km3_bar], VB_d[kVdelgam3d] );
 		B_eps[kBmvgam1eps].Outer( V[km1_bar], VB_eps[kVdelgam1eps] );	
 		B_eps[kBmvgam2eps].Outer( V[km2_bar], VB_eps[kVdelgam2eps] );
 		B_eps[kBmvgam3eps].Outer( V[km3_bar], VB_eps[kVdelgam3eps] );		
}


//##############################################################################################

void APS_kappa_alphaT::Get ( StringT &Name, FEA_dMatrixT &tensor )
{
	if ( Name == "grad_u" )
		tensor = B_gradu[kgrad_u];
	else if ( Name == "grad_gammap" )
		tensor = B_gradgammap[kgrad_gammap];
	else
		cout << " ...ERROR: APS_kappa_alphaT::Get() >> Unknown tensor '"<<Name<<"' requested. \n";
}

//##############################################################################################

void APS_kappa_alphaT::Get ( StringT &Name, FEA_dVectorT &vector )
{
	if ( Name == "gammap" )
		vector = V[kgammap];
	else if ( Name == "out" )
		vector = V_out[kstressstate];	
	else
		cout << " ...ERROR: APS_kappa_alphaT::Get() >> Unknown vector '"<<Name<<"' requested. \n";
}

//##################################################################################

void APS_kappa_alphaT::Get ( StringT &Name, FEA_dScalarT &scalar )
{
	if ( Name == "kappa" )
		scalar = S[kIV_kappa1];
	else
		cout << " ...ERROR: APS_kappa_alphaT::Get() >> Unknown scalar '"<<Name<<"' requested. \n";
}


