// $Id: APS_Bal_EqT.cpp,v 1.5 2003-09-19 00:47:04 raregue Exp $
#include "APS_Bal_EqT.h" 

using namespace Tahoe;

APS_Bal_EqT::APS_Bal_EqT ( FEA_ShapeFunctionT &Shapes, APS_MaterialT *Shear_Matl, APS_VariableT &np1, APS_VariableT &n, 
								int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes,Shear_Matl,np1,n,Integration_Scheme);
}

/* destructor */
//APS_Bal_EqT::~APS_Bal_EqT(void)
//{
//}

//---------------------------------------------------------------------

void APS_Bal_EqT::Construct ( FEA_ShapeFunctionT &Shapes, APS_MaterialT *Shear_Matl, APS_VariableT &np1, APS_VariableT &n, 
			int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Time_Integration_Scheme = Integration_Scheme;

	n_ip 		= np1.fVars[0].IPs(); 
	n_rows		= np1.fVars[0].Rows(); 

//	n_cols		= np1.fVars[0].Cols();
#pragma message("APS_Bal_EqT::Construct: FEA_dVectorT has no Cols() function")

	n_en    	= Shapes.dNdx.Cols();
	n_sd 		= n_rows;	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en = n_sd * n_en;

	if ( fTime_Step != time_step) { 	// New time step
		// Can put flow rule here: Eventually do this with all the kVar_n and make method Next_Step()
		time_step = fTime_Step;
	}

	delta_t = fdelta_t;
	
	Data_Pro.Construct ( Shapes.dNdx );

	Form_C_List		(	Shear_Matl );
	Form_B_List		(	);
	Form_VB_List	(	);
	Form_V_S_List		(	);

	Integral.Construct ( Shapes.j, Shapes.W ); 

}

//---------------------------------------------------------------------

void APS_Bal_EqT::Form_LHS_Keps_Kd	( dMatrixT &Keps, dMatrixT &Kd )  // Untested
{
		Keps 	= Integral.of( B[kB], C[kMu], B[kBgamma] );  
		Keps 	*= -1.0;
	 	Kd  	= Integral.of( B[kB], C[kMu], B[kB] );  	
#pragma message("APS_Bal_EqT::Form_LHS_Keps_Kd: this domain over Gamma_eps")
//		Kd		-= Integral.of( VB[kN], C[kMu], VB[knuB] );
}

//---------------------------------------------------------------------

void APS_Bal_EqT::Form_RHS_F_int ( dArrayT &F_int ) // Untested
{
#pragma message("APS_Bal_EqT::Form_RHS_F_int: how get grad_u and gammap?")
/*		F_int = Integral.of( B[kB], C[kMu], V[kgrad_u] ); 
		F_int -= Integral.of( B[kB], C[kMu], V[kgammap] ); */
#pragma message("APS_Bal_EqT::Form_RHS_F_int: this domain over Gamma_eps")
/*		F_int -= Integral.of( VB[kN], C[kMu], S[knuepsgradu] ); 
		F_int += Integral.of( VB[kN], C[kMu], S[knuepseps] );  */
}

//=== Private =========================================================
	             				
void APS_Bal_EqT::Form_B_List (void)
{
		B.Construct ( kNUM_B_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_en );	
		
		Data_Pro.APS_B(B[kB]);
 		Data_Pro.APS_Ngamma(B[kBgamma]);
}

		
void APS_Bal_EqT::Form_VB_List (void)
{					
		Data_Pro.APS_N(VB[kN]);

 		V[knueps].Dot( B[kB], VB[knuB] );
 		V[knueps].Dot( B[kBgamma], VB[knuNgam] ); 
}


void APS_Bal_EqT::Form_V_S_List (void)
{
#pragma message("APS_Bal_EqT::Form_V_S_List: V[knueps] and V[keps] must be input for BC")
		V[knueps](0) = 1.0;
		V[knueps](1) = 0.0;
		V[keps](0) = 1.0;
		V[keps](1) = 0.0;
#pragma message("APS_Bal_EqT::Form_V_S_List: how get grad_u?")
//		V[knueps].Dot( V[kgrad_u], S[knuepsgradu] );
		V[knueps].Dot( V[keps], S[knuepseps] );
}


void APS_Bal_EqT::Form_C_List (APS_MaterialT *Shear_Matl)
{
		C.Dimension 	( kNUM_C_TERMS );
		C[kMu] 	= Shear_Matl -> Retrieve ( Shear_MatlT::kMu );
}

