// $Id: APS_Bal_EqT.cpp,v 1.27 2006-09-16 15:39:23 regueiro Exp $
#include "APS_Bal_EqT.h" 

using namespace Tahoe;

/* constructors */

APS_Bal_EqT::APS_Bal_EqT ( int& nipsurf, int& nensurf, FEA_ShapeFunctionT &Shapes_displ, FEA_ShapeFunctionT &Shapes_plast, 
							APS_MaterialT *Shear_Matl, 
							APS_MaterialT *APS_Matl, APS_VariableT &np1, APS_VariableT &n, 
							int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Construct (nipsurf, nensurf, Shapes_displ, Shapes_plast, Shear_Matl, APS_Matl, np1, n, Integration_Scheme);
}

/* destructor */
//APS_Bal_EqT::~APS_Bal_EqT(void);


//---------------------------------------------------------------------

void APS_Bal_EqT::Construct ( int& nipsurf, int& nensurf, FEA_ShapeFunctionT &Shapes_displ, FEA_ShapeFunctionT &Shapes_plast, 
							APS_MaterialT *Shear_Matl,
							APS_MaterialT *APS_Matl, APS_VariableT &np1, APS_VariableT &n, 
							int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Time_Integration_Scheme = Integration_Scheme;

	n_ip 		= np1.fVars_vector[0].IPs(); 
	/*n_rows_vector = np1.fVars_vector[0].Rows(); 
	n_rows_matrix = np1.fVars_matrix[0].Rows(); 
	n_cols_matrix = np1.fVars_matrix[0].Cols(); 
	*/
	//#pragma message("APS_Bal_EqT::Construct: FEA_dVectorT has no Cols() function")

	n_en_displ    	= Shapes_displ.dNdx.Cols();
	n_en_plast    	= Shapes_plast.dNdx.Cols();
	n_sd    	= Shapes_displ.dNdx.Rows();
	//n_sd 		= n_rows_vector;	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en_displ = n_sd * n_en_displ;
	n_sd_x_n_en_plast = n_sd * n_en_plast;
	
	n_ip_surf=nipsurf;
	n_en_surf=nensurf;

	delta_t = fdelta_t;
	
	Data_Pro_Displ.Construct ( Shapes_displ.dNdx	);
	Data_Pro_Displ.Insert_N  ( Shapes_displ.N );
	Integral.Construct ( Shapes_displ.j, Shapes_displ.W ); 
	Data_Pro_Plast.Construct ( Shapes_plast.dNdx	);
	Data_Pro_Plast.Insert_N  ( Shapes_plast.N );
	
	Form_C_List		( Shear_Matl, APS_Matl );
	Form_B_List		(  );
	Form_V_S_List	( np1 );
	Form_VB_List	(  );
}

//---------------------------------------------------------------------

void APS_Bal_EqT::Form_LHS_Keps_Kd	( dMatrixT &Keps, dMatrixT &Kd )  
{
		//Keps 	= Integral.of( B_d[kB], C[kMu], B_eps[kBgamma] ); 
		Keps 	= Integral.of( B_d[kB], C[kOne], B_eps[kBgamma] );  
		Keps 	*= -1.0;
	 	//Kd  	= Integral.of( B_d[kB], C[kMu], B_d[kB] );  
	 	Kd  	= Integral.of( B_d[kB], C[kOne], B_d[kB] );  	
}

//---------------------------------------------------------------------

void APS_Bal_EqT::Form_RHS_F_int ( dArrayT &F_int, APS_VariableT &npt ) 
{
		B_gradu[kgrad_u] = npt.Get ( APS::kgrad_u );
		V[kV_Temp2](0)=B_gradu[kgrad_u](0,0);
		V[kV_Temp2](1)=B_gradu[kgrad_u](0,1);
		V[kgammap] = npt.Get ( APS::kgammap );
		/*
		F_int = Integral.of( B_d[kB], C[kMu], V[kV_Temp2] ); 
		F_int -= Integral.of( B_d[kB], C[kMu], V[kgammap] );
		*/
		F_int = Integral.of( B_d[kB], C[kOne], V[kV_Temp2] ); 
		F_int -= Integral.of( B_d[kB], C[kOne], V[kgammap] );
}


//---------------------------------------------------------------------

void APS_Bal_EqT::Form_LHS_Kd_Surf	( dMatrixT &Kd_face, FEA_SurfShapeFunctionT &SurfShapes )  
{		
		Data_Pro_Surf.Construct ( SurfShapes.dNdx	);
		Data_Pro_Surf.Insert_N_surf  ( SurfShapes.N );
		SurfIntegral.Construct ( SurfShapes.j, SurfShapes.W );
		
		Data_Pro_Surf.APS_B_surf(B_d_surf[kB_surf]);
		Data_Pro_Surf.APS_N(VB_d[kN]);

		V_surf[knueps] = SurfShapes.normal;
		//V[knueps](1) = SurfShapes.normal[1];

 		V_surf[knueps].Dot( B_d_surf[kB_surf], VB_d[knuB] ); 
 		
		//Kd_face	= SurfIntegral.of( VB_d[kN], C[kMu], VB_d[knuB] );
		Kd_face	= SurfIntegral.of( VB_d[kN], C[kOne], VB_d[knuB] );
		Kd_face	*= -1.0;
}

//---------------------------------------------------------------------

void APS_Bal_EqT::Form_RHS_F_int_Surf ( dArrayT &F_int_face, APS_VariableT &npt, double &wght  ) 
{
		V_surf[kgammap_surf] = npt.Get ( APS::kgammap_surf );
		V_surf[keps] = V_surf[kgammap_surf];
		V_surf[keps] *= wght;
		/*
		V[keps](0) = C[km1];
		V[keps](0) *= wght;
		V[keps](1) = C[km2];
		V[keps](1) *= wght; 
		*/
		
		B_gradu_surf[kgrad_u_surf] = npt.Get ( APS::kgrad_u_surf );
		V_surf[kV_surf_Temp2](0)=B_gradu_surf[kgrad_u_surf](0,0);
		V_surf[kV_surf_Temp2](1)=B_gradu_surf[kgrad_u_surf](0,1);
		V_surf[knueps].Dot( V_surf[kV_surf_Temp2], S[knuepsgradu] );
		V_surf[knueps].Dot( V_surf[keps], S[knuepseps] );
		
		//F_int_face = SurfIntegral.of( VB_d[kN], C[kMu], S[knuepsgradu] );
		F_int_face = SurfIntegral.of( VB_d[kN], C[kOne], S[knuepsgradu] );
		F_int_face *= -1.0;
		//F_int_face += SurfIntegral.of( VB_d[kN], C[kMu], S[knuepseps] );
		F_int_face += SurfIntegral.of( VB_d[kN], C[kOne], S[knuepseps] );
}



//=== Private =========================================================
	             				
void APS_Bal_EqT::Form_B_List (void)
{
		B_d.Construct 	( kNUM_B_d_TERMS, n_ip, n_sd, n_en_displ);
		B_d_surf.Construct 	( kNUM_B_d_surf_TERMS, n_ip_surf, n_sd, n_en_surf);
		B_eps.Construct ( kNUM_B_eps_TERMS, n_ip, n_sd, n_sd_x_n_en_plast);
		int dum=1;
		B_gradu.Construct ( kNUM_B_gradu_TERMS, n_ip, dum, n_sd);	
		B_gradu_surf.Construct ( kNUM_B_gradu_surf_TERMS, n_ip_surf, dum, n_sd);		
		
		Data_Pro_Displ.APS_B(B_d[kB]);
		
 		Data_Pro_Plast.APS_Ngamma(B_eps[kBgamma]);
}

		
void APS_Bal_EqT::Form_VB_List (void)
{
		VB_d.Construct 	( kNUM_VB_d_TERMS, n_ip_surf, n_en_surf );
		VB_eps.Construct ( kNUM_VB_eps_TERMS, n_ip, n_sd_x_n_en_plast );				
}


void APS_Bal_EqT::Form_V_S_List (APS_VariableT &npt)
{
		S.Construct 	( kNUM_S_TERMS, 	n_ip_surf 		);
		V.Construct 	( kNUM_V_TERMS, 	n_ip, n_sd 	);
		V_surf.Construct ( kNUM_V_surf_TERMS, 	n_ip_surf, n_sd 	);
		int dum=1;
		VS.Construct 	( kNUM_VS_TERMS, 	n_ip, dum 	);
}


void APS_Bal_EqT::Form_C_List (APS_MaterialT *Shear_Matl, APS_MaterialT *APS_Matl)
{
		C.Dimension 	( kNUM_C_TERMS );
		
		C[kMu] 	= Shear_Matl -> Retrieve ( Shear_MatlT::kMu );
		C[km1_x]  = APS_Matl -> Retrieve ( APS_MatlT::km1_x	);
		C[km1_y]  = APS_Matl -> Retrieve ( APS_MatlT::km1_y	);
		C[km2_x]  = APS_Matl -> Retrieve ( APS_MatlT::km2_x	);
		C[km2_y]  = APS_Matl -> Retrieve ( APS_MatlT::km2_y	);
		C[kOne] = 1.0;
}


void APS_Bal_EqT::Get ( StringT &Name, FEA_dMatrixT &tensor )
{
	if ( Name == "grad_u" )
		tensor = B_gradu[kgrad_u];
	else
		cout << " ...ERROR: APS_Bal_EqT::Get() >> Unknown tensor '"<<Name<<"' requested. \n";
}


void APS_Bal_EqT::Get ( StringT &Name, FEA_dVectorT &vector )
{
	if ( Name == "gammap" )
		vector = V[kgammap];
	else
		cout << " ...ERROR: APS_Bal_EqT::Get() >> Unknown vector '"<<Name<<"' requested. \n";
}

void APS_Bal_EqT::Get ( StringT &Name, FEA_dScalarT &scalar )
{
//for now, no scalars to get
/*	if ( Name == "?" )
		scalar = S[?];
	else if ( Name == "?" )
		scalar = S[?];
	else
		cout << " ...ERROR: APS_Bal_EqT::Get() >> Unknown scalar '"<<Name<<"' requested. \n";
		*/
}


