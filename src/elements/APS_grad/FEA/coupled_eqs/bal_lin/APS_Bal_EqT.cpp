// $Id: APS_Bal_EqT.cpp,v 1.3 2003-09-04 15:45:41 paklein Exp $
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

	Integral.Construct ( Shapes.j, Shapes.W ); 

}

//---------------------------------------------------------------------

void APS_Bal_EqT::Form_LHS_Keps_Kd	( dMatrixT &Keps, dMatrixT &Kd )  // Untested
{
		Keps 	= Integral.of( 	B[kB], C[kMu], B[kBgamma] );  	
		Keps 	*= -1.0;
	 	Kd  	= Integral.of( 	B[kB], C[kMu], B[kB] 	);  	
//	Kd		-= Integral.of( VB[kN], C[kMu], VB[knuB] 	); 	
#pragma message("APS_Bal_EqT::Form_LHS_Keps_Kd: there is no (FEA_dVectorT, double FEA_dVectorT) version of this method")
}

//---------------------------------------------------------------------

void APS_Bal_EqT::Form_RHS_F_int ( dArrayT &F_int ) // Untested
{
#pragma message("APS_Bal_EqT::Form_RHS_F_int: tmp, Kepstmp, and tmp2 are not defined")

#if 0
//		F_int = Integral.of	( VB[kN], traction?? );
//		dArrayT tmp1.MultAb (Kd, d???);
		F_int -= tmp;
		Kepstmp = Keps - Integral.of( 	VB[kN], C[kMu], VB[knuNgam] );
//		dArrayT tmp2.MultAb ( Kepstmp, eps??? );
		F_int += tmp2;
#endif
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
		
// 		VB[knuB].Dot( V[knueps]??, B[kB] );
// 		VB[knuNgam].Dot( V[knueps]??, B[kBgamma] );
}
							

void APS_Bal_EqT::Form_C_List (APS_MaterialT *Shear_Matl)
{
		C.Dimension 	( kNUM_C_TERMS );
		C[kMu] 	= Shear_Matl -> Retrieve ( Shear_MatlT::kMu );
}

