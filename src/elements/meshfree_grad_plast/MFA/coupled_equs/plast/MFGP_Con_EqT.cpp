// $Id: MFGP_Con_EqT.cpp
#include "MFGP_Con_EqT.h"

using namespace Tahoe;

/* constructor */
MFGP_Con_EqT::MFGP_Con_EqT(void):
	fData_Pro_Displ(NULL),
	fData_Pro_Plast(NULL)
{

}

/* destructor */
MFGP_Con_EqT::~MFGP_Con_EqT(void)
{
	delete fData_Pro_Displ;
	delete fData_Pro_Plast;
}

/* set dims, derivs, and variables needed */
void MFGP_Con_EqT::Initialize (int &curr_ip, D3MeshFreeShapeFunctionT *Shapes_displ, D3MeshFreeShapeFunctionT *Shapes_plast, 
							GRAD_MRSSKStV *GRAD_MR_Plast_Mat,					
							int &fTime_Step, double fdelta_t) 
{
	//n_ip = np1.fVars_vector[0].IPs(); 
	/*n_rows_vector = np1.fVars_vector[0].Rows(); 
	n_rows_matrix = np1.fVars_matrix[0].Rows(); 
	n_cols_matrix = np1.fVars_matrix[0].Cols(); 
	*/
	//#pragma message("MFGP_Con_EqT::Construct: FEA_dVectorT has no Cols() function")

	n_en_displ = Shapes_displ->Derivatives_U(curr_ip).MajorDim(); //??
	n_en_plast = Shapes_plast->Derivatives_U(curr_ip).MajorDim(); //??
	n_sd = Shapes_displ->Derivatives_U(curr_ip).MinorDim(); //??
	//n_sd = n_rows_vector;	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en_displ = n_sd * n_en_displ;
	n_sd_x_n_en_plast = n_sd * n_en_plast;

	delta_t = fdelta_t;
	
	fData_Pro_Displ = new MFGP_MFA_Data_Processor_DisplT;
	fData_Pro_Displ->Initialize ( Shapes_displ->Derivatives_U(curr_ip), Shapes_displ->DDDerivatives_U(curr_ip) );//??
	fData_Pro_Plast = new MFGP_MFA_Data_Processor_PlastT;
	fData_Pro_Plast->Initialize ( Shapes_displ->IPShapeU(curr_ip), Shapes_displ->DDerivatives_U(curr_ip) );//??

	stress = GRAD_MR_Plast_Mat->s_ij();
	yield = GRAD_MR_Plast_Mat->YieldF();
	moduli = GRAD_MR_Plast_Mat->c_ijkl();

	Form_C_List		( GRAD_MR_Plast_Mat );
	Form_B_List		(  );
	//Form_V_S_List	( np1 );
	//Form_VB_List	(  );
}


//---------------------------------------------------------------------

void MFGP_Con_EqT::Form_LHS_Ku_Klambda(dMatrixT &Ku, dMatrixT &Klambda )  
{
	n_rows = Ku.Rows();
	n_cols = Ku.Cols();
	dMatrixT Kutemp (n_rows,n_cols);
	Ku.MultATBC (phi_lam, Clamu1, B1_d );
	Kutemp.MultATBC ( phi_lam, Clamu2, B3_d );
	Ku += Kutemp;

	n_rows = Klambda.Rows();
	n_cols = Klambda.Cols();
	dMatrixT Klambdatemp (n_rows,n_cols);
	Klambda.MultATBC (phi_lam, Clamlam1, phi_lam  );
	Klambdatemp.MultATBC ( phi_lam, Clamlam2, B4_lam );
	Klambda += Klambdatemp;	
}

//---------------------------------------------------------------------

void MFGP_Con_EqT::Form_RHS_F_int ( dArrayT &F_int) 
{
		//pass column of phi_lam to F_int
		F_int = phi_lam[0];
		F_int *= yield; 
		//updated failure function passed from the constitutive model
		//note: there is no F_int or F_ext in this case, but simply F
		//Reminder: change the constitutive model to pass YieldFunction to the
		//global level in addition to
		//fModuli/Cep and s_ij (Cauchy stress)
}



//=== Private =========================================================
	             				
void MFGP_Con_EqT::Form_B_List (void)
{
		B1_d.Dimension 	( n_sd, n_en_displ);
		B3_d.Dimension 	( n_sd, n_en_displ);
		int dum=1;
		phi_lam.Dimension ( dum, n_en_plast);
		B4_lam.Dimension ( dum, n_sd_x_n_en_plast);
		//NTS: check the allocation of phi_lam and B4
		//B_gradu.Dimension ( dum, n_sd); //B3??	
		
		fData_Pro_Displ->Set_B1(B1_d);
		fData_Pro_Displ->Set_B3(B3_d);
		fData_Pro_Plast->Set_phi(phi_lam);
 		fData_Pro_Plast->Set_B4(B4_lam); //B4 is scalar
}



void MFGP_Con_EqT::Form_C_List (GRAD_MRSSKStV *GRAD_MR_Plast)
{
		//C.Dimension( kNUM_C_TERMS );
		int i,j;
		Cuu1.Dimension(6,6); 
		Cuu2.Dimension(6,6);
		Culam1.Dimension(6,1); 
		Culam2.Dimension(6,1);
		Clamu1.Dimension(1,6); 
		Clamu2.Dimension(1,6);
		//Clamlam1.Dimension(1); 
		//Clamlam2.Dimension(1);
		
		//retrive Cgep/fmoduli, and form the 8 C matrices
		//pass the fModuli from the constitutive model
		for (i=0; i<=7; ++i)
		{
		 for (j=0; j<=14; ++j)
		 {
		  if (i<=5 & j<=5)
		  {
		   Cuu1(i,j)=moduli(i,j);
		  }
		  if (i<=5 & j>5 & j<=11)
		  {
		   Cuu2(i,j-6)=moduli(i,j);
		  }
		  if (i<=5 & j>11 & j<=12)
		  {
		   Culam1(i,j-12)=moduli(i,j);
		  }
		  if (i<=5 & j>12)
		  {
		   Culam2(i,j-13)=moduli(i,j);
		  }
		  if (i>5 & j<=5)
		  {
		   Clamu1(i-6,j)=moduli(i,j);
		  }
		  if (i>5 & j>5 & j<=11)
		  {
		   Clamu2(i-6,j-6)=moduli(i,j);
		  }
		  if (i>5 & j>11 & j<=12)
		  {
		   Clamlam1(i-6,j-12)=moduli(i,j);
		  }
		  if (i>5 & j>12)
		  {
		   Clamlam2(i-6,j-13)=moduli(i,j);
		  }
		 }
		}
}
