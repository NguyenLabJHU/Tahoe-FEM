// $Id: MFGP_Con_EqT.cpp
#include "MFGP_Con_EqT.h"

using namespace Tahoe;

/* constructor */
/*
MFGP_Con_EqT::MFGP_Con_EqT(void):
	fData_Pro_Displ(NULL),
	fData_Pro_Plast(NULL)
	*/
MFGP_Con_EqT::MFGP_Con_EqT(void) { }


/* destructor */
//MFGP_Con_EqT::~MFGP_Con_EqT(void) { }
/*
{
	delete fData_Pro_Displ;
	delete fData_Pro_Plast;
}
*/

/* set dims, derivs, and variables needed */
void MFGP_Con_EqT::Initialize(int& curr_ip, D3MeshFreeShapeFunctionT *Shapes_displ, D3MeshFreeShapeFunctionT *Shapes_plast, 
							GRAD_MRSSKStV *GRAD_MR_Plast_Mat,					
							int& fTime_Step, double fdelta_t) 
{
	//n_ip = np1.fVars_vector[0].IPs(); 
	/*n_rows_vector = np1.fVars_vector[0].Rows(); 
	n_rows_matrix = np1.fVars_matrix[0].Rows(); 
	n_cols_matrix = np1.fVars_matrix[0].Cols(); 
	*/
	//#pragma message("MFGP_Con_EqT::Construct: FEA_dVectorT has no Cols() function")

	n_en_displ = Shapes_displ->Derivatives_U(curr_ip).MinorDim(); // dN:[nsd]x[nnd] 
	n_en_plast = Shapes_plast->Derivatives_U(curr_ip).MinorDim(); 
	n_sd = Shapes_displ->Derivatives_U(curr_ip).MajorDim(); 
	//n_sd = n_rows_vector;
	n_str = dSymMatrixT::NumValues(n_sd);	
	n_sd_x_n_en_displ = n_sd * n_en_displ;
	n_sd_x_n_en_plast = n_sd * n_en_plast;

	delta_t = fdelta_t;
	
	/*
	fData_Pro_Displ = new MFGP_MFA_Data_Processor_DisplT;
	fData_Pro_Displ->Initialize ( Shapes_displ->Derivatives_U(curr_ip), Shapes_displ->DDDerivatives_U(curr_ip) );//??
	fData_Pro_Plast = new MFGP_MFA_Data_Processor_PlastT;
	fData_Pro_Plast->Initialize ( Shapes_displ->IPShapeU(curr_ip), Shapes_displ->DDerivatives_U(curr_ip) );//??
	*/
	Data_Pro_Displ.Initialize(Shapes_displ->Derivatives_U(curr_ip), Shapes_displ->DDDerivatives_U(curr_ip) );//??
	Data_Pro_Plast.Initialize(Shapes_displ->IPShapeU(curr_ip), Shapes_displ->DDerivatives_U(curr_ip) );//??

	stress = GRAD_MR_Plast_Mat->s_ij();
	yield = GRAD_MR_Plast_Mat->YieldF();
	moduli = GRAD_MR_Plast_Mat->c_ijkl();

	Form_C_List(GRAD_MR_Plast_Mat);
	Form_B_List( );
}


//---------------------------------------------------------------------

void MFGP_Con_EqT::Form_LHS_Ku_Klambda(dMatrixT& Ku, dMatrixT& Klambda )  
{
	n_rows = Ku.Rows();
	n_cols = Ku.Cols();
	dMatrixT Kutemp(n_rows,n_cols);
	Ku.MultATBC(phi_lam, Clamu1, B1);
	Kutemp.MultATBC(phi_lam, Clamu2, B3);
	Ku += Kutemp;	// Ku :[nnd]x[nsd*nnd]

	n_rows = Klambda.Rows();
	n_cols = Klambda.Cols();
	dMatrixT Klambdatemp(n_rows,n_cols);
	Klambda.MultATBC(phi_lam, Clamlam1, phi_lam);
	Klambdatemp.MultATBC(phi_lam, Clamlam2, B4);
	Klambda += Klambdatemp;	//Klambda: [nnd]x[nnd]
}

//---------------------------------------------------------------------

void MFGP_Con_EqT::Form_RHS_F_int(dArrayT& F_int) 
{
		//pass column of phi_lam to F_int
		F_int = phi_lam[0];
		F_int *= yield; // F_int: [nnd]x[1]
		//updated failure function passed from the constitutive model
}



//=== Private =========================================================
	             				
void MFGP_Con_EqT::Form_B_List(void)
{
		B1.Dimension(n_str, n_sd_x_n_en_displ);
		B3.Dimension(n_str, n_sd_x_n_en_displ);
		int dum=1;
		phi_lam.Dimension(dum, n_en_plast);
		B4.Dimension(dum, n_en_plast);	
		
		/*
		fData_Pro_Displ->Set_B1(B1);
		fData_Pro_Displ->Set_B3(B3);
		fData_Pro_Plast->Set_psi_lam(phi_lam);
 		fData_Pro_Plast->Set_B4(B4); 
 		*/
 		Data_Pro_Displ.Set_B1(B1);
		Data_Pro_Displ.Set_B3(B3);
		Data_Pro_Plast.Set_psi_lam(phi_lam);
 		Data_Pro_Plast.Set_B4(B4); 
}




void MFGP_Con_EqT::Form_C_List(GRAD_MRSSKStV *GRAD_MR_Plast)
{
		int i,j;
		Cuu1.Dimension(n_str,n_str); 
		Cuu2.Dimension(n_str,n_str);
		Culam1.Dimension(n_str,1); 
		Culam2.Dimension(n_str,1);
		Clamu1.Dimension(1,n_str); 
		Clamu2.Dimension(1,n_str);
		Clamlam1.Dimension(1); 
		Clamlam2.Dimension(1);
		
		/* retrive Cgep/fmoduli: 2D: [4]x[8]; 3D: [7]x[14] 
		*  and form the 8 C matrices
		*  pass the fModuli from the constitutive model */
		for (i=0; i<(n_str+1); ++i)
		{
		 for (j=0; j<(2*n_str+2); ++j)
		 {
		  if (i<n_str & j<n_str)
		  {
		   Cuu1(i,j)=moduli(i,j);
		  }
		  if (i<n_str & j>(n_str-1) & j<2*n_str)
		  {
		   Cuu2(i,j-n_str)=moduli(i,j);
		  }
		  if (i<n_str & j==2*n_str)
		  {
		   Culam1(i,0)=moduli(i,j);
		  }
		  if (i<n_str & j==(2*n_str+1))
		  {
		   Culam2(i,0)=moduli(i,j);
		  }
		  if (i==n_str & j<n_str)
		  {
		   Clamu1(0,j)=moduli(i,j);
		  }
		  if (i==n_str & j>(n_str-1) & j<2*n_str)
		  {
		   Clamu2(0,j-n_str)=moduli(i,j);
		  }
		  if (i==n_str & j==2*n_str)
		  {
		   Clamlam1(0,0)=moduli(i,j);
		  }
		  if (i==n_str & j==(2*n_str+1))
		  {
		   Clamlam2(0,0)=moduli(i,j);
		  }
		 } //inner for loop
		} //outer for loop
		
}
