// $Id: Coupled_EqsT.cpp,v 1.1 2005-01-26 02:13:22 kyonten Exp $
#include "Coupled_EqsT.h" 

using namespace Tahoe;

/* constructor */
Coupled_EqsT::Coupled_EqsT(void) { }


/* destructor */
//Coupled_EqsT::~Coupled_EqsT(void) { }


void Coupled_EqsT::Initialize(int& curr_ip, D3MeshFreeShapeFunctionT* Shapes_displ, D3MeshFreeShapeFunctionT* Shapes_plast, 
							int& fTime_Step, double fdelta_t) 
{
	n_en_displ  = Shapes_displ->Derivatives_U(curr_ip).MinorDim(); //dN: [nsd]x[nnd]
	n_en_plast  = Shapes_plast->Derivatives_U(curr_ip).MinorDim();
	n_sd    	= Shapes_displ->Derivatives_U(curr_ip).MajorDim();
	n_str = dSymMatrixT::NumValues(n_sd);	
	n_sd_x_n_en_displ = n_sd * n_en_displ;
	n_sd_x_n_en_plast = n_sd * n_en_plast;

	delta_t = fdelta_t;
	
	Data_Pro.Initialize(Shapes_displ->IPShapeU(curr_ip), Shapes_displ->Derivatives_U(curr_ip), 
	                    Shapes_displ->DDerivatives_U(curr_ip), Shapes_displ->DDDerivatives_U(curr_ip) ); //??
 	
 	/* select material model */
	if (n_sd == 2)
	{
		stress = Mat2D.s_ij();
		yield = Mat2D.YieldF();
	}
	else if (n_sd == 3)
	{
		stress = Mat3D.s_ij();
		yield = Mat3D.YieldF();
	}
	else
		cout << "\n Coupled_EqsT::Initialize: nsd type not supported: ";
	
	Form_C_List();
	Form_B_List();
}

//---------------------------------------------------------------------

void Coupled_EqsT::Form_KUU_KULam(dMatrixT& Kuu, dMatrixT& Kulam)  
{		
	dMatrixT Ktemp1(Kuu.Rows(), Kuu.Cols());
	Kuu.MultATBC(B1, Cuu1, B1);
	Ktemp1.MultATBC(B1, Cuu2, B3);
	Kuu += Ktemp1; // Kuu: [nsd*nnd]x[nsd*nnd]
	
	dMatrixT Ktemp2(Kulam.Rows(), Kulam.Cols());
	Kulam.MultATBC(B1, Culam1, psi_lam);
	Ktemp2.MultATBC(B1, Culam2, B4);
	Kulam += Ktemp2;	// Kulambda: [nsd*nnd]x[nnd] 
}

void Coupled_EqsT::Form_FU_int(dArrayT& Fu_int) 
{
	
	dMatrixT B1_T(n_sd_x_n_en_displ, n_str);
	dArrayT stress_sym(n_str);
	B1_T.Transpose(B1);
	
	/* take out the symmetric part of the stress */
	for (int ij = 0; ij < n_str; ij++)
	{
		int i, j;
		dSymMatrixT::ExpandIndex(n_sd, ij, i, j);
		stress_sym[ij] = stress(i,j);
	} 
	
	B1_T.MultTx(stress_sym, Fu_int); // Fu_int: [nsd*nnd]x[1]
}


void Coupled_EqsT::Form_KLamU_KLamLam(dMatrixT& Klamu, dMatrixT& Klamlam)  
{
	dMatrixT Ktemp1(Klamlam.Rows(), Klamlam.Cols());
	Klamlam.MultATBC(psi_lam, Clamlam1, psi_lam);
	Ktemp1.MultATBC(psi_lam, Clamlam2, B4);
	Klamlam += Ktemp1;	//Klamlam: [nnd]x[nnd]
	
	dMatrixT Ktemp2(Klamu.Rows(), Klamu.Cols());
	Klamu.MultATBC(psi_lam, Clamu1, B1);
	Ktemp2.MultATBC(psi_lam, Clamu2, B3);
	Klamu += Ktemp2;	// Klamu :[nnd]x[nsd*nnd]
}

//---------------------------------------------------------------------

void Coupled_EqsT::Form_FLambda_int(dArrayT& Flambda_int) 
{
		//pass column of phi_lam to F_int
		Flambda_int = psi_lam[0];
		Flambda_int *= yield; // Flambda_int: [nnd]x[1]
		//updated failure function passed from the constitutive model
}

void Coupled_EqsT::Form_B_List(void)
{
		/* dimension */
		B1.Dimension(n_str, n_sd_x_n_en_displ);
		B3.Dimension(n_str, n_sd_x_n_en_displ);
		int dum=1;
		psi_lam.Dimension(dum, n_en_plast);
		B4.Dimension(dum, n_en_plast);
		
 		Data_Pro.Set_B1(B1);
		Data_Pro.Set_B3(B3);
		Data_Pro.Set_psi_lam(psi_lam);
 		Data_Pro.Set_B4(B4); 
}


void Coupled_EqsT::Form_C_List()
{
		/* dimensions */
		Cuu1.Dimension(n_str,n_str); 
		Cuu2.Dimension(n_str,n_str);
		Culam1.Dimension(n_str,1); 
		Culam2.Dimension(n_str,1);
		Clamu1.Dimension(1,n_str); 
		Clamu2.Dimension(1,n_str);
		Clamlam1.Dimension(1,1); 
		Clamlam2.Dimension(1,1);
		
		/* select material model */
		if (n_sd == 2)
		{
			Cuu1 = Mat2D.c_UU1_ijkl(); 
			Cuu2 = Mat2D.c_UU2_ijkl();
			Culam1 = Mat2D.c_ULam1_ij(); 
			Culam2 = Mat2D.c_ULam2_ij();
			Clamu1 = Mat2D.c_LamU1_ij(); 
			Clamu2 = Mat2D.c_LamU2_ij();
			Clamlam1 = Mat2D.c_LamLam1(); 
			Clamlam2 = Mat2D.c_LamLam2();
		}
		else if (n_sd == 3)
		{
			Cuu1 = Mat3D.c_UU1_ijkl(); 
			Cuu2 = Mat3D.c_UU2_ijkl();
			Culam1 = Mat3D.c_ULam1_ij(); 
			Culam2 = Mat3D.c_ULam2_ij();
			Clamu1 = Mat3D.c_LamU1_ij(); 
			Clamu2 = Mat3D.c_LamU2_ij();
			Clamlam1 = Mat3D.c_LamLam1(); 
			Clamlam2 = Mat3D.c_LamLam2();
		}
		else
			cout << "\n Coupled_EqsT::Initialize: nsd type not supported: ";
}
