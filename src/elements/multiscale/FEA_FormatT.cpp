// $Id: FEA_FormatT.cpp,v 1.14 2003-10-02 19:12:05 raregue Exp $
#include "FEA_FormatT.h"

using namespace Tahoe;

//---------------------------------------------------------------------

void FEA_FormatT::Shapes	(ShapeFunctionT *fShapes, FEA_ShapeFunctionT &FEA_Shapes )
{
	FEA_Shapes.j = fShapes->IPDets(); 		// IPDets() returns double*
	FEA_Shapes.W = fShapes->IPWeights(); 	// IPWeights() returns double*
	
	for	(int l=0; l<fShapes->NumIP(); l++) {
		fShapes->SetIP(l);
		fShapes->GradNa		( FEA_Shapes.dNdx[l] 	); 
	}

}

//---------------------------------------------------------------------

void FEA_FormatT::Na	(int n_en, ShapeFunctionT *fShapes, FEA_ShapeFunctionT &FEA_Shapes )
{
	int a,l,n_ip = fShapes -> NumIP();
	FEA_Shapes.N.FEA_Dimension ( n_ip, n_en );

	for	(l=0; l<n_ip; l++) {
		fShapes->SetIP(l);
		const double *fN = fShapes -> IPShapeU();
		for (a=0; a<n_en; a++) 
			FEA_Shapes.N[l][a] = fN[a]; 
	}
}

//---------------------------------------------------------------------

void FEA_FormatT::Interpolate (	ShapeFunctionT *fShapes,LocalArrayT &gammap_np1,LocalArrayT &gammap_n, 
								FEA_dVectorT &Fgammap_np1, FEA_dVectorT &Fgammap_n)
{
	for	(int l=0; l<fShapes->NumIP(); l++) {
		fShapes->SetIP(l);
		fShapes->InterpolateU	( gammap_n, 	Fgammap_n[l], 	l );
		fShapes->InterpolateU 	( gammap_np1, 	Fgammap_np1[l], l );
	}
}

//not used
void FEA_FormatT::State	( int n_ip, int num_state, dArrayT& fState, FEA_dVectorT& state )
{
	int a,i;
	for (a=0; a<n_ip; a++)
		for (i=0; i<num_state; i++)
			state[a,i] = fState[a*num_state+i];
}

void FEA_FormatT::Copy	( int n_ip, int num_state, dArray2DT& fState, FEA_dVectorT& state )
{
	int a,i;
	for (a=0; a<n_ip; a++)
		for (i=0; i<num_state; i++)
			state[a][i] = fState[a*num_state+i];
}

void FEA_FormatT::Copy	( int n_ip, int num_state, FEA_dVectorT& state, dArray2DT& fState  )
{
	int a,i;
	for (a=0; a<n_ip; a++)
		for (i=0; i<num_state; i++)
			fState[a*num_state+i] = state[a][i];
}




//---------------------------------------------------------------------
/*
void FEA_FormatT::Gradients (	ShapeFunctionT *fShapes,LocalArrayT &u_np1,LocalArrayT &u_n, 
								FEA_dVectorT &GRAD_u_np1, FEA_dVectorT &GRAD_u_n)
{
	for	(int l=0; l<fShapes->NumIP(); l++) {
		fShapes->SetIP(l);
#if 0
		fShapes->GradU	( u_n, 		GRAD_u_n[l], 	l );
		fShapes->GradU 	( u_np1, 	GRAD_u_np1[l], l );
#endif
	}
}
*/

//---------------------------------------------------------------------

void FEA_FormatT::Gradients (	ShapeFunctionT *fShapes,LocalArrayT &u_np1,LocalArrayT &u_n, 
								FEA_dMatrixT &GRAD_u_np1, FEA_dMatrixT &GRAD_u_n)
{
	for	(int l=0; l<fShapes->NumIP(); l++) {
		fShapes->SetIP(l);
		fShapes->GradU	( u_n, 		GRAD_u_n[l], 	l );
		fShapes->GradU 	( u_np1, 	GRAD_u_np1[l], l );
	}
}

//---------------------------------------------------------------------

void FEA_FormatT::Displacements ( LocalArrayT &u_mat, dArrayT &u_vec ) 
{
	int a,i,k=0, n_en=u_mat.NumberOfNodes(), n_dof=u_mat.MinorDim();
	for (a=0; a<n_en; a++)
		for (i=0; i<n_dof; i++)
			u_vec[k++] = u_mat(a,i); 
}

