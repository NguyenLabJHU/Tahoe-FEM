// $Id: FEA_FormatT.cpp,v 1.4 2003-02-03 04:40:20 paklein Exp $
#include "FEA_FormatT.h"

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

void FEA_FormatT::Gradiants (	ShapeFunctionT *fShapes,LocalArrayT &u_np1,LocalArrayT &u_n, 
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
	int i,j,k=0, n_en=u_mat.NumberOfNodes(), n_dof=u_mat.MinorDim();
	for (i=0; i<n_en; i++)
		for (j=0; j<n_dof; j++)
			u_vec[k++] = u_mat(i,j); 
}

