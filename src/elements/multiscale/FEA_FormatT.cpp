// $Id: FEA_FormatT.cpp,v 1.25 2003-10-28 01:52:53 raregue Exp $
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

void FEA_FormatT::Shapes	(ShapeFunctionT &fShapes, FEA_ShapeFunctionT &FEA_Shapes )
{
	FEA_Shapes.j = fShapes.IPDets(); 		// IPDets() returns double*
	FEA_Shapes.W = fShapes.IPWeights(); 	// IPWeights() returns double*
	
	for	(int l=0; l<fShapes.NumIP(); l++) {
		fShapes.SetIP(l);
		fShapes.GradNa		( FEA_Shapes.dNdx[l] 	); 
	}

}


//---------------------------------------------------------------------

void FEA_FormatT::SurfShapeGradient	(int n_en, const ParentDomainT& surf_shapes, 
								FEA_SurfShapeFunctionT &FEA_SurfShapes, 
								LocalArrayT& face_coords, const ParentDomainT& parent,
								LocalArrayT& volume_coords,
								ShapeFunctionT& shapes,
								LocalArrayT &u_np1,LocalArrayT &u_n, 
								FEA_dMatrixT &GRAD_u_np1, FEA_dMatrixT &GRAD_u_n,
								LocalArrayT& face_gamma_p, FEA_dVectorT& fgamma_p_surf, iArrayT& face_nodes )
{
	dMatrixT face_jacobian(2, 1);
	dMatrixT face_Q(2);
	dArrayT fNormal(2);
	dArrayT interp_ip(2), interp_ip_mapped(2);
	bool ismapped;
	dArrayT Na;
	dArray2DT DNa;
	
	int nip_surf = surf_shapes.NumIP();
	
	FEA_SurfShapes.W = surf_shapes.Weight(); 	// IPWeights() returns double*
	
	for	(int l=0; l<nip_surf; l++) 
	{
		surf_shapes.DomainJacobian(face_coords, l, face_jacobian);
		FEA_SurfShapes.j[l] = surf_shapes.SurfaceJacobian(face_jacobian, face_Q); 	// IPDets() returns double*
			
		/* last column is the normal (I think) */
		face_Q.ColumnAlias(face_Q.Cols()-1, fNormal);
		FEA_SurfShapes.normal[l] = fNormal;
		
		surf_shapes.Interpolate(face_coords, interp_ip, l);
		ismapped = parent.MapToParentDomain(volume_coords, interp_ip, interp_ip_mapped);
		shapes.GradU	( u_n, 		GRAD_u_n[l], interp_ip_mapped, Na, DNa);
		shapes.GradU 	( u_np1, 	GRAD_u_np1[l], interp_ip_mapped, Na, DNa );
		
		//use face_nodes?
		for (int a=0; a<n_en; a++) 
		{
			dMatrixT& dN_dx = FEA_SurfShapes.dNdx[l];
			dN_dx(0,a) = DNa(0,face_nodes[a]);
			dN_dx(1,a) = DNa(1,face_nodes[a]);

			dArrayT& N = FEA_SurfShapes.N[l];
			N[a] = Na[face_nodes[a]];
		}	
		
		surf_shapes.Interpolate(face_gamma_p, fgamma_p_surf[l], l);
		
		/*
		const double* DNa = surf_shapes.DShape(l,0);
		//const double* DNa_x = surf_shapes.DShape(l,0);
		//const double* DNa_y = surf_shapes.DShape(l,1);
		const double* Na = surf_shapes.Shape(l); 
		for (int a=0; a<n_en; a++) 
		{
			dMatrixT& dN_dx = FEA_SurfShapes.dNdx[l];
			dN_dx(0,a) = DNa[a];
			//dN_dx(0,a) = DNa_x[a];
			//dN_dx(1,a) = DNa_y[a];

			dArrayT& N = FEA_SurfShapes.N[l];
			N[a] = Na[a];
		}
		*/
		
//		FEA_SurfShapes.N[l] = surf_shapes.Shape(l); 
//		FEA_SurfShapes.dNdx[l].SetRow(0, surf_shapes.DShape(l,0));
//		FEA_SurfShapes.dNdx[l].SetRow(1, surf_shapes.DShape(l,1));
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

void FEA_FormatT::Na	(int n_en, ShapeFunctionT &fShapes, FEA_ShapeFunctionT &FEA_Shapes )
{
	int a,l,n_ip = fShapes.NumIP();
	FEA_Shapes.N.FEA_Dimension ( n_ip, n_en );

	for	(l=0; l<n_ip; l++) {
		fShapes.SetIP(l);
		const double *fN = fShapes.IPShapeU();
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

//---------------------------------------------------------------------

void FEA_FormatT::Interpolate (	ShapeFunctionT &fShapes,LocalArrayT &gammap_np1,LocalArrayT &gammap_n, 
								FEA_dVectorT &Fgammap_np1, FEA_dVectorT &Fgammap_n)
{
	for	(int l=0; l<fShapes.NumIP(); l++) {
		fShapes.SetIP(l);
		fShapes.InterpolateU	( gammap_n, 	Fgammap_n[l], 	l );
		fShapes.InterpolateU 	( gammap_np1, 	Fgammap_np1[l], l );
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

void FEA_FormatT::Gradients (	ShapeFunctionT &fShapes,LocalArrayT &u_np1,LocalArrayT &u_n, 
								FEA_dMatrixT &GRAD_u_np1, FEA_dMatrixT &GRAD_u_n)
{
	for	(int l=0; l<fShapes.NumIP(); l++) {
		fShapes.SetIP(l);
		fShapes.GradU	( u_n, 		GRAD_u_n[l], 	l );
		fShapes.GradU 	( u_np1, 	GRAD_u_np1[l], l );
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

