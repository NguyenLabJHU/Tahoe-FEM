
#ifndef _FEA_FORMATT_
#define _FEA_FORMATT_

#include "FEA.h"
#include "ShapeFunctionT.h"
#include "dArrayT.h"
#include "LocalArrayT.h"

namespace Tahoe {

class	FEA_FormatT {

	public:

		FEA_FormatT (void) { }

		void Shapes				(	ShapeFunctionT *fShapes, FEA_ShapeFunctionT &FEA_Shapes 				);
		void Shapes				(	ShapeFunctionT &fShapes, FEA_ShapeFunctionT &FEA_Shapes 				);
		void SurfShapeGradient	(int n_en, const ParentDomainT& surf_shapes, 
								FEA_SurfShapeFunctionT &FEA_SurfShapes, 
								LocalArrayT& face_coords, const ParentDomainT& parent,
								LocalArrayT& volume_coords, ShapeFunctionT& shapes,
								LocalArrayT &u_np1,LocalArrayT &u_n, 
								FEA_dMatrixT &GRAD_u_np1, FEA_dMatrixT &GRAD_u_n,
								LocalArrayT& face_gamma_p, FEA_dVectorT& fgamma_p_surf, iArrayT& face_nodes   );
		void Na					(	int n_en, ShapeFunctionT *fShapes, FEA_ShapeFunctionT &FEA_Shapes 		);
		void Na					(	int n_en, ShapeFunctionT &fShapes, FEA_ShapeFunctionT &FEA_Shapes 		);
		void Gradients 			(	ShapeFunctionT*,LocalArrayT&,LocalArrayT&,FEA_dMatrixT&,FEA_dMatrixT&	);
		void Gradients 			(	ShapeFunctionT&,LocalArrayT&,LocalArrayT&,FEA_dMatrixT&,FEA_dMatrixT&	);
		//void Gradients 		(	ShapeFunctionT*,LocalArrayT&,LocalArrayT&,FEA_dVectorT&,FEA_dVectorT&	);
		void Displacements 		(	LocalArrayT &u_mat, dArrayT &u_vec 										);
		void Interpolate 		(	ShapeFunctionT*,LocalArrayT&,LocalArrayT&,FEA_dVectorT&,FEA_dVectorT&  	);
		void Interpolate 		(	ShapeFunctionT&,LocalArrayT&,LocalArrayT&,FEA_dVectorT&,FEA_dVectorT&  	);
		void State		 		(	int n_ip, int num_state, dArrayT&, FEA_dVectorT&	);
		void Copy		 		(	int n_ip, int num_state, dArray2DT&, FEA_dVectorT&	);
		void Copy		 		(	int n_ip, int num_state, FEA_dVectorT&, dArray2DT& 	);

};

}

#endif



