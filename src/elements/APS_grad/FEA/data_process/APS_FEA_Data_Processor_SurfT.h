// $Id: APS_FEA_Data_Processor_SurfT.h,v 1.1 2003-10-07 19:40:43 raregue Exp $
#ifndef _APS_FEA_DATAPROCESSORSURFT_H_
#define _APS_FEA_DATAPROCESSORSURFT_H_

#include "FEA.h"

namespace Tahoe {

class APS_FEA_Data_Processor_SurfT  
{
	public:

		enum Spatial_DirectionT { dx1, dx2 };
		enum NumberT { k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 }; // <-- numbers are one less than face value
					
		 APS_FEA_Data_Processor_SurfT 	(); 
		~APS_FEA_Data_Processor_SurfT 	(); 
		APS_FEA_Data_Processor_SurfT 	( int &n_en, iArrayT &face_nodes, FEA_dMatrixT &fdNdx_surf  );
		void Construct 					( int &n_en, iArrayT &face_nodes, FEA_dMatrixT &fdNdx_surf  ); 
		
		void APS_B_surf     		(FEA_dMatrixT &B);
		void APS_N					(FEA_dVectorT &B);
		
		void Insert_N_surf			(FEA_dVectorT &fN_surf) { N_surf = fN_surf; }

    	nMatrixT<int> Map;		
	  	FEA_dMatrixT	dN_surf;	
	  	FEA_dVectorT	N_surf;	
	  	
	  	iArrayT face_nodes;
	  	
		int n_ip, n_en, n_sd;
};

}

// inline routines go here ...

#endif
