// $Id: APS_FEA_Data_Processor_DisplT.h,v 1.1 2004-02-17 19:48:49 raregue Exp $
#ifndef _APS_FEA_DATAPROCESSOR_DISPLT_H_
#define _APS_FEA_DATAPROCESSOR_DISPLT_H_

#include "FEA.h"

namespace Tahoe {

class APS_FEA_Data_Processor_DisplT  
{
	public:

		enum Spatial_DirectionT { dx1, dx2 };
		enum NumberT { k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 }; // <-- numbers are one less than face value
					
		 APS_FEA_Data_Processor_DisplT 	(); 
		~APS_FEA_Data_Processor_DisplT 	(); 
		APS_FEA_Data_Processor_DisplT 	( FEA_dMatrixT &fdNdx  );
		void Construct 					( FEA_dMatrixT &fdNdx  ); 
		
		void APS_B     				(FEA_dMatrixT &B);

		void Insert_N				(FEA_dVectorT &fN) { N = fN; }

    	nMatrixT<int> Map;		
	  	FEA_dMatrixT	dN;	
	  	FEA_dVectorT	N;	
	  	
		int n_ip, n_en, n_sd;
};

}

// inline routines go here ...

#endif
