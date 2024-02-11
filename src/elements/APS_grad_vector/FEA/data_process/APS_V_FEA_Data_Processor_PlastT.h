// $Id: APS_V_FEA_Data_Processor_PlastT.h,v 1.1 2004-03-02 23:46:52 raregue Exp $
#ifndef _APS_V_FEA_DATAPROCESSOR_PLASTT_H_
#define _APS_V_FEA_DATAPROCESSOR_PLASTT_H_

#include "FEA.h"

namespace Tahoe {

class APS_V_FEA_Data_Processor_PlastT  
{
	public:

		enum Spatial_DirectionT { dx1, dx2 };
		enum NumberT { k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 }; // <-- numbers are one less than face value
					
		 APS_V_FEA_Data_Processor_PlastT 	(); 
		~APS_V_FEA_Data_Processor_PlastT 	(); 
		APS_V_FEA_Data_Processor_PlastT 	( FEA_dMatrixT &fdNdx  );
		void Construct 					( FEA_dMatrixT &fdNdx  ); 
		
		void APS_Ngamma     		(FEA_dMatrixT &B);
		void APS_Ngam1d2			(FEA_dVectorT &B);
		void APS_Ngam2d1			(FEA_dVectorT &B);

		void Insert_N				(FEA_dVectorT &fN) { N = fN; }

    	nMatrixT<int> Map;		
	  	FEA_dMatrixT	dN;	
	  	FEA_dVectorT	N;	
	  	
		int n_ip, n_en, n_sd;
};

}

// inline routines go here ...

#endif
