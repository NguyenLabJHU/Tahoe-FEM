// $Id: MeshfreeGradP_MFA_Data_Processor_PlastT.h
#ifndef _MFGP_MFA_DATAPROCESSOR_PLASTT_H_
#define _MFGP_MFA_DATAPROCESSOR_PLASTT_H_

#include "MeshFreeT.h"
#include "MLSSolverGPT.h"

// Pass in nodal info, Laplacian of the plastic multiplier
// shape function to this class. Form B4
namespace Tahoe 
{

class MFGP_MFA_Data_Processor_PlastT  
	{
	public:

		 MFGP_MFA_Data_Processor_PlastT(); 
		~MFGP_MFA_Data_Processor_PlastT(); 
		
		void MFGP_MFA_Data_Processor_PlastT(MLSSolverGPT::SetShapeFunctions
                        (const dArrayT& volume))
		// shape function may be directly passed from
		// MLSSolverGP class??
		void MFGP_MFA_Data_Processor_PlastT(MLSSolverGPT::SetShapeFunctions
                        (const dArrayT& volume))
        void Set_N (double phi);
		void Set_B4 (double B3);
	};

}


#endif /*_MFGP_MFA_DATAPROCESSOR_PLASTT_H_ */
