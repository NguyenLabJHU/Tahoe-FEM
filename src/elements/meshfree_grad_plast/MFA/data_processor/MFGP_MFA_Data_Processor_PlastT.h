// $Id: MeshfreeGradP_MFA_Data_Processor_PlastT.h
#ifndef _MFGP_MFA_DATAPROCESSOR_PLASTT_H_
#define _MFGP_MFA_DATAPROCESSOR_PLASTT_H_

#include "MeshFreeT.h"
#include "MLSSolverT.h"

// Pass in nodal info, Laplacian of the plastic multiplier
// shape function to this class. Form B4
namespace Tahoe 
{

class MFGP_MFA_Data_Processor_PlastT  
	{
	public:

		/* contructor */
		MFGP_MFA_Data_Processor_PlastT(void);

		/* destructor */
		~MFGP_MFA_Data_Processor_PlastT(void); 
		
		void Initialize ( const double *fN, const dArray2DT &fd2Ndx2 );
                        
        void Set_phi (dMatrixT &phi);
		void Set_B4 (dMatrixT &B4);
		
	protected:
		
		const double *N;
		dArray2DT d2N;
	};

}


#endif /*_MFGP_MFA_DATAPROCESSOR_PLASTT_H_ */
