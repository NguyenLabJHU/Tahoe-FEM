// $Id: MeshfreeGradP_MFA_Data_Processor_PlastT.h
#ifndef _MFGP_MFA_DATAPROCESSOR_PLASTT_H_
#define _MFGP_MFA_DATAPROCESSOR_PLASTT_H_

#include "ArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"


/* pass in shape function and it's laplacian of the plastic multiplier
 * to this class. form B4
*/
namespace Tahoe 
{

class MFGP_MFA_Data_Processor_PlastT  
	{
	public:

		/* contructor */
		MFGP_MFA_Data_Processor_PlastT();

		/* destructor */
		//~MFGP_MFA_Data_Processor_PlastT();
		
		void Initialize(const double *fN, const dArray2DT &fd2Ndx2 );
                        
        /* phi: [1]x[nnd] */
        void Set_phi(dMatrixT& phi);
        
        /* B4: [1]x[nnd] */
		void Set_B4(dMatrixT& B4);
		
	protected:
		
		const double *N;
		dArray2DT d2N;
	};

}


#endif /*_MFGP_MFA_DATAPROCESSOR_PLASTT_H_ */
