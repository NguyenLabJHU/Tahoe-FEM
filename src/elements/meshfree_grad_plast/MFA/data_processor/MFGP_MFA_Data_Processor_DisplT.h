// $Id: MFGP_MFA_Data_Processor_DisplT.h
#ifndef _MFGP_MFA_DATAPROCESSOR_DISPLT_H_
#define _MFGP_MFA_DATAPROCESSOR_DISPLT_H_

#include "ArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"

/* pass first and third derivatives of the shape function 
*  of displacement and form B1 and B3 matrices.
*/ 
namespace Tahoe 
{

class MFGP_MFA_Data_Processor_DisplT  
	{
	public:

		/* constructor */
		MFGP_MFA_Data_Processor_DisplT();
		
		/* destructor */
		//~MFGP_MFA_Data_Processor_DisplT();
		
		void Initialize(const dArray2DT& fdNdx, const dArray2DT& fd3Ndx3 );
        
		/* B1: [nstr]x[nsd*nnd] */
		void Set_B1(dMatrixT& B1);
		
		/* B3: [nstr]x[nsd*nnd] */ 
		void Set_B3(dMatrixT& B3);
	
	protected:
		
		dArray2DT dN, d3N;
	};

}

// inline routines go here ...

#endif /* _MFGP_MFA_DATAPROCESSOR_DISPLT_H_ */
