// $Id: MFGP_Data_ProcessorT.h
#ifndef _MFGP_DATAPROCESSORT_H_
#define _MFGP_DATAPROCESSORT_H_

#include "ArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"

/**
 * pass first, second, and third derivatives of the shape functions, 
 * and form B1, B3, psi_lam, and B4.  The first two correspond to 
 * displacement, and the last two to plastic multiplier.  
 * The shape functions for the displacement and plastic multiplier 
 * are the same.
*/ 
namespace Tahoe 
{

class MFGP_Data_ProcessorT  
{
public:

	/* constructor */
	MFGP_Data_ProcessorT();
		
	/* destructor */
	//~MFGP_Data_Processor_DisplT();
		
	void Initialize(const double* fN, const dArray2DT& fDN,
						const dArray2DT& fDDN, const dArray2DT& fDDDN);
        
	/* B1: [nstr]x[nsd*nnd]  */
	void Set_B1(dMatrixT& B1);
		
	/* B3: [nstr]x[nsd*nnd] */ 
	void Set_B3(dMatrixT& B3);
                        
    /* psi_lam: [1]x[nnd] */
    void Set_psi_lam(dMatrixT& psi_lam);
        
    /* B4: [1]x[nnd] */
	void Set_B4(dMatrixT& B4);
	
protected:
		
	const double* N;
	dArray2DT dN, d2N, d3N;
};

}

// inline routines go here ...

#endif /* _MFGP_DATAPROCESSORT_H_ */
