
//DEVELOPMENT

#ifndef _MFGP_BALLINMOM_T_H_ 
#define _MFGP_BALLINMOM_T_H_ 

#include "StringT.h"
#include "GRAD_MRSSKStV.h"
#include "D3MeshFreeShapeFunctionT.h"
#include "D3MeshFreeSupportT.h"

namespace Tahoe 
{

class MFGP_BalLinMomT
{
public:

	enum Eqn_TypeT { kMFGP_Bal_Eq };

	/* constructor */
	MFGP_BalLinMomT(void);
	
	/* destructor */
	virtual ~MFGP_BalLinMomT(void);

	/** Pure virtual functions */
	virtual void Initialize (int&, D3MeshFreeShapeFunctionT*, D3MeshFreeShapeFunctionT*, GRAD_MRSSKStV*, 
							int	&fTime_Step, double fdelta_t = 0.0) =0;
	virtual void Form_LHS_Klambda_Ku ( dMatrixT &Klambda, dMatrixT &Ku ) =0; 
  	virtual void Form_RHS_F_int	( dArrayT  &F_int  ) =0; 
	
};
} // namespace Tahoe 
#endif /* _MFGP_BALLINMOM_T_H_ */

