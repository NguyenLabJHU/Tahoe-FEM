// $Id: MFGP_PlastT.h
#ifndef _MFGP_PLAST_T_H_ 
#define _MFGP_PLAST_T_H_ 

#include "StringT.h"
#include "GRAD_MRSSKStV.h"
#include "D3MeshFreeShapeFunctionT.h"
#include "D3MeshFreeSupportT.h"

namespace Tahoe {

/* forward declaration */
class D3MeshFreeShapeFunctionT; 

class MFGP_PlastT
{
public:

	enum Eqn_TypeT 	{ kMFGP_Con_Eq};

	/* constructor */
	MFGP_PlastT(void);
	
	/* destructor */
	virtual ~MFGP_PlastT(void);

	/** Pure virtual functions */
	virtual void Initialize (  int &curr_ip, D3MeshFreeShapeFunctionT*, D3MeshFreeShapeFunctionT*, GRAD_MRSSKStV*,  
						int	&fTime_Step, double fdelta_t = 0.0 ) = 0; 
	virtual void Form_LHS_Ku_Klambda ( dMatrixT &Ku, dMatrixT &Klambda )	= 0; 
	virtual void Form_RHS_F_int	( dArrayT &F_int ) = 0;

};

} // namespace Tahoe 
#endif /* _MFGP_PLAST_T_H_ */

