/* $Id: FDHookeanMat2DT.h,v 1.1 2004-07-22 21:09:37 paklein Exp $ */
#ifndef _FD_HOOKEAN_MAT_2D_H_
#define _FD_HOOKEAN_MAT_2D_H_

/* base class */
#include "FDHookeanMatT.h"

namespace Tahoe {

/** finite strain 2D Hookean material */
class FDHookeanMat2DT: public FDHookeanMatT
{
public:

	/** constructor */
	FDHookeanMat2DT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/* set inverse of thermal transformation - return true if active */
	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);  			
};

} /* namespace Tahoe */

#endif /* _FD_HOOKEAN_MAT_2D_H_ */
