/* $Id: SSHookeanMat2DT.h,v 1.2 2004-09-10 22:38:57 paklein Exp $ */
#ifndef _SS_HOOKEAN_MAT_2D_H_
#define _SS_HOOKEAN_MAT_2D_H_

/* base classes */
#include "SSHookeanMatT.h"

namespace Tahoe {

/** finite strain 2D Hookean material */
class SSHookeanMat2DT: public SSHookeanMatT
{
public:

	/** constructor */
	SSHookeanMat2DT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/* set the internal thermal strain */
	virtual bool SetThermalStrain(dSymMatrixT& thermal_strain);
};

} /* namespace Tahoe */

#endif /* _SS_HOOKEAN_MAT_2D_H_ */
