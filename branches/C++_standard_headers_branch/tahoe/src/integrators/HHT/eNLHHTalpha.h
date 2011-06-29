/* $Id: eNLHHTalpha.h,v 1.4 2004-07-15 08:30:28 paklein Exp $ */
/* created: paklein (10/17/1996) */
#ifndef _E_NL_HHT_A_H_
#define _E_NL_HHT_A_H_

/* base class */
#include "eLinearHHTalpha.h"

namespace Tahoe {

class eNLHHTalpha: public eLinearHHTalpha
{
public:

	/** constructor */
	eNLHHTalpha(double alpha);

 	/** \name elements of the residual
	 * components of the internal force vector */
	/*@{*/
	virtual int FormMa(double& constMa) const;
	virtual int FormCv(double& constCv) const;
	virtual int FormKd(double& constKd) const;
	/*@}*/

protected:  	
	
	/** recalculate constants */
	virtual void eComputeParameters(void);

private:

	/** \name element residual force coefficients */
	/*@{*/	
	double 	fconstMa;
	double	fconstCv;
	double	fconstKd;
	/*@}*/	
	
};

} // namespace Tahoe 
#endif /* _E_NL_HHT_A_H_ */
