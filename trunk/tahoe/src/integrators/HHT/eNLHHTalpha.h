/* $Id: eNLHHTalpha.h,v 1.3 2002-07-02 19:55:05 cjkimme Exp $ */
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
	eNLHHTalpha(ifstreamT& in, ostream& out, bool auto2ndorder);

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
