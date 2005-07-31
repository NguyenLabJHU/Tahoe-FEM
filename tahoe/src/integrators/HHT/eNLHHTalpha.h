/* $Id: eNLHHTalpha.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/17/1996)                                          */

#ifndef _E_NL_HHT_A_H_
#define _E_NL_HHT_A_H_

/* base class */
#include "eLinearHHTalpha.h"

class eNLHHTalpha: public eLinearHHTalpha
{
public:

	/* constructor */
	eNLHHTalpha(ifstreamT& in, ostream& out, int auto2ndorder = kHHTalphaAuto_O2);

	/* components of the internal force vector */
	virtual int FormMa(double& constMa) const;
	virtual int FormCv(double& constCv) const;
	virtual int FormKd(double& constKd) const;

protected:  	
	
	/* recalculate constants */
	virtual void eComputeParameters(void);

private:
	
	/* element residual force coefficients */
	double 	fconstMa;
	double	fconstCv;
	double	fconstKd;
	
};

#endif /* _E_NL_HHT_A_H_ */
