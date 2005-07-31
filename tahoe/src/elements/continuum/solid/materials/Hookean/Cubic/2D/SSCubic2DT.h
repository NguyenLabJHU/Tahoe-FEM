/* $Id: SSCubic2DT.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/11/97)                                            */

#ifndef _SS_CUBIC_2D_T_H_
#define _SS_CUBIC_2D_T_H_

/* base classes */
#include "SSHookeanMatT.h"
#include "Cubic2DT.h"

class SSCubic2DT: public SSHookeanMatT, public Cubic2DT
{
public:

	/* constructor */
	SSCubic2DT(ifstreamT& in, const ElasticT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

private:

	/* set the internal thermal strain */
	virtual bool SetThermalStrain(dSymMatrixT& thermal_strain);
};

#endif /* _SS_CUBIC_2D_T_H_ */
