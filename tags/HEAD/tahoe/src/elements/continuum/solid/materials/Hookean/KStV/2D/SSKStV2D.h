/* $Id: SSKStV2D.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/10/97)                                            */

#ifndef _SS_KSTV_2D_H_
#define _SS_KSTV_2D_H_

/* base classes */
#include "SSHookeanMatT.h"
#include "KStV2D.h"

class SSKStV2D: public SSHookeanMatT, public KStV2D
{
public:

	/* constructor */
	SSKStV2D(ifstreamT& in, const ElasticT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables */
	virtual int NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

private:

	/* set the internal thermal strain */
	virtual bool SetThermalStrain(dSymMatrixT& thermal_strain);
};

#endif /* _SS_KSTV_2D_H_ */
