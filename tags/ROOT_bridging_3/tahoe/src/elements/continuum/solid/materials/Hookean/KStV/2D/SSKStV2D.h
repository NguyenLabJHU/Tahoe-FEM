/* $Id: SSKStV2D.h,v 1.4 2002-11-14 17:06:06 paklein Exp $ */
/* created: paklein (06/10/97) */
#ifndef _SS_KSTV_2D_H_
#define _SS_KSTV_2D_H_

/* base classes */
#include "SSKStV.h"
#include "Material2DT.h"

namespace Tahoe {

class SSKStV2D: public SSKStV, public Material2DT
{
public:

	/* constructor */
	SSKStV2D(ifstreamT& in, const SSMatSupportT& support);

	/* print parameters */
	virtual void Print(ostream& out) const;

	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables */
	virtual int NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

private:

	/* set the internal thermal strain */
	virtual bool SetThermalStrain(dSymMatrixT& thermal_strain);
};

} // namespace Tahoe 
#endif /* _SS_KSTV_2D_H_ */
