/* $Id: SSKStV2D.h,v 1.4.32.2 2004-03-02 17:46:15 paklein Exp $ */
/* created: paklein (06/10/97) */
#ifndef _SS_KSTV_2D_H_
#define _SS_KSTV_2D_H_

/* base classes */
#include "SSKStV.h"

namespace Tahoe {

class SSKStV2D: public SSKStV
{
public:

	/** constructor */
	SSKStV2D(ifstreamT& in, const SSMatSupportT& support);
	SSKStV2D(void);
	
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
