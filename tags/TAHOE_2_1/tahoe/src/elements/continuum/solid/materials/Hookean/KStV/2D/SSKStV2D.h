/* $Id: SSKStV2D.h,v 1.5 2004-07-15 08:27:18 paklein Exp $ */
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
	SSKStV2D(void);
	
	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables */
	virtual int NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/*@}*/

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

private:

	/* set the internal thermal strain */
	virtual bool SetThermalStrain(dSymMatrixT& thermal_strain);
};

} // namespace Tahoe 
#endif /* _SS_KSTV_2D_H_ */
