/* $Id: J2SSKStV.h,v 1.8.34.1 2005-04-05 23:29:02 thao Exp $ */
/* created: paklein (06/18/1997) */
#ifndef _J2_SS_KSTV_H_
#define _J2_SS_KSTV_H_

/* base classes */
#include "SSSolidMatT.h"
#include "IsotropicT.h"
#include "HookeanMatT.h"
#include "J2SSC0HardeningT.h"

namespace Tahoe {

class J2SSKStV: public SSSolidMatT,
				public IsotropicT,
				public HookeanMatT,
				public J2SSC0HardeningT
{
public:

	/* constructor */
	J2SSKStV(ifstreamT& in, const SSMatSupportT& support);

	/* initialization */
	virtual void Initialize(void);

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fStress.Trace()/3.0; };
	/*@}*/

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables */
	virtual int NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

	/*internal variable accessors*/
	virtual const iArrayT& InternalDOF(void) const;
	virtual const dArrayT& InternalStrainVars(void);
	virtual const dArrayT& InternalStressVars(void);
	 	 	
private:

	/* return values */
	dSymMatrixT	fStress;
	dMatrixT	fModulus;
};

} // namespace Tahoe 
#endif /* _J2_SS_KSTV_H_ */
