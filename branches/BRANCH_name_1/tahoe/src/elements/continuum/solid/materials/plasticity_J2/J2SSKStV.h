/* $Id: J2SSKStV.h,v 1.3.4.1 2002-06-27 18:03:35 cjkimme Exp $ */
/* created: paklein (06/18/1997) */

#ifndef _J2_SS_KSTV_H_
#define _J2_SS_KSTV_H_

/* base classes */
#include "SSStructMatT.h"
#include "IsotropicT.h"
#include "HookeanMatT.h"
//#include "J2SSLinHardT.h"
#include "J2SSC0HardeningT.h"


namespace Tahoe {

class J2SSKStV: public SSStructMatT,
				public IsotropicT,
				public HookeanMatT,
//				public J2SSLinHardT
				public J2SSC0HardeningT
{
public:

	/* constructor */
	J2SSKStV(ifstreamT& in, const SmallStrainT& element);

	/* initialization */
	virtual void Initialize(void);

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

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
	 	 	
private:

	/* return values */
	dSymMatrixT	fStress;
	dMatrixT	fModulus;
};

} // namespace Tahoe 
#endif /* _J2_SS_KSTV_H_ */
