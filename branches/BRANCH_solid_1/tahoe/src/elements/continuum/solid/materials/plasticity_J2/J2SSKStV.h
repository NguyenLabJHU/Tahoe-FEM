/* $Id: J2SSKStV.h,v 1.1.1.1.2.1 2001-06-06 16:27:56 paklein Exp $ */
/* created: paklein (06/18/1997)                                          */

#ifndef _J2_SS_KSTV_H_
#define _J2_SS_KSTV_H_

/* base classes */
#include "SSStructMatT.h"
#include "IsotropicT.h"
#include "HookeanMatT.h"
#include "J2SSLinHardT.h"

class J2SSKStV: public SSStructMatT,
				public IsotropicT,
				public HookeanMatT,
				public J2SSLinHardT
{
public:

	/* constructor */
	J2SSKStV(ifstreamT& in, const ElasticT& element);

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

#endif /* _J2_SS_KSTV_H_ */
