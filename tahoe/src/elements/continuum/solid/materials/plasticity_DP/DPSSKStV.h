/* $Id: DPSSKStV.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: myip (06/01/1999)                                             */

#ifndef _DP_SS_KSTV_H_
#define _DP_SS_KSTV_H_

/* base classes */
#include "SSStructMatT.h"
#include "KStV.h"
#include "HookeanMatT.h"
#include "DPSSLinHardT.h"

class DPSSKStV: public SSStructMatT,
				public KStV,
				public HookeanMatT,
				public DPSSLinHardT
{
public:

	/* constructor */
	DPSSKStV(ifstreamT& in, const ElasticT& element);

	/* form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;

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

private:

	/* return values */
	dSymMatrixT	fStress;
	dMatrixT	fModulus;

	/* elastic modulus */
	dMatrixT    fElasticModulus;
};

#endif /* _DP_SS_KSTV_H_ */
