/* $Id: ThreeBodyT.h,v 1.1.1.1 2001-01-29 08:20:26 paklein Exp $ */
/* created: paklein (10/11/1997)                                          */
/* Base class for the 3 body contribution to the strain energy density    */

#ifndef _THREE_BODY_T_H_
#define _THREE_BODY_T_H_

/* direct members */
#include "dArrayT.h"
#include "dArray2DT.h"

/* forward declaration */
class iArray2DT;
class dMatrixT;
class ThermalDilatationT;

class ThreeBodyT
{
public:

	/* constructor */
	ThreeBodyT(const dArrayT& lengths, const dArrayT& angles,
		const iArray2DT& bondpairs, const ThermalDilatationT* thermal);

	/* triggers recomputation */
	virtual void Set(void) = 0;

	/* Accessors */
	const dArrayT& Phi(void) const;
	const dArray2DT& dPhi(void) const;
	const dArray2DT& ddPhi(void) const;

protected:

	/* 3-body bond pairs */
	const dArrayT&		fLengths;
	const dArrayT&		fAngles;
	const iArray2DT& 	fPairs;
	
	/* potential function values and derivatives */
	dArrayT		fPhi;
	dArray2DT	fdPhi;
	dArray2DT	fddPhi;
	
	/* thermal dilatation LTf */
	const ThermalDilatationT* fThermal;
	
};

#endif /* _THREE_BODY_T_H_ */
