/* $Id: TwoBodyT.h,v 1.1.1.1.10.1 2002-06-27 18:03:30 cjkimme Exp $ */
/* created: paklein (10/11/1997)                                          */
/* Base class for the 2 body contribution to the strain energy density    */

#ifndef _TWO_BODY_T_H_
#define _TWO_BODY_T_H_

/* direct members */
#include "dArrayT.h"

/* forward declarations */

namespace Tahoe {

class ThermalDilatationT;

class TwoBodyT
{
public:

	/* constructor */
	TwoBodyT(const dArrayT& lengths, const ThermalDilatationT* thermal);

	/* set free dof - triggers recomputation */
	virtual void Set(void) = 0;

	/* Accessors */
	const dArrayT& Phi(void) const;
	const dArrayT& dPhi(void) const;
	const dArrayT& ddPhi(void) const;

protected:

	/* geometry */
	const dArrayT& fLengths;

	/* potential function values and derivatives */
	dArrayT	fPhi;
	dArrayT	fdPhi;
	dArrayT	fddPhi;
	
	/* thermal dilatation LTf */
	const ThermalDilatationT* fThermal;
	
};

} // namespace Tahoe 
#endif /* _TWO_BODY_T_H_ */