/* $Id: TotalLagrangianT.h,v 1.6 2002-07-05 22:28:03 paklein Exp $ */
/* created: paklein (07/03/1996) */

#ifndef _TOTAL_LAGRANGRIAN_T_H_
#define _TOTAL_LAGRANGRIAN_T_H_

/* base class */
#include "FiniteStrainT.h"

/* direct members */
#include "dMatrixT.h"


namespace Tahoe {

class TotalLagrangianT: public FiniteStrainT
{
public:

	/* constructors */
	TotalLagrangianT(const ElementSupportT& support, const FieldT& field);

	/* data initialization */
	virtual void Initialize(void);
		
protected:

	/* form the element stiffness matrix */
	virtual void FormStiffness(double constK);
	
//DEV - Rayleigh damping should be added to the constitutive level
#if 0
	/* compute the effective acceleration and velocities based
	 * on the algorithmic flags formXx and the given constants
	 * constXx.
	 *
	 *		acc_eff  = constMa acc  + constCv a vel
	 *      disp_eff = constKd disp
	 *
	 * where a and b are the Rayleigh damping coefficients.
	 *
	 *        ***The effective displacement does not include
	 *           velocity since the internal force is a nonlinear
	 *           function of the displacements
*/
	virtual void ComputeEffectiveDVA(int formBody,
		int formMa, double constMa, int formCv, double constCv,
		int formKd, double constKd);
#endif

	/* calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

protected:

	/* workspace */
	dMatrixT fStressMat;   // space for a stress tensor
	dMatrixT fStressStiff; // compact stress stiffness contribution
	dMatrixT fGradNa;      // shape function gradients matrix
	
	dArrayT   fTemp2;
	dMatrixT  fTempMat1, fTempMat2;
	dArray2DT fDNa_x;
};

} // namespace Tahoe 
#endif /* _TOTAL_LAGRANGRIAN_T_H_ */
