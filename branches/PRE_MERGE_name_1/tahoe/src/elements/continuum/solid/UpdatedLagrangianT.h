/* $Id: UpdatedLagrangianT.h,v 1.4 2002-06-08 20:20:22 paklein Exp $ */
/* created: paklein (07/03/1996) */

#ifndef _UPDATED_LAGRANGIAN_T_H_
#define _UPDATED_LAGRANGIAN_T_H_

/* base class */
#include "FiniteStrainT.h"

/* direct members */
#include "dMatrixT.h"

/** update Lagrangian, finite strain solid */
class UpdatedLagrangianT: public FiniteStrainT
{
public:

	/* constructors */
	UpdatedLagrangianT(const ElementSupportT& support, const FieldT& field);

	/* destructors */
	virtual ~UpdatedLagrangianT(void);

	/* data initialization */
	virtual void Initialize(void);
		
protected:

	/* initialization functions */
	virtual void SetLocalArrays(void);
	virtual void SetShape(void);

	/* form shape functions and derivatives */
	virtual void SetGlobalShape(void);

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
	 *           function of the displacements */
	virtual void ComputeEffectiveDVA(int formBody,
		int formMa, double constMa, int formCv, double constCv,
		int formKd, double constKd);
#endif

	/* calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

protected:

	/* shape functions (wrt current config) */
	ShapeFunctionT* fCurrShapes;

	dMatrixT fCauchyStress;	// symmetric Cauchy stress tensor
	dMatrixT fStressStiff;	// compact stress stiffness contribution
	dMatrixT fGradNa;       // shape function gradients matrix
	
	/* arrays with local ordering */
	LocalArrayT fLocCurrCoords;	// current coords with local ordering
		
	/* workspace */
	dArrayT fTemp2;
};

#endif /* _UPDATED_LAGRANGIAN_T_H_ */
