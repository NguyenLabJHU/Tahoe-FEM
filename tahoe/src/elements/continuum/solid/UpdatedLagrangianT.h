/* $Id: UpdatedLagrangianT.h,v 1.2 2001-07-03 01:34:53 paklein Exp $ */
/* created: paklein (07/03/1996)                                          */

#ifndef _UPDATED_LAGRANGIAN_T_H_
#define _UPDATED_LAGRANGIAN_T_H_

/* base class */
#include "FiniteStrainT.h"

/* direct members */
#include "dMatrixT.h"

class UpdatedLagrangianT: public FiniteStrainT
{
public:

	/* constructors */
	UpdatedLagrangianT(FEManagerT& fe_manager);

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
	
	/* return a pointer to a new material list */
	virtual MaterialListT* NewMaterialList(int size) const;

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

	/* calculate the damping force contribution ("-c*v") */
	virtual void FormCv(double constC);

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
