/* $Id: UpdatedLagrangianAxiT.h,v 1.1 2004-02-03 08:24:57 paklein Exp $ */
#ifndef _UPDATED_LAGRANGIAN_AXI_T_H_
#define _UPDATED_LAGRANGIAN_AXI_T_H_

/* base class */
#include "FiniteStrainAxiT.h"

namespace Tahoe {

/** update Lagrangian, finite strain, axisymmatrix solid */
class UpdatedLagrangianAxiT: public FiniteStrainAxiT
{
public:

	/* constructors */
	UpdatedLagrangianAxiT(const ElementSupportT& support, const FieldT& field);

	/* destructors */
	virtual ~UpdatedLagrangianAxiT(void);

	/* data initialization */
	virtual void Initialize(void);
		
protected:

	/* initialization functions */
	virtual void SetShape(void);

	/* form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/* form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/* calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

protected:

  	/** array of shape functions */
  	dArrayT fIPShape;

  	/** 2D-axis stress */
  	dSymMatrixT fStress2D_axi;

	/** \name work space */
	/*@{*/
	dMatrixT fStressMat;    /**< space for a stress 3D tensor */
	dMatrixT fStressStiff;	/**< "compact" stress stiffness contribution: [nen] x [nen] */
	dMatrixT fGradNa;       /**< shape function gradients matrix: [nsd] x [nen] */
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _UPDATED_LAGRANGIAN_AXI_T_H_ */
