/* $Id: UpdatedLagrangianT.h,v 1.8 2002-10-10 01:38:22 paklein Exp $ */
/* created: paklein (07/03/1996) */

#ifndef _UPDATED_LAGRANGIAN_T_H_
#define _UPDATED_LAGRANGIAN_T_H_

/* base class */
#include "FiniteStrainT.h"

/* direct members */
#include "dMatrixT.h"

namespace Tahoe {

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

	/* calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

protected:

	/** \name work space */
	/*@{*/
	dMatrixT fCauchyStress;	/**< matrix for Cauchy stress tensor: [nsd] x [nsd] */
	dMatrixT fStressStiff;	/**< "compact" stress stiffness contribution: [nen] x [nen] */
	dMatrixT fGradNa;       /**< shape function gradients matrix: [nsd] x [nen] */
	/*@}*/
	
	/** current coords with local ordering */
	LocalArrayT fLocCurrCoords;	
};

} // namespace Tahoe 
#endif /* _UPDATED_LAGRANGIAN_T_H_ */
