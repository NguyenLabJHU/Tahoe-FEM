/* $Id: UpdatedLagrangianT.h,v 1.7 2002-07-17 00:02:10 paklein Exp $ */
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

	dMatrixT fCauchyStress;	// symmetric Cauchy stress tensor
	dMatrixT fStressStiff;	// compact stress stiffness contribution
	dMatrixT fGradNa;       // shape function gradients matrix
	
	/* arrays with local ordering */
	LocalArrayT fLocCurrCoords;	// current coords with local ordering
		
	/* workspace */
	dArrayT fTemp2;
};

} // namespace Tahoe 
#endif /* _UPDATED_LAGRANGIAN_T_H_ */
