/* $Id: UpdatedLangrianMF.h,v 1.1 2003-03-19 18:46:05 thao Exp $ */

#ifndef _UPDATED_LAGRANGIAN_MF_H_
#define _UPDATED_LAGRANGIAN_MF_H_

/* base class */
#include "FiniteStrainMF.h"

/* direct members */
#include "dMatrixT.h"

namespace Tahoe {

/** update Lagrangian, finite strain solid with material force calc*/
class UpdatedLagrangianMF: public FiniteStrainMF
{
public:

	/* constructors */
	UpdatedLagrangianMF(const ElementSupportT& support, const FieldT& field);

	/* destructors */
	virtual ~UpdatedLagrangianMF(void);

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
#endif /* _UPDATED_LAGRANGIAN_MF_H_ */
