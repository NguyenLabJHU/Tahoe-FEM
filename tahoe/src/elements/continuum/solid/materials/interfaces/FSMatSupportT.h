/* $Id: FSMatSupportT.h,v 1.1.2.1 2002-10-28 06:49:15 paklein Exp $ */
#ifndef _FD_MAT_SUPPORT_T_H_
#define _FD_MAT_SUPPORT_T_H_

/* base class */
#include "MaterialSupportT.h"

/* direct members */
#include "dArrayT.h"
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class FiniteStrainT;

/** support for the finite strain Tahoe materials classes */
class FSMatSupportT: public MaterialSupportT
{
public:

	/** constructor */
	FSMatSupportT(int nsd, int ndof, int nip);

	/** destructor */
	~FSMatSupportT(void);

	/** \name deformation gradients */
	/*@{*/
	/** total deformation gradient at the current integration point */
	const dMatrixT& DeformationGradient(void) const;

	/** total deformation gradient at the specified integration point */
	const dMatrixT& DeformationGradient(int ip) const;

	/** total strain from the end of the previous time step at the current integration point */
	const dMatrixT& DeformationGradient_last(void) const;

	/** total strain from the end of the previous time step at the specified integration point */
	const dMatrixT& DeformationGradient_last(int ip) const;
	/*@}*/

	/** \name field gradients */
	/*@{*/
	/** compute field gradients with respect to current coordinates at the current integration point */
	void ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const;

	/** compute field gradients with respect to current coordinates at the specified integration point */
	void ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;

	/** compute field gradients with respect to reference coordinates at the current integration point */
	void ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u) const;

	/** compute field gradients with respect to reference coordinates at the specified integration point */
	void ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;
	/*@}*/

	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const FiniteStrainT* FiniteStrain(void) const { return fFiniteStrain; };

	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);
	
	/** return a pointer the specified local array, or NULL if the array is not
	 * available. During calls the materials routines these will contain the
	 * values for the current element. */
	virtual const LocalArrayT* LocalArray(LocalArrayT::TypeT t) const;
	/*@}*/

  private:

  	/** \name work space  */
  	/*@{*/
  	ArrayT<dMatrixT> fF_List;      /**< deformation gradient */
  	dArrayT          fF_all;       /**< grouped memory for all deformation gradients */
  	ArrayT<dMatrixT> fF_last_List; /**< last deformation gradient */
  	dArrayT          fF_last_all;  /**< grouped memory for all last deformation gradients */
  	/*@}*/

  	/** pointer to the finite strain element */
	const FiniteStrainT* fFiniteStrain;	
};

} /* namespace Tahoe */
#endif /* _FD_MAT_SUPPORT_T_H_ */
