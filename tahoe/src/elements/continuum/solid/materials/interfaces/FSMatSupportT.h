/* $Id: FSMatSupportT.h,v 1.1.2.2 2002-10-30 09:18:11 paklein Exp $ */
#ifndef _FD_MAT_SUPPORT_T_H_
#define _FD_MAT_SUPPORT_T_H_

/* base class */
#include "MaterialSupportT.h"

/* direct members */
#include "ArrayT.h"
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

	/** set source for the deformation gradient */
	void SetDeformationGradient(const ArrayT<dMatrixT>* F_List);

	/** set source for the deformation gradient from the end of the previous time step */
	void SetDeformationGradient_last(const ArrayT<dMatrixT>* F_last_List);
	/*@}*/

	/** \name field gradients */
	/*@{*/
	/** compute field gradients with respect to current coordinates at the current integration point.
	 * Returns true of the operation is supported, false otherwise. */
	bool ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const;

	/** compute field gradients with respect to current coordinates at the specified integration point.
	 * Returns true of the operation is supported, false otherwise. */
	bool ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;

	/** compute field gradients with respect to reference coordinates at the current integration point.
	 * Returns true of the operation is supported, false otherwise. */	
	bool ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u) const;

	/** compute field gradients with respect to reference coordinates at the specified integration point.
	 * Returns true of the operation is supported, false otherwise. */
	bool ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;
	/*@}*/

	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const FiniteStrainT* FiniteStrain(void) const { return fFiniteStrain; };

	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);
	
	/** nodal temperatures. Returns NULL if not available */
	const LocalArrayT* Temperatures(void) const;

	/** nodal temperatures from the last time step. Returns NULL if 
	 * not available */
	const LocalArrayT* LastTemperatures(void) const;

	/** return a pointer the specified local array, or NULL if the array is not
	 * available. During calls the materials routines these will contain the
	 * values for the current element. */
	virtual const LocalArrayT* LocalArray(LocalArrayT::TypeT t) const;
	/*@}*/

  private:

  	/** \name sources for deformation gradients */
  	/*@{*/
  	const ArrayT<dMatrixT>* fF_List; /**< deformation gradient */
  	const ArrayT<dMatrixT>* fF_last_List; /**< last deformation gradient */
  	/*@}*/

  	/** pointer to the finite strain element */
	const FiniteStrainT* fFiniteStrain;	
};

/* inlines */
const dMatrixT& FSMatSupportT::DeformationGradient(void) const
{
	if (!fF_List) throw ExceptionT::kGeneralFail;
	return (*fF_List)[CurrIP()]; 
}

const dMatrixT& FSMatSupportT::DeformationGradient(int ip) const 
{ 
	if (!fF_List) throw ExceptionT::kGeneralFail;
	return (*fF_List)[ip]; 
}

const dMatrixT& FSMatSupportT::DeformationGradient_last(void) const 
{
	if (!fF_last_List) throw ExceptionT::kGeneralFail;
	return (*fF_last_List)[CurrIP()]; 
}

const dMatrixT& FSMatSupportT::DeformationGradient_last(int ip) const 
{
	if (!fF_last_List) throw ExceptionT::kGeneralFail;
	return (*fF_last_List)[ip]; 
}

} /* namespace Tahoe */
#endif /* _FD_MAT_SUPPORT_T_H_ */
