/* $Id: DiffusionMatSupportT.h,v 1.3 2003-01-29 07:34:33 paklein Exp $ */
#ifndef _DIFF_MAT_SUPPORT_T_H_
#define _DIFF_MAT_SUPPORT_T_H_

/* base class */
#include "MaterialSupportT.h"

/* direct members */
#include "dArrayT.h"
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class DiffusionElementT;

/** support for the finite strain Tahoe materials classes */
class DiffusionMatSupportT: public MaterialSupportT
{
public:

	/** constructor */
	DiffusionMatSupportT(int nsd, int ndof, int nip);

	/** \name field gradients.
	 * Field gradients can only be access after the source for the
	 * gradient information is set using DiffusionMatSupportT::SetGradient. */
	/*@{*/
	/** field gradient at the current integration point */
	const dArrayT& Gradient(void) const;

	/** field gradient at the specified integration point */
	const dArrayT& Gradient(int ip) const;

	/** set the source for the gradient information */
	void SetGradient(const ArrayT<dArrayT>* gradient_list);
	/*@}*/

	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const DiffusionElementT* Diffusion(void) const { return fDiffusion; };

	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);
	/*@}*/
	
  private:

	/** field gradient. Pointer to the array that always contains the
	 * current values of the field gradient over the element being
	 * calculated: [nip] x [nsd] */
	const ArrayT<dArrayT>* fGradient_list;

  	/** pointer to the diffusion element */
	const DiffusionElementT* fDiffusion;
};

/* inlines */
inline const dArrayT& DiffusionMatSupportT::Gradient(int ip) const
{
	if (!fGradient_list) throw ExceptionT::kGeneralFail;
	return (*fGradient_list)[ip];
}

inline const dArrayT& DiffusionMatSupportT::Gradient(void) const
{
	if (!fGradient_list) throw ExceptionT::kGeneralFail;
	return (*fGradient_list)[CurrIP()];
}


} /* namespace Tahoe */
#endif /* _DIFF_MAT_SUPPORT_T_H_ */
