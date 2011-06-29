/* $Id: DiffusionMatSupportT.h,v 1.5 2004-07-15 08:26:22 paklein Exp $ */
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
	DiffusionMatSupportT(int ndof, int nip);

	/** \name field values at the integration points
	 * The field values can only be access after the source for the
	 * field source is set using DiffusionMatSupportT::SetField. */
	/*@{*/
	/** field value at the current integration point */
	double Field(void) const;

	/** field value at the specified integration point */
	double Field(int ip) const;

	/** set the source for the gradient information */
	void SetField(const dArrayT* field_list);
	/*@}*/

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

	/** field values. Pointer to the array that always contains the
	 * current values of the field over the element being calculated: [nip] */
	const dArrayT* fField_list;

	/** field gradient. Pointer to the array that always contains the
	 * current values of the field gradient over the element being
	 * calculated: [nip] x [nsd] */
	const ArrayT<dArrayT>* fGradient_list;

  	/** pointer to the diffusion element */
	const DiffusionElementT* fDiffusion;
};

/* inlines */
inline double DiffusionMatSupportT::Field(int ip) const
{
#if __option(extended_errorcheck)
	if (!fField_list) ExceptionT::GeneralFail("DiffusionMatSupportT::Field");
#endif
	return (*fField_list)[ip];
}

inline double DiffusionMatSupportT::Field(void) const
{
#if __option(extended_errorcheck)
	if (!fField_list) ExceptionT::GeneralFail("DiffusionMatSupportT::Field");
#endif
	return (*fField_list)[CurrIP()];
}

inline const dArrayT& DiffusionMatSupportT::Gradient(int ip) const
{
#if __option(extended_errorcheck)
	if (!fGradient_list) ExceptionT::GeneralFail("DiffusionMatSupportT::Gradient");
#endif
	return (*fGradient_list)[ip];
}

inline const dArrayT& DiffusionMatSupportT::Gradient(void) const
{
#if __option(extended_errorcheck)
	if (!fGradient_list) ExceptionT::GeneralFail("DiffusionMatSupportT::Gradient");
#endif
	return (*fGradient_list)[CurrIP()];
}

} /* namespace Tahoe */

#endif /* _DIFF_MAT_SUPPORT_T_H_ */
