/* $Id: FluidMatSupportT.h,v 1.1 2006-07-13 17:57:28 a-kopacz Exp $ */
#ifndef _FLUID_SUPPORT_T_H_
#define _FLUID_SUPPORT_T_H_

/* base class */
#include "MaterialSupportT.h"

/* direct members */
#include "dArrayT.h"
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class FluidElementT;

/** support for the finite strain Tahoe materials classes */
class FluidMatSupportT: public MaterialSupportT
{
public:

	/** constructor */
	FluidMatSupportT(int ndof, int nip);

	/** \name field values at the integration points
	 * The field values can only be access after the source for the
	 * field source is set using FluidMatSupportT::SetField. */
	/*@{*/
	/** current velocity value at the specified integration point */
	dArrayT Velocity(int ip) const;

	/** current velocity value at the current integration point */
	dArrayT Velocity(void) const;

	/** current pressure value at the specified integration point */
	double Pressure(int ip) const;

	/** current Pressure value at the current integration point */
	double Pressure(void) const;

	/** set the source for the gradient information */
	void SetField(const ArrayT<dArrayT>* fvel_list, const dArrayT* fpres_list);
	/*@}*/

	/** \name field gradients.
	 * Field gradients can only be access after the source for the
	 * gradient information is set using FluidMatSupportT::SetGradient. */
	/*@{*/
	/** field gradient at the current integration point */
	const dMatrixT& VelGrad(void) const;

	/** field gradient at the current integration point */
	const dArrayT& PresGrad(void) const;

	/** field gradient at the specified integration point */
	const dMatrixT& VelGrad(int ip) const;

	/** field gradient at the specified integration point */
	const dArrayT& PresGrad(int ip) const;

	/** set the source for the gradient information */
	void SetGradient(const ArrayT<dMatrixT>* gradient_vel_list, const ArrayT<dArrayT>* gradient_pres_list);
	/*@}*/

	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const FluidElementT* Fluid(void) const { return fFluid; };

	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);
	/*@}*/
	
  private:

	/** field values. Pointer to the array that always contains the
	 * current values of the field over the element being calculated: [nip] */
	const ArrayT<dArrayT>* fVel_list;
	const dArrayT* fPres_list;
	
	/** field gradient. Pointer to the array that always contains the
	 * current values of the field gradient over the element being
	 * calculated: [nip] x [nsd] */
	const ArrayT<dMatrixT>* fGradVel_list;
	const ArrayT<dArrayT>* fGradPres_list;

  	/** pointer to the diffusion element */
	const FluidElementT* fFluid;
};

/* inlines */
inline dArrayT FluidMatSupportT::Velocity(int ip) const
{
#if __option(extended_errorcheck)
	if (!fVel_list) ExceptionT::GeneralFail("FluidMatSupportT::Field");
#endif
	return (*fVel_list)[ip];
}

inline dArrayT FluidMatSupportT::Velocity(void) const
{
#if __option(extended_errorcheck)
	if (!fVel_list) ExceptionT::GeneralFail("FluidMatSupportT::Field");
#endif
	return (*fVel_list)[CurrIP()];
}

inline double FluidMatSupportT::Pressure(int ip) const
{
#if __option(extended_errorcheck)
	if (!fPres_list) ExceptionT::GeneralFail("FluidMatSupportT::Field");
#endif
	return (*fPres_list)[ip];
}

inline double FluidMatSupportT::Pressure(void) const
{
#if __option(extended_errorcheck)
	if (!fPres_list) ExceptionT::GeneralFail("FluidMatSupportT::Field");
#endif
	return (*fPres_list)[CurrIP()];
}

inline const dMatrixT& FluidMatSupportT::VelGrad(int ip) const
{
#if __option(extended_errorcheck)
	if (!fGradVel_list) ExceptionT::GeneralFail("FluidMatSupportT::Gradient");
#endif
	return (*fGradVel_list)[ip];
}

inline const dMatrixT& FluidMatSupportT::VelGrad(void) const
{
#if __option(extended_errorcheck)
	if (!fGradVel_list) ExceptionT::GeneralFail("FluidMatSupportT::Gradient");
#endif
	return (*fGradVel_list)[CurrIP()];
}

inline const dArrayT& FluidMatSupportT::PresGrad(int ip) const
{
#if __option(extended_errorcheck)
	if (!fGradPres_list) ExceptionT::GeneralFail("FluidMatSupportT::Gradient");
#endif
	return (*fGradPres_list)[ip];
}

inline const dArrayT& FluidMatSupportT::PresGrad(void) const
{
#if __option(extended_errorcheck)
	if (!fGradPres_list) ExceptionT::GeneralFail("FluidMatSupportT::Gradient");
#endif
	return (*fGradPres_list)[CurrIP()];
}

} /* namespace Tahoe */

#endif /* _FLUID_SUPPORT_T_H_ */
