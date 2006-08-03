/* $Id: FSFiberMatSupportT.h,v 1.1 2006-08-03 01:10:41 thao Exp $ */
#ifndef _FD__FIBER_MAT_SUPPORT_T_H_
#define _FD__FIBER_MAT_SUPPORT_T_H_

/* base class */
#include "FSMatSupportT.h"

/* direct members */
#include "ArrayT.h"
#include "dMatrixT.h"
#include "ParameterListT.h"

namespace Tahoe {

/* forward declarations */
class UpLagFiberCompT;

/** support for the finite strain Tahoe materials classes */
class FSFiberMatSupportT: public FSMatSupportT
{
public:

	/** constructor */
	FSFiberMatSupportT(int ndof, int nip);

	/** \name retrieve orientation vectors P_vec of of fiber families for current element */
	/*@{*/
	const dArray2DT& Fiber_Vec(void) const;

	/** \name returns the number of fiber families for current element*/
	const int NumFibers(void) const;

	/** \name retrieve orientation vectors P_vec of of fiber families for given element */
	/*@{*/
	const dArray2DT& Fiber_Vec(const int elem) const;

	/** \name returns the number of fiber families for given element*/
	const int NumFibers(const int elem) const;
	
	/** set source for the deformation gradient */
	void SetFibers(const ArrayT<dArray2DT>* Fiber_list);

	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const UpLagFiberCompT* FiberElement(void) const { return fFiberElement; };

	virtual void SetContinuumElement(const ContinuumElementT* p);
	
  private:

  	/** \name sources for deformation gradients */
  	/*@{*/
  	const ArrayT<dArray2DT>* fFiber_list; /**< fiber orientation vectors */
  	/*@}*/

  	/** pointer to the finite strain element */
	const UpLagFiberCompT* fFiberElement;	
};

/* inlines */
inline const dArray2DT& FSFiberMatSupportT::Fiber_Vec(void) const
{
	if (!fFiber_list) throw ExceptionT::kGeneralFail;
	return (*fFiber_list)[CurrElementNumber()]; 
}

inline const dArray2DT& FSFiberMatSupportT::Fiber_Vec(const int elem) const
{
	if (!fFiber_list) throw ExceptionT::kGeneralFail;
	return (*fFiber_list)[elem]; 
}

inline const int FSFiberMatSupportT::NumFibers(void) const
{
	if (!fFiber_list) throw ExceptionT::kGeneralFail;
	return (*fFiber_list)[CurrElementNumber()].MajorDim(); 
}

inline const int FSFiberMatSupportT::NumFibers(const int elem) const
{
	if (!fFiber_list) throw ExceptionT::kGeneralFail;
	return (*fFiber_list)[elem].MajorDim(); 
}

} /* namespace Tahoe */
#endif /* _FD__FIBER_MAT_SUPPORT_T_H_ */
