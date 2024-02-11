/* $Id: GradSSMatSupportT.h,v 1.10 2004-09-02 18:25:04 rdorgan Exp $ */
#ifndef _GRAD_SS_MAT_SUPPORT_T_H_
#define _GRAD_SS_MAT_SUPPORT_T_H_

/* base class */
#include "SSMatSupportT.h"

/* direct members */
#include "dArrayT.h"
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class GradSmallStrainT;

/** support for the small strain Tahoe materials classes */
class GradSSMatSupportT: public SSMatSupportT
{
public:

	/** constructor */
	GradSSMatSupportT(int ndof_disp, int ndof_field, int nip_disp, int nip_field);
	
	/** destructor */
	~GradSSMatSupportT(void);
	
	/** \name field */
	/*@{*/
	const dMatrixT& LinearPMultiplier(void) const;
	const dMatrixT& LinearPMultiplier(int ip) const;
	const dMatrixT& LinearPMultiplier_last(void) const;
	const dMatrixT& LinearPMultiplier_last(int ip) const;
	/*@}*/
	
	/** \name gradient field */
	/*@{*/
	const dMatrixT& LinearGradPMultiplier(void) const;
	const dMatrixT& LinearGradPMultiplier(int ip) const;
	const dMatrixT& LinearGradPMultiplier_last(void) const;
	const dMatrixT& LinearGradPMultiplier_last(int ip) const;
	/*@}*/
	
	/** \name Laplacian field */
	/*@{*/
	const dMatrixT& LinearLapPMultiplier(void) const;
	const dMatrixT& LinearLapPMultiplier(int ip) const;
	const dMatrixT& LinearLapPMultiplier_last(void) const;
	const dMatrixT& LinearLapPMultiplier_last(int ip) const;
	/*@}*/

	/** set source for the field */
	void SetLinearPMultiplier(const ArrayT<dMatrixT>* pmultiplier_List);
	void SetLinearPMultiplier_last(const ArrayT<dMatrixT>* pmultiplier_last_List);
	
	/** set source for the gradient of field */
	void SetLinearGradPMultiplier(const ArrayT<dMatrixT>* gradpmultiplier_List);
	void SetLinearGradPMultiplier_last(const ArrayT<dMatrixT>* gradpmultiplier_last_List);
	
	/** set source for the Laplacian of field */
	void SetLinearLapPMultiplier(const ArrayT<dMatrixT>* lappmultiplier_List);
	void SetLinearLapPMultiplier_last(const ArrayT<dMatrixT>* lappmultiplier_last_List);
	
	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const GradSmallStrainT* GradSmallStrain(void) const { return fGradSmallStrain; };

	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);
	/*@}*/
	
private:
	
	/** \name return values */
	/*@{*/
	const ArrayT<dMatrixT>* fPMultiplier_List;
	const ArrayT<dMatrixT>* fPMultiplier_last_List;
	
	const ArrayT<dMatrixT>* fGradPMultiplier_List;
	const ArrayT<dMatrixT>* fGradPMultiplier_last_List;

	const ArrayT<dMatrixT>* fLapPMultiplier_List;
	const ArrayT<dMatrixT>* fLapPMultiplier_last_List;
	/*@}*/
	
  	/** pointer to the small strain element */
	const GradSmallStrainT* fGradSmallStrain;	
};

/* inlines */
inline const dMatrixT& GradSSMatSupportT::LinearPMultiplier(void) const
{
	if (!fPMultiplier_List) throw ExceptionT::kGeneralFail;
	return (*fPMultiplier_List)[CurrIP()]; 
}

inline const dMatrixT& GradSSMatSupportT::LinearPMultiplier(int ip) const
{
	if (!fPMultiplier_List) throw ExceptionT::kGeneralFail;
	return (*fPMultiplier_List)[ip]; 
}

inline const dMatrixT& GradSSMatSupportT::LinearPMultiplier_last(void) const
{
	if (!fPMultiplier_last_List) throw ExceptionT::kGeneralFail;
	return (*fPMultiplier_last_List)[CurrIP()]; 
}

inline const dMatrixT& GradSSMatSupportT::LinearPMultiplier_last(int ip) const
{
	if (!fPMultiplier_last_List) throw ExceptionT::kGeneralFail;
	return (*fPMultiplier_last_List)[ip]; 
}

inline const dMatrixT& GradSSMatSupportT::LinearGradPMultiplier(void) const
{
	if (!fGradPMultiplier_List) throw ExceptionT::kGeneralFail;
	return (*fGradPMultiplier_List)[CurrIP()]; 
}

inline const dMatrixT& GradSSMatSupportT::LinearGradPMultiplier(int ip) const
{
	if (!fGradPMultiplier_List) throw ExceptionT::kGeneralFail;
	return (*fGradPMultiplier_List)[ip]; 
}

inline const dMatrixT& GradSSMatSupportT::LinearGradPMultiplier_last(void) const
{
	if (!fGradPMultiplier_last_List) throw ExceptionT::kGeneralFail;
	return (*fGradPMultiplier_last_List)[CurrIP()]; 
}

inline const dMatrixT& GradSSMatSupportT::LinearGradPMultiplier_last(int ip) const
{
	if (!fGradPMultiplier_last_List) throw ExceptionT::kGeneralFail;
	return (*fGradPMultiplier_last_List)[ip]; 
}

inline const dMatrixT& GradSSMatSupportT::LinearLapPMultiplier(void) const
{
	if (!fLapPMultiplier_List) throw ExceptionT::kGeneralFail;
	return (*fLapPMultiplier_List)[CurrIP()]; 
}

inline const dMatrixT& GradSSMatSupportT::LinearLapPMultiplier(int ip) const
{
	if (!fLapPMultiplier_List) throw ExceptionT::kGeneralFail;
	return (*fLapPMultiplier_List)[ip]; 
}

inline const dMatrixT& GradSSMatSupportT::LinearLapPMultiplier_last(void) const
{
	if (!fLapPMultiplier_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapPMultiplier_last_List)[CurrIP()]; 
}

inline const dMatrixT& GradSSMatSupportT::LinearLapPMultiplier_last(int ip) const
{
	if (!fLapPMultiplier_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapPMultiplier_last_List)[ip]; 
}

} /* namespace Tahoe */
#endif /* _GRAD_SS_MAT_SUPPORT_T_H_ */
