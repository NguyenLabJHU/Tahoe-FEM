/* $Id: GradSSMatSupportT.h,v 1.9 2004-07-20 23:16:50 rdorgan Exp $ */
#ifndef _GRAD_SS_MAT_SUPPORT_T_H_
#define _GRAD_SS_MAT_SUPPORT_T_H_

/* base class */
#include "SSMatSupportT.h"

/* direct members */
#include "dArrayT.h"
#include "dSymMatrixT.h"

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
	const double& LinearPMultiplier(void) const;
	const double& LinearPMultiplier(int ip) const;
	/*@}*/
	
	/** \name field from the end of the previous time step */
	/*@{*/
	const double& LinearPMultiplier_last(void) const;
	const double& LinearPMultiplier_last(int ip) const;
	/*@}*/
	
	/** \name gradient field */
	/*@{*/
	const double& LinearGradPMultiplier(void) const;
	const double& LinearGradPMultiplier(int ip) const;
	/*@}*/
	
	/** \name gradient field from the end of the previous time step */
	/*@{*/
	const double& LinearGradPMultiplier_last(void) const;
	const double& LinearGradPMultiplier_last(int ip) const;
	/*@}*/
	
	/** \name Laplacian field */
	/*@{*/
	const double& LinearLapPMultiplier(void) const;
	const double& LinearLapPMultiplier(int ip) const;
	/*@}*/
	
	/** \name Laplacian field from the end of the previous time step */
	/*@{*/
	const double& LinearLapPMultiplier_last(void) const;
	const double& LinearLapPMultiplier_last(int ip) const;
	/*@}*/
	
	/** set source for the field */
	void SetLinearPMultiplier(const dArrayT* pmultiplier_List);
	
	/** set source for the field from the end of the previous time step */
	void SetLinearPMultiplier_last(const dArrayT* pmultiplier_last_List);
	
	/** set source for the gradient of field */
	void SetLinearGradPMultiplier(const dArrayT* gradpmultiplier_List);
	
	/** set source for the gradient of field from the end of the previous time step */
	void SetLinearGradPMultiplier_last(const dArrayT* gradpmultiplier_last_List);
	
	/** set source for the Laplacian of field */
	void SetLinearLapPMultiplier(const dArrayT* lappmultiplier_List);
	
	/** set source for the Laplacian of field from the end of the previous time step */
	void SetLinearLapPMultiplier_last(const dArrayT* lappmultiplier_last_List);
	
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
	const dArrayT* fPMultiplier_List;
	const dArrayT* fPMultiplier_last_List;
	
	const dArrayT* fGradPMultiplier_List;
	const dArrayT* fGradPMultiplier_last_List;

	const dArrayT* fLapPMultiplier_List;
	const dArrayT* fLapPMultiplier_last_List;
	/*@}*/
	
  	/** pointer to the small strain element */
	const GradSmallStrainT* fGradSmallStrain;	
};

/* inlines */
inline const double& GradSSMatSupportT::LinearPMultiplier(void) const
{
	if (!fPMultiplier_List) throw ExceptionT::kGeneralFail;
	return (*fPMultiplier_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearPMultiplier(int ip) const
{
	if (!fPMultiplier_List) throw ExceptionT::kGeneralFail;
	return (*fPMultiplier_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearPMultiplier_last(void) const
{
	if (!fPMultiplier_last_List) throw ExceptionT::kGeneralFail;
	return (*fPMultiplier_last_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearPMultiplier_last(int ip) const
{
	if (!fPMultiplier_last_List) throw ExceptionT::kGeneralFail;
	return (*fPMultiplier_last_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearGradPMultiplier(void) const
{
	if (!fGradPMultiplier_List) throw ExceptionT::kGeneralFail;
	return (*fGradPMultiplier_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearGradPMultiplier(int ip) const
{
	if (!fGradPMultiplier_List) throw ExceptionT::kGeneralFail;
	return (*fGradPMultiplier_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearGradPMultiplier_last(void) const
{
	if (!fGradPMultiplier_last_List) throw ExceptionT::kGeneralFail;
	return (*fGradPMultiplier_last_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearGradPMultiplier_last(int ip) const
{
	if (!fGradPMultiplier_last_List) throw ExceptionT::kGeneralFail;
	return (*fGradPMultiplier_last_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearLapPMultiplier(void) const
{
	if (!fLapPMultiplier_List) throw ExceptionT::kGeneralFail;
	return (*fLapPMultiplier_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearLapPMultiplier(int ip) const
{
	if (!fLapPMultiplier_List) throw ExceptionT::kGeneralFail;
	return (*fLapPMultiplier_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearLapPMultiplier_last(void) const
{
	if (!fLapPMultiplier_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapPMultiplier_last_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearLapPMultiplier_last(int ip) const
{
	if (!fLapPMultiplier_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapPMultiplier_last_List)[ip]; 
}

} /* namespace Tahoe */
#endif /* _GRAD_SS_MAT_SUPPORT_T_H_ */
