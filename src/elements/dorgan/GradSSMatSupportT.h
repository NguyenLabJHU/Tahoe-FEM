/* $Id: GradSSMatSupportT.h,v 1.8 2004-06-09 00:25:53 rdorgan Exp $ */
#ifndef _GRAD_SS_MAT_SUPPORT_T_H_
#define _GRAD_SS_MAT_SUPPORT_T_H_

/* base class */
#include "SSMatSupportT.h"

/* direct members */
#include "dArrayT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/** support for the small strain Tahoe materials classes */
class GradSSMatSupportT: public SSMatSupportT
{
public:

	/** constructor */
	GradSSMatSupportT(int nsd, int ndof_disp, int ndor_field, int nip_disp, int nip_field);
	
	/** destructor */
	~GradSSMatSupportT(void);
	
	/** \name field */
	/*@{*/
	const double& LinearField(void) const;
	const double& LinearField(int ip) const;
	/*@}*/
	
	/** \name field from the end of the previous time step */
	/*@{*/
	const double& LinearField_last(void) const;
	const double& LinearField_last(int ip) const;
	/*@}*/
	
	/** \name gradient field */
	/*@{*/
	const double& LinearGradField(void) const;
	const double& LinearGradField(int ip) const;
	/*@}*/
	
	/** \name gradient field from the end of the previous time step */
	/*@{*/
	const double& LinearGradField_last(void) const;
	const double& LinearGradField_last(int ip) const;
	/*@}*/
	
	/** \name Laplacian field */
	/*@{*/
	const double& LinearLapField(void) const;
	const double& LinearLapField(int ip) const;
	/*@}*/
	
	/** \name Laplacian field from the end of the previous time step */
	/*@{*/
	const double& LinearLapField_last(void) const;
	const double& LinearLapField_last(int ip) const;
	/*@}*/
	
	/** set source for the field */
	void SetLinearField(const dArrayT* field_List);
	
	/** set source for the field from the end of the previous time step */
	void SetLinearField_last(const dArrayT* field_last_List);
	
	/** set source for the gradient of field */
	void SetLinearGradField(const dArrayT* gradfield_List);
	
	/** set source for the gradient of field from the end of the previous time step */
	void SetLinearGradField_last(const dArrayT* gradfield_last_List);
	
	/** set source for the Laplacian of field */
	void SetLinearLapField(const dArrayT* lapfield_List);
	
	/** set source for the Laplacian of field from the end of the previous time step */
	void SetLinearLapField_last(const dArrayT* lapfield_last_List);
	
	/** \name dimensions */
	/*@{*/
	/** number of degrees of freedom of field (per node) */
	int NumDOF_Field(void) const { return fNumDOF_Field; };
	
	/** total number of degrees of freedom (per node) */
	int NumDOF_Total(void) const { return fNumDOF_Total; };
	
	/** stress evaluation points per element for field */
	int NumIP_Field(void) const { return fNumIP_Field; };
	/*@}*/
	
	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);
	/*@}*/
	
private:
	
	/** \name return values */
	/*@{*/
	const dArrayT* fField_List;
	const dArrayT* fField_last_List;
	
	const dArrayT* fGradField_List;
	const dArrayT* fGradField_last_List;

	const dArrayT* fLapField_List;
	const dArrayT* fLapField_last_List;
	/*@}*/
	
	/** \name dimensions */
	/*@{*/
	/** number of degrees of freedom for field */
	int fNumDOF_Field;
	
	/** total number of degrees of freedom */
	int fNumDOF_Total;
	
	/** number of integration points for field */
	int fNumIP_Field;
	/*@}*/
};

/* inlines */
inline const double& GradSSMatSupportT::LinearField(void) const
{
	if (!fField_List) throw ExceptionT::kGeneralFail;
	return (*fField_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearField(int ip) const
{
	if (!fField_List) throw ExceptionT::kGeneralFail;
	return (*fField_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearField_last(void) const
{
	if (!fField_last_List) throw ExceptionT::kGeneralFail;
	return (*fField_last_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearField_last(int ip) const
{
	if (!fField_last_List) throw ExceptionT::kGeneralFail;
	return (*fField_last_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearGradField(void) const
{
	if (!fGradField_List) throw ExceptionT::kGeneralFail;
	return (*fGradField_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearGradField(int ip) const
{
	if (!fGradField_List) throw ExceptionT::kGeneralFail;
	return (*fGradField_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearGradField_last(void) const
{
	if (!fGradField_last_List) throw ExceptionT::kGeneralFail;
	return (*fGradField_last_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearGradField_last(int ip) const
{
	if (!fGradField_last_List) throw ExceptionT::kGeneralFail;
	return (*fGradField_last_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearLapField(void) const
{
	if (!fLapField_List) throw ExceptionT::kGeneralFail;
	return (*fLapField_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearLapField(int ip) const
{
	if (!fLapField_List) throw ExceptionT::kGeneralFail;
	return (*fLapField_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearLapField_last(void) const
{
	if (!fLapField_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapField_last_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearLapField_last(int ip) const
{
	if (!fLapField_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapField_last_List)[ip]; 
}

} /* namespace Tahoe */
#endif /* _GRAD_SS_MAT_SUPPORT_T_H_ */
