#ifndef _GRAD_SS_MAT_SUPPORT_T_H_
#define _GRAD_SS_MAT_SUPPORT_T_H_

/* base class */
#include "SSMatSupportT.h"

/* direct members */
#include "dArrayT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class DorganVoyiadjisMarin;

/** support for the small strain Tahoe materials classes */
class GradSSMatSupportT: public SSMatSupportT
{
public:

	/** constructor */
	GradSSMatSupportT(int nsd, int ndof, int nip);

	/** destructor */
	~GradSSMatSupportT(void);

	/** \name isotropic hardening */
	/*@{*/
	double& LinearR(void) const;
	double& LinearR(int ip) const;
	/*@}*/

	/** \name isotropic hardening from the end of the previous time step */
	/*@{*/
	double& LinearR_last(void) const;
	double& LinearR_last(int ip) const;
	/*@}*/

	/** \name Laplacian isotropic hardening */
	/*@{*/
	double& LinearLaplacianR(void) const;
	double& LinearLaplacianR(int ip) const;
	/*@}*/

	/** \name Laplacian isotropic hardening from the end of the previous time step */
	/*@{*/
	double& LinearLaplacianR_last(void) const;
	double& LinearLaplacianR_last(int ip) const;
	/*@}*/

	/** set source for the isotropic hardening */
	void SetLinearR(const dArrayT* r_List);

	/** set source for the isotropic hardening from the end of the previous time step */
	void SetLinearR_last(const dArrayT* r_last_List);

	/** set source for the Laplacian of isotropic hardening */
	void SetLinearLaplacianR(const dArrayT* lapr_List);

	/** set source for the Laplacian of isotropic hardening from the end of the previous time step */
	void SetLinearLaplacianR_last(const dArrayT* lapr_last_List);

	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const DorganVoyiadjisMarin* DVM(void) const { return fDorganVoyiadjisMarin; };

	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);
	/*@}*/

  private:

  	/** \name return values */
	/*@{*/
 	const dArrayT* fR_List;
  	const dArrayT* fR_last_List;

  	const dArrayT* fLapR_List;
  	const dArrayT* fLapR_last_List;
	/*@}*/

  	/** pointer to the small strain element */
	const DorganVoyiadjisMarin* fDorganVoyiadjisMarin;	
};

/* inlines */
inline double& GradSSMatSupportT::LinearR(void) const
{
	if (!fR_List) throw ExceptionT::kGeneralFail;
	return (*fR_List)[CurrIP()]; 
}

inline double& GradSSMatSupportT::LinearR(int ip) const
{
	if (!fR_List) throw ExceptionT::kGeneralFail;
	return (*fR_List)[ip]; 
}

inline double& GradSSMatSupportT::LinearR_last(void) const
{
	if (!fR_last_List) throw ExceptionT::kGeneralFail;
	return (*fR_last_List)[CurrIP()]; 
}

inline double& GradSSMatSupportT::LinearR_last(int ip) const
{
	if (!fR_last_List) throw ExceptionT::kGeneralFail;
	return (*fR_last_List)[ip]; 
}

inline double& GradSSMatSupportT::LinearLaplacianR(void) const
{
	if (!fLapR_List) throw ExceptionT::kGeneralFail;
	return (*fLapR_List)[CurrIP()]; 
}

inline double& GradSSMatSupportT::LinearLaplacianR(int ip) const
{
	if (!fLapR_List) throw ExceptionT::kGeneralFail;
	return (*fLapR_List)[ip]; 
}

inline double& GradSSMatSupportT::LinearLaplacianR_last(void) const
{
	if (!fLapR_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapR_last_List)[CurrIP()]; 
}

inline double& GradSSMatSupportT::LinearLaplacianR_last(int ip) const
{
	if (!fLapR_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapR_last_List)[ip]; 
}

} /* namespace Tahoe */
#endif /* _GRAD_SS_MAT_SUPPORT_T_H_ */
