/* $Id: GradSSMatSupportT.h,v 1.5 2003-11-21 22:54:37 paklein Exp $ */
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
        GradSSMatSupportT(int nsd, int ndof_disp, int ndor_r, int nip_disp, int nip_r);
        
        /** destructor */
        ~GradSSMatSupportT(void);
        
        /** \name isotropic hardening */
        /*@{*/
        const double& LinearR(void) const;
        const double& LinearR(int ip) const;
        /*@}*/
        
        /** \name isotropic hardening from the end of the previous time step */
        /*@{*/
        const double& LinearR_last(void) const;
        const double& LinearR_last(int ip) const;
        /*@}*/
        
        /** \name Laplacian isotropic hardening */
        /*@{*/
        const double& LinearLaplacianR(void) const;
        const double& LinearLaplacianR(int ip) const;
        /*@}*/
        
        /** \name Laplacian isotropic hardening from the end of the previous time step */
        /*@{*/
        const double& LinearLaplacianR_last(void) const;
        const double& LinearLaplacianR_last(int ip) const;
        /*@}*/
        
        /** set source for the isotropic hardening */
        void SetLinearR(const dArrayT* r_List);
        
        /** set source for the isotropic hardening from the end of the previous time step */
        void SetLinearR_last(const dArrayT* r_last_List);
        
        /** set source for the Laplacian of isotropic hardening */
        void SetLinearLaplacianR(const dArrayT* lapr_List);
        
        /** set source for the Laplacian of isotropic hardening from the end of the previous time step */
        void SetLinearLaplacianR_last(const dArrayT* lapr_last_List);
        
        /** \name dimensions */
        /*@{*/
        /** number of degrees of freedom of iso_hard field (per node) */
        int NumDOF_R(void) const { return fNumDOF_R; };
        
        /** total number of degrees of freedom (per node) */
        int NumDOF_Total(void) const { return fNumDOF_Total; };
        
        /** stress evaluation points per element for iso_hard field */
        int NumIP_R(void) const { return fNumIP_R; };
        /*@}*/
        
        /** \name host code information */
        /*@{*/
        /** return a pointer to the host element. Returns NULL if no
         * no element information in available. The ContinuumElementT
         * pointer is set using MaterialSupportT::SetContinuumElement. */
        const GradSmallStrainT* DVM(void) const { return fGradSmallStrainT; };
        
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
        
        /** \name dimensions */
        /*@{*/
        /** number of degrees of freedom for iso_hard field */
        int fNumDOF_R;
	
        /** total number of degrees of freedom */
        int fNumDOF_Total;
	
        /** number of integration points for iso_hard field */
        int fNumIP_R;
        /*@}*/
  	
        /** pointer to the small strain element */
        const GradSmallStrainT* fGradSmallStrainT;	
};

/* inlines */
inline const double& GradSSMatSupportT::LinearR(void) const
{
        if (!fR_List) throw ExceptionT::kGeneralFail;
        return (*fR_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearR(int ip) const
{
        if (!fR_List) throw ExceptionT::kGeneralFail;
        return (*fR_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearR_last(void) const
{
        if (!fR_last_List) throw ExceptionT::kGeneralFail;
        return (*fR_last_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearR_last(int ip) const
{
        if (!fR_last_List) throw ExceptionT::kGeneralFail;
        return (*fR_last_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearLaplacianR(void) const
{
        if (!fLapR_List) throw ExceptionT::kGeneralFail;
        return (*fLapR_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearLaplacianR(int ip) const
{
        if (!fLapR_List) throw ExceptionT::kGeneralFail;
        return (*fLapR_List)[ip]; 
}

inline const double& GradSSMatSupportT::LinearLaplacianR_last(void) const
{
        if (!fLapR_last_List) throw ExceptionT::kGeneralFail;
        return (*fLapR_last_List)[CurrIP()]; 
}

inline const double& GradSSMatSupportT::LinearLaplacianR_last(int ip) const
{
        if (!fLapR_last_List) throw ExceptionT::kGeneralFail;
        return (*fLapR_last_List)[ip]; 
}

} /* namespace Tahoe */
#endif /* _GRAD_SS_MAT_SUPPORT_T_H_ */
