/* $Id: GradSSMatSupportT.cpp,v 1.5 2004-01-14 19:31:16 rdorgan Exp $ */ 
#include "GradSSMatSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "GradSmallStrainT.h"
#endif

using namespace Tahoe;

/* constructor */
GradSSMatSupportT::GradSSMatSupportT(int nsd, int ndof_disp, int ndof_r, int nip_disp, int nip_r):
        SSMatSupportT(nsd, ndof_disp, nip_disp),
        fGradSmallStrainT(NULL),

        fR_List(NULL),
        fR_last_List(NULL),
        fLapR_List(NULL),
        fLapR_last_List(NULL),

        fNumDOF_R(ndof_r),
        fNumDOF_Total(ndof_disp + ndof_r),

        fNumIP_R(nip_r)
{

}
 
/* destructor */
GradSSMatSupportT::~GradSSMatSupportT(void)
{

}

/* set source for the isotropic hadening */
void GradSSMatSupportT::SetLinearR(const dArrayT* r_List)
{
        /* keep pointer */
        fR_List = r_List;
}

/** set source for the isotropic hardening from the end of the previous time step */
void GradSSMatSupportT::SetLinearR_last(const dArrayT* r_last_List)
{
        /* keep pointer */
        fR_last_List = r_last_List;
}

/* set source for the laplacian of isotropic hadening */
void GradSSMatSupportT::SetLinearLaplacianR(const dArrayT* lapr_List)
{
        /* keep pointer */
        fLapR_List = lapr_List;
}

/** set source for the laplacian of isotropic hardening from the end of the previous time step */
void GradSSMatSupportT::SetLinearLaplacianR_last(const dArrayT* lapr_last_List)
{
        /* keep pointer */
        fLapR_last_List = lapr_last_List;
}

/* set the element group pointer */
void GradSSMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
        /* inherited */
        SSMatSupportT::SetContinuumElement(p);

#ifdef CONTINUUM_ELEMENT
        /* cast to GradSmallStrainT pointer */
        fGradSmallStrainT = TB_DYNAMIC_CAST(const GradSmallStrainT*, p);
//        fGradSmallStrainT = dynamic_cast<const GradSmallStrainT*>(p);
#endif
}
