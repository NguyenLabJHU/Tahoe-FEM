/* $Id: SMRPrimitiveT.h,v 1.1 2006-07-27 13:20:08 kyonten Exp $ */
/* created: Majid T. Manzari (04/16/2003)                */
/*
/* Base class for a nonassociative, small strain,        */
/* pressure dependent plasticity model with nonlinear    */
/* isotropic hardening/softening.                        */
/* The model is consistent with the traction sepration   */
/* cohesive surface models MR2DT and MR_RP2D.            */

/* The model is simplified version of MR model */

#ifndef _SMR_PRIMITIVET_H_
#define _SMR_PRIMITIVET_H_

/* base class */
#include "ParameterInterfaceT.h"

namespace Tahoe {

/* forward declarations */
class dSymMatrixT;

class SMRPrimitiveT: public ParameterInterfaceT
{
  public:

    /* constructor */
    SMRPrimitiveT(void);

    /* destructor */
    virtual ~SMRPrimitiveT(void);

    /** \name implementation of the ParameterInterfaceT interface */
    /*@{*/
    /** describe the parameters needed by the interface */
    virtual void DefineParameters(ParameterListT& list) const;

    /** accept parameter list */
    virtual void TakeParameterList(const ParameterListT& list);
    /*@}*/
    
    double YieldCondition(const dSymMatrixT& devstress, 
            const double meanstress, double ftan_phi, double fc) const;

  protected:
    
    //double fGf_I;    /* Mode_I Fracture Energy */
    double fGf_II;   /* Mode_II Fracture Energy */

    /* Inelastic response parameters */
    //double fchi_p; /* peak tensile strength*/  
    //double fchi_r; /* residual tensile strength */
    //double fc_p;   /* peak cohesion */
    //double fc_r;   /* residual cohesion */
    double fc;   /* coefficient of cohesion */
    double fphi_p; /* peak friction angle */
    double fphi_r; /* critical state friction angle */
    double fpsi_p; /* peak dilation angle */
    //double falpha_chi; /* Coefficient of chi degredation */
    //double falpha_c; /* Coefficient of c degredation */
    double falpha_phi; /*  Coefficient of phi degredation */
    double falpha_psi; /*  Coefficient of psi degredatione */
    double fTol_1;    /*  Tolerance for Yield Function */
    double fTol_2; /*  Tolerance for Residuals */
};

} // namespace Tahoe 
#endif /* _SMR_PRIMITIVET_H_ */
