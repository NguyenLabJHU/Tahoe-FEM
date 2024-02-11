/* $Id: GRAD_MRPrimitiveT.h,v 1.5 2006-08-23 21:57:03 kyonten Exp $ */
/* created: Karma Yonten (03/04/2004)                   
   Gradient Enhanced MR Model
*/
   
/* base class for a nonassociative, small strain,        */
/* pressure dependent gradient plasticity model          */
/* with nonlinear isotropic hardening/softening.         */
/* The model is consistent with the traction separation  */
/* cohesive surface models MR2DT and MR_RP2D.            */

#ifndef _GRAD_MR_PRIMITIVET_H_ 
#define _GRAD_MR_PRIMITIVET_H_

/* base class */
#include "ParameterInterfaceT.h"

namespace Tahoe 
{

/* forward declarations */
class dSymMatrixT;

class GRAD_MRPrimitiveT: public ParameterInterfaceT
{
  public:

	/* constructor */
	GRAD_MRPrimitiveT(void);

	/* destructor */
	virtual ~GRAD_MRPrimitiveT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
	double YieldCondition(const dSymMatrixT& devstress, 
            const double meanstress, double fchi,
            double fc, double ftan_phi) const;

  protected:
	
	double fGf_I;    /* Mode_I Fracture Energy */
	double fGf_II;   /* Mode_II Fracture Energy */

	/* length scale parameters */
	double flse_v; /* pore space length scale (elastic) */
	double flse_s; /* grain size length scale (elastic) */
	double flsp_v; /* pore space length scale (plastic) */
	double flsp_s; /* grain size length scale (plastic) */
	
	/* inelastic response parameters */
	double fchi_p; /* peak tensile strength*/  
	double fchi_r; /* residual tensile strength */
	double fc_p;   /* peak cohesion */
	double fc_r;   /* residual cohesion */
	double fphi_p; /* peak friction angle */
	double fphi_r; /* critical state friction angle */
	double fpsi_p; /* peak dilation angle */
	double falpha_chi; /* Coefficient of chi degredation */
	double falpha_c; /* Coefficient of c degredation */
	double falpha_phi; /*  Coefficient of phi degredation */
	double falpha_psi; /*  Coefficient of psi degredatione */
    double fTol_1;    /*  Tolerance for Yield Function */
    double fchi, fc, fphi, fpsi;
};

} // namespace Tahoe 
#endif /* _GRAD_MR_PRIMITIVET_H_ */