/* created: Majid T. Manzari (04/16/2003)                */
/*
/* Base class for a nonassociative, small strain,        */
/* pressure dependent gradient plasticity model          */
/* with nonlinear isotropic hardening/softening.         */
/* The model is consistent with the traction sepration   */
/* cohesive surface models MR2DT and MR_RP2D.            */

#ifndef _GRAD_MR_PRIMITIVET_H_ 
#define _GRAD_MR_PRIMITIVET_H_

/* project headers */
#include "Environment.h"

#include "ios_fwd_decl.h"

namespace Tahoe 
{

/* forward declarations */
class ifstreamT;
class dSymMatrixT;

class GRAD_MRPrimitiveT
{
  public:

	/* constructor */
	GRAD_MRPrimitiveT(ifstreamT& in);

	/* destructor */
	virtual ~GRAD_MRPrimitiveT(void);

	/* write parameters to stream */
   	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	double YieldCondition(const dSymMatrixT& devstress, 
			const double meanstress);

  protected:
	
	double fE;       /* elastic modulus */
	double fnu;      /* Poisson's ratio */
	double fGf_I;    /* Mode_I Fracture Energy */
	double fGf_II;   /* Mode_II Fracture Energy */

	/* Length scale parameters */
	double flse_v; /* pore space length scale (elastic) */
	double flse_s; /* grain size length scale (elastic) */
	double flsp_v; /* pore space length scale (plastic) */
	double flsp_s; /* grain size length scale (plastic) */
	
	/* Inelastic response parameters */
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