/* created: Majid T. Manzari (04/16/2003)                */

/* Base class for a nonassociative, small strain,        */
/* pressure dependent gradient plasticity model          */
/* with nonlinear isotropic hardening/softening.         */
/* The model is consistent with the traction sepration   */
/* cohesive surface models MR2DT and MR_RP2D.            */


#include "GRAD_MRPrimitiveT.h"

#include <iostream.h>
#include <math.h>

#include "ifstreamT.h"
#include "dSymMatrixT.h"


using namespace Tahoe;


/* constructor */
GRAD_MRPrimitiveT::GRAD_MRPrimitiveT(ifstreamT& in)

{
	/* read parameters */
	
	in >> fE;     if (fE < 0) throw ExceptionT::kBadInputValue;
	in >> fnu;    if (fnu < 0) throw ExceptionT::kBadInputValue;
	in >> fGf_I;  if (fGf_I < 0) throw ExceptionT::kBadInputValue;
	in >> fGf_II; if (fGf_II < 0) throw ExceptionT::kBadInputValue;
	
	/* length scale parameters */
	in >> flse_v;     if (flse_v < 0) throw ExceptionT::kBadInputValue;
    in >> flse_s;     if (flse_s < 0) throw ExceptionT::kBadInputValue;
    in >> flsp_v;     if (flsp_v < 0) throw ExceptionT::kBadInputValue;
    in >> flsp_s;     if (flsp_s < 0) throw ExceptionT::kBadInputValue;
    
	/* Inelastic Response initiation parameters */
	in >> fchi_p; if (fchi_p < 0) throw ExceptionT::kBadInputValue;
	in >> fchi_r; if (fchi_r < 0) throw ExceptionT::kBadInputValue;
	in >> fc_p; if (fc_p < 0) throw ExceptionT::kBadInputValue;
	in >> fc_r; if (fc_r < 0) throw ExceptionT::kBadInputValue;
    in >> fphi_p; if (fphi_p < 0) throw ExceptionT::kBadInputValue;
	in >> fphi_r; if (fphi_r < 0) throw ExceptionT::kBadInputValue;
	in >> fpsi_p; if (fpsi_p < 0) throw ExceptionT::kBadInputValue;
	in >> falpha_chi; if (falpha_chi < 0) throw ExceptionT::kBadInputValue;
	in >> falpha_c; if (falpha_c < 0) throw ExceptionT::kBadInputValue;
	in >> falpha_phi; if (falpha_phi < 0) throw ExceptionT::kBadInputValue;
	in >> falpha_psi; if (falpha_psi < 0) throw ExceptionT::kBadInputValue;
	in >> fTol_1; if (fTol_1 < 0) throw ExceptionT::kBadInputValue;
}

/* destructor */
GRAD_MRPrimitiveT::~GRAD_MRPrimitiveT(void) { }

/* write parameters */
void GRAD_MRPrimitiveT::Print(ostream& out) const
{
    out << " Elastic tangential stiffness. . . . . . . . . . = " << fE << '\n';
	out << " Elastic Normal stiffness . . .. . . . . . . . . = " << fnu     << '\n';
	out << " Mode_I Fracture Energy            . . . . . . . = " << fGf_I     << '\n';
	out << " Mode_II Fracture Energy            . . .  . . . = " << fGf_II << '\n';
	out << " Pore space length scale (elastic)  . . .  . . . = " << flse_v << '\n';
	out << " Grain size length scale (elastic)  . . .  . . . = " << flse_s << '\n';
	out << " Pore space length scale (plastic)  . . .  . . . = " << flsp_v << '\n';
	out << " Grain size length scale (plastic)  . . .  . . . = " << flsp_s << '\n';
	out << " Peak Cohesion                 . . . . . . . . . = " << fc_p    << '\n';
	out << " Residual Cohesion             . . . . . . . . . = " << fc_r    << '\n';
	out << " Peak Tensile Strength   . . . . . . . . . . . . = " << fchi_p << '\n';
	out << " Residual Tensile Strength . . . . . . . . . . . = " << fchi_r << '\n';
	out << " Peak Friction Angle         . . . . . . . . . . = " << fphi_p   << '\n';
	out << " Critical State Friction Angle       . . . . . . = " << fphi_r  << '\n';
	out << " Peak Dilation Angle.. . . . . . . . . . . . . . = " << fpsi_p   << '\n';
	out << " Coefficient of Tensile Strength Degradation ..  = " << falpha_chi   << '\n';
	out << " Coefficient of Cohesion Degradation. .. . . . . = " << falpha_c   << '\n';
	out << " Coefficient for Frictional Angle Degradation .  = " << falpha_phi   << '\n';
	out << " Coefficient for Dilation Angle Degradation  . . = " << falpha_psi   << '\n';
	out << " Error Tolerance for Yield Function. . . . . . . = " << fTol_1 << '\n';
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void GRAD_MRPrimitiveT::PrintName(ostream& out) const
{
  out << "    Nonassociative Manzari-Regueiro, Pressure-Dependent, Small Strain, \n";
  out << "    Gradient Plasticity Model with Nonlinear Isotropic Hardening/Softening\n";
}

/*
 * Returns the value of the yield function given the
 * stress vector and state variables, where alpha
 * represents isotropic hardening.
 */
double GRAD_MRPrimitiveT::YieldCondition(const dSymMatrixT& devstress, 
			const double meanstress)
{
  double kTemp1, kTemp2, kTemp3, kTemp4;
  double fc, fchi, ffriction, ff, ftan_phi, fpress;

  fpress  = meanstress;
  double enp  = 0.;
  double esp  = 0.;
  fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
  fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
  ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
  ffriction = ftan_phi;
  ff   = (devstress.ScalarProduct())/2.0;
  kTemp2  = (fc - ffriction*fpress);
  kTemp1  = kTemp2;
  kTemp1 *= kTemp2;
  ff  -= kTemp1;
  kTemp3  = (fc - ffriction*fchi);
  kTemp4  = kTemp3;
  kTemp4 *= kTemp3;
  ff  += kTemp4;
  return  ff;
}