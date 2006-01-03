/* $Id: MRPrimitiveT.cpp,v 1.7 2006-01-03 22:18:34 kyonten Exp $ */
/* created: Majid T. Manzari (04/16/2003)                */

/* Base class for a nonassociative, small strain,        */
/* pressure dependent plasticity model with nonlinear    */
/* isotropic hardening/softening.                        */
/* The model is consistent with the traction sepration   */
/* cohesive surface models MR2DT and MR_RP2D.            */

#include "MRPrimitiveT.h"

#include "dSymMatrixT.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
MRPrimitiveT::MRPrimitiveT(void):
	ParameterInterfaceT("MP_primitive")
{

}

/* destructor */
MRPrimitiveT::~MRPrimitiveT(void) { }

/* describe the parameters needed by the interface */
void MRPrimitiveT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT Gf_I(fGf_I, "Gf_I");
	Gf_I.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Gf_I);
	
	ParameterT Gf_II(fGf_II, "Gf_II");
	Gf_II.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Gf_II);
	
	ParameterT chi_p(fchi_p, "chi_p");
	chi_p.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(chi_p);
	
	ParameterT chi_r(fchi_r, "chi_r");
	chi_r.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(chi_r);
	
	ParameterT c_p(fc_p, "c_p");
	c_p.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(c_p);
	
	ParameterT c_r(fc_r, "c_r");
	c_r.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(c_r);
	
	ParameterT phi_p(fphi_p, "phi_p");
	phi_p.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(phi_p);
	
	ParameterT phi_r(fphi_r, "phi_r");
	phi_r.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(phi_r);
	
	ParameterT psi_p(fpsi_p, "psi_p");
	psi_p.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(psi_p);
	
	ParameterT alpha_chi(falpha_chi, "alpha_chi");
	alpha_chi.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_chi);
	
	ParameterT alpha_c(falpha_c, "alpha_c");
	alpha_c.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_c);
	
	ParameterT alpha_phi(falpha_phi, "alpha_phi");
	alpha_phi.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_phi);
	
	ParameterT alpha_psi(falpha_psi, "alpha_psi");
	alpha_psi.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_psi);
	
	ParameterT Tol_1(fTol_1, "Tol_1");
	Tol_1.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Tol_1);
	
	ParameterT Tol_2(fTol_2, "Tol_2");
	Tol_2.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Tol_2);
}

/* accept parameter list */
void MRPrimitiveT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	fGf_I = list.GetParameter("Gf_I");
	fGf_II = list.GetParameter("Gf_II");
	fc_p = list.GetParameter("c_p");
	fc_r = list.GetParameter("c_r");
	fchi_p = list.GetParameter("chi_p");
	fchi_r = list.GetParameter("chi_r");
	fphi_p = list.GetParameter("phi_p");
	fphi_r = list.GetParameter("phi_r");
	fpsi_p = list.GetParameter("psi_p");
	falpha_chi = list.GetParameter("alpha_chi");
	falpha_c = list.GetParameter("alpha_c");
	falpha_phi = list.GetParameter("alpha_phi");
	falpha_psi = list.GetParameter("alpha_psi");
	fTol_1 = list.GetParameter("Tol_1");
	fTol_2 = list.GetParameter("Tol_2");
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/*
 * Returns the value of the yield function given the
 * stress vector and state variables
 */
double MRPrimitiveT::YieldCondition(const dSymMatrixT& devstress, 
			const double meanstress) const
{
  double fpress  = meanstress;
  double enp  = 0.;
  double esp  = 0.;
  double fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
  double fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
  double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
  double ffriction = ftan_phi;
  double temp  = (devstress.ScalarProduct())/2.0;
  temp += pow((fc - ffriction*fchi), 2);
  double ff = sqrt(temp); 
  ff -= (fc - ffriction*fpress);
  return  ff;
}