/* $Id: LennardJonesPairT.cpp,v 1.10 2003-10-31 20:51:11 paklein Exp $ */
#include "LennardJonesPairT.h"
#include "toolboxConstants.h"
#include <iostream.h>

using namespace Tahoe;

/* initialize static parameters */
double LennardJonesPairT::s_eps = 1.0;
double LennardJonesPairT::s_sigma = 1.0;
double LennardJonesPairT::s_alpha =-1.0;
double LennardJonesPairT::s_phi_rc = 0.0;
double LennardJonesPairT::s_dphi_rc = 0.0;

/* constructor */
LennardJonesPairT::LennardJonesPairT(double mass, double eps, double sigma, double alpha):
	f_eps(eps),
	f_sigma(sigma),
	f_alpha(alpha),
	f_dphi_rc(0.0),
	f_phi_rc(0.0)
{
	SetName("Lennard_Jones");
	SetRange(f_sigma*f_alpha);
	SetMass(mass);
	
	/* evaluate unmodified force at the cut-off */
	if (f_alpha > kSmall) {
		s_eps = f_eps;
		s_sigma = f_sigma;
		s_alpha = -1.0;
		s_phi_rc = 0.0;
		s_dphi_rc = 0.0;
		
		f_phi_rc = Energy(f_sigma*f_alpha, NULL, NULL);
		f_dphi_rc = Force(f_sigma*f_alpha, NULL, NULL);
	}
	else
		f_alpha = 0.0;
}

LennardJonesPairT::LennardJonesPairT(void):
	f_eps(0.0),
	f_sigma(0.0),
	f_alpha(0.0),
	f_dphi_rc(0.0),
	f_phi_rc(0.0)
{
	SetName("Lennard_Jones");
}

/* return a pointer to the energy function */
PairPropertyT::EnergyFunction LennardJonesPairT::getEnergyFunction(void)
{
	/* copy my data to static */
	s_eps = f_eps;
	s_sigma = f_sigma;
	s_alpha = f_alpha;
	s_phi_rc = f_phi_rc;
	s_dphi_rc = f_dphi_rc;

	/* return function pointer */
	return LennardJonesPairT::Energy;
}

PairPropertyT::ForceFunction LennardJonesPairT::getForceFunction(void)
{
	/* copy my data to static */
	s_eps = f_eps;
	s_sigma = f_sigma;
	s_alpha = f_alpha;
	s_phi_rc = f_phi_rc;
	s_dphi_rc = f_dphi_rc;

	/* return function pointer */
	return LennardJonesPairT::Force;
}

PairPropertyT::StiffnessFunction LennardJonesPairT::getStiffnessFunction(void)
{
	/* copy my data to static */
	s_eps = f_eps;
	s_sigma = f_sigma;
	s_alpha = f_alpha;
	s_phi_rc = f_phi_rc;
	s_dphi_rc = f_dphi_rc;

	/* return function pointer */
	return LennardJonesPairT::Stiffness;
}

/* write properties to output */
void LennardJonesPairT::Write(ostream& out) const
{
	/* inherited */
	PairPropertyT::Write(out);
	out << " Energy scaling (epsilon). . . . . . . . . . . . = " << f_eps << '\n';
	out << " Length scaling (sigma). . . . . . . . . . . . . = " << f_sigma << '\n';
	out << " Cut-off parameter (alpha) . . . . . . . . . . . = " << f_alpha << '\n';
}

/* describe the parameters needed by the interface */
void LennardJonesPairT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PairPropertyT::DefineParameters(list);

	ParameterT eps(f_eps, "energy_scaling");
	eps.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(eps);

	ParameterT sigma(f_sigma, "length_scaling");
	sigma.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(sigma);

	ParameterT alpha(f_alpha, "cut_off_distance");
	alpha.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(alpha, ParameterListT::ZeroOrOnce);
}

/* accept parameter list */
void LennardJonesPairT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PairPropertyT::TakeParameterList(list);

	f_eps = list.GetParameter("energy_scaling");
	f_sigma = list.GetParameter("length_scaling");
	
	/* optional cut-off distance */
	const ParameterT* alpha = list.Parameter("cut_off_distance");
	if (alpha)
		f_alpha = *alpha;
	else
		f_alpha = 0.0;
	SetRange(f_sigma*f_alpha);

	/* evaluate unmodified force at the cut-off */
	if (f_alpha > kSmall) {
		s_eps = f_eps;
		s_sigma = f_sigma;
		s_alpha = -1.0;
		s_phi_rc = 0.0;
		s_dphi_rc = 0.0;
		
		f_phi_rc = Energy(f_sigma*f_alpha, NULL, NULL);
		f_dphi_rc = Force(f_sigma*f_alpha, NULL, NULL);
	}
}

/***********************************************************************
 * Private
 ***********************************************************************/

double LennardJonesPairT::Energy(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)
	double r_c = s_sigma*s_alpha;
	if (s_alpha > kSmall && r_ab > r_c)
		return 0.0;
	else
	{
		double r = s_sigma/r_ab;
		double r_6 = r*r*r*r*r*r;
		double r_12 = r_6*r_6;
		double sigma_6 = s_sigma*s_sigma*s_sigma*s_sigma*s_sigma*s_sigma;
		double sigma_12 = sigma_6*sigma_6;
		return 4.0*s_eps*(r_12 - r_6) - s_phi_rc - (r_ab - r_c)*s_dphi_rc;
	}
}

double LennardJonesPairT::Force(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	double r_c = s_sigma*s_alpha;
	if (s_alpha > kSmall && r_ab > r_c)
		return 0.0;
	else
	{
		double r = s_sigma/r_ab;
	
		double r_6 = r*r*r*r*r*r;
		double r_7 = r_6*r;
		double r_13 = r_6*r_7;
	
		return 4.0*s_eps*(-12.0*r_13 + 6.0*r_7)/s_sigma - s_dphi_rc;
	}
}

double LennardJonesPairT::Stiffness(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	double r_c = s_sigma*s_alpha;
	if (s_alpha > kSmall && r_ab > r_c)
		return 0.0;
	else
	{
		double r = s_sigma/r_ab;

		double r_6 = r*r*r*r*r*r;
		double r_8 = r_6*r*r;
		double r_14 = r_6*r_8;
	
		return 4.0*s_eps*(156.0*r_14 - 42.0*r_8)/s_sigma/s_sigma;
	}
}
