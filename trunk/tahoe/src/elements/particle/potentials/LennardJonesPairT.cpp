/* $Id: LennardJonesPairT.cpp,v 1.2 2002-11-26 01:55:37 paklein Exp $ */
#include "LennardJonesPairT.h"

using namespace Tahoe;

/* initialize static parameters */
double LennardJonesPairT::s_eps = 1.0;
double LennardJonesPairT::s_sigma = 1.0;
double LennardJonesPairT::s_cut = 1.0;

/* constructor */
LennardJonesPairT::LennardJonesPairT(double mass, double eps, double sigma, double cut_off):
	f_eps(eps),
	f_sigma(sigma),
	f_cut(cut_off)
{
	SetRange(f_cut);
	SetMass(mass);
}

/* return a pointer to the energy function */
PairPropertyT::EnergyFunction LennardJonesPairT::getEnergyFunction(void)
{
	/* copy my data to static */
	s_eps = f_eps;
	s_sigma = f_sigma;
	s_cut = f_cut;

	/* return function pointer */
	return LennardJonesPairT::Energy;
}

PairPropertyT::ForceFunction LennardJonesPairT::getForceFunction(void)
{
	/* copy my data to static */
	s_eps = f_eps;
	s_sigma = f_sigma;
	s_cut = f_cut;

	/* return function pointer */
	return LennardJonesPairT::Force;
}

PairPropertyT::StiffnessFunction LennardJonesPairT::getStiffnessFunction(void)
{
	/* copy my data to static */
	s_eps = f_eps;
	s_sigma = f_sigma;
	s_cut = f_cut;

	/* return function pointer */
	return LennardJonesPairT::Stiffness;
}

/***********************************************************************
 * Private
 ***********************************************************************/

double LennardJonesPairT::Energy(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	double r = s_sigma/r_ab;
	
	double r_6 = r*r*r*r*r*r;
	double r_12 = r_6*r_6;
	
	return s_eps*(0.5*r_12 - r_6);
}

double LennardJonesPairT::Force(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	double r = s_sigma/r_ab;
	
	double r_6 = r*r*r*r*r*r;
	double r_7 = r_6*r;
	double r_13 = r_6*r_7;
	
	return 6.0*s_eps*(-r_13 + r_7);
}

double LennardJonesPairT::Stiffness(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	double r = s_sigma/r_ab;

	double r_6 = r*r*r*r*r*r;
	double r_8 = r_6*r*r;
	double r_14 = r_6*r_8;
	
	return s_eps*(78.0*r_14 - 42.0*r_8);
}
