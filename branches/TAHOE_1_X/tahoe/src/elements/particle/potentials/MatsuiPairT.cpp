/* $Id: MatsuiPairT.cpp,v 1.2 2003-10-28 23:31:51 paklein Exp $ */

#include "MatsuiPairT.h"
#include "toolboxConstants.h"
#include <iostream.h>
#include <math.h>

using namespace Tahoe;

/* initialize static parameters */
double MatsuiPairT::s_sqr_q = 1.0;
double MatsuiPairT::s_two_A = 1.0;
double MatsuiPairT::s_two_B = 1.0;
double MatsuiPairT::s_sqr_C = 1.0;
double MatsuiPairT::s_f = 1.0;
double MatsuiPairT::s_rc = 1.0;
double MatsuiPairT::s_phi_rc = 0.0;
double MatsuiPairT::s_dphi_rc = 0.0;

/* constructor */
MatsuiPairT::MatsuiPairT(double mass, double sqr_q, double two_A, double two_B, double sqr_C, double f, double rc):
	f_sqr_q(sqr_q),
	f_two_A(two_A),
	f_two_B(two_B),
	f_sqr_C(sqr_C),
	f_f(f),
	f_rc(rc),
	f_dphi_rc(0.0),
	f_phi_rc(0.0)
{
	SetName("Matsui");
	SetRange(f_rc);
	SetMass(mass);
	
	/* evaluate unmodified force at the cut-off */
	if (f_rc > kSmall)
	{
		s_f = f_f;
		s_rc = f_rc;
		s_phi_rc = 0.0;
		s_dphi_rc = 0.0;
		s_sqr_q = f_sqr_q;
		s_two_A = f_two_A;
		s_two_B = f_two_B;
		s_sqr_C = f_sqr_C;

		f_phi_rc = Energy(f_rc, NULL, NULL);
		f_dphi_rc = Force(f_rc, NULL, NULL);
	}
	else
		f_rc = 0.0;
}

MatsuiPairT::MatsuiPairT(void):
	f_sqr_q(0.0),
	f_two_A(0.0),
	f_two_B(0.0),
	f_sqr_C(0.0),
	f_f(0.0),
	f_rc(0.0),
	f_dphi_rc(0.0),
	f_phi_rc(0.0)
{
	SetName("Matsui");
}

/* return a pointer to the energy function */
PairPropertyT::EnergyFunction MatsuiPairT::getEnergyFunction(void)
{
	/* copy my data to static */
	s_sqr_q = f_sqr_q;
	s_two_A = f_two_A;
	s_two_B = f_two_B;
	s_sqr_C = f_sqr_C;
	s_f = f_f;
	s_rc = f_rc;
	s_phi_rc = f_phi_rc;
	s_dphi_rc = f_dphi_rc;

	/* return function pointer */
	return MatsuiPairT::Energy;
}

PairPropertyT::ForceFunction MatsuiPairT::getForceFunction(void)
{
	/* copy my data to static */
	s_sqr_q = f_sqr_q;
	s_two_A = f_two_A;
	s_two_B = f_two_B;
	s_sqr_C = f_sqr_C;
	s_f = f_f;
	s_rc = f_rc;
	s_phi_rc = f_phi_rc;
	s_dphi_rc = f_dphi_rc;

	/* return function pointer */
	return MatsuiPairT::Force;
}

PairPropertyT::StiffnessFunction MatsuiPairT::getStiffnessFunction(void)
{
	/* copy my data to static */
	s_sqr_q = f_sqr_q;
	s_two_A = f_two_A;
	s_two_B = f_two_B;
	s_sqr_C = f_sqr_C;
	s_f = f_f;
	s_rc = f_rc;
	s_phi_rc = f_phi_rc;
	s_dphi_rc = f_dphi_rc;

	/* return function pointer */
	return MatsuiPairT::Stiffness;
}

/* write properties to output */
void MatsuiPairT::Write(ostream& out) const
{
	/* inherited */
	PairPropertyT::Write(out);
	out << " Charge (sqr_q). . . . . . . . . . . . . . . . . = " << f_sqr_q << '\n';
	out << " A1 + A2 (two_A) . . . . . . . . . . . . . . . . = " << f_two_A << '\n';
	out << " B1 + B2 (two_B) . . . . . . . . . . . . . . . . = " << f_two_B << '\n';
	out << " C1 * C2 (sqr_C) . . . . . . . . . . . . . . . . = " << f_sqr_C << '\n';
	out << " Standard Force (f). . . . . . . . . . . . . . . = " << f_f << '\n';
	out << " Cut-off radius (rc) . . . . . . . . . . . . . . = " << f_rc << '\n';
}

/* describe the parameters needed by the interface */
void MatsuiPairT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PairPropertyT::DefineParameters(list);

	list.AddParameter(f_sqr_q, "sqr_q");
	list.AddParameter(f_two_A, "two_A");
	list.AddParameter(f_two_B, "two_B");
	list.AddParameter(f_sqr_C, "sqr_C");
	list.AddParameter(f_f    , "f");
	list.AddParameter(f_rc   , "rc", ParameterListT::ZeroOrOnce);
}

/* accept parameter list */
void MatsuiPairT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PairPropertyT::TakeParameterList(list);

	f_sqr_q = list.GetParameter("sqr_q");
	f_two_A = list.GetParameter("two_A");
	f_two_B = list.GetParameter("two_B");
	f_sqr_C = list.GetParameter("sqr_C");
	f_f     = list.GetParameter("f");

	/* optional cut-off */
	const ParameterT* rc = list.Parameter("rc");
	if (rc)
		f_rc = *rc;
	else
		f_rc = 0.0;
	SetRange(f_rc);
	
	/* evaluate unmodified force at the cut-off */
	if (f_rc > kSmall)
	{
		s_f = f_f;
		s_rc = f_rc;
		s_phi_rc = 0.0;
		s_dphi_rc = 0.0;
		s_sqr_q = f_sqr_q;
		s_two_A = f_two_A;
		s_two_B = f_two_B;
		s_sqr_C = f_sqr_C;

		f_phi_rc = Energy(f_rc, NULL, NULL);
		f_dphi_rc = Force(f_rc, NULL, NULL);
	}
}

/***********************************************************************
 * Private
 ***********************************************************************/

double MatsuiPairT::Energy(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	if (s_rc > kSmall && r_ab > s_rc)
		return 0.0;
	else
	{
	        double r_ab6 = r_ab*r_ab*r_ab*r_ab*r_ab*r_ab;

		return s_sqr_q/r_ab - s_sqr_C/r_ab6 + s_f*s_two_B*exp((s_two_A-r_ab)/s_two_B);
	}
}

double MatsuiPairT::Force(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	if (s_rc > kSmall && r_ab > s_rc)
		return 0.0;
	else
	{
		double r_ab2 = r_ab*r_ab;
		double r_ab7 = r_ab*r_ab*r_ab*r_ab*r_ab*r_ab*r_ab;

		return -s_sqr_q/r_ab2 + 6*s_sqr_C/r_ab7 - s_f*exp((s_two_A-r_ab)/s_two_B);
	}
}

double MatsuiPairT::Stiffness(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	if (s_rc > kSmall && r_ab > s_rc)
		return 0.0;
	else
	{
		double r_ab3 = r_ab*r_ab*r_ab;
		double r_ab8 = r_ab*r_ab*r_ab*r_ab*r_ab*r_ab*r_ab*r_ab;

		return 2*s_sqr_q/r_ab3 - 42*s_sqr_C/r_ab8 + s_f/s_two_B*exp((s_two_A-r_ab)/s_two_B);
	}
}
