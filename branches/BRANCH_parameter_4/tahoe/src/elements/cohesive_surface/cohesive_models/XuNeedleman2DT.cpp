/* $Id: XuNeedleman2DT.cpp,v 1.16.2.3 2004-07-13 16:42:27 paklein Exp $ */
/* created: paklein (11/14/1997) */
#include "XuNeedleman2DT.h"

#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"

using namespace Tahoe;

/* class parameters */
const int    knumDOF = 2;
const double kExpMax = 100;

/* constructor */
XuNeedleman2DT::XuNeedleman2DT(void): 
	SurfacePotentialT(knumDOF),
	q(0.0),
	r(0.0),
	d_n(0.0),
	d_t(0.0),
	phi_n(0.0),
	r_fail(0.0),
	fKratio(0.0),
	fK(0.0)
{
	SetName("Xu-Needleman_2D");
}

/* surface potential */
double XuNeedleman2DT::FractureEnergy(const ArrayT<double>& state)
{
#pragma unused(state)
	return phi_n; 
}

double XuNeedleman2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8;

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	z1 = 1./d_n;
	z2 = 1./(d_t*d_t);
	z3 = -q;
	z4 = -1. + r;
	z5 = -r;
	z6 = u_t*u_t;
	z7 = -u_n*z1;
	z1 = u_n*z1;
	z8 = 1. + z3;
	z3 = r + z3;
	z4 = 1./z4;
	z2 = -z2*z6;
	z2 = exp(z2);
	z6 = exp(z7);
	z3 = z1*z3*z4;
	z1 = 1. + z1 + z5;
	z3 = q + z3;
	z1 = z1*z4*z8;
	z2 = -z2*z3;
	z1 = z1 + z2;
	z1 = phi_n*z1*z6;
	// phi_n + z1

	/* penetration */
	if (u_n < 0.0) z1 += 0.5*u_n*u_n*fK;

	return phi_n + z1;
}
	
/* traction vector given displacement jump vector */	
const dArrayT& XuNeedleman2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(state)
#pragma unused(sigma)
#pragma unused(qIntegrate)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11;

	double u_t = jump_u[0];
	double u_n = jump_u[1];


	z1 = 1./d_n;
	z2 = 1./(d_t*d_t);
	z3 = -q;
	z4 = -1. + r;
	z5 = -r;
	z6 = u_t*u_t;
	z7 = -u_n*z1;
	z8 = -z7;
	z9 = 1. + z3;
	z3 = r + z3;
	z4 = 1./z4;
	z6 = -z2*z6;
	z10 = exp(z6); //don't limit shear opening
	// limit compressive deformation
	if (z7 > kExpMax)
		ExceptionT::BadJacobianDet("XuNeedleman2DT::Traction", "exp(x): x = %g > kExpMax", z7);
	z11 = exp(z7);
	z6 = z6 + z7; // since (z6 < 0), (z6' < z7) and z7 is checked above
	z7 = z3*z4*z8;
	z5 = 1. + z5 + z8;
	z8 = z1*z4*z9;
	z6 = exp(z6);
	z3 = -z1*z10*z3*z4;
	z7 = q + z7;
	z4 = z4*z5*z9;
	z3 = z3 + z8;
	z5 = -z10*z7;
	z2 = 2.*phi_n*u_t*z2*z6*z7;
	z3 = phi_n*z11*z3;
	z4 = z4 + z5;
	z1 = -phi_n*z1*z11*z4;
	z1 = z1 + z3;
	//z1 = List(z2,z1);

	fTraction[0] = z2;
	fTraction[1] = z1;

	/* penetration */
	if (u_n < 0.0) fTraction[1] += u_n*fK;

	return fTraction;
}

/* potential stiffness */
const dMatrixT& XuNeedleman2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(state)
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16;

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	z1 = 1./(d_n*d_n);
	z2 = 1./d_n;
	z3 = pow(d_t,-4.);
	z4 = 1./(d_t*d_t);
	z5 = -q;
	z6 = -1. + r;
	z7 = -r;
	z8 = u_t*u_t;
	z9 = -u_n*z2;
	z10 = u_n*z2;
	z11 = 1. + z5;
	z5 = r + z5;
	z6 = 1./z6;
	z12 = -z4*z8;
	z13 = exp(z12);
	z14 = exp(z9);
	z15 = z10*z5*z6;
	z16 = z11*z2*z6;
	z7 = 1. + z10 + z7;
	z9 = z12 + z9;
	z9 = exp(z9);
	z10 = q + z15;
	z7 = z11*z6*z7;
	z11 = -z13*z2*z5*z6;
	z12 = -z10*z13;
	z11 = z11 + z16;
	z5 = 2.*phi_n*u_t*z2*z4*z5*z6*z9;
	z6 = 2.*phi_n*z10*z4*z9;
	z4 = -2.*phi_n*u_t*z10*z2*z4*z9;
	z3 = -4.*phi_n*z10*z3*z8*z9;
	z7 = z12 + z7;
	z2 = -2.*phi_n*z11*z14*z2;
	z4 = z4 + z5;
	z3 = z3 + z6;
	z1 = phi_n*z1*z14*z7;

	// {{z3, z4}, {z4, z1 + z2}}

	fStiffness[0] = z3;
	fStiffness[1] = z4;
	fStiffness[2] = z4;
	fStiffness[3] = z1 + z2;

	/* penetration */
	if (u_n < 0.0) fStiffness[3] += fK;

	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT XuNeedleman2DT::Status(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	double u_t1 = jump_u[0];
	double u_t  = sqrt(u_t1*u_t1);
	double u_n  = jump_u[1];
	
	/* square box for now */
	if (u_t > r_fail*d_t || u_n > r_fail*d_n)
		return Failed;
	else if (u_t > d_t || u_n > d_n)
		return Critical;
	else
		return Precritical;
}

#if 0
/* print parameters to the output stream */
void XuNeedleman2DT::Print(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << " Surface energy ratio (phi_t/phi_n). . . . . . . = " << q       << '\n';
	out << " Critical opening ratio (delta_n* /d_n). . . . . = " << r       << '\n';
	out << " Characteristic normal opening to failure. . . . = " << d_n     << '\n';
	out << " Characteristic tangential opening to failure. . = " << d_t     << '\n';
	out << " Mode I work to fracture (phi_n) . . . . . . . . = " << phi_n   << '\n';
	out << " Failure ratio (d_n/delta_n or d_t/delta_t). . . = " << r_fail   << '\n';
	out << " Penetration stiffness multiplier. . . . . . . . = " << fKratio << '\n';
#endif
}
#endif

/* describe the parameters  */
void XuNeedleman2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SurfacePotentialT::DefineParameters(list);

	ParameterT q_(q, "q");
	q_.AddLimit(0.0, LimitT::Lower);
	q_.AddLimit(1.0, LimitT::UpperInclusive);
	q_.SetDefault(1.0);
	list.AddParameter(q_);

	ParameterT r_(r, "r");
	r_.AddLimit(0.0, LimitT::LowerInclusive);
	r_.AddLimit(1.0, LimitT::Upper);
	r_.SetDefault(0.0);
	list.AddParameter(r_);

	ParameterT d_n_(d_n, "d_n");
	d_n_.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(d_n_);

	ParameterT d_t_(d_t, "d_t");
	d_t_.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(d_t_);

	ParameterT phi_n_(d_t, "phi_n");
	phi_n_.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(phi_n_);

	ParameterT r_fail_(r_fail, "r_fail");
	r_fail_.AddLimit(1.0, LimitT::LowerInclusive);
	list.AddParameter(r_fail_);

	ParameterT Kratio(fKratio, "K_ratio");
	Kratio.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Kratio);
}

/* accept parameter list */
void XuNeedleman2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SurfacePotentialT::TakeParameterList(list);

	q = list.GetParameter("q");
	r = list.GetParameter("r");
	d_n = list.GetParameter("d_n");
	d_t = list.GetParameter("d_t");
	phi_n = list.GetParameter("phi_n");
	r_fail = list.GetParameter("r_fail");
	fKratio = list.GetParameter("K_ratio");

	/* penetration stiffness */
	fK = fKratio*phi_n/(d_n*d_n);
}
