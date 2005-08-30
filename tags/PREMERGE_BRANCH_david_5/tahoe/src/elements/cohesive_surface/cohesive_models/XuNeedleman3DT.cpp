/* $Id: XuNeedleman3DT.cpp,v 1.21 2004-07-15 08:26:03 paklein Exp $ */
/* created: paklein (06/23/1999)*/
#include "XuNeedleman3DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"
using namespace Tahoe;

/* class parameters */
const int    knumDOF = 3;
const double kExpMax = 20;

XuNeedleman3DT::XuNeedleman3DT(dArrayT& params): SurfacePotentialT(knumDOF)
{
	SetName("Xu-Needleman_3D");

	q = params[0];; // phi_t/phi_n
	r = params[1]; // delta_n* /d_n
	if (q < 0.0 || r < 0.0) throw ExceptionT::kBadInputValue;
	
	d_n = params[2]; // characteristic normal opening
	d_t = params[3]; // characteristic tangent opening
	if (d_n < 0.0 || d_t < 0.0) throw ExceptionT::kBadInputValue;
	
	phi_n = params[4]; // mode I work to fracture
	if (phi_n < 0.0) throw ExceptionT::kBadInputValue;

	r_fail = params[5]; // d/d_(n/t) for which surface is considered failed
	if (r_fail < 1.0) throw ExceptionT::kBadInputValue;

	fKratio = params[6]; // stiffening ratio
	if (fKratio < 0.0) throw ExceptionT::kBadInputValue;
	
	fK = fKratio*phi_n/(d_n*d_n);
}

XuNeedleman3DT::XuNeedleman3DT(void): 
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
	SetName("Xu-Needleman_3D");
}

/* surface potential */
double XuNeedleman3DT::FractureEnergy(const ArrayT<double>& state) 
{
#pragma unused(state)
	return phi_n; 
}

double XuNeedleman3DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;	
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11;

	double u_t1 = jump_u[0];
	double u_t2 = jump_u[1];
	double u_n  = jump_u[2];

	z1 = 1./d_n;
	z2 = -d_n;
	z3 = 1./(d_t*d_t);
	z4 = -q;
	z5 = -1. + r;
	z6 = -r;
	z7 = d_n*r;
	z8 = u_t1*u_t1;
	z9 = u_t2*u_t2;
	z10 = -u_n*z1;
	z1 = u_n*z1;
	z11 = 1. + z4;
	z4 = r + z4;
	z5 = 1./z5;
	z2 = z2 + z7;
	z7 = z8 + z9;
	z8 = exp(z10);
	z2 = 1./z2;
	z3 = -z3*z7;
	z1 = 1. + z1 + z6;
	z3 = exp(z3);
	z2 = u_n*z2*z4;
	z1 = z1*z11*z5;
	z2 = q + z2;
	z2 = -z2*z3;
	z1 = z1 + z2;
	z1 = phi_n*z1*z8;
	// phi_n + z1
	
	/* penetration */
	if (u_n < 0.0) z1 += 0.5*u_n*u_n*fK;

	return phi_n + z1;
}

/* traction vector given displacement jump vector */	
const dArrayT& XuNeedleman3DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(state)
#pragma unused(sigma)
#pragma unused(qIntegrate)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;

	double u_t1 = jump_u[0];
	double u_t2 = jump_u[1];
	double u_n  = jump_u[2];
	
	z1 = 1./d_n;
	z2 = 1./(d_t*d_t);
	z3 = r-q;
	z4 = 1./(-1. + r);
	z5 = 1 - r;
	z6 = -z2*(u_t1*u_t1 + u_t2*u_t2);
	z8 = -u_n*z1;
	z9 = -z8;
	z10 = 1. - q;
	
	// limit compressive deformation
	if (z8 > kExpMax)
		ExceptionT::BadJacobianDet("XuNeedleman2DT::Traction", "exp(x): x = %g > kExpMax", z8);

	z7 = exp(z8);
	z11 = z1*z10*z4;
	z12 = z3*z4*z9 + q;
	z5 += z9;
	z9 = exp(z6); //don't limit shear opening
	z5 *= z10*z4;
	z6 += z8; // since (z6 < 0), (z6' < z8) and z8 is checked above
	z6 = exp(z6);
	z8 = -z9;
	z3 *= z1*z4*z8;
	z4 = z12*z8 + z5;
	z2 *= 2.*phi_n*z12*z6;
	z3 += z11;
	z5 = u_t1*z2;
	z2 *= u_t2;
	z6 = phi_n*z7;
	z3 *= z6;
	z1 *= -z4*z6;
	z1 += z3;
	
	fTraction[0] = z5;
	fTraction[1] = z2;
	fTraction[2] = z1;
	
	/* penetration */
	if (u_n < 0.0) fTraction[2] += u_n*fK;

	return fTraction;
}


/* potential stiffness */
const dMatrixT& XuNeedleman3DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#pragma unused(state)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18;

	double u_t1 = jump_u[0];
	double u_t2 = jump_u[1];
	double u_n  = jump_u[2];
	
	z1 = 1./(d_n*d_n);
	z2 = 1./d_n;
	z3 = pow(d_t,-4.);
	z4 = 1./(d_t*d_t);
	z5 = -q;
	z6 = -1. + r;
	z7 = -r;
	z8 = u_t1*u_t1;
	z9 = u_t2*u_t2;
	z10 = -u_n*z2;
	z11 = u_n*z2;
	z12 = 1. + z5;
	z5 = r + z5;
	z6 = 1./z6;
	z13 = z8 + z9;
	z14 = exp(z10);
	z15 = z11*z5*z6;
	z16 = z12*z2*z6;
	z13 = -z13*z4;
	z7 = 1. + z11 + z7;
	z11 = exp(z13);
	z15 = q + z15;
	z10 = z10 + z13;
	z7 = z12*z6*z7;
	z10 = exp(z10);
	z11 = -z11;
	z12 = z11*z2*z5*z6;
	z11 = z11*z15;
	z10 = phi_n*z10;
	z12 = z12 + z16;
	z7 = z11 + z7;
	z11 = u_t1*z10;
	z13 = u_t2*z10;
	z16 = z10*z15;
	z17 = z11*z15;
	z15 = z13*z15;
	z18 = -4.*z16*z3;
	z3 = -4.*u_t2*z17*z3;
	z10 = 2.*z10*z4;
	z10 = z11*z2*z4;
	z11 = 2.*z13*z2*z4*z5*z6;
	z13 = 2.*z16*z4;
	z16 = -2.*z17*z2*z4;
	z4 = -2.*z15*z2*z4;
	z5 = 2.*z10*z5*z6;
	z6 = z18*z8;
	z8 = z18*z9;
	z9 = phi_n*z14;
	z4 = z11 + z4;
	z5 = z16 + z5;
	z6 = z13 + z6;
	z8 = z13 + z8;
	z2 = -2.*z12*z2*z9;
	z1 = z1*z7*z9;

	// {{z6, z3, z5},
	//  {z3, z8, z4},
	//  {z5, z4, z1 + z2}}

	fStiffness[0] = z6;
	fStiffness[1] = z3;
	fStiffness[2] = z5;

	fStiffness[3] = z3;
	fStiffness[4] = z8;
	fStiffness[5] = z4;

	fStiffness[6] = z5;
	fStiffness[7] = z4;
	fStiffness[8] = z1 + z2;

	/* penetration */
	if (u_n < 0.0) fStiffness[8] += fK;

	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT XuNeedleman3DT::Status(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	double u_t1 = jump_u[0];
	double u_t2 = jump_u[1];
	double u_t  = sqrt(u_t2*u_t2 + u_t1*u_t1);
	double u_n  = jump_u[2];
	
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
void XuNeedleman3DT::Print(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	out << " Surface energy ratio (phi_t/phi_n). . . . . . . = " << q       << '\n';
	out << " Critical opening ratio (delta_n* /d_n). . . . . = " << r       << '\n';
	out << " Characteristic normal opening to failure. . . . = " << d_n     << '\n';
	out << " Characteristic tangential opening to failure. . = " << d_t     << '\n';
	out << " Mode I work to fracture (phi_n) . . . . . . . . = " << phi_n   << '\n';
	out << " Penetration stiffness multiplier. . . . . . . . = " << fKratio << '\n';
#else
#pragma unused(out)
#endif
}
#endif

/* describe the parameters  */
void XuNeedleman3DT::DefineParameters(ParameterListT& list) const
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
void XuNeedleman3DT::TakeParameterList(const ParameterListT& list)
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
