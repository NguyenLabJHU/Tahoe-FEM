/* $Id: ViscTvergHutch2DT.cpp,v 1.15.2.3 2004-07-13 16:42:27 paklein Exp $ */
/* created: paklein (02/05/2000) */
#include "ViscTvergHutch2DT.h"

#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"
#include "StringT.h"

using namespace Tahoe;

/* class parameters */
const int knumDOF = 2;

/* constructor */
ViscTvergHutch2DT::ViscTvergHutch2DT(void): 
	SurfacePotentialT(knumDOF),
	fTimeStep(NULL),
	fsigma_max(0.0),
	fd_c_n(0.0),
	fd_c_t(0.0),
	fL_1(0.0),
	fL_2(0.0),
	fL_fail(0.0),
	fbeta(0.0),
	feta0(0.0),
	fpenalty(0.0),
	fK(0.0)
{
	SetName("viscous_Tvergaard-Hutchinson_2D");
}

/* return the number of state variables needed by the model */
int ViscTvergHutch2DT::NumStateVariables(void) const { return 2*knumDOF + 1; }

/* incremental heat */
double ViscTvergHutch2DT::IncrementalHeat(const dArrayT& jump, const ArrayT<double>& state)
{
#pragma unused(jump)
	return state[kIncHeat];
}

/* surface potential */
double ViscTvergHutch2DT::FractureEnergy(const ArrayT<double>& state)
{
#pragma unused(state)

	return fd_c_n*fsigma_max*0.5*(1 - fL_1 + fL_2); 
}

double ViscTvergHutch2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/fd_c_t;
	double r_n = u_n/fd_c_n;
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)

	double u;
	if (L < fL_1)
		u = fsigma_max*0.5*(L/fL_1)*L;
	else if (L < fL_2)
		u = fsigma_max*(0.5*fL_1 + (L - fL_1));
	else if (L < 1)
	{
		double z1 = (1.0 - fL_2);
		double z2 = (1.0 - L);
		u = fsigma_max*(0.5*fL_1 + (fL_2 - fL_1) + 0.5*(z1 - (z2/z1)*z2));
	}
	else
		u = fsigma_max*0.5*(1 - fL_1 + fL_2);

	/* penetration */
	if (u_n < 0) u += 0.5*u_n*u_n*fK/fd_c_n;

	return u*fd_c_n; // (1.2)
}
	
/* traction vector given displacement jump vector */	
const dArrayT& ViscTvergHutch2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
	const char caller[] = "ViscTvergHutch2DT::Traction";
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) ExceptionT::SizeMismatch(caller);
	if (state.Length() != NumStateVariables()) ExceptionT::SizeMismatch(caller);
	if (*fTimeStep < 0.0)
		ExceptionT::BadInputValue(caller, "expecting non-negative time increment: %g", fTimeStep);
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/fd_c_t;
	double r_n = u_n/fd_c_n;
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)

	double sigbyL;
	if (L < fL_1)
		sigbyL = fsigma_max/fL_1;
	else if (L < fL_2)
		sigbyL = fsigma_max/L;
	else if (L < 1)
		sigbyL = fsigma_max*(1 - (L - fL_2)/(1 - fL_2))/L;
	else
		sigbyL = 0.0;	

	fTraction[0] = sigbyL*(u_t/fd_c_t)*(fd_c_n/fd_c_t);
	fTraction[1] = sigbyL*(u_n/fd_c_n);

	/* penetration */
	if (u_n < 0) fTraction[1] += fK*u_n;

	/* incremental opening */
	double dd_t = u_t - state[qIntegrate ? 0 : 2];
	double dd_n = u_n - state[qIntegrate ? 1 : 3];

	/* viscous part */
	double T_visc_t = 0.0;
	double T_visc_n = 0.0;
	if (L < 1)
	{
		/* allow dt -> 0 */
		double eta_dt = 0.0;
		if (fabs(*fTimeStep) > kSmall) eta_dt = feta0*(1 - L)/(*fTimeStep);

		T_visc_t = eta_dt*dd_t;
		T_visc_n = eta_dt*dd_n;
	}
	
	/* update state variables */	
	if (qIntegrate)
	{
		/* compute heat generation */
		double d_heat = T_visc_t*dd_t + T_visc_n*dd_n;
		if (dd_n > 0) /* only heat on opening */
			d_heat += fTraction[1]*dd_n;
		if (u_t*dd_t > 0)
			d_heat += fabs(fTraction[0]*dd_t); /* too lazy to figure out the correct sign ;) */
		state[kIncHeat] = 0.9*d_heat; /* work to heat conversion factor */

		/* integrate rest of state variables */
		state[2] = state[0];
		state[3] = state[1];
		state[0] = u_t;
		state[1] = u_n;
	}
	
	/* add viscous stress */
	fTraction[0] += T_visc_t;
	fTraction[1] += T_visc_n;
	return fTraction;
}

/* potential stiffness */
const dMatrixT& ViscTvergHutch2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double z1, z2, z3, z4, z5, z6;	

	/* effective opening */
	z1 = u_n*u_n;
	z2 = 1./(fd_c_n*fd_c_n);
	z3 = u_t*u_t;
	z4 = 1./(fd_c_t*fd_c_t);
	z5 = z1*z2;
	z6 = z3*z4;
	z5 = z5 + z6;
	double L = sqrt(z5);

	if (L < fL_1) // K1
	{
		fStiffness[0] = (fd_c_n/fd_c_t)*fsigma_max/(fL_1*fd_c_t);
		fStiffness[1] = 0.0;
		fStiffness[2] = 0.0;
		fStiffness[3] = fsigma_max/(fL_1*fd_c_n);
	}
	else if (L < fL_2) // K2
	{
		z5 = 1./fd_c_n;
		z2 = z1*z2;
		z3 = z3*z4;
		z2 = z2 + z3;
		z2 = pow(z2,-1.5);
		z2 = fsigma_max*z2*z5;
		z3 = z2*z3;
		z2 = z2*z4;
		z4 = -u_n*u_t*z2;
		z1 = z1*z2;

		//{{z1, z4},
		// {z4, z3}}
		//fStiffness[0] = z1;
		fStiffness[1] = z4;
		fStiffness[2] = z4;
		//fStiffness[3] = z3;

		/* don't allow zero tangent stiffness */
		if (fabs(z1) < kSmall)
			fStiffness[0] = (fd_c_n/fd_c_t)*fsigma_max/(L*fd_c_t); // secant stiffness
		else
			fStiffness[0] = z1;

		if (fabs(z3) < kSmall)
			fStiffness[3] = fsigma_max/(L*fd_c_n); // secant stiffness
		else
			fStiffness[3] = z3;
	}
	else if (L < 1) // K3
	{
		double z7, z8, z9, z10, z11, z12, z13, z14, z15, z16;

		z5 = -1. + fL_2;
		z6 = -fL_2;
		z7 = pow(fd_c_n,-3.);
		z8 = 1./fd_c_n;
		z9 = fd_c_n*fd_c_n;
		z10 = pow(fd_c_n,3.);
		z11 = pow(fd_c_t,-4.);
		z12 = fd_c_t*fd_c_t;
		z2 = z1*z2;
		z13 = z3*z4;
		z6 = 1. + z6;
		z12 = z1*z12;
		z2 = z13 + z2;
		z9 = z3*z9;
		z5 = 1./z5;
		z6 = 1./z6;
		z14 = pow(z2,-1.5);
		z15 = 1./z2;
		z16 = 1./sqrt(z2);
		z2 = sqrt(z2);
		z9 = z12 + z9;
		z12 = -fsigma_max*z1*z15*z6*z7;
		z2 = -z2;
		z9 = 1./z9;
		z2 = fL_2 + z2;
		z5 = fsigma_max*z5*z9;
		z9 = u_n*fd_c_n*u_t*z5;
		z5 = z10*z13*z5;
		z2 = z2*z6;
		z2 = 1. + z2;
		z2 = fsigma_max*z2;
		z1 = -z1*z14*z2*z7;
		z6 = z16*z2*z8;
		z3 = -fd_c_n*z11*z14*z2*z3;
		z7 = -u_n*u_t*z14*z2*z4*z8;
		z2 = fd_c_n*z16*z2*z4;
		z1 = z1 + z12 + z6;
		z4 = z7 + z9;
		z2 = z2 + z3 + z5;
	
		//{{z2, z4},
		// {z4, z1}}
		fStiffness[0] = z2;
		fStiffness[1] = z4;
		fStiffness[2] = z4;
		fStiffness[3] = z1;
	}
	else
	{
		fStiffness[0] = 0.0;
		fStiffness[1] = 0.0;
		fStiffness[2] = 0.0;
		fStiffness[3] = 0.0;	
	}

	/* penetration */
	if (u_n < 0) fStiffness[3] += fK;

	/* viscous part */
	if (L < 1)
	{
		/* well-behaved viscous part */
		double eta_dt = 0.0;
		if (fabs(*fTimeStep) > kSmall) /* allow dt -> 0 */
			eta_dt = feta0*(1 - L)/(*fTimeStep);
		fStiffness.PlusIdentity(eta_dt);

		/* viscous part that needs special treatment near small openings */
		double v0 = 0.0;
		double v1 = 0.0;
		if (fabs(*fTimeStep) > kSmall) { /* allow dt -> 0 */
			v0 = -feta0*(jump_u[0] - state[2])/(*fTimeStep);
			v1 = -feta0*(jump_u[1] - state[3])/(*fTimeStep);
		}
		if (L < kSmall) /* small openings */
		{
			fStiffness[0] += v0/fd_c_t;
			fStiffness[1] += v1/fd_c_t;
			fStiffness[2] += v0/fd_c_n;
			fStiffness[3] += v1/fd_c_n;
		}
		else
		{
			double u0 = jump_u[0]/fd_c_t/fd_c_t/L;
			double u1 = jump_u[1]/fd_c_n/fd_c_n/L;
	
			fStiffness[0] += v0*u0;
			fStiffness[1] += v1*u0;
			fStiffness[2] += v0*u1;
			fStiffness[3] += v1*u1;
		}
	}
	
	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT ViscTvergHutch2DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/fd_c_t;
	double r_n = u_n/fd_c_n;
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)
	
	if (L > fL_fail)
		return Failed;
	else if (L > fL_1)
		return Critical;
	else
		return Precritical;
}

#if 0
/* print parameters to the output stream */
void ViscTvergHutch2DT::Print(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	out << " Cohesive stress . . . . . . . . . . . . . . . . = " << fsigma_max << '\n';
	out << " Normal opening to failure . . . . . . . . . . . = " << fd_c_n     << '\n';
	out << " Tangential opening to failure . . . . . . . . . = " << fd_c_t     << '\n';
	out << " Non-dimensional opening to peak traction. . . . = " << fL_1       << '\n';
	out << " Non-dimensional opening to declining traction . = " << fL_2       << '\n';
	out << " Non-dimensional opening to failure. . . . . . . = " << fL_fail    << '\n';
	out << " Damping parameter . . . . . . . . . . . . . . . = " << feta0      << '\n';
	out << " Penetration stiffness multiplier. . . . . . . . = " << fpenalty   << '\n';
#endif
}
#endif

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int ViscTvergHutch2DT::NumOutputVariables(void) const { return 2; }
void ViscTvergHutch2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(2);
	labels[0] = "lambda";
	labels[1] = "dw_visc";
}

void ViscTvergHutch2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/fd_c_t;
	double r_n = u_n/fd_c_n;
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)
	output[0] = L;
	
	/* approximate incremental viscous dissipation, assuming a constant
	 * average viscosity over the increment */
	if (L < 1)
	{
		/* increment displacement */
		double d_t = u_t - state[2];
		double d_n = u_n - state[3];

		/* previous lambda */
		r_t = state[2]/fd_c_t;
		r_n = state[3]/fd_c_n;
		double L_last = sqrt(r_t*r_t + r_n*r_n); // (1.1)
		
		/* average viscosity */
		double eta = feta0*(1.0 - 0.5*(L + L_last));
		
		/* approximate incremental dissipation */
		if (fabs(*fTimeStep) > kSmall)
			output[1] = 0.5*eta*(d_t*d_t + d_n*d_n)/(*fTimeStep);
		else /* for dt -> 0 */
			output[1] = 0.0;
	}
	else
		output[1] = 0.0;
}

bool ViscTvergHutch2DT::CompatibleOutput(const SurfacePotentialT& potential) const
{
#ifdef __NO_RTTI__
#pragma unused(potential)
	return false;
#else
	const ViscTvergHutch2DT* pTH = dynamic_cast<const ViscTvergHutch2DT*>(&potential);
	return pTH != NULL;
#endif
}

/* describe the parameters  */
void ViscTvergHutch2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SurfacePotentialT::DefineParameters(list);

	ParameterT sigma_max(fsigma_max, "sigma_max");
	sigma_max.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(sigma_max);

	ParameterT d_c_n(fd_c_n, "d_c_n");
	d_c_n.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(d_c_n);

	ParameterT d_c_t(fd_c_t, "d_c_t");
	d_c_t.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(d_c_t);

	ParameterT L_1(fL_1, "L_1");
	L_1.AddLimit(0.0, LimitT::Lower);
	L_1.AddLimit(1.0, LimitT::Upper);
	list.AddParameter(L_1);

	ParameterT L_2(fL_2, "L_2");
	L_2.AddLimit(0.0, LimitT::Lower);
	L_2.AddLimit(1.0, LimitT::Upper);
	list.AddParameter(L_2);

	ParameterT beta(fbeta, "beta");
	beta.AddLimit(0.0, LimitT::LowerInclusive);
	beta.SetDefault(0.9);
	list.AddParameter(beta);

	ParameterT eta_0(feta0, "eta_0");
	eta_0.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(eta_0);

	ParameterT L_fail(fL_fail, "L_fail");
	L_fail.AddLimit(1.0, LimitT::LowerInclusive);
	list.AddParameter(L_fail);

	ParameterT penalty(fpenalty, "penalty");
	penalty.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(penalty);
}

/* accept parameter list */
void ViscTvergHutch2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SurfacePotentialT::TakeParameterList(list);

	fsigma_max = list.GetParameter("sigma_max");
	fd_c_n = list.GetParameter("d_c_n");
	fd_c_t = list.GetParameter("d_c_t");

	fL_1 = list.GetParameter("L_1");
	fL_2 = list.GetParameter("L_2");
	if (fL_2 < fL_1) ExceptionT::BadInputValue("ViscTvergHutch2DT::TakeParameterList",
		"L2 < L1: %g < %g", fL_2, fL_1);

	fbeta = list.GetParameter("beta");
	feta0 = list.GetParameter("eta_0");

	fL_fail = list.GetParameter("L_fail");
	fpenalty = list.GetParameter("penalty");

	/* penetration stiffness */
	fK = fpenalty*fsigma_max/(fL_1*fd_c_n);
}
