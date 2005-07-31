/* $Id: TvergHutch2DT.cpp,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (02/05/2000)                                          */
/* cohesive potential from Tvergaard and Hutchinson,                      */
/* JMPS v41, n6, 1995, 1119-1135.                                         */

#include "TvergHutch2DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionCodes.h"
#include "fstreamT.h"
#include "StringT.h"

/* class parameters */
const int knumDOF = 2;

/* constructor */
TvergHutch2DT::TvergHutch2DT(ifstreamT& in): SurfacePotentialT(knumDOF)
{
	/* traction potential parameters */
	in >> fsigma_max; if (fsigma_max < 0) throw eBadInputValue;
	in >> fd_c_n; if (fd_c_n < 0) throw eBadInputValue;
	in >> fd_c_t; if (fd_c_t < 0) throw eBadInputValue;
	
	/* non-dimensional opening parameters */
	in >> fL_1; if (fL_1 < 0 || fL_1 > 1) throw eBadInputValue;
	in >> fL_2; if (fL_2 < fL_1 || fL_2 > 1) throw eBadInputValue;
	in >> fL_fail;
	if (fL_fail < 1.0) fL_fail = 1.0;

	/* stiffness multiplier */
	in >> fpenalty; if (fpenalty < 0) throw eBadInputValue;

	/* penetration stiffness */
	fK = fpenalty*fsigma_max/(fL_1*fd_c_n);
}

/* surface potential */
double TvergHutch2DT::Potential(const dArrayT& jump_u)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw eSizeMismatch;
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
const dArrayT& TvergHutch2DT::Traction(const dArrayT& jump_u)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw eSizeMismatch;
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

	return fTraction;
}

/* potential stiffness */
const dMatrixT& TvergHutch2DT::Stiffness(const dArrayT& jump_u)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw eSizeMismatch;
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

	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT TvergHutch2DT::Status(const dArrayT& jump_u)
{
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

void TvergHutch2DT::PrintName(ostream& out) const
{
	out << "    Tvergaard-Hutchinson 2D\n";
}

/* print parameters to the output stream */
void TvergHutch2DT::Print(ostream& out) const
{
	out << " Cohesive stress . . . . . . . . . . . . . . . . = " << fsigma_max << '\n';
	out << " Normal opening to failure . . . . . . . . . . . = " << fd_c_n     << '\n';
	out << " Tangential opening to failure . . . . . . . . . = " << fd_c_t     << '\n';
	out << " Non-dimensional opening to peak traction. . . . = " << fL_1       << '\n';
	out << " Non-dimensional opening to declining traction . = " << fL_2       << '\n';
	out << " Non-dimensional opening to failure. . . . . . . = " << fL_fail    << '\n';
	out << " Penetration stiffness multiplier. . . . . . . . = " << fpenalty   << '\n';
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int TvergHutch2DT::NumOutputVariables(void) const { return 1; }
void TvergHutch2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Allocate(1);
	labels[0] = "lambda";
}

void TvergHutch2DT::ComputeOutput(const dArrayT& jump_u, dArrayT& output)
{
	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/fd_c_t;
	double r_n = u_n/fd_c_n;
	output[0]  = sqrt(r_t*r_t + r_n*r_n); // (1.1)
}

bool TvergHutch2DT::CompatibleOutput(const SurfacePotentialT& potential) const
{
#ifdef __NO_RTTI__
#pragma unused(potential)
	return false;
#else
	const TvergHutch2DT* pTH = dynamic_cast<const TvergHutch2DT*>(&potential);
	return pTH != NULL;
#endif
}
