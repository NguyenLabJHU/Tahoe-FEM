/* $Id: TvergHutchRigid2DT.cpp,v 1.3 2004-07-15 08:26:02 paklein Exp $ */
#include "TvergHutchRigid2DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"

#include "StringT.h"

using namespace Tahoe;

/* class parameters */
const int knumDOF = 2;

/* constructor */
TvergHutchRigid2DT::TvergHutchRigid2DT(ifstreamT& in): 
	SurfacePotentialT(knumDOF)
{
ExceptionT::GeneralFail("TvergHutchRigid2DT::TvergHutchRigid2DT", "out of date");
#pragma unused(in)
#if 0
	/* traction potential parameters */
	in >> fsigma_max; if (fsigma_max < 0) throw ExceptionT::kBadInputValue;
	in >> fd_c_n; if (fd_c_n < 0) throw ExceptionT::kBadInputValue;
	in >> fd_c_t; if (fd_c_t < 0) throw ExceptionT::kBadInputValue;
	
	/* non-dimensional opening parameters */
	in >> fL_1; if (fL_1 < 0 || fL_1 > 1) throw ExceptionT::kBadInputValue;
	in >> fL_2; if (fL_2 < 0 || fL_2 > 1) throw ExceptionT::kBadInputValue;
	in >> fT_2; if (fT_2 < 0.0) throw ExceptionT::kBadInputValue;
	in >> fL_fail; if (fL_fail < 1.0) fL_fail = 1.0;

	/* stiffness multiplier */
	in >> fpenalty; if (fpenalty < 0) throw ExceptionT::kBadInputValue;

	/* penetration stiffness */
	fK = fpenalty*fsigma_max/(fL_1*fd_c_n);
#endif
}

/* initialize the state variable array */
void TvergHutchRigid2DT::InitStateVariables(ArrayT<double>& state)
{
	/* scale initial traction to the cohesive stress */
	dArrayT tmp;
	tmp.Alias(state);
	tmp.UnitVector();
}

/* surface potential */
double TvergHutchRigid2DT::FractureEnergy(const ArrayT<double>& state)
{
#pragma unused(state)
	return fd_c_n*fsigma_max*0.5*(fT_2 + fL_2);
}

double TvergHutchRigid2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	/* effective opening */
	double u_t = jump_u[0];
	double u_n = jump_u[1];
	double r_t = u_t/fd_c_t;
	double r_n = u_n/fd_c_n;
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)

	double u;
	if (L < fL_2)
		u = fsigma_max*(L + (fT_2 - 1.0)/(2.0*fL_2)*L*L);
	else if (L < 1)
	{
		/* first branch */
		u = 0.5*fsigma_max*(1.0 + fT_2)*fL_2;
		
		/* second branch */		
		u += fsigma_max*fT_2*(L*(1.0 - 0.5*L) + fL_2*(1.0 - 0.5*fL_2))/(1 - fL_2);
	}
	else
		u = fsigma_max*0.5*(fT_2 + fL_2);

	/* penetration */
	if (u_n < 0) u += 0.5*u_n*u_n*fK/fd_c_n;

	return u*fd_c_n; // (1.2)
}
	
/* traction vector given displacement jump vector */	
const dArrayT& TvergHutchRigid2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(state)
#pragma unused(sigma)
#pragma unused(qIntegrate)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	/* effective opening */
	double u_t = jump_u[0];
	double u_n = jump_u[1];
	double r_t = u_t/fd_c_t;
	double r_n = u_n/fd_c_n;
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)

	/* no opening */
	if (fabs(L) < kSmall) {
		fTraction[0] = state[0];
		fTraction[1] = state[1];
	}
	else /* opening */
	{
		double sigbyL;
		if (L < fL_2)
			sigbyL = fsigma_max*(1.0/L + (fT_2 - 1.0)/fL_2);
		else if (L < 1)
			sigbyL = fsigma_max*fT_2*(1.0 - L)/(1.0 - fL_2)/L;
		else
			sigbyL = 0.0;	
	
		fTraction[0] = sigbyL*r_t*(fd_c_n/fd_c_t);
		fTraction[1] = sigbyL*r_n;

		/* penetration */
		if (u_n < 0) fTraction[1] += fK*u_n;
	}

	return fTraction;

}

/* potential stiffness */
const dMatrixT& TvergHutchRigid2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(state)
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

//TEMP
ExceptionT::GeneralFail("TvergHutchRigid2DT::Stiffness", "not implemented");
//see TvergHutchRigid2DT::Stiffness for a start

	fStiffness = 0.0;
	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT TvergHutchRigid2DT::Status(const dArrayT& jump_u, 
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
void TvergHutchRigid2DT::Print(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	out << " Cohesive stress . . . . . . . . . . . . . . . . = " << fsigma_max << '\n';
	out << " Normal opening to failure . . . . . . . . . . . = " << fd_c_n     << '\n';
	out << " Tangential opening to failure . . . . . . . . . = " << fd_c_t     << '\n';
	out << " Non-dimensional compressive reference opening . = " << fL_1       << '\n';
	out << " Non-dimensional shape factor opening. . . . . . = " << fL_2       << '\n';
	out << " Non-dimensional shape factor traction . . . . . = " << fT_2       << '\n';
	out << " Non-dimensional opening to failure. . . . . . . = " << fL_fail    << '\n';
	out << " Penetration stiffness multiplier. . . . . . . . = " << fpenalty   << '\n';
#endif
}
#endif

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int TvergHutchRigid2DT::NumOutputVariables(void) const { return 1; }
void TvergHutchRigid2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(1);
	labels[0] = "lambda";
}

void TvergHutchRigid2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/fd_c_t;
	double r_n = u_n/fd_c_n;
	output[0]  = sqrt(r_t*r_t + r_n*r_n); // (1.1)
}

bool TvergHutchRigid2DT::CompatibleOutput(const SurfacePotentialT& potential) const
{
#ifdef __NO_RTTI__
#pragma unused(potential)
	return false;
#else
	const TvergHutchRigid2DT* pTH = dynamic_cast<const TvergHutchRigid2DT*>(&potential);
	return pTH != NULL;
#endif
}
