/* $Id: TvergHutch2DT.cpp,v 1.20 2004-06-17 07:13:28 paklein Exp $ */
/* created: paklein (02/05/2000) */

#include "TvergHutch2DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"
#include "ifstreamT.h"
#include "StringT.h"

using namespace Tahoe;

/* class parameters */
const int knumDOF = 2;

/* constructor */
TvergHutch2DT::TvergHutch2DT(ifstreamT& in): SurfacePotentialT(knumDOF)
{
	/* traction potential parameters */
	in >> fsigma_max; if (fsigma_max < 0) throw ExceptionT::kBadInputValue;
	in >> fd_c_n; if (fd_c_n < 0) throw ExceptionT::kBadInputValue;
	in >> fd_c_t; if (fd_c_t < 0) throw ExceptionT::kBadInputValue;
	
	/* non-dimensional opening parameters */
	in >> fL_1; if (fL_1 < 0 || fL_1 > 1) throw ExceptionT::kBadInputValue;
	in >> fL_2; if (fL_2 < fL_1 || fL_2 > 1) throw ExceptionT::kBadInputValue;
	in >> fL_fail; if (fL_fail < 1.0) fL_fail = 1.0;

	/* stiffness multiplier */
	in >> fpenalty; if (fpenalty < 0) throw ExceptionT::kBadInputValue;

	/* penetration stiffness */
	fK = fpenalty*fsigma_max/(fL_1*fd_c_n);
}

/* surface potential */
double TvergHutch2DT::FractureEnergy(const ArrayT<double>& state)
{
#pragma unused(state)

	return fd_c_n*fsigma_max*0.5*(1 - fL_1 + fL_2); 
}

double TvergHutch2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
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
const dArrayT& TvergHutch2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(state)
#pragma unused(sigma)
#pragma unused(qIntegrate)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
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
	else if (L < 1.)
		sigbyL = fsigma_max*(1. - L)/(1. - fL_2)/L;
	else
		sigbyL = 0.0;	

	fTraction[0] = sigbyL*r_t*(fd_c_n/fd_c_t);
	fTraction[1] = sigbyL*r_n;

	/* penetration */
	if (u_n < 0) fTraction[1] += fK*u_n;

	return fTraction;

}

/* potential stiffness */
const dMatrixT& TvergHutch2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(state)
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif
	
	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double dtm2 = 1./fd_c_t/fd_c_t;
	double dnm2 = 1./fd_c_n/fd_c_n;
	double L = sqrt(u_t*u_t*dtm2 + u_n*u_n*dnm2);
	
	if (L < fL_1) // K1
	{
		fStiffness[0] = (fd_c_n/fd_c_t)*fsigma_max/(fL_1*fd_c_t);
		fStiffness[1] = 0.0;
		fStiffness[2] = 0.0;
		fStiffness[3] = fsigma_max/(fL_1*fd_c_n);
	}
	else 
	{
		double lt_0 = u_t*dtm2;
		double lt_2 = u_n*dnm2;
		
		if (L < fL_2) // K2
		{
			double dijTerm = fsigma_max/L*fd_c_n;
			
			fStiffness[0] = dijTerm*dtm2;
			fStiffness[3] = dijTerm*dnm2;
			dijTerm /= -L*L;
			fStiffness[0] += dijTerm*lt_0*lt_0;
			fStiffness[2] = fStiffness[1] = dijTerm*lt_0*lt_2;
			fStiffness[3] += dijTerm*lt_2*lt_2;
		}
		else 
			if (L < 1.) // K3
			{
				double dijTerm = fsigma_max*(1./L-1.)/(1.-fL_2)*fd_c_n;
			
				fStiffness[0] = dijTerm*dtm2;
				fStiffness[3] = dijTerm*dnm2;
				dijTerm = -fsigma_max/(1.-fL_2)*fd_c_n/L/L/L;
				fStiffness[0] += dijTerm*lt_0*lt_0;
				fStiffness[3] += dijTerm*lt_2*lt_2;
				fStiffness[1] = fStiffness[2] = dijTerm*lt_0*lt_2;
			}
			else
			{
				fStiffness[0] = 0.0;
				fStiffness[1] = 0.0;
				fStiffness[2] = 0.0;
				fStiffness[3] = 0.0;	
			}
	}

	/* penetration */
	if (u_n < 0) fStiffness[3] += fK;
	
	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT TvergHutch2DT::Status(const dArrayT& jump_u, 
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

void TvergHutch2DT::PrintName(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	out << "    Tvergaard-Hutchinson 2D\n";
#endif
}

/* print parameters to the output stream */
void TvergHutch2DT::Print(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	out << " Cohesive stress . . . . . . . . . . . . . . . . = " << fsigma_max << '\n';
	out << " Normal opening to failure . . . . . . . . . . . = " << fd_c_n     << '\n';
	out << " Tangential opening to failure . . . . . . . . . = " << fd_c_t     << '\n';
	out << " Non-dimensional opening to peak traction. . . . = " << fL_1       << '\n';
	out << " Non-dimensional opening to declining traction . = " << fL_2       << '\n';
	out << " Non-dimensional opening to failure. . . . . . . = " << fL_fail    << '\n';
	out << " Penetration stiffness multiplier. . . . . . . . = " << fpenalty   << '\n';
#endif
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int TvergHutch2DT::NumOutputVariables(void) const { return 1; }
void TvergHutch2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(1);
	labels[0] = "lambda";
}

void TvergHutch2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
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
