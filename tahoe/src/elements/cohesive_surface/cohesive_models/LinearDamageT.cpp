/* $Id: LinearDamageT.cpp,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (08/21/2000)                                          */

#include "LinearDamageT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionCodes.h"
#include "fstreamT.h"

/* map to internal variables */
const int   kMaxOpening = 0;
const int kTrialOpening = 1;
const int kInitTraction = 2;

/* constructor */
LinearDamageT::LinearDamageT(ifstreamT& in, const dArrayT& init_traction,
	iArrayT& i_store, dArrayT& d_store):
	DecohesionLawT(init_traction.Length(), i_store, d_store),
	fInitTraction(init_traction)
{
	/* traction potential parameters */
	in >> fd_c_n; if (fd_c_n < 0.0) throw eBadInputValue;
	in >> fd_c_t; if (fd_c_t < 0.0) throw eBadInputValue;
	
	/* stiffness multiplier */
	in >> fpenalty; if (fpenalty < 0.0) throw eBadInputValue;

	/* penetration stiffness */
	fK = 0.0;
//TEMP: decide on penalty stiffness
}

/* surface potential */
double LinearDamageT::Potential(const dArrayT& jump_u)
{
#pragma unused(jump_u)
	
	/* not meaningful to define this quantity */
	return 0.0;
}
	
/* traction vector given displacement jump vector */	
const dArrayT& LinearDamageT::Traction(const dArrayT& jump_u)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != fTraction.Length()) throw eSizeMismatch;
#endif

	double u_n = jump_u.Last();
	
	/* opening parameter */
	double r_n = (u_n < 0.0) ? 0.0 : u_n/fd_c_n; // no damage evolution in compression
	double r_t2;
	if (jump_u.Length() == 2)
		r_t2 = jump_u[0]*jump_u[0]/fd_c_t/fd_c_t;
	else
		r_t2 = (jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1])/fd_c_t/fd_c_t;
	double L = fd_store[kTrialOpening] = sqrt(r_t2 + r_n*r_n);

	if (L > kSmall)
	{
		/* damage law */
		double L_max = fd_store[kMaxOpening];
		double f;
		if (L > 1.0)
			f = 0.0;     // failed
		else if (L > L_max || L_max < kSmall)
			f = 1.0 - L; // loading
		else
			f = L*(1.0 - L_max)/L_max; // unloading to origin

		/* traction */
		double* init_traction = fd_store.Pointer(kInitTraction);
		fTraction[0] = f*(*init_traction++);
		fTraction[1] = f*(*init_traction++);
		if (jump_u.Length() == 3) fTraction[2] = f*(*init_traction);
	}
	else
		fTraction = 0.0;

	/* penetration */
	if (u_n < 0) fTraction.Last() += fK*u_n;

	return fTraction;
}

/* potential stiffness */
const dMatrixT& LinearDamageT::Stiffness(const dArrayT& jump_u)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != fTraction.Length()) throw eSizeMismatch;
#endif

	int nsd = jump_u.Length();
	double u_n = jump_u.Last();
	
	/* opening parameter */
	double r_n = (u_n < 0.0) ? 0.0 : u_n/fd_c_n; // no damage evolution in compression
	double r_t2;
	if (nsd == 2)
		r_t2 = jump_u[0]*jump_u[0]/fd_c_t/fd_c_t;
	else
		r_t2 = (jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1])/fd_c_t/fd_c_t;
	double L = fd_store[kTrialOpening] = sqrt(r_t2 + r_n*r_n);

	if (L > kSmall)
	{
		/* damage law */
		double L_max = fd_store[kMaxOpening];
		double Df;
		if (L > 1.0)
			Df = 0.0;     // failed
		else if (L > L_max || L_max < kSmall)
			Df = -1.0; // loading
		else
			Df = (1.0 - L_max)/L_max; // unloading
			
		double* init_traction = fd_store.Pointer(kInitTraction);
		double DfbyL = Df/L;
		if (nsd == 2)
		{
			double r_t = DfbyL*jump_u[0]/fd_c_t/fd_c_t;
			double r_n = DfbyL*jump_u[1]/fd_c_n/fd_c_n;
			fStiffness(0,0) = init_traction[0]*r_t;
			fStiffness(0,1) = init_traction[0]*r_n;
			fStiffness(1,0) = init_traction[1]*r_t;
			fStiffness(1,1) = init_traction[1]*r_n;
		}
		else
		{
			double r_t0 = DfbyL*jump_u[0]/fd_c_t/fd_c_t;
			double r_t1 = DfbyL*jump_u[1]/fd_c_t/fd_c_t;
			double r_n  = DfbyL*jump_u[2]/fd_c_n/fd_c_n;
			fStiffness(0,0) = init_traction[0]*r_t0;
			fStiffness(0,1) = init_traction[0]*r_t1;
			fStiffness(0,2) = init_traction[0]*r_n;
			fStiffness(1,0) = init_traction[1]*r_t0;
			fStiffness(1,1) = init_traction[1]*r_t1;
			fStiffness(1,2) = init_traction[1]*r_n;
			fStiffness(2,0) = init_traction[2]*r_t0;
			fStiffness(2,1) = init_traction[2]*r_t1;
			fStiffness(2,2) = init_traction[2]*r_n;
		}
	}
	else
		fStiffness = 0.0;

	/* penetration */
	int i_n = nsd - 1;
	if (u_n < 0) fStiffness(i_n, i_n) += fK;

	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT LinearDamageT::Status(const dArrayT& jump_u)
{
	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/fd_c_t;
	double r_n = u_n/fd_c_n;
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)
	
	if (L > 1.0)
		return Failed;
	else if (L > 0.0)
		return Critical;
	else
		return Precritical;
}

void LinearDamageT::PrintName(ostream& out) const
{
	out << "    Linear Damage\n";
}

/* print parameters to the output stream */
void LinearDamageT::Print(ostream& out) const
{
	out << " Normal opening to failure . . . . . . . . . . . = " << fd_c_n     << '\n';
	out << " Tangential opening to failure . . . . . . . . . = " << fd_c_t     << '\n';
	out << " Penetration stiffness multiplier. . . . . . . . = " << fpenalty   << '\n';
}

/* storage dimensions */
int LinearDamageT::IntegerStorage(void) const
{
	return 0;
}
int LinearDamageT::DoubleStorage(void) const
{
	return fInitTraction.Length() + // initiation traction
	       1 +                      // max opening
	       1;                       // trial max opening
}

/* initialize surface */
void LinearDamageT::InitializeFacet(void) // facet at a time
{
	/* initialization traction */
	double* ptraction = fd_store.Pointer(kInitTraction);
	for (int i = 0; i < fInitTraction.Length(); i++)
		*ptraction++ = fInitTraction[i];

	fd_store[  kMaxOpening] = 0.0; // max opening
	fd_store[kTrialOpening] = 0.0; // trial opening	
}

/* update/reset internal variables */
void LinearDamageT::UpdateHistory(void) // facet at a time
{
	if (fd_store[kTrialOpening] > fd_store[kMaxOpening])
		fd_store[kMaxOpening] = fd_store[kTrialOpening];
}

void LinearDamageT::ResetHistory(void)  // facet at a time
{
	// nothing to do
}
