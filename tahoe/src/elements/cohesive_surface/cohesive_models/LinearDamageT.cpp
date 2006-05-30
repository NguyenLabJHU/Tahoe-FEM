/* $Id: LinearDamageT.cpp,v 1.19 2006-05-30 23:22:04 tdnguye Exp $ */
/* created: paklein (08/21/2000) */
#include "LinearDamageT.h"

#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"

using namespace Tahoe;

/* map to internal variables */
const int   kMaxOpening = 0;
const int kTrialOpening = 1;
const int kInitTraction = 2;

/* constructor */
LinearDamageT::LinearDamageT(ifstreamT& in, const dArrayT& init_traction):
	SurfacePotentialT(init_traction.Length()),
	fInitTraction(init_traction)
{
ExceptionT::GeneralFail("LinearDamageT::LinearDamageT", "out of date");
#if 0
	/* traction potential parameters */
	in >> fd_c_n; if (fd_c_n < 0.0) throw ExceptionT::kBadInputValue;
	in >> fd_c_t; if (fd_c_t < 0.0) throw ExceptionT::kBadInputValue;
#endif
}

/* describe the parameters  */
void LinearDamageT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SurfacePotentialT::DefineParameters(list);

	ParameterT del_c_n(fd_c_n, "delta_crit_n");
	del_c_n.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(del_c_n);

	ParameterT del_c_t(fd_c_t, "delta_crit_t");
	del_c_t.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(del_c_t);

	ParameterT L_max(fL_max, "max opening ratio");
	L_max.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(L_max);

	ParameterT K(fK, "penetration stiffness");
	K.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(K);
}

/* accept parameter list */
void LinearDamageT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SurfacePotentialT::TakeParameterList(list);

	fd_c_n = list.GetParameter("delta_crit_n");
	fd_c_t = list.GetParameter("delta_crit_t");
	fL_max = list.GetParameter("max opening ratio");
	fK = list.GetParameter("penetration stiffness");
}


/* return the number of state variables */
int LinearDamageT::NumStateVariables(void) const
{
	return 2*fInitTraction.Length() + // initiation traction
	       1 +                      // max opening
  	       1;                       // trial max opening
}

/* initialize the state variable array */
void LinearDamageT::InitStateVariables(ArrayT<double>& state)
{
	/* initialization traction */
	double* ptraction = state.Pointer(kInitTraction);
	for (int i = 0; i < fInitTraction.Length(); i++)
		*ptraction++ = fInitTraction[i];

	state[  kMaxOpening] = 0.0; // max opening
	state[kTrialOpening] = 0.0; // trial opening	
}


/* surface potential */
double LinearDamageT::FractureEnergy(const ArrayT<double>& state)
{
#pragma unused(state)

	return 0.5*fInitTraction.Magnitude()*fd_c_n;
}

double LinearDamageT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(jump_u)
#pragma unused(state)
	
	/* not meaningful to define this quantity */
	return 0.0;
}
	
/* traction vector given displacement jump vector */	
const dArrayT& LinearDamageT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != fTraction.Length()) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	double u_n = jump_u.Last();
	
	/* opening parameter */
	double r_n = (u_n < 0.0) ? 0.0 : u_n/fd_c_n;
	double r_t2;
	if (jump_u.Length() == 2)
		r_t2 = jump_u[0]*jump_u[0]/fd_c_t/fd_c_t;
	else
		r_t2 = (jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1])/fd_c_t/fd_c_t;
		
	double L  = sqrt(r_t2 + r_n*r_n);

	if (L > kSmall)
	{
		/* damage law */
		double L_max = state[kMaxOpening];
		double f;
		if (L > fL_max)
			f = 0.0;     // failed
		else
			f = L/fL_max*(fL_max/L - 1.0); // unloading to origin

		/* traction */
		double* init_traction = state.Pointer(kInitTraction);
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
const dMatrixT& LinearDamageT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != fTraction.Length()) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
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
	double L = sqrt(r_t2 + r_n*r_n);

	if (L > kSmall && (fL_max - L) > kSmall)
	{
		/* damage law */
		const double* init_traction = state.Pointer(kInitTraction);
		if (nsd == 2)
		{
			double k22 = (fL_max/L - 1.0) - (u_n*u_n*fL_max)/(L*L*L);
			double k11 = (fL_max/L - 1.0) - (jump_u[0]*jump_u[0]*fL_max)/(L*L*L);
			double k12 = - u_n*jump_u[0]*fL_max/(L*L*L);
			
			fStiffness(0,0) = init_traction[0]*k11;
			fStiffness(0,1) = init_traction[0]*k12;
			fStiffness(1,0) = init_traction[1]*k12;
			fStiffness(1,1) = init_traction[1]*k22;
		}
		else
		{
			double k33 = (fL_max/L - 1.0) - (u_n*u_n*fL_max)/(L*L*L);
			double k22 = (fL_max/L - 1.0) - (jump_u[1]*jump_u[1]*fL_max)/(L*L*L);
			double k11 = (fL_max/L - 1.0) - (jump_u[0]*jump_u[0]*fL_max)/(L*L*L);
			double k13 = - u_n*jump_u[0]*fL_max/(L*L*L);
			double k23 = - u_n*jump_u[1]*fL_max/(L*L*L);
			double k12 = - jump_u[1]*jump_u[0]*fL_max/(L*L*L);
			
			fStiffness(0,0) = init_traction[0]*k11;
			fStiffness(0,1) = init_traction[0]*k12;
			fStiffness(0,2) = init_traction[0]*k23;
			fStiffness(1,0) = init_traction[1]*k12;
			fStiffness(1,1) = init_traction[1]*k22;
			fStiffness(1,2) = init_traction[1]*k23;
			fStiffness(2,0) = init_traction[2]*k13;
			fStiffness(2,1) = init_traction[2]*k23;
			fStiffness(2,2) = init_traction[2]*k33;
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
SurfacePotentialT::StatusT LinearDamageT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(state)

//	double u_t = jump_u[0];
//	double u_n = jump_u[1];

	int nsd = jump_u.Length();
	double u_n = jump_u.Last();
	
	/* opening parameter */
	double r_n = (u_n < 0.0) ? 0.0 : u_n/fd_c_n; // no damage evolution in compression
	double r_t2;
	if (nsd == 2)
		r_t2 = jump_u[0]*jump_u[0]/fd_c_t/fd_c_t;
	else
		r_t2 = (jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1])/fd_c_t/fd_c_t;
	double L = sqrt(r_t2 + r_n*r_n);

	if (L > fL_max)
		return Failed;
	else if (L > 0.0)
		return Critical;
	else
		return Precritical;
}
