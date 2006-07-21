/* $Id: LinearDamage2DT.cpp,v 1.2 2006-06-18 01:05:57 tdnguye Exp $ */
/* created: paklein (08/21/2000) */
#include "LinearDamage2DT.h"

#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"

using namespace Tahoe;

/* map to internal variables */
const int kOpening = 2;
const int kTraction = 0;

const int knumdof = 2;

/* constructor */
LinearDamage2DT::LinearDamage2DT(void):
	SurfacePotentialT(knumdof),
//	fd_c_n(0.0),
//	fd_c_t(0.0),
	fK(0.0),
	fdel_max(1.0),
	fsigma_max(0.0)
{
	SetName("Linear_Damage_2D");
}

/* describe the parameters  */
void LinearDamage2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SurfacePotentialT::DefineParameters(list);
/*
	ParameterT del_c_n(fd_c_n, "delta_crit_n");
	del_c_n.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(del_c_n);

	ParameterT del_c_t(fd_c_t, "delta_crit_t");
	del_c_t.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(del_c_t);
*/
	ParameterT sigma_max(fsigma_max, "effective_cohesive_strength");
	sigma_max.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(sigma_max);

	ParameterT del_max(fdel_max, "effective_critical_crack_opening");
	del_max.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(del_max);

	ParameterT K(fK, "penetration_stiffness");
	K.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(K);
}

/* accept parameter list */
void LinearDamage2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SurfacePotentialT::TakeParameterList(list);

//	fd_c_n = list.GetParameter("delta_crit_n");
//	fd_c_t = list.GetParameter("delta_crit_t");
	fdel_max = list.GetParameter("effective_critical_crack_opening");
	fsigma_max = list.GetParameter("effective_cohesive_strength");
	fK = list.GetParameter("penetration_stiffness");
	fTraction.Dimension(knumdof);
}
/* return the number of state variables */
int LinearDamage2DT::NumStateVariables(void) const
{
	return knumdof + // tractions previous time step
		   2*knumdof;  // crack opening previous time step
}

/* initialize the state variable array */
void LinearDamage2DT::InitStateVariables(ArrayT<double>& state)
{
	/* initialization traction */

	for (int i = 0; i<state.Length(); i++)
		state[i] = 0.0;
}


/* surface potential */
double LinearDamage2DT::FractureEnergy(const ArrayT<double>& state)
{
#pragma unused(state)

	return 0.5*fTraction.Magnitude()*fdel_max;
}

double LinearDamage2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(jump_u)
#pragma unused(state)
	
	/* not meaningful to define this quantity */
	return 0.0;
}
	
/* traction vector given displacement jump vector */	
const dArrayT& LinearDamage2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != fTraction.Length()) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif
	double k11=-33.3;
	double k22=-33.3;
	double k12=0.0;
	
	double* p = state.Pointer(kTraction);
	fTraction_Last.Set(knumdof, p);
	
	p = state.Pointer(kOpening);
	fOpening_Last.Set(knumdof, p);
	
	if (!qIntegrate)
	{
		fTraction[0] = fTraction_Last[0];
		fTraction[1] = fTraction_Last[1];
	}
	else
	{
		double u_n = jump_u[1];
		if (u_n < kSmall) u_n = 0.0;
		
		/* opening parameter */
		double u_t = jump_u[0];
		double del  = sqrt(u_t*u_t + u_n*u_n);
		double k0 = fsigma_max/fdel_max;
	
		if (del > kSmall && del < fdel_max)
		{
			double rn2 = u_n*u_n/(del*del);
			double rt2 = u_t*u_t/(del*del);
			double du_n = u_n - fOpening_Last[1];
			double du_t = u_t - fOpening_Last[0];
			
			k22 = -(1.0-fdel_max/del*rt2)*k0;
			k11 = -(1.0 - fdel_max/del*rn2)*k0;
			k12 = -(fdel_max/del * u_t*u_n/(del*del))*k0;


			fTraction[0] = fTraction_Last[0] + k12*du_n + k11*du_t;
			fTraction[1] = fTraction_Last[1] + k12*du_t + k22*du_n;
		}
		else
			fTraction = 0.0;
	
		/* penetration */
		if (u_n < 0) fTraction.Last() += fK*u_n;
//		cout << "\nfTraction_Last: "<<fTraction_Last;
//		cout << "\nfOpening_Last: "<<fOpening_Last;
//		cout << "\nfTraction: "<<fTraction;
		fTraction_Last= fTraction;
		
		state[4] = fOpening_Last[0];
		state[5] = fOpening_Last[1];
		fOpening_Last = jump_u;
	}
	return fTraction;
}

/* potential stiffness */
const dMatrixT& LinearDamage2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != fTraction.Length()) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	double k11=-33.3;
	double k22=-33.3;
	double k12=0.0;

	double u_n = jump_u[1];
	if (u_n < kSmall) u_n = 0.0;
	
	/* opening parameter */
	double u_t = jump_u[0];
	double u_t0 = state[4];
	double u_n0 = state[5];
	
	double del  = sqrt(u_t*u_t + u_n*u_n);
	double k0 = fsigma_max/fdel_max;
	if (del > kSmall && del < fdel_max)
	{
		double rn = u_n/del;
		double rt = u_t/del;
		double rn0 = u_n0/del;
		double rt0 = u_t0/del;
		double rmax = fdel_max/del;
			
		fStiffness(1,1) = -k0*(1.0-rmax*(
			(3.0*rn*rn0+rt*rt0)*rt*rt - 2.0*rn*rn*rt*rt0));
		fStiffness(0,1) = k0*rmax*(
			(-2.0*rt*rn0+rn*rt0)*rn*rn +(-2.0*rn*rt0+rt*rn0)*rt*rt); 
		fStiffness(1,0) = k0*rmax*(
			(-2.0*rt*rn0+rn*rt0)*rn*rn +(-2.0*rn*rt0+rt*rn0)*rt*rt); 
		fStiffness(0,0) = -k0*(1.0-rmax*(
			(3.0*rt*rt0+rn*rn0)*rn*rn - 2.0*rt*rt*rn*rn0));
	}
	else
		fStiffness = 0.0;

	/* penetration */
	if (u_n < 0) fStiffness(1,1) += fK;

	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT LinearDamage2DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
	double u_n = jump_u[1];
	if (u_n < kSmall) u_n = 0.0;
	
	/* opening parameter */
	double u_t = jump_u[0];
	double del  = sqrt(u_t*u_t + u_n*u_n);

	if (del > fdel_max)
		return Failed;
	else if (del > 0.0)
		return Critical;
	else
		return Precritical;
}
