/* $Id: LinearDamage2DT.cpp,v 1.1 2006-06-03 16:26:41 tdnguye Exp $ */
/* created: paklein (08/21/2000) */
#include "LinearDamage2DT.h"

#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"

using namespace Tahoe;

/* map to internal variables */
const int   kMaxOpening = 3;
const int kTrialOpening = 2;
const int kInitTraction = 0;

const int knumdof = 2;

/* constructor */
LinearDamage2DT::LinearDamage2DT(void):
	SurfacePotentialT(knumdof),
//	fd_c_n(0.0),
//	fd_c_t(0.0),
	fK(0.0),
	fL_max(1.0)
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
	ParameterT G_max(fG_max, "fracture_energy");
	G_max.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(G_max);

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
	fG_max = list.GetParameter("fracture_energy");
	fK = list.GetParameter("penetration_stiffness");
	
	//TEMP
//	fG_max = 1.67;
}
/* return the number of state variables */
int LinearDamage2DT::NumStateVariables(void) const
{
	return knumdof + // initiation traction
	       1 +       // max opening
  	       1;                       // trial max opening
}

/* initialize the state variable array */
void LinearDamage2DT::InitStateVariables(ArrayT<double>& state)
{
	/* initialization traction */
	state[  kMaxOpening] = 0.0; // max opening
	state[kTrialOpening] = 0.0; // trial opening	

	double* ptraction = state.Pointer(kInitTraction);
	fInitTraction.Set(knumdof, ptraction);
	fInitTraction = 0.0;
}


/* surface potential */
double LinearDamage2DT::FractureEnergy(const ArrayT<double>& state)
{
#pragma unused(state)

//	return 0.5*fInitTraction.Magnitude()*fd_c_n;
	return fG_max;
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

	double u_n = jump_u[1];
	
	/* opening parameter */
	double u_t = jump_u[0];
	double L  = sqrt(u_t*u_t + u_n*u_n);
	
	const double* init_traction = state.Pointer(kInitTraction);
	double sigma_max = sqrt(init_traction[0]*init_traction[0] + init_traction[1]*init_traction[1]);
	double L_max = 2.0*fG_max/sigma_max;
	if (L > kSmall)
	{
		double f;
		if (L > L_max)
			f = 0.0;     // failed
		else {
			f = (L_max/L - 1)/L_max; // unloading to origin
		}
		/* traction */
//		cout <<"\nL "<<L<<"\nsigma_max: "<<sigma_max;		
		fTraction[0] = init_traction[0]*f*u_t;
		fTraction[1] = init_traction[1]*f*u_n;
	}
	else
		fTraction = 0.0;
	
/*	for (int j = 0; j < state.Length(); j++)
		cout << "\nState: "<<state[j]; */

	/* penetration */
	if (u_n < 0) fTraction.Last() += fK*u_n;

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

	double u_n = jump_u[1];
	
	/* opening parameter */
	double u_t = jump_u[0];
	double L  = sqrt(u_t*u_t + u_n*u_n);

	const double* init_traction = state.Pointer(kInitTraction);
	double sigma_max = sqrt(init_traction[0]*init_traction[0] + init_traction[1]*init_traction[1]);
	double L_max = 2.0*fG_max/sigma_max;

	if (L > kSmall)
	{
		if (L > L_max)
			fStiffness = 0.0;     // failed
		else {
//			cout << "\nLHS sigma_max: "<<sigma_max;
			double r1 = (L_max/L - 1.0)/L_max;
			double r2 = 1.0/(L*L*L);

			fStiffness(0,0) = init_traction[0]*(r1 - r2*u_t*u_t);
			fStiffness(0,1) = -init_traction[0]*r2*u_t*u_n;
			fStiffness(1,0) = -init_traction[1]*r2*u_t*u_n;
			fStiffness(1,1) = init_traction[1]*(r1 - r2*u_n*u_n);
		}
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
	const double* init_traction = state.Pointer(kInitTraction);
	double sigma_max = sqrt(init_traction[0]*init_traction[0] + init_traction[1]*init_traction[1]);
	double L_max = 2.0*fG_max/sigma_max;

	int nsd = jump_u.Length();
	double u_n = jump_u.Last();	
	
	/* opening parameter */
	double u_t = jump_u[0];
	double L  = sqrt(u_t*u_t + u_n*u_n);

	if (L > L_max)
		return Failed;
	else if (L > 0.0)
		return Critical;
	else
		return Precritical;
}
