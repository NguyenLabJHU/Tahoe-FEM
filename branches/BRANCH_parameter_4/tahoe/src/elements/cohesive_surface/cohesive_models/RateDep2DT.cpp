/* $Id: RateDep2DT.cpp,v 1.16.2.3 2004-07-13 16:42:27 paklein Exp $  */
/* created: cjkimme (10/23/2001) */
#include "RateDep2DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"
#include "StringT.h"
#include "SecantMethodT.h"

using namespace Tahoe;

/* class parameters */
const int knumDOF = 2;
const int kL_1 = 0;
const int kL_2 = 2;
const int ksigma_max = 4;
const int kd_c_n = 6;
const int kd_c_t = 8;
const int kDelta = 10;

/* constructor */
RateDep2DT::RateDep2DT(void): 
	SurfacePotentialT(knumDOF),
	fNeedsStress(false),
	fGroup(-1),
	fsigma_max(0.0),
	fd_c_n(0.0),
	fd_c_t(0.0),
	fL_1(0.0),
	fL_2(0.0),
	fL_fail(0.0),
	fpenalty(0.0),
	fK(0.0),
	fTimeStep(NULL),
	L_2_b(0.0),
	L_2_m(0.0),
	fslope(0.0)
{
	SetName("Tvergaard-Hutchinson_rate_dep_2D");
}

/*initialize state variables with values from the rate-independent model */
void RateDep2DT::InitStateVariables(ArrayT<double>& state)
{
  	int num_state = NumStateVariables();
	if (state.Length() != num_state) 
		ExceptionT::SizeMismatch("RateDep2DT::InitStateVariables", 
			"expecting %d not %d state variables", num_state, state.Length());

	/* clear */
	if (num_state > 0) state = 0.0;

	/* state variables below are with their corresponding parameters 
         * in TvergHutch2DT. state[n] has a flag in state[n+1] that
         * can be used to allow the traction routine to change the
         * parameter to reflect local opening rates. For instance,
         * state[6] represents the critical normal length scale and
         * state[7] = 0. indicates state[6] has not been modified
         * from its initialized (rate-independent) value while
         * state[7] = 1. indicates that the length-scale has
         * been modified as a function of rate. */
	state[kL_1] = fL_1;
	state[kL_2] = fL_2;
	state[ksigma_max] = fsigma_max;
	state[kd_c_n] = fd_c_n;
	state[kd_c_t] = fd_c_t;
	/* state[10] = previous timestep's u_t
	 * state[11] = previous timestep's u_n
	 * state[12] = previous timestep's u_t (after integration)
	 * state[13] = previous timestep's u_n (after integration)
	 */
  	
}

/* return the number of state variables needed by the model */
int RateDep2DT::NumStateVariables(void) const { return 10+2*knumDOF; }

/* surface potential */ 
double RateDep2DT::FractureEnergy(const ArrayT<double>& state) 
{
   	return state[kd_c_n]*state[ksigma_max]*0.5*(1 - state[1] + state[kL_2]); 
}

double RateDep2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/state[kd_c_t];
	double r_n = u_n/state[kd_c_n];
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)

	double u;
	if (L < state[kL_1])
		u = state[ksigma_max]*0.5*(L/state[kL_1])*L;
	else if (L < state[kL_2])
		u = state[ksigma_max]*(0.5*fL_1 + (L - state[kL_1]));
	else if (L < 1)
	{
		double z1 = (1.0 - state[kL_2]);
		double z2 = (1.0 - L);
		u = state[ksigma_max]*(0.5*state[kL_1] + (state[kL_2]- state[kL_1]) + 0.5*(z1 - (z2/z1)*z2));
	}
	else
		u = state[ksigma_max]*0.5*(1 - state[kL_1] + state[kL_2]);

	/* penetration */
	if (u_n < 0) u += 0.5*u_n*u_n*fK/state[kd_c_n];

	return u*state[kd_c_n]; // (1.2)
}
	
/* traction vector given displacement jump vector */	
const dArrayT& RateDep2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
	if (*fTimeStep < 0.0) 
		ExceptionT::BadInputValue("RateDep2DT::Traction", "expecting non-negative time increment %g", fTimeStep);
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/state[kd_c_t];
	double r_n = u_n/state[kd_c_n];
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)

	double sigbyL;
	if (L < state[kL_1]) 
		sigbyL = state[ksigma_max]/state[kL_1];
	else if (L < state[kL_2]) /* we're at or beyond the plateau stress */
	{ 
		if (state[kd_c_n + 1] == 0.) 
		{
			double u_n_dot = 0.0;
			if (fabs(*fTimeStep) > kSmall) 
				u_n_dot = (u_n-state[qIntegrate ? kDelta : kDelta + 2])/(*fTimeStep);
			
			if (u_n_dot > kSmall)
			{
				if (qIntegrate) 
				{
					state[kd_c_n + 1] = 1.;
					state[kd_c_n] = L_2_b + L_2_m * log(u_n_dot);
					/* make sure new length scale is greater than current
		 		 	 * gap vector 
                 	 */
					state[kL_1] *= fd_c_n/state[kd_c_n];
	        	}
	        	r_n = u_n/state[kd_c_n];
	        	L = sqrt(r_t*r_t+r_n*r_n);
	        	sigbyL = state[ksigma_max]*(1+fslope*(L-state[kL_1]))/L;
				if (state[kd_c_n] < u_n || state[kd_c_n] < fd_c_n*fL_1)
				{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
//		  			cout <<  "\n RateDep2DT::Traction: rate-dependent length scale " << state[kd_c_n] << " is incompatible with rate-independent one. Check input parameters. \n ";
#endif	
					if (qIntegrate)
					{	
						state[kd_c_n] = fd_c_n;
						state[kL_2] = state[kL_1]; /* start unloading now */
	          		}
	          		r_n = u_n/state[kd_c_n];
	          		L = sqrt(r_t*r_t+r_n*r_n);
	          		sigbyL = state[ksigma_max]*(1-L)/(1-state[kL_2])/L;
				}
			}
		}
		else 
	    	sigbyL = state[ksigma_max]*(1+fslope*(L-state[kL_1]))/L;
	}
	else if (L < 1)
		sigbyL = state[ksigma_max]*(1+fslope*(state[kL_2]-state[kL_1]))*(1 - L)/(1 - state[kL_2])/L;
	else
		sigbyL = 0.0;	

	fTraction[0] = sigbyL*(u_t/state[kd_c_t])*(state[kd_c_n]/state[kd_c_t]);
	fTraction[1] = sigbyL*(u_n/state[kd_c_n]);

	/* penetration */
	if (u_n < 0) fTraction[1] += fK*u_n;

	if (qIntegrate)
	{
		state[kDelta + 2] = state[kDelta];
		state[kDelta + 3] = state[kDelta + 1];
		state[kDelta] = u_t;
		state[kDelta+1] = u_n;
	}

	return fTraction;
	
}

/* potential stiffness */
const dMatrixT& RateDep2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif	

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double dtm2 = 1./state[kd_c_t]/state[kd_c_t];
	double dnm2 = 1./state[kd_c_n]/state[kd_c_n];
	double L = sqrt(u_t*u_t*dtm2 + u_n*u_n*dnm2);
	
	if (L < state[kL_1]) // K1
	{
		fStiffness[0] = (state[kd_c_n]/state[kd_c_t])*state[ksigma_max]/(state[kL_1]*state[kd_c_t]);
		fStiffness[1] = 0.0;
		fStiffness[2] = 0.0;
		fStiffness[3] = state[ksigma_max]/(state[kL_1]*state[kd_c_n]);
	}
	else 
	{
		double lt_0 = u_t*dtm2;
		double lt_2 = u_n*dnm2;
		
		if (L < state[kL_2]) // K2
		{
			double dijTerm = state[ksigma_max]/L*state[kd_c_n]*(1+fslope*(L-state[kL_1]));
			
			fStiffness[0] = dijTerm*dtm2;
			fStiffness[3] = dijTerm*dnm2;
			
			dijTerm = -state[ksigma_max]/L/L/L*state[kd_c_n]*(1-fslope*state[kL_1]);
			
			fStiffness[0] += dijTerm*lt_0*lt_0;
			fStiffness[2] = fStiffness[1] = dijTerm*lt_0*lt_2;
			fStiffness[3] += dijTerm*lt_2*lt_2;
		}
		else 
			if (L < 1 && state[kL_2] < 1.) // K3
			{
				double dijTerm = state[ksigma_max]*(1.+fslope*(state[kL_2]-state[kL_1]))*(1./L-1.)/(1.-state[kL_2])*state[kd_c_n];
				
				fStiffness[0] = dijTerm*dtm2;
				fStiffness[3] = dijTerm*dnm2;
				
				dijTerm = -state[ksigma_max]*(1.+fslope*(state[kL_2]-state[kL_1]))/(1.-state[kL_2])*state[kd_c_n]/L/L/L;
				
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
SurfacePotentialT::StatusT RateDep2DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif
       
	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/state[kd_c_t];
	double r_n = u_n/state[kd_c_n];
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)
	
	if (L > fL_fail)
		return Failed;
	else if (L > state[kL_1])
		return Critical;
	else
		return Precritical;

}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int RateDep2DT::NumOutputVariables(void) const { return 4; }

void RateDep2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(4);
	labels[0] = "lambda";
	labels[1] = "D_t_dot";
	labels[2] = "D_n_dot";
  	labels[3] = "fd_c_n";
}

void RateDep2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif	
	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/state[kd_c_t];
	double r_n = u_n/state[kd_c_n];
	output[0]  = sqrt(r_t*r_t + r_n*r_n); // (1.1)
	if (fabs(*fTimeStep) > kSmall)
	{
		output[1] = (u_t-state[kDelta + 2])/(*fTimeStep);
		output[2] = (u_n-state[kDelta + 3])/(*fTimeStep);
	}
	else
	{
		output[1] = 0.0;
		output[2] = 0.0;
	}
	output[3] = state[kd_c_n];
}

/*
double RateDep2DT::ComputeNodalValue(const dArrayT& nodalRow) 
{
        return (nodalRow[0]+nodalRow[1])/3;
}

void RateDep2DT::UpdateStateVariables(const dArrayT& IPdata, ArrayT<double>& state)
{
  //state[7] = IPdata[0];
}
*/

bool RateDep2DT::CompatibleOutput(const SurfacePotentialT& potential) const
{
#ifdef __NO_RTTI__
#pragma unused(potential)
	return false;
#else
	const RateDep2DT* pTH = dynamic_cast<const RateDep2DT*>(&potential);
	return pTH != NULL;
#endif
}

/* describe the parameters  */
void RateDep2DT::DefineParameters(ParameterListT& list) const
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

	ParameterT L_fail(fL_fail, "L_fail");
	L_fail.AddLimit(1.0, LimitT::LowerInclusive);
	list.AddParameter(L_fail);

	ParameterT penalty(fpenalty, "penalty");
	penalty.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(penalty);

	/* rate-dependent parameters */
	list.AddParameter(L_2_b, "L_2_b");
	list.AddParameter(L_2_m, "L_2_m");
	list.AddParameter(fslope, "fslope");	
}

/* accept parameter list */
void RateDep2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SurfacePotentialT::TakeParameterList(list);

	fsigma_max = list.GetParameter("sigma_max");
	fd_c_n = list.GetParameter("d_c_n");
	fd_c_t = list.GetParameter("d_c_t");

	fL_1 = list.GetParameter("L_1");
	fL_2 = list.GetParameter("L_2");
	if (fL_2 < fL_1) ExceptionT::BadInputValue("TvergHutch2DT::TakeParameterList",
		"L2 < L1: %g < %g", fL_2, fL_1);

	fL_fail = list.GetParameter("L_fail");
	fpenalty = list.GetParameter("penalty");

	/* penetration stiffness */
	fK = fpenalty*fsigma_max/(fL_1*fd_c_n);

	/* rate-dependent parameters */
	L_2_b = list.GetParameter("L_2_b");
	L_2_m = list.GetParameter("L_2_m");
	fslope = list.GetParameter("fslope");	
}
