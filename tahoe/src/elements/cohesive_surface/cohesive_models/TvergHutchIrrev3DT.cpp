/* $Id: TvergHutchIrrev3DT.cpp,v 1.1 2007-03-12 02:20:07 cjkimme Exp $ */
/* created: paklein (02/05/2000) */
#include "TvergHutchIrrev3DT.h"

#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"
#include "StringT.h"

using namespace Tahoe;

/* class parameters */
const int knumDOF = 3;

/* clarity for state variables arrays */
const int kIrrev = 0;
const int kLambda = 1;
const int kSigma = 2;
const int kGap = 3;
const int kKill = 5;

/* constructor */
TvergHutchIrrev3DT::TvergHutchIrrev3DT(dArrayT& params): 
	SurfacePotentialT(knumDOF),
	fSecantStiffness(false),
	fTimeStep(NULL)
{
	SetName("Tvergaard-Hutchinson_Irreversible_3D");
	
	/* traction potential parameters */
	fsigma_max = params[0]; if (fsigma_max < 0) throw ExceptionT::kBadInputValue;
	fd_c_n = params[1]; if (fd_c_n < 0) throw ExceptionT::kBadInputValue;
	fd_c_t = params[2]; if (fd_c_t < 0) throw ExceptionT::kBadInputValue;
	/* non-dimensional opening parameters */
	fL_1 = params[3]; if (fL_1 < 0 || fL_1 > 1) throw ExceptionT::kBadInputValue;
	fL_2 = params[4]; if (fL_2 < fL_1 || fL_2 > 1) throw ExceptionT::kBadInputValue;
	fL_fail = params[5]; if (fL_fail < 1.0) fL_fail = 1.0;

	/* stiffness multiplier */
	fpenalty = params[6]; if (fpenalty < 0) throw ExceptionT::kBadInputValue;

	/* penetration stiffness */
	fK = fpenalty*fsigma_max/(fL_1*fd_c_n);
}

TvergHutchIrrev3DT::TvergHutchIrrev3DT(void): 
	SurfacePotentialT(knumDOF),
	fsigma_max(0.0),
	fd_c_n(0.0),
	fd_c_t(0.0),
	fL_1(0.0),
	fL_2(0.0),
	fL_fail(0.0),
	fpenalty(0.0),
	fK(0.0),
	fSecantStiffness(false),
	fTimeStep(NULL)
{
	SetName("Tvergaard-Hutchinson_Irreversible_3D");
}

/* return the number of state variables needed by the model */
int TvergHutchIrrev3DT::NumStateVariables(void) const 
{ 
  // state[0] is flag to swith irrev model, state[1] is \lambda_{peak}, state[2] is T_{peak}, 
  // state[3-5] is gap vector for rate calculation
  // state[6-8] is gap vector for rate calculation
  // state[9] tells us the IP has been killed 
	return kKill + 1; 
}

void TvergHutchIrrev3DT::InitStateVariables(ArrayT<double>& state)
{
  state[kIrrev] = 0.; // IP is initially following TvergHutch cohesive zone model
  state[kLambda] = fL_1; // Initially following TH CZM
  state[kSigma] = fsigma_max; // Initially following TH CZM
  state[kGap] = state[kGap+1] = 0.; // Zero initial displacments
  state[kKill] = 0.; // IP is initially active
}

/* surface potential */
double TvergHutchIrrev3DT::FractureEnergy(const ArrayT<double>& state)
{
#pragma unused(state)

	return fd_c_n*fsigma_max*0.5*(1 - fL_1 + fL_2); 
}

double TvergHutchIrrev3DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	// Just Leave Potential alone for now since its unused code

	double u_t1 = jump_u[0];
	double u_t2 = jump_u[1];
	double u_n = jump_u[2];

	double r_t = 1./fd_c_t/fd_c_t;
	r_t *= u_t1*u_t1 + u_t2*u_t2;
	double r_n = u_n/fd_c_n;
	double L = sqrt(r_t + r_n*r_n); // (1.1)

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
const dArrayT& TvergHutchIrrev3DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(state)
#pragma unused(sigma)
#pragma unused(qIntegrate)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	double u_t1 = jump_u[0];
	double u_t2 = jump_u[1];
	double u_n = jump_u[2];
	
	double r_t1 = u_t1/fd_c_t;
	double r_t2 = u_t2/fd_c_t;
	double r_n = u_n/fd_c_n;
	double L = sqrt(r_t1*r_t1 + r_t2*r_t2 + r_n*r_n); // (1.1)

	double delta_dot = 0.0;
	double L_old = state[qIntegrate ? kGap : kGap + 1];
	if (fabs(*fTimeStep) > kSmall) 
	  delta_dot = (L - L_old)/(*fTimeStep);
	else
	  delta_dot = 0.;
	  //ExceptionT::GeneralFail("TvergHutchIrrev3DT::Traction::Time Step too small :  \n");

	if (delta_dot < kSmall || state[kIrrev] == 1.) { // we are unloading
	  if ( state[kIrrev] == 0. && L > fL_1 ) {
	    if (qIntegrate) { // set state variables for next timestep	      
	      state[kIrrev] = 1.;
	      state[kLambda] = L_old;
	      if (state[kLambda] < fL_2)
		state[kSigma] = fsigma_max;//*state[kLambda]/fL_1; 
	      else 
		if (state[kLambda] < 1.) 
		  state[kSigma] = fsigma_max*(1. - state[kLambda])/(1 - fL_2); // fsigma_max*state[kLambda]/fL_2;
	        else {
		  state[kSigma] = 0.; // model is finished and only interpenetration should be enforced
		  state[kKill] = 1.;
		}
	    }
	  } else {
	    if ( state[kIrrev] == 1. && L > state[kLambda] && qIntegrate ) { // put us back on original TH CZM
	      //cout << " Passed test L = " << L << " ";
	      state[kIrrev] = 0.;
	    }
	  }
	  // cout << "Integrated : L " << L << " kIrrev " << state[kIrrev] << " kSigma " << state[kSigma] << " kLambda " << state[kLambda] << "\n";
	} 

	double sigbyL(0.);
	if (state[kKill] == 0.) {
	  if (L < state[kLambda]) // original test: fL_1)
	    sigbyL = state[kSigma]/state[kLambda]; // original : fsigma_max/fL_1;
	  else {
	      if (L < fL_2)
		sigbyL = fsigma_max/L;
	      else if (L < 1.)
		sigbyL = fsigma_max*(1. - L)/(1. - fL_2)/L;
	      else {
		state[kKill] = 1.0;
		sigbyL = 0.0;	
	      }
	  }
	}

	fTraction[0] = sigbyL*r_t1*(fd_c_n/fd_c_t);
	fTraction[1] = sigbyL*r_t2*(fd_c_n/fd_c_t);
	fTraction[2] = sigbyL*r_n;

	/* penetration */
	if (u_n < 0) fTraction[2] += fK*u_n;

	if (qIntegrate) { // Store previous time step's Lambda
	  state[kGap + 1] = state[kGap];
	  //state[kDelta + 4] = state[kDelta + 1];
	  //state[kDelta + 5] = state[kDelta + 2];
	  state[kGap] = L; //u_t1;
	  //state[kDelta + 1] = u_t2;
	  //state[kDelta + 2] = u_n;
        }

	return fTraction;

}

/* potential stiffness */
const dMatrixT& TvergHutchIrrev3DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(state)
#pragma unused(sigma)
#if __option(extended_errorcheck)
	const char caller[] = "TvergHutchIrrev3DT::Stiffness";
	if (jump_u.Length() != knumDOF) ExceptionT::SizeMismatch(caller);
	if (state.Length() != NumStateVariables()) ExceptionT::GeneralFail(caller);
#endif

	double u_t1 = jump_u[0];
	double u_t2 = jump_u[1];
	double u_n = jump_u[2];
	
	double dtm1 = 1./fd_c_t/fd_c_t;
	double dnm1 = 1./fd_c_n/fd_c_n;
	double L = sqrt((u_t1*u_t1 + u_t2*u_t2)*dtm1 + u_n*u_n*dnm1);

	if (state[kKill] == 1.)
	  fStiffness = 0.;
	else {
	if (fSecantStiffness) /* positive-definite approximation */
	{
		/* slope */
		double sigbyL;
		if (L < state[kLambda]) // orig: fL_1)
		  sigbyL = state[kSigma]/state[kLambda]; // orig: fsigma_max/fL_1;
		else if (L < fL_2)
			sigbyL = fsigma_max/L;
		else if (L < 1.)
			sigbyL = fsigma_max*(1. - L)/(1. - fL_2)/L;
		else
			sigbyL = 0.0;	
	
		/* stiffness */
		fStiffness[0] = fStiffness[4] = fd_c_n/fd_c_t*sigbyL/fd_c_t;
		fStiffness[1] = fStiffness[2] = fStiffness[3] = fStiffness[5] = 0.0;
		fStiffness[6] = fStiffness[7] = 0.0;
		fStiffness[8] = sigbyL/fd_c_n;
	}
	else /* tangent stiffness */
	{
	        if (L < state[kLambda]) // orig: fL_1) // K1
		{
			fStiffness[0] = fStiffness[4] = fd_c_n/fd_c_t*state[kSigma]/(state[kLambda]*fd_c_t);
			fStiffness[1] = fStiffness[2] = fStiffness[3] = fStiffness[5] = 0.0;
			fStiffness[6] = fStiffness[7] = 0.0;
			fStiffness[8] = state[kSigma]/(fL_1*fd_c_n);
		}
		else 
		{
			double lt_0 = jump_u[0]*dtm1;
			double lt_1 = jump_u[1]*dtm1;
			double lt_2 = jump_u[2]*dnm1;
			
			if (L < fL_2) // K2
			{
				double dijTerm = fsigma_max/L*fd_c_n;
			
				fStiffness[0] = fStiffness[4] = dijTerm*dtm1;
				fStiffness[8] = dijTerm*dnm1;
				dijTerm /= -L*L;
				fStiffness[0] += dijTerm*lt_0*lt_0;
				fStiffness[1] = fStiffness[3] = dijTerm*lt_0*lt_1;
				fStiffness[2] = fStiffness[6] = dijTerm*lt_0*lt_2;
				fStiffness[4] += dijTerm*lt_1*lt_1;
				fStiffness[5] = fStiffness[7] = dijTerm*lt_1*lt_2;
				fStiffness[8] += dijTerm*lt_2*lt_2;
			} 
			else 
			{
				if (L < 1.) // K3
				{
					double dijTerm = fsigma_max*(1./L-1.)/(1.-fL_2)*fd_c_n;

					fStiffness[0] = fStiffness[4] = dijTerm*dtm1;
					fStiffness[8] = dijTerm*dnm1;
					dijTerm = -fsigma_max/(1.-fL_2)*fd_c_n/L/L/L;
					fStiffness[0] += dijTerm*lt_0*lt_0;
					fStiffness[4] += dijTerm*lt_1*lt_1;
					fStiffness[8] += dijTerm*lt_2*lt_2;
					fStiffness[1] = fStiffness[3] = dijTerm*lt_0*lt_1;
					fStiffness[2] = fStiffness[6] = dijTerm*lt_0*lt_2;
					fStiffness[5] = fStiffness[7] = dijTerm*lt_1*lt_2;
				}
				else /*Failure*/
				{
					fStiffness[0] = fStiffness[1] = fStiffness[2] = fStiffness[3] = 0.0;
					fStiffness[4] = fStiffness[5] = fStiffness[6] = fStiffness[7] = 0.;
					fStiffness[8] = 0.0;
				}
			}
		}
	}
	}

	/* penetration */
	if (u_n < 0) fStiffness[3] += fK;

	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT TvergHutchIrrev3DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	double r_t = 1./fd_c_t/fd_c_t;
	r_t *= jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1];
	double r_n = jump_u[2]/fd_c_n;
	double L = sqrt(r_t + r_n*r_n); // (1.1)
	
	if (L > fL_fail)
		return Failed;
	else if (L > fL_1)
		return Critical;
	else
		return Precritical;
}

/* describe the parameters  */
void TvergHutchIrrev3DT::DefineParameters(ParameterListT& list) const
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

	ParameterT secant(fSecantStiffness, "secant_stiffness");
	secant.SetDefault(fSecantStiffness);
	list.AddParameter(secant);
}

/* accept parameter list */
void TvergHutchIrrev3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SurfacePotentialT::TakeParameterList(list);

	fsigma_max = list.GetParameter("sigma_max");
	fd_c_n = list.GetParameter("d_c_n");
	fd_c_t = list.GetParameter("d_c_t");

	fL_1 = list.GetParameter("L_1");
	fL_2 = list.GetParameter("L_2");
	if (fL_2 < fL_1) ExceptionT::BadInputValue("TvergHutchIrrev3DT::TakeParameterList",
		"L2 < L1: %g < %g", fL_2, fL_1);

	fL_fail = list.GetParameter("L_fail");
	fpenalty = list.GetParameter("penalty");

	/* penetration stiffness */
	fK = fpenalty*fsigma_max/(fL_1*fd_c_n);

	/* secant stiffness */
	fSecantStiffness = list.GetParameter("secant_stiffness");
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int TvergHutchIrrev3DT::NumOutputVariables(void) const { return 1; }
void TvergHutchIrrev3DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(1);
	labels[0] = "lambda";
}

void TvergHutchIrrev3DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	double r_t = 1./fd_c_t/fd_c_t;
	r_t *= jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1];
	double r_n = jump_u[2]/fd_c_n;
	
	//output[0]  = sqrt(r_t + r_n*r_n); // (1.1)
	output[0] = state[kLambda];
}

bool TvergHutchIrrev3DT::CompatibleOutput(const SurfacePotentialT& potential) const
{
#ifdef __NO_RTTI__
#pragma unused(potential)
	return false;
#else
	const TvergHutchIrrev3DT* pTH = dynamic_cast<const TvergHutchIrrev3DT*>(&potential);
	return pTH != NULL;
#endif
}
