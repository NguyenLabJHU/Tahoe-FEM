/* $Id: YoonAllen3DT.cpp,v 1.15 2004-07-15 08:26:03 paklein Exp $ */
#include "YoonAllen3DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"

#include "StringT.h"
#include "YoonAllen2DT.h" /* needed for YoonAllen3DT::NewSub */

using namespace Tahoe;

/* class parameters */
const int knumDOF = 3;

YoonAllen3DT::YoonAllen3DT(dArrayT& fparams, iArrayT& iparams, const double& time_step): 
	SurfacePotentialT(knumDOF),
	fTimeStep(&time_step)
{
	SetName("Yoon-Allen_3D");

	/* traction potential parameters */
	fsigma_0 = fparams[0]; if (fsigma_0 < 0) throw ExceptionT::kBadInputValue;
	fd_c_n = fparams[1]; if (fd_c_n < 0) throw ExceptionT::kBadInputValue;
	fd_c_t = fparams[2]; if (fd_c_t < 0) throw ExceptionT::kBadInputValue;
	
	/* moduli and time constants */
	fE_infty = fparams[3]; if (fE_infty < 0) throw ExceptionT::kBadInputValue;
	int n_prony = iparams[0]; if (n_prony < 0) throw ExceptionT::kBadInputValue;
	
	ftau.Allocate(n_prony);
	fE_t.Allocate(n_prony);
	fexp_tau.Allocate(n_prony);
	
	int i;
	for (i = 0;i < n_prony; i++)
	{
		fE_t[i] = fparams[i+4]; if (fE_t[i] < 0) throw ExceptionT::kBadInputValue;
		ftau[i] = fparams[i+5]; if (ftau[i] < 0) throw ExceptionT::kBadInputValue;
	
		fexp_tau[i] = exp(-(*fTimeStep)/ftau[i]) - 1.;
		
		/* scale ftau by fE_t to reduce multiplications in traction
		 * and stiffness routines
		 */
		ftau[i] *= fE_t[i];
	}
	i += 4;

	idamage = iparams[1]; if (idamage < 1 && idamage > 3) throw ExceptionT::kBadInputValue;
	/* damage evolution law parameters */
	falpha_exp = fparams[i++]; //if (falpha_exp < 1.) throw ExceptionT::kBadInputValue;
	falpha_0 = fparams[i++]; if (falpha_0 <= kSmall) throw ExceptionT::kBadInputValue;
	if (idamage == 2) 
	{
		flambda_exp = fparams[i++]; if (flambda_exp > -1.) throw ExceptionT::kBadInputValue;	
		flambda_0 = fparams[i++]; if (flambda_0 < 1.) throw ExceptionT::kBadInputValue;
	}
	
	/* stiffness multiplier */
	fpenalty = fparams[i++]; if (fpenalty < 0) throw ExceptionT::kBadInputValue;

	
	/* penetration stiffness */
	fK = fpenalty*fsigma_0/fd_c_n;

}

YoonAllen3DT::YoonAllen3DT(void): 
	SurfacePotentialT(knumDOF),
	fTimeStep(NULL),
	
	/* traction potential parameters */
	fsigma_0(0.0),
	fd_c_n(0.0),
	fd_c_t(0.0),
	
	/* moduli and time constants */
	fE_infty(0.0),
	fpenalty(0.0),
	
	/* damage flag */
	idamage(-1),
	fCurrentTimeStep(-1.0)
{
	SetName("Yoon-Allen_3D");
}

/*initialize state variables with values from the rate-independent model */
void YoonAllen3DT::InitStateVariables(ArrayT<double>& state)
{
 	int num_state = NumStateVariables();
	if (state.Length() != num_state) 
		ExceptionT::SizeMismatch("YoonAllen3DT::InitStateVariables", 
			"expecting %d not %d state variables", num_state, state.Length());

	/* clear */
	if (num_state > 0) state = 0.0;

	/* 	The first n_prony slots in state are the previous
	 *  steps hereditary integral (more or less). 0..n_prony-1
	 *	The next kNumDOF slots are the previous step's components of
	 *  the gap vector.   n_prony..n_prony+knumDOF-1
	 *	The next kNumDOF slots are the ith components of the previous
	 *	step's tractions. n_prony+knumDOF..n_prony+2knumDOF-1
	 *  The next slot is the rate of change of the previous step's
	 *  lambda. n_prony+2knumDOF
	 *	The next slot is the previous step's value of alpha, the damage
	 *	coefficient.  n_prony+2knumDOF+1
	 *	The final slot holds the integral of T dot dDelta. n_prony+2knumDOF+2
	 */ 

}

/* return the number of state variables needed by the model */
int YoonAllen3DT::NumStateVariables(void) const 
{ 
	return 4*knumDOF + ftau.Length() + 5; 
}

/* surface potential */ 
double YoonAllen3DT::FractureEnergy(const ArrayT<double>& state) 
{
   	return state[4*knumDOF + ftau.Length() + 4]; 
}

double YoonAllen3DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(jump_u)
#pragma unused(state)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	ExceptionT::GeneralFail("YoonAllen3DT::Potential", "not implemented");
	return 0.;
}
	
/* traction vector given displacement jump vector */	
const dArrayT& YoonAllen3DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
	const char caller[] = "YoonAllen3DT::Traction";
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
	if (*fTimeStep < 0.0) 
		ExceptionT::BadInputValue(caller, "expecting non-negative time increment %g", fTimeStep);
#endif

	int n_prony = ftau.Length();
	if (!qIntegrate)
	{
		fTraction[0] = state[n_prony+knumDOF];
		fTraction[1] = state[n_prony+knumDOF+1];
		fTraction[2] = state[n_prony+knumDOF+2];
		return fTraction;
	}
	else
	{
		/* update decay functions */
		if (fabs(fCurrentTimeStep - *fTimeStep) > kSmall)
		{
			/* decay function */
			for (int i = 0; i < fexp_tau.Length(); i++)
				fexp_tau[i] = exp(-(*fTimeStep)*fE_t[i]/ftau[i]) - 1.0;

			/* step size */
			fCurrentTimeStep = *fTimeStep;
		}	

		double u_t0 = jump_u[0];
		double u_t1 = jump_u[1];
		double u_n = jump_u[2];
		
		double *state2 = state.Pointer(n_prony);

		/* These calculations don't get used. They're useless! */
		/*double u_t0_dot = (u_t0 - state2[0])/fTimeStep;
		double u_t1_dot = (u_t1 - state2[1])/fTimeStep;
		double u_n_dot = (u_n - state2[2])/fTimeStep;
		*/
			
		double l_0 = u_t0/fd_c_t;
		double l_1 = u_t1/fd_c_t;
		double l_2 = u_n/fd_c_n;
		double l = sqrt(l_0*l_0+l_1*l_1+l_2*l_2); // l stands for lambda 
		
		fTraction = 0.;
		
		double l_0_old = state2[0]/fd_c_t;
		double l_1_old = state2[1]/fd_c_t;
		double l_2_old = state2[2]/fd_c_n;
		double l_old = sqrt(l_0_old*l_0_old+l_1_old*l_1_old+l_2_old*l_2_old);
		double prefactold = 1./(1-state2[2*knumDOF+1]);
		double l_dot = 0.0;
		if (fabs(*fTimeStep) > kSmall) l_dot = (l-l_old)/(*fTimeStep); /* for dt -> 0 */

		/* do the bulk of the computation now */
		if (l_old > kSmall)
		{
			prefactold *= l_old;
			if (fabs(l_0_old) > kSmall)
				fTraction[0] = prefactold/l_0_old*state2[knumDOF]; 
			if (fabs(l_1_old) > kSmall)
				fTraction[1] = prefactold/l_1_old*state2[knumDOF+1];  
			if (fabs(l_2_old) > kSmall)
				fTraction[2] = prefactold/l_2_old*state2[knumDOF+2];
		}
		else
		{
			fTraction[0] = prefactold*state2[knumDOF];
			fTraction[1] = prefactold*state2[knumDOF+1];
			fTraction[2] = prefactold*state2[knumDOF+2];
		}	
		
		double tmpSum = fE_infty*(*fTimeStep);
		for (int i = 0; i < n_prony; i++)
		  tmpSum -= fexp_tau[i]*ftau[i];
		tmpSum *= l_dot;
		
		fTraction += tmpSum;
		
		/* update the S coefficients */
		for (int i = 0; i < n_prony; i++)
		{
			state[i] += (state[i]-state2[2*knumDOF]*ftau[i])*fexp_tau[i];
			fTraction += state[i]*fexp_tau[i];	
		}
		
		/* evolve the damage parameter */
		double alpha = state2[2*knumDOF+1];
		if (l_dot > kSmall)
		{
			switch (idamage) 
			{
				case 1:
				{
					alpha += (*fTimeStep)*falpha_0*pow(l,falpha_exp);
					break;
				}
				case 2:
				{
					alpha += (*fTimeStep)*falpha_0*pow(1.-alpha,falpha_exp)*pow(1.-flambda_0*l,flambda_exp);
					break;
				}
				case 3:
				{
					alpha += (*fTimeStep)*falpha_0*pow(l_dot,falpha_exp);
					break;
				}
			}
		}
		 	
		if (alpha >= 1.)
		{
			fTraction = 0.;
			return fTraction;
		}
		
		/* scale the final tractions */
		fTraction *= 1.-alpha;
		if (l > kSmall)
		{
			fTraction[0] *= l_0/l;
			fTraction[1] *= l_1/l;
			fTraction[2] *= l_2/l;
		}
		else
		{
			if (fabs(l_0 - l_0_old) < kSmall)
				fTraction[0] = 0.;
			if (fabs(l_1 - l_1_old) < kSmall)
				fTraction[1] = 0.;
			if (fabs(l_2 - l_2_old) < kSmall)
				fTraction[2] = 0.; 
		}
		
		/* handle penetration */
	//	if (u_n < 0) fTraction[1] += fK*u_n;
		
		/* update the rest of the state variables */
		state2[8] = state2[0];
		state2[9] = state2[1];
		state2[10] = state2[2];
		state2[11] = state2[3];
		state2[12] = state2[4];
		state2[13] = state2[5];
		state2[14] = state2[6];
		state2[15] = state2[7];
		state2[0] = u_t0;
		state2[1] = u_t1;
		state2[2] = u_n;
		state2[3] = fTraction[0];
		state2[4] = fTraction[1];
		state2[5] = fTraction[2];
		state2[6] = l_dot;
		state2[7] = alpha;
		state2[16] += fTraction[0]*(u_t0 - state2[8]) + fTraction[1]*(u_t1 - state2[9])+ fTraction[2]*(u_n - state2[10]);
		
	}
	
	return fTraction;
	
}

/* potential stiffness */
const dMatrixT& YoonAllen3DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif	

	double u_t0 = jump_u[0];
	double u_t1 = jump_u[1];
	double u_n = jump_u[2];

	/* compute the current tractions first */
	int n_prony = ftau.Length();
	const double *state2 = state.Pointer(n_prony);

	/*double u_t0_dot = (u_t0 - state2[0])/fTimeStep;
	double u_t1_dot = (u_t1 - state2[1])/fTimeStep;
	double u_n_dot = (u_n - state2[2])/fTimeStep;
	*/
	
	double l_0 = u_t0/fd_c_t;
	double l_1 = u_t1/fd_c_t;
	double l_2 = u_n/fd_c_n;
	double l = sqrt(l_0*l_0+l_1*l_1+l_2*l_2); // l stands for lambda 
	
	fStiffness = 0.;
	
	/* take care of the first time through */
	if (l < kSmall)
	{
		return fStiffness;
	}
	
//	double l_dot = (l_0*u_t_dot/fd_c_t+l_1*u_n_dot/fd_c_n)/l;
	
	double l_0_old = state2[8]/fd_c_t;
	double l_1_old = state2[9]/fd_c_t;
	double l_2_old = state2[10]/fd_c_n;
	double l_old = sqrt(l_0_old*l_0_old+l_1_old*l_1_old+l_2_old*l_2_old);
	double prefactold = 1./(1-state2[2*knumDOF+9]);
	double l_dot = 0.0;
	if (fabs(*fTimeStep) > kSmall) l_dot = (l-l_old)/(*fTimeStep); /* for dt -> 0 */

	dArrayT currTraction(knumDOF);
	currTraction = 0.;
	
	/* do the bulk of the computation now */
	if (l_old > kSmall)
	{
		prefactold *= l_old;
		if (fabs(l_0_old) > kSmall)
			currTraction[0] = prefactold/l_0_old*state2[knumDOF+8]; 
		if (fabs(l_1_old) > kSmall)
			currTraction[1] = prefactold/l_1_old*state2[knumDOF+9];  
		if (fabs(l_2_old) > kSmall)
			currTraction[2] = prefactold/l_2_old*state2[knumDOF+10];
	}
	else
	{
		currTraction[0] = prefactold*state2[knumDOF+8];
		currTraction[1] = prefactold*state2[knumDOF+9];
		currTraction[2] = prefactold*state2[knumDOF+10];
	}	
		
	double tmpSum = fE_infty*(*fTimeStep);
	for (int i = 0; i < n_prony; i++)
	  tmpSum -= fexp_tau[i]*ftau[i];
	tmpSum *= l_dot;
	
	currTraction += tmpSum;

	/* update the S coefficients */
	for (int i = 0; i < n_prony; i++)
	{
		currTraction += state[i]*fexp_tau[i];	
	}
	
	/* grab the damage parameter */
	double alpha = state[n_prony+2*knumDOF+1];
	
	if (alpha >= 1.)
	{
		fStiffness = 0.;
		return fStiffness;
	}
	
	if (fabs(l_dot) > kSmall)
	{
		/* Now tackle some stiffnesses */
		fStiffness = 0.0;
		if (fabs(*fTimeStep) > kSmall) fStiffness = tmpSum/l_dot*(1-alpha)/l/l/(*fTimeStep);
		if (l < kSmall)
		{
			fStiffness[0] *= l_0*l_0;
			fStiffness[1] *= l_0*l_1;
			fStiffness[2] *= l_0*l_2;
			fStiffness[3] *= l_1*l_0;
			fStiffness[4] *= l_1*l_1;
			fStiffness[5] *= l_1*l_2;
			fStiffness[6] *= l_2*l_0;
			fStiffness[7] *= l_2*l_1;
			fStiffness[8] *= l_2*l_2;
		}
	}
	
	/* scale the final tractions up to components of the gap vector */
	if (l < kSmall)
	{
		currTraction *= (1.-alpha)/l;

		if (fabs(l_0 - l_0_old) < kSmall)
			currTraction[0] = 0.;
		if (fabs(l_1 - l_1_old) < kSmall)
			currTraction[1] = 0.;
		if (fabs(l_2 - l_2_old) < kSmall)
			currTraction[2] = 0.;
	}
	else
	{
		currTraction *= (1-alpha)/l;
		
		/* delta_ij terms added to stiffness */
		fStiffness[0] += currTraction[0];
		fStiffness[4] += currTraction[1];
		fStiffness[8] += currTraction[2];
	
		/* now scale currTraction to actually be the traction */
		currTraction[0] *= l_0;
		currTraction[1] *= l_1;
		currTraction[2] *= l_2;
	}
	
	double scratch;
	if (l > kSmall)
		scratch = 1./l/l;
	else
		scratch = 0;
	if (l_dot > kSmall)
	{
		double ftmp;
		switch (idamage)
		{
			case 1:
			{
				ftmp = falpha_exp*pow(l,falpha_exp-1.)*falpha_0*(*fTimeStep)/(1.-alpha);
				break;
			}
			case 2:
			{
				ftmp = -flambda_0*(*fTimeStep)*flambda_exp*pow(1.-flambda_0*l,flambda_exp-1.)/(1-alpha);
				break;
			}
			case 3:
			{
				ftmp = falpha_exp*pow(l_dot,falpha_exp-1)*falpha_0/(1.-alpha);
				break;
			}
		}
		if (l > kSmall)
			ftmp /= l;
		scratch += ftmp;
	}
		
	if (l > kSmall)
	{
		l_0 *= scratch;
		l_1 *= scratch;
		l_2 *= scratch;
	}
	else
	{
		l_0 = scratch;
		l_1 = scratch;
		l_2 = scratch;
	}
			
	
	fStiffness[0] -= currTraction[0]*l_0;
	fStiffness[1] -= currTraction[0]*l_1;
	fStiffness[2] -= currTraction[0]*l_2;
	fStiffness[3] -= currTraction[1]*l_0;
	fStiffness[4] -= currTraction[1]*l_1;
	fStiffness[5] -= currTraction[1]*l_2;
	fStiffness[6] -= currTraction[2]*l_0;
	fStiffness[7] -= currTraction[2]*l_1;
	fStiffness[8] -= currTraction[2]*l_2;
	
	/*scale stiffnesses by critical lengths */
	fStiffness[0] /= fd_c_t;
	fStiffness[1] /= fd_c_t;
	fStiffness[2] /= fd_c_n;
	fStiffness[3] /= fd_c_t;
	fStiffness[4] /= fd_c_t;
	fStiffness[5] /= fd_c_n;
	fStiffness[6] /= fd_c_t;
	fStiffness[7] /= fd_c_t;
	fStiffness[8] /= fd_c_n; 

	/* penetration */
//	if (u_n < 0) fStiffness[3] += fK;

	return fStiffness;

}

/* surface status */
SurfacePotentialT::StatusT YoonAllen3DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif
	
	int n_prony = ftau.Length();
	if (state[n_prony+2*knumDOF+1]+kSmall > 1.)
		return Failed;
	else if (state[n_prony+2*knumDOF+1] > kSmall)
		return Critical;
	else
		return Precritical;

}

#if 0
/* print parameters to the output stream */
void YoonAllen3DT::Print(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	out << " Cohesive stress . . . . . . . . . . . . . . . . = " << fsigma_0   << '\n';
	out << " Normal length scale . . . . . . . . . . . . . . = " << fd_c_n     << '\n';
	out << " Tangential length scale . . . . . . . . . . . . = " << fd_c_t     << '\n';
	out << " Long-time modulus . . . . . . . . . . . . . . . = " << fE_infty   << '\n';
	out << " Number of terms in prony series . . . . . . . . = " << n_prony << "\n";
	for (int i = 0; i < n_prony;i++)
	{
		out << " Transient modulus for mode "<<i<<" . . . . . . . . .  = " << fE_t[i]    << '\n';
		out << " Time constant for mode "<<i<<" . . . . . . . . . . .  = " << ftau[i]/fE_t[i] << '\n';
	}
	out << " Damage evolution law code . . . . . . . . . . . = " << idamage << "\n";
	if (idamage == 2) 
	{
		out << " Damage evolution law exponent . . . . . . . . . = " << falpha_exp << "\n";
		out << " Damage evolution law constant . . . . . . . . . = " << falpha_0 << "\n";
	}
	out << " Damage evolution law lambda exponent. . . . . . = " << flambda_exp << "\n";
	out << " Damage evolution law lambda prefactor . . . . . = " << flambda_0 << "\n";
	out << " Penetration stiffness multiplier. . . . . . . . = " << fpenalty   << '\n';
#else
#pragma unused(out)
#endif
}
#endif

/* returns the number of variables computed for nodal extrapolation
 * during for element output, ie. internal variables. Returns 0
 * by default 
 */
int YoonAllen3DT::NumOutputVariables(void) const { return 3; }

void YoonAllen3DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(3);
	labels[0] = "lambda";
	labels[1] = "lambda_dot";
	labels[2] = "alpha";
}

void YoonAllen3DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif	

	double u_t0 = jump_u[0];
	double u_t1 = jump_u[1];
	double u_n = jump_u[2];
	
	double l_t0 = u_t0/fd_c_t;
	double l_t1 = u_t1/fd_c_t;
	double l_n = u_n/fd_c_n;

	int n_prony = ftau.Length();
	output[0] = sqrt(l_t0*l_t0 + l_t1*l_t1 + l_n*l_n); 
	output[1] = state[n_prony+2*knumDOF+8];
	output[2] = state[n_prony+2*knumDOF+9];
	
	double u_t0_dot = 0.0;
	double u_t1_dot = 0.0;
	double u_n_dot = 0.0;
	if (fabs(*fTimeStep) > kSmall) /* allow dt -> 0 */
	{
		u_t0_dot = (u_t0 - state[n_prony+8])/(*fTimeStep);
		u_t1_dot = (u_t1 - state[n_prony+9])/(*fTimeStep);
		u_n_dot = (u_n - state[n_prony+10])/(*fTimeStep);	
	}

	double l_dot = (l_t0*u_t0_dot/fd_c_t+l_t1*u_t1_dot/fd_c_t+l_n*u_n_dot/fd_c_n)/output[0];
	if (l_dot > kSmall)
		output[2] += falpha_0*pow(output[0],falpha_exp);
	
}

bool YoonAllen3DT::CompatibleOutput(const SurfacePotentialT& potential) const
{
#ifdef __NO_RTTI__
#pragma unused(potential)
	return false;
#else
	const YoonAllen3DT* pTH = dynamic_cast<const YoonAllen3DT*>(&potential);
	return pTH != NULL;
#endif
}

void YoonAllen3DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SurfacePotentialT::DefineParameters(list);

	/* common limit */
	LimitT non_negative(0, LimitT::LowerInclusive);

	/* traction potential parameters */
	ParameterT sigma_0(fsigma_0, "sigma_0");
	sigma_0.AddLimit(non_negative);
	ParameterT d_c_n(fd_c_n, "d_c_n");
	d_c_n.AddLimit(non_negative);
	ParameterT d_c_t(fd_c_t, "d_c_t");
	d_c_t.AddLimit(non_negative);
	list.AddParameter(sigma_0);
	list.AddParameter(d_c_n);
	list.AddParameter(d_c_t);

	 /* asymptotic modulus of cohesive zone */
	ParameterT E_infty(fE_infty, "E_infty");
	E_infty.AddLimit(non_negative);
	list.AddParameter(E_infty);

	/* penalty stiffness */
	ParameterT penalty(fpenalty, "penalty");
	penalty.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(penalty);
}

/* information about subordinate parameter lists */
void YoonAllen3DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SurfacePotentialT::DefineSubs(sub_list);

	/* prony series */
	sub_list.AddSub("Prony_series");

	/* damage evolution */
	sub_list.AddSub("Yoon-Allen_damage_choice", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* YoonAllen3DT::NewSub(const StringT& name) const
{
	if (name == "Prony_series" || name == "Yoon-Allen_damage_choice")
	{
		/* use definitions from 2D model */
		YoonAllen2DT YH2D;
		return YH2D.NewSub(name);
	}
	else /* inherited */
		return SurfacePotentialT::NewSub(name);
}
	
/* accept parameter list */
void YoonAllen3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SurfacePotentialT::TakeParameterList(list);

	/* traction potential parameters */
	fsigma_0 = list.GetParameter("sigma_0");
	fd_c_n = list.GetParameter("d_c_n");
	fd_c_t = list.GetParameter("d_c_t");

	 /* asymptotic modulus of cohesive zone */
	fE_infty = list.GetParameter("E_infty");

	/* penetration stiffness */
	fpenalty = list.GetParameter("penalty");
	fK = fpenalty*fsigma_0/fd_c_n;

	/* prony series */
	const ParameterListT& prony_series = list.GetList("Prony_series");
	int n_terms = prony_series.NumLists("Prony_pair");
	ftau.Dimension(n_terms);
	fE_t.Dimension(n_terms);
	fexp_tau.Dimension(n_terms);
	for (int i = 0; i < n_terms; i++)
	{
		const ParameterListT& prony_pair = prony_series.GetList("Prony_pair", i);
		fE_t[i] = prony_pair.GetParameter("transient_modulus");
		ftau[i] = prony_pair.GetParameter("time_constant");
	
		/* decay function */
		fexp_tau[i] = exp(-(*fTimeStep)/ftau[i]) - 1.0;
		
		/* scale ftau by fE_t to reduce multiplications in traction
		 * and stiffness routines */
		ftau[i] *= fE_t[i];
	}
	fCurrentTimeStep = *fTimeStep;

	/* damage evolution */
	const ParameterListT& damage = list.GetListChoice(*this, "Yoon-Allen_damage_choice");
	if (damage.Name() == "YA_damage_1")
	{
		idamage = 1;
		falpha_exp = damage.GetParameter("alpha_exp");
		falpha_0 = damage.GetParameter("alpha_0");
	}
	else if (damage.Name() == "YA_damage_2")
	{
		idamage = 2;
		falpha_exp = damage.GetParameter("alpha_exp");
		falpha_0 = damage.GetParameter("alpha_0");
		flambda_exp = damage.GetParameter("lambda_exp");
		flambda_0 = damage.GetParameter("lambda_0");
	}
	else if (damage.Name() == "YA_damage_3")
	{
		idamage = 3;
		falpha_exp = damage.GetParameter("alpha_exp");
		falpha_0 = damage.GetParameter("alpha_0");	
	}
	else
		ExceptionT::GeneralFail("YoonAllen3DT::TakeParameterList", "unrecognized damage option \"%s\"",
			damage.Name().Pointer());
}
