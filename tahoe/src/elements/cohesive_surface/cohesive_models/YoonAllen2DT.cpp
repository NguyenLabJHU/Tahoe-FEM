/* $Id: YoonAllen2DT.cpp,v 1.3 2002-07-02 19:55:17 cjkimme Exp $ */

#include "YoonAllen2DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionCodes.h"
#include "fstreamT.h"
#include "StringT.h"
#include "SecantMethodT.h"

/* class parameters */

using namespace Tahoe;

const int knumDOF = 2;

/* constructor */
YoonAllen2DT::YoonAllen2DT(ifstreamT& in, const double& time_step): 
	SurfacePotentialT(knumDOF),
	fTimeStep(time_step)
{
	/* traction potential parameters */
	in >> fsigma_0; if (fsigma_0 < 0) throw eBadInputValue;
	in >> fd_c_n; if (fd_c_n < 0) throw eBadInputValue;
	in >> fd_c_t; if (fd_c_t < 0) throw eBadInputValue;
	
	/* moduli and time constants */
	in >> fE_infty; if (fE_infty < 0) throw eBadInputValue;
	in >> fNumRelaxTimes; if (fNumRelaxTimes < 0) throw eBadInputValue;
	ftau.Allocate(fNumRelaxTimes);
	fE_t.Allocate(fNumRelaxTimes);
	fexp_tau.Allocate(fNumRelaxTimes);
	for (int i = 0;i < fNumRelaxTimes; i++)
	{
		in >> fE_t[i]; if (fE_t[i] < 0) throw eBadInputValue;
		in >> ftau[i]; if (ftau[i] < 0) throw eBadInputValue;
		fexp_tau[i] = exp(-fTimeStep/ftau[i]) - 1.;
		
		/* scale ftau by fE_t to reduce multiplications in traction
		 * and stiffness routines
		 */
		ftau[i] *= fE_t[i];
	}

	/* damage evolution law parameters */
	in >> falpha_exp; //if (falpha_exp < 1.) throw eBadInputValue;
	in >> falpha_0; if (falpha_0 <= kSmall) throw eBadInputValue;
	in >> flambda_exp; if (flambda_exp > -1.) throw eBadInputValue;
	in >> flambda_0; if (flambda_0 < 1.) throw eBadInputValue;

	/* stiffness multiplier */
	in >> fpenalty; if (fpenalty < 0) throw eBadInputValue;
	
	/* penetration stiffness */
	fK = fpenalty*fsigma_0/fd_c_n;

}

/*initialize state variables with values from the rate-independent model */
void YoonAllen2DT::InitStateVariables(ArrayT<double>& state)
{
 	int num_state = NumStateVariables();
	if (state.Length() != num_state) 
	{
	  	cout << "\n SurfacePotentialT::InitStateVariables: expecting state variable array\n"
		     <<   "     length " << num_state << ", found length " << state.Length() << endl;
		throw eSizeMismatch;
	}

	/* clear */
	if (num_state > 0) state = 0.0;

	/* 	The first fNumRelaxTimes slots in state are the previous
	 *  steps hereditary integral (more or less). 0..fNumRelaxTimes-1
	 *	The next kNumDOF slots are the previous step's components of
	 *  the gap vector.   fNumRelaxTimes..fNumRelaxTimes+knumDOF-1
	 *	The next kNumDOF slots are the ith components of the previous
	 *	step's tractions. fNumRelaxTimes+knumDOF..fNumRelaxTimes+2knumDOF-1
	 *  The next slot is the rate of change of the previous step's
	 *  lambda. fNumRelaxTimes+2knumDOF
	 *	The next slot is the previous step's value of alpha, the damage
	 *	coefficient.  fNumRelaxTimes+2knumDOF+1
	 *	The final slot holds the integral of T dot dDelta. fNumRelaxTimes+2knumDOF+2
	 */ 

}

/* return the number of state variables needed by the model */
int YoonAllen2DT::NumStateVariables(void) const 
{ 
	return 2*knumDOF + fNumRelaxTimes + 3; 
}

/* surface potential */ 
double YoonAllen2DT::FractureEnergy(const ArrayT<double>& state) 
{
   	return state[2*knumDOF + fNumRelaxTimes + 2]; 
}

double YoonAllen2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(jump_u)
#pragma unused(state)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw eSizeMismatch;
	if (state.Length() != NumStateVariables()) throw eSizeMismatch;
#endif

	cout << "YoonAllen2DT::Potential is not implemented. It's viscoelastic \n";
	return 0.;
}
	
/* traction vector given displacement jump vector */	
const dArrayT& YoonAllen2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw eSizeMismatch;
	if (state.Length() != NumStateVariables()) throw eSizeMismatch;
	if (fTimeStep <= 0.0) {
		cout << "\n YoonAllen2DT::Traction: expecting positive time increment: "
		     << fTimeStep << endl;
		throw eBadInputValue;
	}
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];
	
	double *state2 = state.Pointer(fNumRelaxTimes);

	double u_t_dot = (u_t - state2[0])/fTimeStep;
	double u_n_dot = (u_n - state2[1])/fTimeStep;
	
	double l_0 = u_t/fd_c_t;
	double l_1 = u_n/fd_c_n;
	double l = sqrt(l_0*l_0+l_1*l_1); // l stands for lambda 
	
	fTraction = 0.;
	
	/* handle the first time through */
	if (l < kSmall)
	{
		return fTraction;
	}
	
//	double l_dot = (l_0*u_t_dot/fd_c_t+l_1*u_n_dot/fd_c_n)/l;
	
	double l_0_old = state2[0]/fd_c_t;
	double l_1_old = state2[1]/fd_c_n;
	double l_old = sqrt(l_0_old*l_0_old+l_1_old*l_1_old);
	double prefactold = l_old/(1-state2[2*knumDOF+1]);
//	cout << " l_dot comp : " << l_dot << " " << (l-l_old)/fTimeStep<<"\n";
	double l_dot = (l-l_old)/fTimeStep;

//	cout << "l_0_old " << l_0_old;

	if (prefactold > kSmall)
	{
		/* do the bulk of the computation now */
		if (fabs(l_0_old) > kSmall)
			fTraction[0] = prefactold/l_0_old*state2[knumDOF]; 
		if (fabs(l_1_old) > kSmall)
			fTraction[1] = prefactold/l_1_old*state2[knumDOF+1];  
	}
	
//	cout << "fTraction (1st) : " << fTraction[0] << " " << fTraction[1] << "\n";
	
	double tmpSum = fE_infty*fTimeStep;
	for (int i = 0; i < fNumRelaxTimes; i++)
	  tmpSum -= fexp_tau[i]*ftau[i];
	tmpSum *= l_dot;
	
	fTraction += tmpSum;
	
//	cout << "fTraction (2nd) : " << fTraction[0] << " " << fTraction[1] << "\n";
	
	/* update the S coefficients */
	for (int i = 0; i < fNumRelaxTimes; i++)
	{
		state[i] += (state[i]-state2[2*knumDOF]*ftau[i])*fexp_tau[i];
		fTraction += state[i]*fexp_tau[i];	
	}
	
	/* evolve the damage parameter */
	double alpha = state2[2*knumDOF+1];
	if (l_dot > kSmall)
		alpha += fTimeStep*falpha_0*pow(l,falpha_exp);
//		alpha += fTimeStep*falpha_0*pow(1.-alpha,falpha_exp)*pow(1.-flambda_0*l,flambda_exp);
//		alpha += fTimeStep*falpha_0*pow(l_dot,falpha_exp);

	 	
	if (alpha >= 1.)
	{
		fTraction = 0.;
		return fTraction;
	}
	
	/* scale the final tractions */
	fTraction[0] *= l_0/l*(1.-alpha);
	fTraction[1] *= l_1/l*(1.-alpha);

	/* handle penetration */
//	if (u_n < 0) fTraction[1] += fK*u_n;
	
	/* overwrite rates of change with gap vector increments in order to 
	 * compute the work done in this step
	 */
	u_t_dot = u_t - state2[0];
	u_n_dot = u_n - state2[1];
	
	/* update the rest of the state variables */
	state2[0] = u_t;
	state2[1] = u_n;
	state2[2] = fTraction[0];
	state2[3] = fTraction[1];
	state2[4] = l_dot;
	state2[5] = alpha;
	state2[6] += fTraction[0]*u_t_dot + fTraction[1]*u_n_dot;
	
//	cout << "fTraction : " << fTraction[0] << " " << fTraction[1] << " alpha " << alpha <<"\n";

	return fTraction;
	
}

/* potential stiffness */
const dMatrixT& YoonAllen2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw eSizeMismatch;
	if (state.Length() != NumStateVariables()) throw eGeneralFail;
#endif	

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	/* fStiffness = {dT_0/dd_0,dT_0/dd_1,dT_1,dd_0,dT_1/dd_1} */
	/* compute the current tractions first */
	double *state2 = state.Pointer(fNumRelaxTimes);

	double u_t_dot = (u_t - state2[0])/fTimeStep;
	double u_n_dot = (u_n - state2[1])/fTimeStep;
	
	double l_0 = u_t/fd_c_t;
	double l_1 = u_n/fd_c_n;
	double l = sqrt(l_0*l_0+l_1*l_1); // l stands for lambda 
	
	fStiffness = 0.;
	
	/* take care of the first time through */
	if (l < kSmall)
	{
		return fStiffness;
	}
	
//	double l_dot = (l_0*u_t_dot/fd_c_t+l_1*u_n_dot/fd_c_n)/l;
	
	double l_0_old = state2[0]/fd_c_t;
	double l_1_old = state2[1]/fd_c_n;
	double l_old = sqrt(l_0_old*l_0_old+l_1_old*l_1_old);
	double prefactold = l_old/(1-state2[2*knumDOF+1]);
	double l_dot = (l-l_old)/fTimeStep;

	dArrayT currTraction(2);
	currTraction = 0.;
	
	if (prefactold > kSmall) 
	{
		/* do the bulk of the computation now */
		if (fabs(l_0_old) > kSmall)
			currTraction[0] = prefactold/l_0_old*state2[knumDOF]; 
		if (fabs(l_1_old) > kSmall)
			currTraction[1] = prefactold/l_1_old*state2[knumDOF+1];  
	}
		
	double tmpSum = fE_infty*fTimeStep;
	for (int i = 0; i < fNumRelaxTimes; i++)
	  tmpSum -= fexp_tau[i]*ftau[i];
	tmpSum *= l_dot;
	
	currTraction += tmpSum;

	/* update the S coefficients */
	double ftmp;
	for (int i = 0; i < fNumRelaxTimes; i++)
	{
		ftmp = state[i] + (state[i]-state2[2*knumDOF]*ftau[i])*fexp_tau[i];
		currTraction += ftmp*fexp_tau[i];	
	}
	
	/* evolve the damage parameter */
	double alpha = state[fNumRelaxTimes+2*knumDOF+1];
	if (l_dot > kSmall)
	{
		alpha += fTimeStep*falpha_0*pow(l,falpha_exp);
//		alpha += fTimeStep*falpha_0*pow(1.-alpha,falpha_exp)*pow(1.-flambda_0*l,flambda_exp);
//		alpha += fTimeStep*falpha_0*pow(l_dot,falpha_exp);
	}
		
	if (alpha >= 1.)
	{
		fStiffness = 0.;
		return fStiffness;
	}
	
	if (fabs(l_dot) > kSmall)
	{
		/* Now tackle some stiffnesses */
		fStiffness = tmpSum/l_dot*(1-alpha)/l/l/fTimeStep;
		fStiffness[0] *= l_0*l_0;
		fStiffness[1] *= l_0*l_1;
		fStiffness[2] *= l_1*l_0;
		fStiffness[3] *= l_1*l_1;
	}
	
	/* scale the final tractions up to components of the gap vector */
	currTraction *= (1.-alpha)/l;

	/* delta_ij terms added to stiffness */
	fStiffness[0] += currTraction[0];
	fStiffness[3] += currTraction[1];
	
	/* now scale currTraction to actually be the traction */
	currTraction[0] *= l_0;
	currTraction[1] *= l_1;
	
	double scratch;
	scratch = 1./l/l;
	if (l_dot > kSmall)
		scratch += falpha_exp*pow(l,falpha_exp-2.)*falpha_0*fTimeStep/(1.-alpha);
//	    scratch -= flambda_0*fTimeStep*flambda_exp*pow(1.-flambda_0*l,flambda_exp-1.)/(1-alpha);
//		scratch += falpha_exp*pow(l_dot,falpha_exp-1)*falpha_0*fTimeStep/(1.-alpha)/l;
	l_0 *= scratch;
	l_1 *= scratch;
	
	fStiffness[0] -= currTraction[0]*l_0;
	fStiffness[1] -= currTraction[0]*l_1;
	fStiffness[2] -= currTraction[1]*l_0;
	fStiffness[3] -= currTraction[1]*l_1;
	
	/*scale stiffnesses by critical lengths */
	fStiffness[0] /= fd_c_t;
	fStiffness[1] /= fd_c_n;
	fStiffness[2] /= fd_c_t;
	fStiffness[3] /= fd_c_n; 

	/* penetration */
//	if (u_n < 0) fStiffness[3] += fK;

	return fStiffness;

}

/* surface status */
SurfacePotentialT::StatusT YoonAllen2DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw eSizeMismatch;
#endif
	
	if (state[fNumRelaxTimes+2*knumDOF+1]+kSmall > 1.)
		return Failed;
	else if (state[fNumRelaxTimes+2*knumDOF+1] > kSmall)
		return Critical;
	else
		return Precritical;

}

void YoonAllen2DT::PrintName(ostream& out) const
{
	out << " Yoon-Allen 2D \n";
}

/* print parameters to the output stream */
void YoonAllen2DT::Print(ostream& out) const
{
	out << " Cohesive stress . . . . . . . . . . . . . . . . = " << fsigma_0   << '\n';
	out << " Normal length scale . . . . . . . . . . . . . . = " << fd_c_n     << '\n';
	out << " Tangential length scale . . . . . . . . . . . . = " << fd_c_t     << '\n';
	out << " Long-time modulus . . . . . . . . . . . . . . . = " << fE_infty   << '\n';
	out << " Number of terms in prony series . . . . . . . . = " << fNumRelaxTimes << "\n";
	for (int i = 0; i < fNumRelaxTimes;i++)
	{
		out << " Transient modulus for mode "<<i<<" . . . . . . . . .  = " << fE_t[i]    << '\n';
		out << " Time constant for mode "<<i<<" . . . . . . . . . . .  = " << ftau[i]/fE_t[i] << '\n';
	}
	out << " Damage evolution law exponent . . . . . . . . . = " << falpha_exp << "\n";
	out << " Damage evolution law constant . . . . . . . . . = " << falpha_0 << "\n";
	out << " Damage evolution law lambda exponent. . . . . . = " << flambda_exp << "\n";
	out << " Damage evolution law lambda prefactor . . . . . = " << flambda_0 << "\n";
	out << " Penetration stiffness multiplier. . . . . . . . = " << fpenalty   << '\n';
}

/* returns the number of variables computed for nodal extrapolation
 * during for element output, ie. internal variables. Returns 0
 * by default 
 */
int YoonAllen2DT::NumOutputVariables(void) const { return 3; }

void YoonAllen2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Allocate(3);
	labels[0] = "lambda";
	labels[1] = "lambda_dot";
	labels[2] = "alpha";
}

void YoonAllen2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw eGeneralFail;
#endif	
	double u_t = jump_u[0];
	double u_n = jump_u[1];
	
	double l_t = u_t/fd_c_t;
	double l_n = u_n/fd_c_n;

	output[0] = sqrt(l_t*l_t + l_n*l_n); 
	output[1] = state[fNumRelaxTimes+4];
	output[2] = state[fNumRelaxTimes+2*knumDOF+1];
	
	double u_t_dot = (u_t - state[fNumRelaxTimes])/fTimeStep;
	double u_n_dot = (u_n - state[fNumRelaxTimes+1])/fTimeStep;
	double l_dot = (l_t*u_t_dot/fd_c_t+l_t*u_n_dot/fd_c_n)/output[0];
	if (l_dot > kSmall)
		output[2] += falpha_0*pow(output[0],falpha_exp);
	
}

bool YoonAllen2DT::NeedsNodalInfo(void) { return false; }

int YoonAllen2DT::NodalQuantityNeeded(void) 
{ 
        return -1; 
}

int YoonAllen2DT::ElementGroupNeeded(void) 
{
	return -1;
}

bool YoonAllen2DT::CompatibleOutput(const SurfacePotentialT& potential) const
{
#ifdef __NO_RTTI__
#pragma unused(potential)
	return false;
#else
	const YoonAllen2DT* pTH = dynamic_cast<const YoonAllen2DT*>(&potential);
	return pTH != NULL;
#endif
}


