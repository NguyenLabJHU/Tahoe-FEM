/* $Id: YoonAllen2DT.cpp,v 1.14 2004-06-17 07:13:28 paklein Exp $ */
#include "YoonAllen2DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"
#include "ifstreamT.h"
#include "StringT.h"

using namespace Tahoe;

/* class parameters */
const int knumDOF = 2;

/* constructor */
YoonAllen2DT::YoonAllen2DT(ifstreamT& in, const double& time_step): 
	SurfacePotentialT(knumDOF),
	fTimeStep(time_step)
{
	/* traction potential parameters */
	in >> fsigma_0; if (fsigma_0 < 0) throw ExceptionT::kBadInputValue;
	in >> fd_c_n; if (fd_c_n < 0) throw ExceptionT::kBadInputValue;
	in >> fd_c_t; if (fd_c_t < 0) throw ExceptionT::kBadInputValue;
	
	/* moduli and time constants */
	in >> fE_infty; if (fE_infty < 0) throw ExceptionT::kBadInputValue;
	in >> iNumRelaxTimes; if (iNumRelaxTimes < 0) throw ExceptionT::kBadInputValue;

	ftau.Dimension(iNumRelaxTimes);
	fE_t.Dimension(iNumRelaxTimes);
	fexp_tau.Dimension(iNumRelaxTimes);

	for (int i = 0;i < iNumRelaxTimes; i++)
	{
		in >> fE_t[i]; if (fE_t[i] < 0) throw ExceptionT::kBadInputValue;
		in >> ftau[i]; if (ftau[i] < 0) throw ExceptionT::kBadInputValue;
	
		fexp_tau[i] = exp(-fTimeStep/ftau[i]) - 1.;
		
		/* scale ftau by fE_t to reduce multiplications in traction
		 * and stiffness routines
		 */
		ftau[i] *= fE_t[i];
	}

	in >> idamage; if (idamage < 1 && idamage > 3) throw ExceptionT::kBadInputValue;
	/* damage evolution law parameters */
	in >> falpha_exp; //if (falpha_exp < 1.) throw ExceptionT::kBadInputValue;
	in >> falpha_0; if (falpha_0 <= kSmall) throw ExceptionT::kBadInputValue;
	if (idamage == 2) 
	{
		in >> flambda_exp; if (flambda_exp > -1.) throw ExceptionT::kBadInputValue;	
		in >> flambda_0; if (flambda_0 < 1.) throw ExceptionT::kBadInputValue;
	}

	/* stiffness multiplier */
	in >> fpenalty; if (fpenalty < 0) throw ExceptionT::kBadInputValue;

	
	/* penetration stiffness */
	fK = fpenalty*fsigma_0/fd_c_n;

}

/*initialize state variables with values from the rate-independent model */
void YoonAllen2DT::InitStateVariables(ArrayT<double>& state)
{
 	int num_state = NumStateVariables();
	if (state.Length() != num_state) 
	{
#ifndef _SIERRA_TEST_	
	  	cout << "\n SurfacePotentialT::InitStateVariables: expecting state variable array\n"
		     <<   "     length " << num_state << ", found length " << state.Length() << endl;
#endif
		throw ExceptionT::kSizeMismatch;	
	}

	/* clear */
	if (num_state > 0) state = 0.0;

	/* 	The first iNumRelaxTimes slots in state are the previous
	 *  steps hereditary integral (more or less). 0..iNumRelaxTimes-1
	 *	The next kNumDOF slots are the previous step's components of
	 *  the gap vector.   iNumRelaxTimes..iNumRelaxTimes+knumDOF-1
	 *	The next kNumDOF slots are the ith components of the previous
	 *	step's tractions. iNumRelaxTimes+knumDOF..iNumRelaxTimes+2 knumDOF-1
	 *  The next slot is the previous step's
	 *  lambda. iNumRelaxTimes+2knumDOF
	 *	The next slot is the previous step's value of alpha, the damage
	 *	coefficient.  iNumRelaxTimes+2knumDOF+1
	 *	The final slot holds the integral of T dot dDelta. iNumRelaxTimes+2knumDOF+2
	 */ 

}

/* return the number of state variables needed by the model */
int YoonAllen2DT::NumStateVariables(void) const 
{ 
	return 4*knumDOF + iNumRelaxTimes + 5; 
}

/* surface potential */ 
double YoonAllen2DT::FractureEnergy(const ArrayT<double>& state) 
{
   	return state[4*knumDOF + iNumRelaxTimes + 4]; 
}

double YoonAllen2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(jump_u)
#pragma unused(state)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

#ifndef _SIERRA_TEST_
	cout << "YoonAllen2DT::Potential is not implemented. It's viscoelastic \n";
#endif
	return 0.;
}
	
/* traction vector given displacement jump vector */	
const dArrayT& YoonAllen2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
	const char caller[] = "YoonAllen2DT::Traction";
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) ExceptionT::SizeMismatch(caller);
	if (state.Length() != NumStateVariables()) ExceptionT::SizeMismatch(caller);
	if (fTimeStep < 0.0) {
#ifndef _SIERRA_TEST_
		ExceptionT::BadInputValue(caller, "expecting non-negative time increment: %g", fTimeStep);
#endif		     
		throw ExceptionT::kBadInputValue;
	}
#endif

	if (!qIntegrate)
	{
		fTraction[0] = state[iNumRelaxTimes+knumDOF];
		fTraction[1] = state[iNumRelaxTimes+knumDOF+1];
		return fTraction;
	}
	else
	{
		double u_t = jump_u[0];
		double u_n = jump_u[1];
	
		double *state2 = state.Pointer(iNumRelaxTimes);
		
		// optional way to calculate rates
		/*
		double u_t_dot = (u_t - state2[0])/fTimeStep;
		double u_n_dot = (u_n - state2[1])/fTimeStep;
		*/
		
		double l_0 = u_t/fd_c_t;
		double l_1 = u_n/fd_c_n;
		double l = sqrt(l_0*l_0+l_1*l_1); // l stands for lambda 
		
		fTraction = 0.;
		
		double l_0_old = state2[0]/fd_c_t;
		double l_1_old = state2[1]/fd_c_n;
		double l_old = sqrt(l_0_old*l_0_old+l_1_old*l_1_old);
		double prefactold = 1./(1-state2[2*knumDOF+1]);
		double l_dot = 0.0;
		if (fabs(fTimeStep) > kSmall) l_dot = (l-l_old)/fTimeStep; /* allow for dt -> 0 */

		/* do the bulk of the computation now */
		if (l_old > kSmall) 
		{
			prefactold *= l_old;
			if (fabs(l_0_old) > kSmall)
				fTraction[0] = prefactold/l_0_old*state2[knumDOF]; 
			if (fabs(l_1_old) > kSmall)
				fTraction[1] = prefactold/l_1_old*state2[knumDOF+1];  
		}
		else
		{
			fTraction[0] = prefactold*state2[knumDOF];
			fTraction[1] = prefactold*state2[knumDOF+1];
		}
			
		
		double tmpSum = fE_infty*fTimeStep;
		for (int i = 0; i < iNumRelaxTimes; i++)
		  	tmpSum -= fexp_tau[i]*ftau[i];
		tmpSum *= l_dot;
		
		fTraction += tmpSum;
		
		/* update the S coefficients */
		for (int i = 0; i < iNumRelaxTimes; i++)
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
					alpha += fTimeStep*falpha_0*pow(l,falpha_exp);
					break;
				}
				case 2:
				{
					alpha += fTimeStep*falpha_0*pow(1.-alpha,falpha_exp)*pow(1.-flambda_0*l,flambda_exp);
					break;
				}
				case 3:
				{
					alpha += fTimeStep*falpha_0*pow(l_dot,falpha_exp);
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
		}
		else
		{
			if (fabs(l_0 - l_0_old) < kSmall)
				fTraction[0] = 0.;
			if (fabs(l_1 - l_1_old) < kSmall)
				fTraction[1] = 0.;
		}

		/* handle penetration */
	//	if (u_n < 0) fTraction[1] += fK*u_n;
		
		/* update the rest of the state variables */
		state2[6] = state2[0];
		state2[7] = state2[1];
		state2[8] = state2[2];
		state2[9] = state2[3];
		state2[10] = state2[4];
		state2[11] = state2[5];
		state2[0] = u_t;
		state2[1] = u_n;
		state2[2] = fTraction[0];
		state2[3] = fTraction[1];
		state2[4] = l_dot;
		state2[5] = alpha;
		state2[12] += fTraction[0]*(u_t - state2[6]) + fTraction[1]*(u_n - state2[7]);
		
		return fTraction;
	}
	
}

/* potential stiffness */
const dMatrixT& YoonAllen2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif	

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	/* compute the current tractions first */
	const double *state2 = state.Pointer(iNumRelaxTimes);

	/*double u_t_dot = (u_t - state2[0])/fTimeStep;
	double u_n_dot = (u_n - state2[1])/fTimeStep;
	*/
	
	double l_0 = u_t/fd_c_t;
	double l_1 = u_n/fd_c_n;
	double l = sqrt(l_0*l_0+l_1*l_1); // l stands for lambda 
	
	fStiffness = 0.;
	
	/* take care of the first time through */
//	if (l < kSmall)
//	{
//		return fStiffness;
//	}
	
//	double l_dot = (l_0*u_t_dot/fd_c_t+l_1*u_n_dot/fd_c_n)/l;
	
	double l_0_old = state2[6]/fd_c_t; 
	double l_1_old = state2[7]/fd_c_n; 
	double l_old = sqrt(l_0_old*l_0_old+l_1_old*l_1_old);
	double prefactold = 1./(1-state2[2*knumDOF+7]);
	double l_dot = 0.0;
	if (fabs(fTimeStep) > kSmall) l_dot = (l-l_old)/fTimeStep; /* allow for dt -> 0 */

	dArrayT currTraction(2);
	currTraction = 0.;
	
	/* do the bulk of the computation now */
	if (l_old > kSmall) 
	{
		prefactold *= l_old;
		if (fabs(l_0_old) > kSmall)
			currTraction[0] = prefactold/l_0_old*state2[knumDOF+6]; 
		if (fabs(l_1_old) > kSmall)
			currTraction[1] = prefactold/l_1_old*state2[knumDOF+7];  
	}
	else
	{
		currTraction[0] = prefactold*state2[knumDOF+6];
		currTraction[1] = prefactold*state2[knumDOF+7];
	}
		
	double tmpSum = fE_infty*fTimeStep;
	for (int i = 0; i < iNumRelaxTimes; i++)
	  tmpSum -= fexp_tau[i]*ftau[i];
	tmpSum *= l_dot;
	
	currTraction += tmpSum;

	/* update the S coefficients */
	for (int i = 0; i < iNumRelaxTimes; i++)
	{
		currTraction += state[i]*fexp_tau[i];	
	}
	
	/* grab the damage parameter */
	double alpha = state[iNumRelaxTimes+2*knumDOF+1];
		
	if (alpha >= 1.)
	{
		fStiffness = 0.;
		return fStiffness;
	}
	
	if (fabs(l_dot) > kSmall)
	{
		/* Now tackle some stiffnesses */
		fStiffness = 0.0;
		if (fabs(fTimeStep) > kSmall) fStiffness = tmpSum/l_dot*(1-alpha)/fTimeStep;
		if (l > kSmall) 
		{
			fStiffness /= l*l;
			fStiffness[0] *= l_0*l_0;
			fStiffness[1] *= l_0*l_1;
			fStiffness[2] *= l_1*l_0;
			fStiffness[3] *= l_1*l_1;
		}
	}
	
	/* scale the final tractions up to lambda-dependence */
	if (l < kSmall)
	{
		currTraction *= (1.-alpha);
		
		if (fabs(l_0 - l_0_old) < kSmall)
			currTraction[0] = 0.;
		if (fabs(l_1 - l_1_old) < kSmall)
			currTraction[1] = 0.;
	}
	else
	{
		currTraction *= (1.-alpha)/l;
	
		/* delta_ij terms added to stiffness */
		fStiffness[0] += currTraction[0];
		fStiffness[3] += currTraction[1];
	
		/* now scale currTraction to actually be the traction */
		currTraction[0] *= l_0;
		currTraction[1] *= l_1;
	}
	
	double scratch;
	if (l > kSmall)
		scratch = 1./l/l;
	else 
		scratch = 0.;

	if (l_dot > kSmall)
	{
		double ftmp;
		switch (idamage)
		{
			case 1:
			{
				ftmp = falpha_exp*pow(l,falpha_exp-1.)*falpha_0*fTimeStep/(1.-alpha);
				break;
			}
			case 2:
			{
			    ftmp = -flambda_0*fTimeStep*flambda_exp*pow(1.-flambda_0*l,flambda_exp-1.)/(1-alpha);
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
	}
	else
	{
		l_0 = scratch;
		l_1 = scratch;
	}
	
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
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif
	
	if (state[iNumRelaxTimes+2*knumDOF+1]+kSmall > 1.)
		return Failed;
	else if (state[iNumRelaxTimes+2*knumDOF+1] > kSmall)
		return Critical;
	else
		return Precritical;

}

void YoonAllen2DT::PrintName(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	out << " Yoon-Allen 2D \n";
#endif
}

/* print parameters to the output stream */
void YoonAllen2DT::Print(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	out << " Cohesive stress . . . . . . . . . . . . . . . . = " << fsigma_0   << '\n';
	out << " Normal length scale . . . . . . . . . . . . . . = " << fd_c_n     << '\n';
	out << " Tangential length scale . . . . . . . . . . . . = " << fd_c_t     << '\n';
	out << " Long-time modulus . . . . . . . . . . . . . . . = " << fE_infty   << '\n';
	out << " Number of terms in prony series . . . . . . . . = " << iNumRelaxTimes << "\n";
	for (int i = 0; i < iNumRelaxTimes;i++)
	{
		out << " Transient modulus for mode "<<i<<" . . . . . . . . .  = " << fE_t[i]    << '\n';
		out << " Time constant for mode "<<i<<" . . . . . . . . . . .  = " << ftau[i]/fE_t[i] << '\n';
	}
	out << " Damage evolution law code . . . . . . . . . . . = " << idamage << "\n";
	out << " Damage evolution law exponent . . . . . . . . . = " << falpha_exp << "\n";
	out << " Damage evolution law constant . . . . . . . . . = " << falpha_0 << "\n";
	if (idamage == 2)
	{
		out << " Damage evolution law lambda exponent. . . . . . = " << flambda_exp << "\n";
		out << " Damage evolution law lambda prefactor . . . . . = " << flambda_0 << "\n";
	}
	out << " Penetration stiffness multiplier. . . . . . . . = " << fpenalty   << '\n';
#endif
}

/* returns the number of variables computed for nodal extrapolation
 * during for element output, ie. internal variables. Returns 0
 * by default 
 */
int YoonAllen2DT::NumOutputVariables(void) const { return 3; }

void YoonAllen2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(3);
	labels[0] = "lambda";
	labels[1] = "lambda_dot";
	labels[2] = "alpha";
}

void YoonAllen2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#ifndef _FRACTURE_INTERFACE_LIBRARY_
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif	
#endif
	double u_t = jump_u[0];
	double u_n = jump_u[1];
	
	double l_t = u_t/fd_c_t;
	double l_n = u_n/fd_c_n;

	output[0] = sqrt(l_t*l_t + l_n*l_n); 
	output[1] = state[iNumRelaxTimes+2*knumDOF+6];
	output[2] = state[iNumRelaxTimes+2*knumDOF+7];
	
	double u_t_dot = (u_t - state[iNumRelaxTimes+6])/fTimeStep;
	double u_n_dot = (u_n - state[iNumRelaxTimes+7])/fTimeStep;
	double l_dot = (l_t*u_t_dot/fd_c_t+l_t*u_n_dot/fd_c_n)/output[0];
	if (l_dot > kSmall)
		output[2] += falpha_0*pow(output[0],falpha_exp);
	
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


