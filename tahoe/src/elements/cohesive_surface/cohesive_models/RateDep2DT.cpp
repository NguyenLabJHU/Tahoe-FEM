/* $Id: RateDep2DT.cpp,v 1.11 2002-10-23 00:18:03 cjkimme Exp $  */
/* created: cjkimme (10/23/2001) */

#include "RateDep2DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"
#include "fstreamT.h"
#include "StringT.h"
#include "SecantMethodT.h"

/* class parameters */

using namespace Tahoe;

const int knumDOF = 2;

/* constructor */
RateDep2DT::RateDep2DT(ifstreamT& in, const double& time_step): 
	SurfacePotentialT(knumDOF),
	fTimeStep(time_step)
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
	in >> L_2_b;
	in >> L_2_m;
  	in >> fslope;	
	
  	/* penetration stiffness */
	fK = fpenalty*fsigma_max/(fL_1*fd_c_n);

}

/*initialize state variables with values from the rate-independent model */
void RateDep2DT::InitStateVariables(ArrayT<double>& state)
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

	/* state variables below are with their corresponding parameters 
         * in TvergHutch2DT. state[n] has a flag in state[n+1] that
         * can be used to allow the traction routine to change the
         * parameter to reflect local opening rates. For instance,
         * state[6] represents the critical normal length scale and
         * state[7] = 0. indicates state[6] has not been modified
         * from its initialized (rate-independent) value while
         * state[7] = 1. indicates that the length-scale has
         * been modified as a function of rate. */
	state[0] = fL_1;
	state[2] = fL_2;
	state[4] = fsigma_max;
	state[6] = fd_c_n;
	state[7] = 0.;
	/* state[8] = previous timestep's u_t
	 * state[9] = previous timestep's u_n
	 */
  	state[10] = fd_c_t;
}

/* return the number of state variables needed by the model */
int RateDep2DT::NumStateVariables(void) const { return 11; }

/* surface potential */ 
double RateDep2DT::FractureEnergy(const ArrayT<double>& state) 
{
   	return state[6]*state[4]*0.5*(1 - state[1] + state[2]); 
}

double RateDep2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/state[10];
	double r_n = u_n/state[6];
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)

	double u;
	if (L < state[0])
		u = state[4]*0.5*(L/state[0])*L;
	else if (L < state[2])
		u = state[4]*(0.5*fL_1 + (L - state[0]));
	else if (L < 1)
	{
		double z1 = (1.0 - state[2]);
		double z2 = (1.0 - L);
		u = state[4]*(0.5*state[0] + (state[2]- state[0]) + 0.5*(z1 - (z2/z1)*z2));
	}
	else
		u = state[4]*0.5*(1 - state[0] + state[2]);

	/* penetration */
	if (u_n < 0) u += 0.5*u_n*u_n*fK/state[6];

	return u*state[6]; // (1.2)
}
	
/* traction vector given displacement jump vector */	
const dArrayT& RateDep2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
	if (fTimeStep <= 0.0) {
#ifndef _SIERRA_TEST_	
		cout << "\n RateDep2DT::Traction: expecting positive time increment: "
		     << fTimeStep << endl;
#endif
		throw ExceptionT::kBadInputValue;
	}
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/state[10];
	double r_n = u_n/state[6];
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)

	double sigbyL;
	if (L < state[0])
		sigbyL = state[4]/state[0];
	else if (L < state[2]) /*L > state[0] means we're at the plateau stress */
	{ 
	  	if (state[7] == 0.) 
	  	{ 
	    	double u_n_dot = (u_n-state[9])/fTimeStep;
	      	if (u_n_dot > kSmall) 
	      	{
				state[7] = 1.;
				state[6] = L_2_b + L_2_m * log(u_n_dot);
				/* make sure new length scale is greater than current
		 		 * gap vector 
                 */
	        	state[0] *= fd_c_n/state[6];
	        	r_n = u_n/state[6];
	        	L = sqrt(r_t*r_t+r_n*r_n);
	        	sigbyL = state[4]*(1+fslope*(L-state[0]))/L;
				if (state[6] < u_n || state[6] < fd_c_n*fL_1)
				{
#ifndef _SIERRA_TEST_
		  			cout <<  "\n RateDep2DT::Traction: rate-dependent length scale " << state[6] << " is incompatible with rate-independent one. Check input parameters. \n ";
#endif
	          		state[6] = fd_c_n;
	          		state[2] = state[0]; /* start unloading now */
	          		r_n = u_n/state[6];
	          		L = sqrt(r_t*r_t+r_n*r_n);
                  	sigbyL = state[4]*(1-L)/(1-state[2])/L;
				}
	      	}
	  	}
	  	else 
	    	sigbyL = state[4]*(1+fslope*(L-state[0]))/L;
	}
	else if (L < 1)
		sigbyL = state[4]*(1+fslope*(state[2]-state[0]))*(1 - L)/(1 - state[2])/L;
	else
		sigbyL = 0.0;	

	fTraction[0] = sigbyL*(u_t/state[10])*(state[6]/state[10]);
	fTraction[1] = sigbyL*(u_n/state[6]);

	/* penetration */
	if (u_n < 0) fTraction[1] += fK*u_n;

	state[8] = u_t;
	state[9] = u_n;

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

//	double z1, z2, z3, z4, z5, z6;	

	double dtm2 = 1./state[10]/state[10];
	double dnm2 = 1./state[6]/state[6];
	double L = sqrt(u_t*u_t*dtm2 + u_n*u_n*dnm2);
	
	/* effective opening */
/*	z1 = u_n*u_n;
	z2 = 1./(state[6]*state[6]);
	z3 = u_t*u_t;
	z4 = 1./(state[10]*state[10]);
	z5 = z1*z2;
	z6 = z3*z4;
	z5 += z6;
	double L = sqrt(z5);
*/
	if (L < state[0]) // K1
	{
		fStiffness[0] = (state[6]/state[10])*state[4]/(state[0]*state[10]);
		fStiffness[1] = 0.0;
		fStiffness[2] = 0.0;
		fStiffness[3] = state[4]/(state[0]*state[6]);
	}
	else 
	{
		double lt_0 = u_t*dtm2;
		double lt_2 = u_n*dnm2;
		
		if (L < state[2]) // K2
		{
			double dijTerm = state[4]/L*state[6]*(1+fslope*(L-state[0]));
			
			fStiffness[0] = dijTerm*dtm2;
			fStiffness[3] = dijTerm*dnm2;
			
			dijTerm = -state[4]/L/L/L*state[6]*(1-fslope*state[0]);
			
			fStiffness[0] += dijTerm*lt_0*lt_0;
			fStiffness[2] = fStiffness[1] = dijTerm*lt_0*lt_2;
			fStiffness[3] += dijTerm*lt_2*lt_2;
				
			/*z5 = 1./state[6];
			z2 *= z1;
			z3 *= z4;
			z2 += z3;
			z2 = pow(z2,-1.5);
			z2 *= state[4]*(1+fslope*(L-state[0]))*z5;
			z3 *= z2;
			z2 *= z4;
			z4 = -u_n*u_t*z2;
			z1 *= z2;

			//cout << "compare: " << z1 << " " << fStiffness[0] <<" "<<z3 <<  " " <<fStiffness[3]<<"\n";
			//cout << "line 2 : " << z4 << " " << fStiffness[2] << " \n";
			//{{z1, z4},
			// {z4, z3}}
			fStiffness[0] = z1;
			fStiffness[1] = z4;
			fStiffness[2] = z4;
			fStiffness[3] = z3;*/
		}
		else 
			if (L < 1 && state[2] < 1.) // K3
			{
				/*double z7, z8, z9, z10, z11, z12, z13, z14, z15, z16;

				z5 = -1. + state[2];
				z6 = -state[2];
				z7 = pow(state[6],-3.);
				z8 = 1./state[6];
				z9 = state[6]*state[6];
				z10 = pow(state[6],3.);
				z11 = pow(state[10],-4.);
				z12 = state[10]*state[10];
				z2 = z1*z2;
				z13 = z3*z4;
				z6 += 1.;
				z12 *= z1;
				z2 += z13;
				z9 *= z3;
				z5 = 1./z5;
				z6 = 1./z6;
				z14 = pow(z2,-1.5);
				z15 = 1./z2;
				z16 = 1./sqrt(z2);
				z2 = sqrt(z2);
				z9 += z12;
				z12 = -state[4]*(1+fslope*(state[2]-state[0]))*z1*z15*z6*z7;
				z2 = -z2;
				z9 = 1./z9;
				z2 += state[2];
				z5 *= state[4]*(1+fslope*(state[2]-state[0]))*z9;
				z9 = u_n*state[6]*u_t*z5;
				z5 *= z10*z13;
				z2 *= z6;
				z2 += 1.;
				z2 *= state[4]*(1+fslope*(state[2]-state[0]));
				z1 *= -z14*z2*z7;
				z6 = z16*z2*z8;
				z3 *= -state[6]*z11*z14*z2;
				z7 = -u_n*u_t*z14*z2*z4*z8;
				z2 *= state[6]*z16*z4;
				z1 += z12 + z6;
				z4 = z7 + z9;
				z2 += z3 + z5;
			
				//{{z2, z4},
				// {z4, z1}}
				fStiffness[0] = z2;
				fStiffness[1] = z4;
				fStiffness[2] = z4;
				fStiffness[3] = z1;
				*/
				double dijTerm = state[4]*(1.+fslope*(state[2]-state[0]))*(1./L-1.)/(1.-state[2])*state[6];
				fStiffness[0] = dijTerm*dtm2;
				fStiffness[3] = dijTerm*dnm2;
				
				dijTerm = -state[4]*(1.+fslope*(state[2]-state[0]))/(1.-state[2])*state[6]/L/L/L;
				
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

	double r_t = u_t/state[10];
	double r_n = u_n/state[6];
	double L = sqrt(r_t*r_t + r_n*r_n); // (1.1)
	
	if (L > fL_fail)
		return Failed;
	else if (L > state[0])
		return Critical;
	else
		return Precritical;

}

void RateDep2DT::PrintName(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << "    RateDep 2D \n";
#endif
}

/* print parameters to the output stream */
void RateDep2DT::Print(ostream& out) const
{
#ifndef _SIERRA_TEST_
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

	double r_t = u_t/state[10];
	double r_n = u_n/state[6];
	output[0]  = sqrt(r_t*r_t + r_n*r_n); // (1.1)
	output[1] = (u_t-state[8])/fTimeStep;
	output[2] = (u_n-state[9])/fTimeStep;
	output[3] = state[6];

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


