/* $Id: Tijssens2DT.cpp,v 1.1 2001-10-25 22:22:12 cjkimme Exp $  */
/* created: cjkimme (10/23/2001) */

#include "Tijssens2DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionCodes.h"
#include "fstreamT.h"
#include "StringT.h"

/* class parameters */
const int knumDOF = 2;

/* constructor */
Tijssens2DT::Tijssens2DT(ifstreamT& in, const double& time_step): 
	SurfacePotentialT(knumDOF),
	fTimeStep(time_step)
{
	/* traction rate parameters */
        in >> fk_t0; if (fk_t0 < 0) throw eBadInputValue;
	in >> fk_n; if (fk_n < 0) throw eBadInputValue;
	in >> fc_1; if (fc_1 < 0) throw eBadInputValue;
	in >> fDelta_n_ccr; if (fDelta_n_ccr < 0) throw eBadInputValue;

	/* craze initiation parameters */
	in >> fA_0; if (fA < 0) throw eBadInputValue;
	in >> fB_0; if (fB < 0) throw eBadInputValue;
	in >> fQ_A; if (fQ_A < 0) throw eBadInputValue;
	in >> fQ_B; if (fQ_B < 0) throw eBadInputValue;

        /* crazing state variables' parameters */
	in >> fLambda_0; if (fLambda_0 < 0) throw eBadInputValue;
	in >> fDelta_0; if (fDelta_0 < 0) throw eBadInputValue;
	in >> fsigma_c; if (fsigma_c < 0) throw eBadInputValue;
        in >> ftau_c; if (ftau_c < 0) throw eBadInputValue;
	in >> fastar; if (fastar < 0) throw eBadInputValue;
        in >> ftemp; if (ftemp < 0) throw eBadInputValue;

	fA = fA_0/2.*exp(fQ_A/ftemp);
	fB = fB_0/6.*exp(fQ_B/ftemp);
        fc_1 /= fDelta_n_ccr;

}

/* return the number of state variables needed by the model */
int Tijssens2DT::NumStateVariables(void) const { return 3*knumDOF+1; }

/* surface potential */
double Tijssens2DT::FractureEnergy(void) 
{
  return 0.0;
  //	return state[6]; 
}

double Tijssens2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
        return state[6];
}
	
/* traction vector given displacement jump vector */	
const dArrayT& Tijssens2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw eSizeMismatch;
	if (state.Length() != NumStateVariables()) throw eSizeMismatch;
	if (fTimeStep <= 0.0) {
		cout << "\n Tijssens2DT::Traction: expecting positive time increment: "
		     << fTimeStep << endl;
		throw eBadInputValue;
	}
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];
	double sigma_m = state[1];
	// I need the trace of the stress tensor
	/* see if crazing has been initiated */
	double delta_tc, delta_nc;
	if ((1.5*sigma_m - fA + fB/sigma_m - state[1]) > 0.)
        {
	        delta_tc = 0.;
		delta_nc = 0.;
	} 
	else
	{
	        if (state[5] <= fDelta_n_ccr) {
		  delta_tc *= fLambda_0;
		  if (state[0] <= ftau_c)
		    delta_tc *= exp(-fastar*ftau_c/ftemp*(1-state[0]/ftau_c))-exp(-fastar*ftau_c/ftemp*(1+state[0]/ftau_c));
		  delta_nc = fDelta_0;
		  if (state[1] <= fsigma_c)
		    delta_nc *= exp(-fastar*fsigma_c/ftemp*(1-state[1]/fsigma_c));
		}
		else
		{
		  delta_tc = 0.;
		  delta_nc = 0.;
		}
	}

	/* calculate gap vector rates */
        double du_t = u_t - state[2];
        double du_n = u_n - state[3];
	double u_t_dot = (du_t)/fTimeStep;
        double u_n_dot = (du_n)/fTimeStep;

	/* calculate traction rate and evolve traction*/

        if (state[5] <= fDelta_n_ccr) {
          state[0] += fk_t0*exp(-fc_1*state[5])* (u_t_dot - delta_tc) * fTimeStep;
	  /* update craze widths */
	  state[4] += delta_tc*fTimeStep;
	  state[5] += delta_nc*fTimeStep;
	}
	state[1] += fk_n * (u_n_dot - delta_nc) * fTimeStep;

	state[2] = u_t;
	state[3] = u_n;
	/* penetration */ //What do I do here?
	//	if (u_n < 0) fTraction[1] += fK*u_n;

	fTraction[0] = state[0];
	fTraction[1] = state[1];

	/* add to integrated energy */
        state[6] += fTraction[0]*du_t + fTraction[1]*du_n;
	
	return fTraction;
}

/* potential stiffness */
const dMatrixT& Tijssens2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw eSizeMismatch;
	if (state.Length() != NumStateVariables()) throw eGeneralFail;
#endif

	double u_t = jump_u[0];
	double u_n = jump_u[1];

	fStiffness[1] = fStiffness[2] = 0.;
	if (state[5] <= fDelta_n_ccr) 
	  fStiffness[0] = fk_t0*exp(-fc_1*state[4]);
	else 
	  fStiffness[0] = 0.;
	fStiffness[3] = fk_n;
	
	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT Tijssens2DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw eSizeMismatch;
#endif
	double sigma_m = state[1];

	if ((1.5*sigma_m - fA + fB/sigma_m - state[1]) > 0.)
	         return Precritical;
	else 
	         if (state[5] > fDelta_n_ccr) 
		        return Failed;
		 else
		        return Critical;

}

void Tijssens2DT::PrintName(ostream& out) const
{
	out << "    Tijssens 2D \n";
}

/* print parameters to the output stream */
void Tijssens2DT::Print(ostream& out) const
{
	out << " Initial tangential stiffness. . . . . . . . . . = " << fk_t0 << '\n';
	out << " Normal stiffness . . . .  . . . . . . . . . . . = " << fk_n     << '\n';
	out << " Tangential stiffness rate constant. . . . . . . = " << fDelta_n_ccr*fc_1     << '\n';
	out << " Critical crazing width. . . . . . . . . . . . . = " << fDelta_n_ccr << '\n';
	out << " Crazing initiation parameter A. . . . . . . . . = " << fA_0    << '\n';
	out << " Crazing initiation parameter B. . . . . . . . . = " << fB_0    << '\n';
	out << " Thermal activation for A. . . . . . . . . . . . = " << fQ_A << '\n';
	out << " Thermal activation for B. . . . . . . . . . . . = " << fQ_B << '\n';
	out << " Tangential crazing rate constant. . . . . . . . = " << fLambda_0      << '\n';
	out << " Normal crazing rate constant. . . . . . . . . . = " << fDelta_0   << '\n';
	out << " Critical normal traction for crazing. . . . . . = " << fsigma_c  << '\n';
	out << " Critical tangential traction for crazing. . . . = " << ftau_c   << '\n';
	out << " Material parameter. . . . . . . . . . . . . . . = " << fastar   << '\n';
	out << " Temperature. . . . . . . . . . .. . . . . . . . = " << ftemp   << '\n';
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int Tijssens2DT::NumOutputVariables(void) const { return 2; }

void Tijssens2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Allocate(2);
	labels[0] = "lambda";
	labels[1] = "dw_visc";
}

void Tijssens2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw eGeneralFail;
#endif
/*
	double u_t = jump_u[0];
	double u_n = jump_u[1];

	double r_t = u_t/fd_c_t;
	double r_n = u_n/fd_c_n;
	double   L = sqrt(r_t*r_t + r_n*r_n); // (1.1)
	output[0] = L;
*/	
	/* approximate incremental viscous dissipation, assuming a constant
	 * average viscosity over the increment */
/*	if (L < 1)
	{*/
		/* increment displacement */
	/*	double d_t = u_t - state[0];
		double d_n = u_n - state[1];
*/
		/* previous lambda */
/*		r_t = state[0]/fd_c_t;
		r_n = state[1]/fd_c_n;
		double L_last = sqrt(r_t*r_t + r_n*r_n); // (1.1)
*/		
		/* average viscosity */
/*		double eta = feta0*(1.0 - 0.5*(L + L_last));
		
*/		/* approximate incremental dissipation */
/*		output[1] = 0.5*eta*(d_t*d_t + d_n*d_n)/fTimeStep;
	}
	else
		output[1] = 0.0;*/
	output[1] = 0.0;
}

bool Tijssens2DT::CompatibleOutput(const SurfacePotentialT& potential) const
{
#ifdef __NO_RTTI__
#pragma unused(potential)
	return false;
#else
	const Tijssens2DT* pTH = dynamic_cast<const Tijssens2DT*>(&potential);
	return pTH != NULL;
#endif
}




