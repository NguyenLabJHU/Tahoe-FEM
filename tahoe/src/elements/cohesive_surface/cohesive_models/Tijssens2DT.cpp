/* $Id: Tijssens2DT.cpp,v 1.6 2001-12-17 00:15:51 paklein Exp $  */
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
	//	in >> fGamma_0; if (fGamma_0 < 0) throw eBadInputValue;
	in >> fDelta_0; if (fDelta_0 < 0) throw eBadInputValue;
	in >> fsigma_c; if (fsigma_c < 0) throw eBadInputValue;
	//        in >> ftau_c; if (ftau_c < 0) throw eBadInputValue;
	in >> fastar; if (fastar < 0) throw eBadInputValue;
        in >> ftemp; if (ftemp < 0) throw eBadInputValue;
	in >> fGroup; if (fGroup <= 0) throw eBadInputValue;
        //in >> fY; if (fY < 0) throw eBadInputValue;

	fA = fA_0/2.*exp(fQ_A/ftemp);
	fB = fB_0/6.*exp(fQ_B/ftemp);
        fc_1 /= fDelta_n_ccr;
	double root3 = sqrt(3.);
	ftau_c = fsigma_c/root3;
	fGamma_0 = fDelta_0*root3;
	fastar /= ftemp;

}

/* return the number of state variables needed by the model */
int Tijssens2DT::NumStateVariables(void) const { return 4*knumDOF; }

/* surface potential */ 
double Tijssens2DT::FractureEnergy(const ArrayT<double>& state) 
{
   	return state[6]; 
}

double Tijssens2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(jump_u)

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

	/* see if crazing has been initiated */
	double delta_tc, delta_nc, ddelta_tc;
//	if ((1.5*state[7] - fA + fB/state[7] - state[1]) > kSmall)
	if (state[1]*state[1]+fA*state[1] < -fA*state[0]+fB+state[0]*state[0])
        {
	        delta_tc = 0.;
		delta_nc = 0.;
		ddelta_tc = 0.;
	} 
	else
	{
	        if (state[5] <= fDelta_n_ccr) {
		  delta_tc = fGamma_0;
		  ddelta_tc = fGamma_0*fastar;
		  if (state[0] <= ftau_c)
		  {
		    delta_tc *= exp(-fastar*(ftau_c-state[0]))-exp(-fastar*(ftau_c+state[0]));
		    ddelta_tc *= exp(-fastar*(ftau_c-state[0]))+exp(-fastar*(ftau_c+state[0]));
		  }
		  delta_nc = fDelta_0;
		  if (state[1] <= fsigma_c)
		    delta_nc *= exp(-fastar*(fsigma_c-state[1]));
		}
		else
		{
		  delta_tc = 0.;
		  delta_nc = 0.;
		  ddelta_tc = 0.;
		}
	}

	/* calculate gap vector rates */
        double du_t = u_t - state[2];
        double du_n = u_n - state[3];
	double u_t_dot = (du_t)/fTimeStep;
        double u_n_dot = (du_n)/fTimeStep;

	/* calculate traction rate and evolve traction*/

        if (state[5] <= fDelta_n_ccr) {
	  double k_t = fk_t0*exp(-fc_1*state[5]);
          state[0] += k_t/(1+k_t*ddelta_tc*fTimeStep) * (u_t_dot - delta_tc) * fTimeStep;
	  state[1] += fk_n/(1+fk_n*delta_nc*fastar*fTimeStep) * (u_n_dot - delta_nc) * fTimeStep;
	  
	  /* update craze widths */
	  state[4] += delta_tc*fTimeStep;
	  state[5] += delta_nc*fTimeStep;
	}

	state[2] = u_t;
	state[3] = u_n;

	fTraction[0] = state[0];
	fTraction[1] = state[1];

	/* add to integrated energy */
        state[6] += fTraction[0]*du_t + fTraction[1]*du_n;
	state[7] = (state[0]+state[1])/3.;
	
	return fTraction;
}

/* potential stiffness */
const dMatrixT& Tijssens2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw eSizeMismatch;
	if (state.Length() != NumStateVariables()) throw eGeneralFail;
#endif

	double delta_nc, ddelta_tc;

	if ((1.5*state[7] - fA + fB/state[7] - state[1]) > kSmall)
        {
		delta_nc = 0.;
		ddelta_tc = 0.;
	} 
	else
	{
	        if (state[5] <= fDelta_n_ccr) {
		  ddelta_tc = fGamma_0*fastar;
		  if (state[0] <= ftau_c)
		    ddelta_tc *= exp(-fastar*(ftau_c-state[0]))+exp(-fastar*(ftau_c+state[0]));
		  delta_nc = fDelta_0;
		  if (state[1] <= fsigma_c)
		    delta_nc *= exp(-fastar*(fsigma_c-state[1]));
		}
		else
		{
		  delta_nc = 0.;
		  ddelta_tc = 0.;
		}
	}

	fStiffness[1] = fStiffness[2] = 0.;
	if (state[5] <= fDelta_n_ccr) 
	{
	  double k_t = fk_t0*exp(-fc_1*state[5]);
	  fStiffness[0] = k_t/(1+k_t*ddelta_tc*fTimeStep);
	}	
	else 
	  fStiffness[0] = fk_t0*exp(-fc_1);
	fStiffness[3] = fk_n/(1+fk_n*fastar*delta_nc*fTimeStep);
	
	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT Tijssens2DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw eSizeMismatch;
#endif

	if ((1.5*state[7] - fA + fB/state[7] - state[1]) > kSmall)
	         return Precritical;
	else 
	  //	         if (state[5] > fDelta_n_ccr) 
	  //            return Failed;
	  //     else
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
	//	out << " Tangential crazing rate constant. . . . . . . . = " << fGamma_0      << '\n';
	out << " Normal crazing rate constant. . . . . . . . . . = " << fDelta_0   << '\n';
	out << " Critical normal traction for crazing. . . . . . = " << fsigma_c  << '\n';
	//	out << " Critical tangential traction for crazing. . . . = " << ftau_c   << '\n';
	out << " Material parameter. . . . . . . . . . . . . . . = " << fastar*ftemp   << '\n';
	out << " Temperature . . . . . . . . . . . . . . . . . . = " << ftemp   << '\n';
	//	out << " Bulk Yield Stress . . . . . . . . . . . . . . . = " << fY << '\n';
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int Tijssens2DT::NumOutputVariables(void) const { return 1; }

void Tijssens2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Allocate(1);
	labels[0] = "sigma_m";
}

void Tijssens2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw eGeneralFail;
#endif	
	output[0] = state[7];
}

bool Tijssens2DT::NeedsNodalInfo(void) { return true; }

int Tijssens2DT::NodalQuantityNeeded(void) 
{ 
        return 2; 
}

double Tijssens2DT::ComputeNodalValue(const dArrayT& nodalRow) 
{
        return (nodalRow[0]+nodalRow[1])/3;
}

void Tijssens2DT::UpdateStateVariables(const dArrayT& IPdata, ArrayT<double>& state)
{
        state[7] = IPdata[0];
}

int Tijssens2DT::ElementGroupNeeded(void) 
{
	return fGroup-1;
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




