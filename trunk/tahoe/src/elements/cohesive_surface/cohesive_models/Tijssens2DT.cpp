/* $Id: Tijssens2DT.cpp,v 1.23 2004-06-17 07:13:28 paklein Exp $  */
/* created: cjkimme (10/23/2001) */

#include "Tijssens2DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"
#include "ifstreamT.h"
#include "StringT.h"
#include "SecantMethodT.h"

using namespace Tahoe;

/* class parameters */
const int knumDOF = 2;
const int ku_c = 5;
const int kIntShift = 6;

/* constructor */
Tijssens2DT::Tijssens2DT(ifstreamT& in, const double& time_step): 
	SurfacePotentialT(knumDOF),
	fTimeStep(time_step)
{
	/* traction rate parameters */
	in >> fk_t0; if (fk_t0 < 0) throw ExceptionT::kBadInputValue;
	in >> fk_n; if (fk_n < 0) throw ExceptionT::kBadInputValue;
	in >> fc_1; if (fc_1 < 0) throw ExceptionT::kBadInputValue;
	in >> fDelta_n_ccr; if (fDelta_n_ccr < 0) throw ExceptionT::kBadInputValue;

	/* craze initiation parameters */
	in >> fA_0; if (fA_0 < 0) throw ExceptionT::kBadInputValue;
	in >> fB_0; if (fB_0 < 0) throw ExceptionT::kBadInputValue;
	in >> fQ_A; if (fQ_A < 0) throw ExceptionT::kBadInputValue;
	in >> fQ_B; if (fQ_B < 0) throw ExceptionT::kBadInputValue;
	
	/* crazing state variables' parameters */
	in >> fDelta_0; if (fDelta_0 < 0) throw ExceptionT::kBadInputValue;
	in >> fsigma_c; if (fsigma_c < 0) throw ExceptionT::kBadInputValue;
	in >> fastar; if (fastar < 0) throw ExceptionT::kBadInputValue;
	in >> ftemp; if (ftemp < 0) throw ExceptionT::kBadInputValue;
	in >> fGroup; if (fGroup <= 0) throw ExceptionT::kBadInputValue;
	in >> fSteps; if (fSteps < 0) throw ExceptionT::kBadInputValue;

	fA = fA_0/2.*exp(fQ_A/ftemp);
	fB = fB_0/6.*exp(fQ_B/ftemp);
 	fc_1 /= fDelta_n_ccr;
	double root3 = sqrt(3.);
	ftau_c = fsigma_c/root3;
	fGamma_0 = fDelta_0*root3;
	fastar /= ftemp;
}

/* return the number of state variables needed by the model */
int Tijssens2DT::NumStateVariables(void) const { return 6*knumDOF+1; }
/* State variable array: 
 * 0..knumDOF - 1 are components of previous step's traction vector
 * knumDOF .. 2 knumDOF - 1 are components of previous step's gap vector
 * 2 knumDOF .. 3 knumDOF - 1 are components of previous step's craze vector
 * Then, they all repeat again to store the previous timestep's after
 * integration.
 * Final one is integral of T dot dDelta
 */

/* surface potential */ 
double Tijssens2DT::FractureEnergy(const ArrayT<double>& state) 
{
	return state[6*knumDOF]; 
}

double Tijssens2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(jump_u)

	/* There isn't a potential to return. Just return the work done so far */
	return state[6*knumDOF];
}
	
/* traction vector given displacement jump vector */	
const dArrayT& Tijssens2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
	if (fTimeStep < 0.0) {
#ifndef _FRACTURE_INTERFACE_LIBRARY_
		cout << "\n Tijssens2DT::Traction: expecting non-negative time increment: "
		     << fTimeStep << endl;
#endif
		throw ExceptionT::kBadInputValue;
	}
#endif

	if (!qIntegrate)
	{
		fTraction[0] = state[0];
		fTraction[1] = state[1];
		
		return fTraction;
	}
	else
	{
		double du_t = jump_u[0]-state[2];
		double du_n = jump_u[1]-state[3];
		
		/* Store previous step's data */
		/* (You can practically see the for loop writing itself here) */
		state[kIntShift] = state[0];
		state[kIntShift+1] = state[1];
		state[kIntShift+2] = state[2];
		state[kIntShift+3] = state[3];
		state[kIntShift+4] = state[4];
		state[kIntShift+5] = state[5];

		/* see if crazing has been initiated */
		//if (state[7] < kSmall || (state[8] <= kSmall && (1.5*state[7] - fA + fB/state[7] - state[1]) > kSmall) || (state[5] >= fDelta_n_ccr && state[9] <= kSmall))
		if (jump_u[1] < 1.01*fsigma_c/fk_n) 
		{
			state[1] += fk_n*du_n;
			state[0] += fk_t0*du_t;

			/* interpenetration */
			if (jump_u[1] < kSmall)
				state[1] += 2.*fk_n*du_n;
		}
		else 
			if (state[5] > fDelta_n_ccr)
			{
				if (state[1] - fk_n/fSteps*du_n > kSmall)
				{
					state[1] -= fk_n/fSteps*du_n;
					state[0] -= fk_t0*exp(-fc_1)*du_t;
				}
				else
					state[1] = state[0] = 0.;
		  	}
		  	else
		  	{ 
		    	/*NormalTraction*/
				if (state[1] < kSmall) 
					state[1] = 1.1*fsigma_c;
				double Tnp1 = state[1];
				SecantMethodT secant(20);
				double du_nd = 0.0;
				if (fabs(fTimeStep) > kSmall) du_nd = du_n/fTimeStep;

				secant.Reset(fsigma_c,-state[1]-fk_n*fTimeStep*(du_nd-fDelta_0),1.5*state[1],.5*state[1]-fk_n*fTimeStep*(du_nd-fDelta_0*exp(-fastar*(fsigma_c-Tnp1))));
				Tnp1 = secant.NextGuess();
				while (!secant.NextPoint(Tnp1,Tnp1-state[1]-(fk_n*fTimeStep*(du_nd-fDelta_0*exp(-fastar*(fsigma_c-Tnp1))))))
					Tnp1 = secant.NextGuess();

				double du_c = fTimeStep*fDelta_0*exp(-fastar*(fsigma_c-Tnp1));
				state[1] += fk_n*(du_nd*fTimeStep-du_c);

				state[5] += du_c;
				state[6*knumDOF] += state[1]*du_n;

			    /* Tangential traction *//*
			    Tnp1 = state[0];
		    	double fk_t = fk_t0*exp(-fc_1*state[4]/fDelta_n_ccr);
		    	secant.Reset(0,-state[0]-fk_t*(du_t-fTimeStep*fGamma_0*(exp(-fastar*(ftau_c-Tnp1))-exp(-fastar*(ftau_c+Tnp1)))),1.5*state[0],.5*state[0]-fk_t*(du_t-fTimeStep*fGamma_0*(exp(-fastar*(ftau_c-Tnp1))-exp(-fastar*(ftau_c+Tnp1)))));
		    	Tnp1 = secant.NextGuess();
		    	while (!secant.NextPoint(Tnp1,Tnp1-state[0]-fk_t*(du_t-fTimeStep*fGamma_0*(exp(-fastar*(ftau_c-Tnp1))-exp(-fastar*(ftau_c+Tnp1))))))
		      		Tnp1 = secant.NextGuess();

		      	state[0] += fk_t*(du_t-fTimeStep*fGamma_0*(exp(-fastar*(ftau_c-Tnp1))-exp(-fastar*(ftau_c+Tnp1))));*/

		    	/*if a craze will fail, request a smaller timestep*/
				/*not implemented */
			
				double dw_t = fk_t0*exp(-fc_1*state[5]/fDelta_n_ccr)*du_t;
				state[0] += dw_t;
				state[6*knumDOF] += state[0]*du_t;
			}

		fTraction[0] = state[0];
		fTraction[1] = state[1];
		state[2] = jump_u[0];
		state[3] = jump_u[1];

		return fTraction;
	}
	
}

/* potential stiffness */
const dMatrixT& Tijssens2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	fStiffness[1] = fStiffness[2] = 0.;
	if (jump_u[1] <= 1.01*fsigma_c/fk_n) 
	{
		fStiffness[3] = fk_n;
		fStiffness[0] = fk_t0;

		/* interpenetration */
		if (jump_u[1] < kSmall)
			fStiffness[3] += 2.*fk_n;
	}
	else 
	{
		const double *state_old = state.Pointer(kIntShift);
		if (state_old[5] > fDelta_n_ccr) 
		{
			if (state_old[1] - fk_n/fSteps*(jump_u[1]-state_old[3]) > kSmall)
			{
				fStiffness[3] = -fk_n/fSteps;
				fStiffness[0] = -fk_t0*exp(-fc_1);
			}
			else 
				fStiffness[3] = fStiffness[0] = 0.;
		}
		else
		{
			/* Normal stiffness */	
			double du_n = jump_u[1]-state_old[3];
		
			//if (state_old[1] < kSmall)
			//	state_old[1] = 1.1 * fsigma_c;
				
			/* Grab current traction */
			double Tnp1 = state[1];
	    	
			if (state_old[1] < 1.01*fsigma_c && du_n > kSmall)
			{ 
				fStiffness[3] = (Tnp1-state_old[1])/(jump_u[1]-state_old[3]);
			}
			else
			{
				fStiffness[3] = fk_n/(1.+fk_n*fTimeStep*fastar*fDelta_0*exp(-fastar*(fsigma_c-Tnp1)));
			}
	      	
			double fk_t = fk_t0*exp(-fc_1*state[5+kIntShift]/fDelta_n_ccr);
			fStiffness[0] = fk_t;
		}
	}
	  
	return fStiffness;

}

/* surface status */
SurfacePotentialT::StatusT Tijssens2DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif
       
	//	if ((1.5*state[7] - fA + fB/state[7] - state[1]) > kSmall)
	if (state[1+kIntShift] < 1.1*fsigma_c)
	         return Precritical;
	else 
	  	if (state[5+kIntShift] >= fDelta_n_ccr) 
	    	return Failed;
	  	else
	    	return Critical;

}

void Tijssens2DT::PrintName(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	out << "    Tijssens 2D \n";
#endif
}

/* print parameters to the output stream */
void Tijssens2DT::Print(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	out << " Initial tangential stiffness. . . . . . . . . . = " << fk_t0 << '\n';
	out << " Normal stiffness . . . .  . . . . . . . . . . . = " << fk_n     << '\n';
	out << " Tangential stiffness rate constant. . . . . . . = " << fDelta_n_ccr*fc_1     << '\n';
	out << " Critical crazing width. . . . . . . . . . . . . = " << fDelta_n_ccr << '\n';
	out << " Crazing initiation parameter A. . . . . . . . . = " << fA_0    << '\n';
	out << " Crazing initiation parameter B. . . . . . . . . = " << fB_0    << '\n';
	out << " Thermal activation for A. . . . . . . . . . . . = " << fQ_A << '\n';
	out << " Thermal activation for B. . . . . . . . . . . . = " << fQ_B << '\n';
	out << " Normal crazing rate constant. . . . . . . . . . = " << fDelta_0/exp(-fastar*fsigma_c)   << '\n';
	out << " Critical normal traction for crazing. . . . . . = " << fsigma_c  << '\n';
	out << " Material parameter. . . . . . . . . . . . . . . = " << fastar*ftemp   << '\n';
	out << " Temperature . . . . . . . . . . . . . . . . . . = " << ftemp   << '\n';
#endif
}

int Tijssens2DT::NumOutputVariables(void) const { return 4; }

void Tijssens2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(4);
	labels[0] = "sigma_m";
	labels[1] = "Delta_t_c";
	labels[2] = "Delta_n_c";
	labels[3] = "f(sigma_m)";
}

void Tijssens2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif	
	output[0] = (jump_u[1]-state[9])/fTimeStep;
	output[1] = state[10];
	output[2] = state[11];
	output[3] = 0;
}

bool Tijssens2DT::NeedsNodalInfo(void) { return false; }

int Tijssens2DT::NodalQuantityNeeded(void) 
{ 
        return 2; 
}

/*double Tijssens2DT::ComputeNodalValue(const dArrayT& nodalRow) 
{
        return (nodalRow[0]+nodalRow[1])/3;
}

void Tijssens2DT::UpdateStateVariables(const dArrayT& IPdata, ArrayT<double>& state)
{
        state[7] = IPdata[0];
}*/

void Tijssens2DT::SetElementGroupsNeeded(iArrayT& iGroups) 
{	
	iGroups[0] = 1;
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

