/* $Id: TiedPotentialT.cpp,v 1.4 2002-08-05 19:27:55 cjkimme Exp $  */
/* created: cjkimme (10/23/2001) */

#include "TiedPotentialT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionCodes.h"
#include "fstreamT.h"
#include "StringT.h"
#include "SecantMethodT.h"

/* class parameters */

using namespace Tahoe;

const int    knumDOF = 2;
const double kExpMax = 100;
double TiedPotentialT::fsigma_critical = 0.;

/* constructor */
TiedPotentialT::TiedPotentialT(ifstreamT& in, const double& time_step): 
	SurfacePotentialT(knumDOF),
	fTimeStep(time_step)
{
#pragma unused(time_step)

	in >> q; /*Traction plateau*/
	in >> r; /*Critical Length*/

	d_n = r; 
	fsigma_critical = q;
}

/*initialize state variables with values from the rate-independent model */
void TiedPotentialT::InitStateVariables(ArrayT<double>& state)
{
	state = 0.;
}

/* return the number of state variables needed by the model */
int TiedPotentialT::NumStateVariables(void) const { return 3; }

/* surface potential */ 
double TiedPotentialT::FractureEnergy(const ArrayT<double>& state) 
{
#pragma unused(state)

   	return .5*q*r;
}

double TiedPotentialT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw eSizeMismatch;
	if (state.Length() != NumStateVariables()) throw eSizeMismatch;
#endif

	return 0.;
}
	
/* traction vector given displacement jump vector */	
const dArrayT& TiedPotentialT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw eSizeMismatch;
	if (state.Length() != NumStateVariables()) throw eSizeMismatch;
	if (fTimeStep <= 0.0) {
#ifndef _TAHOE_FRACTURE_INTERFACE_
		cout << "\n TiedPotentialT::Traction: expecting positive time increment: "
		     << fTimeStep << endl;
#endif
		throw eBadInputValue;
	}
#endif

	if (state[0] != 1./* && state[1] < kSmall*/)
	{
		fTraction = 0.;
		if (state[0] == -10.)
		{
			state[0] = 1.;
			fTraction[1] = fsigma_critical*(1. - jump_u[1]/d_n);
		}
		return fTraction;
	}
	else
	{
		fTraction[0] = 0.;
		if (jump_u[1] < d_n)
			fTraction[1] = fsigma_critical*(1. - jump_u[1]/d_n); //Dugdale model
		else
			fTraction[1] = 0.;
		return fTraction;
	}

	double u_t = jump_u[0] + state[0];
	double u_n = jump_u[1] + state[1];
	
}

/* potential stiffness */
const dMatrixT& TiedPotentialT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw eSizeMismatch;
	if (state.Length() != NumStateVariables()) throw eGeneralFail;
#endif

	double u_t = jump_u[0] + state[0];
	double u_n = jump_u[1] + state[1]; 
	
	if (state[0] != 1.)
	{
		fStiffness = 0.;
		if (state[0] == -10.)
			fStiffness[3] = -fsigma_critical/d_n;
		return fStiffness;
	}
	else
	{
		fStiffness = 0.;
		fStiffness[3] = -fsigma_critical/d_n;
		return fStiffness;
	}

}

/* surface status */
SurfacePotentialT::StatusT TiedPotentialT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw eSizeMismatch;
#endif
       
	double u_t1 = jump_u[0];
	double u_t  = sqrt(u_t1*u_t1);
	double u_n  = jump_u[1];
	
	/* square box for now */
	if (u_n > d_n)
		return Failed;
	else //if (u_t > d_t || u_n > d_n)
		if (u_n > kSmall)
		return Critical;
	else
		return Precritical;
}

void TiedPotentialT::PrintName(ostream& out) const
{
#ifndef _TAHOE_FRACTURE_INTERFACE_
	out << "    TiedPotentialT (modified Xu-Needleman) 2D \n";
#endif
}

/* print parameters to the output stream */
void TiedPotentialT::Print(ostream& out) const
{
#ifndef _TAHOE_FRACTURE_INTERFACE_
	out << " Surface energy ratio (phi_t/phi_n). . . . . . . = " << q       << '\n';
	out << " Critical opening ratio (delta_n* /d_n). . . . . = " << r       << '\n';
	out << " Characteristic normal opening to failure. . . . = " << d_n     << '\n';
	out << " Characteristic tangential opening to failure. . = " << d_t     << '\n';
	out << " Mode I work to fracture (phi_n) . . . . . . . . = " << phi_n   << '\n';
	out << " Failure ratio (d_n/delta_n or d_t/delta_t). . . = " << r_fail   << '\n';
	out << " Penetration stiffness multiplier. . . . . . . . = " << fKratio << '\n';
#endif
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int TiedPotentialT::NumOutputVariables(void) const { return 1; }

void TiedPotentialT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Allocate(1);
	labels[0] = "state[1]";
}

void TiedPotentialT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw eGeneralFail;
#endif	

	output[0] = state[1];
}

bool TiedPotentialT::NeedsNodalInfo(void) { return true; }

int TiedPotentialT::NodalQuantityNeeded(void) 
{ 
        return kAverageCode; /*get stress tensor from bulk */ 
}

/*double TiedPotentialT::ComputeNodalValue(const dArrayT& nodalRow) 
{
       return (nodalRow[0]+nodalRow[1])/3;
}

void TiedPotentialT::UpdateStateVariables(const dArrayT& IPdata, ArrayT<double>& state)
{
  //state[7] = IPdata[0];
}
*/
int TiedPotentialT::ElementGroupNeeded(void) 
{
	return 0;
}

bool TiedPotentialT::InitiationQ(const double* sigma) 
{
#pragma unused(sigma)
	return sigma[1] >= fsigma_critical;
}

/*void TiedPotentialT::AllocateSpace(int MajorDim, int MinorDim) 
{
  	
  	nodal_stresses.Allocate(MajorDim,MinorDim);
  	
}*/

bool TiedPotentialT::CompatibleOutput(const SurfacePotentialT& potential) const
{
#ifdef __NO_RTTI__
#pragma unused(potential)
	return false;
#else
	const TiedPotentialT* pTH = dynamic_cast<const TiedPotentialT*>(&potential);
	return pTH != NULL;
#endif
}

