/* $Id: From2Dto3DT.cpp,v 1.4 2004-06-17 07:13:28 paklein Exp $ */
/* created: paklein (06/23/1999)*/
#include "From2Dto3DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"
#include "ifstreamT.h"

#include "XuNeedleman2DT.h"
#include "TvergHutch2DT.h"
#include "ViscTvergHutch2DT.h"
#include "YoonAllen2DT.h"

using namespace Tahoe;

/* class parameters */
const int    knumDOF = 3;

#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* constructor */
From2Dto3DT::From2Dto3DT(ifstreamT& in, int code, const double& time_step): SurfacePotentialT(knumDOF)
{
	switch(code)
	{
		case kXuNeedleman:
		{
			f2DModel = new XuNeedleman2DT(in);
			break;
		}
		case kTvergaardHutchinson:
		{
			f2DModel = new TvergHutch2DT(in);
			break;
		}
		case kViscTvergaardHutchinson:
		{
			f2DModel = new ViscTvergHutch2DT(in, time_step);
			break;
		}
		case kYoonAllen:
		{
			f2DModel = new YoonAllen2DT(in, time_step);
			break;
		}
		default:
		{
			ExceptionT::BadInputValue("From2Dto3DT::From2Dto3DT","Cohesive model code not supported\n");
		}
	}
}
#endif

From2Dto3DT::From2Dto3DT(dArrayT& params): SurfacePotentialT(knumDOF)
{
#pragma unused(params)
	ExceptionT::GeneralFail("From2Dto3DT::From2Dto3DT","Constructor not implemented yet");
//	f2DModel = new XuNeedleman2DT(params);
}

From2Dto3DT::~From2Dto3DT(void)
{
	delete f2DModel;
}

int From2Dto3DT::NumStateVariables(void) const
{
	return f2DModel->NumStateVariables();
}

void From2Dto3DT::InitStateVariables(ArrayT<double>& state)
{	
	f2DModel->InitStateVariables(state);
}

double From2Dto3DT::FractureEnergy(const ArrayT<double>& state) 
{
	return f2DModel->FractureEnergy(state); 
}

double From2Dto3DT::IncrementalHeat(const dArrayT& jump, const ArrayT<double>& state)
{
	dArrayT jump_2D(2);
	jump_2D[0] = sqrt(jump[0]*jump[0] + jump[1]*jump[1]);
	jump_2D[1] = jump[2];

	return f2DModel->IncrementalHeat(jump_2D, state);
}

double From2Dto3DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;	
#endif

	dArrayT new_jump_u(2);
	new_jump_u[0] = sqrt(jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1]);
	new_jump_u[1] = jump_u[2];
	
	return f2DModel->Potential(new_jump_u,state);
}

/* traction vector given displacement jump vector */	
const dArrayT& From2Dto3DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(state)
#pragma unused(sigma)
#pragma unused(qIntegrate)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif
	
	dArrayT new_jump_u(2);
	new_jump_u[0] = sqrt(jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1]);
	new_jump_u[1] = jump_u[2];
	
	dArrayT new_traction;
	new_traction =  f2DModel->Traction(new_jump_u, state, sigma, qIntegrate);
	
	if (fabs(jump_u[0]) > kSmall)
		fTraction[0] = jump_u[0]/new_jump_u[0]*new_traction[0];
	else
		fTraction[0] = 0.;
	if (fabs(jump_u[1]) > kSmall)
		fTraction[1] = jump_u[1]/new_jump_u[0]*new_traction[0];
	else
		fTraction[1] = 0.;
	
	fTraction[2] = new_traction[1];
	
	return fTraction;
}


/* potential stiffness */
const dMatrixT& From2Dto3DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#pragma unused(state)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif
	
	dArrayT new_jump_u(2);
	new_jump_u[0] = sqrt(jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1]);
	new_jump_u[1] = jump_u[2];
	
	dMatrixT new_stiffness;
	new_stiffness =  f2DModel->Stiffness(new_jump_u, state, sigma);
	
	double f0, f1;
	if (fabs(jump_u[0]) > kSmall)
		f0 = jump_u[0]/new_jump_u[0];
	else
		f0 = 0.;
	if (fabs(jump_u[1]) > kSmall)
		f1 = jump_u[1]/new_jump_u[0];
	else
		f1 = 0.;
		
	if (fabs(new_jump_u[0]) < kSmall)
		f0 = f1 = 1/sqrt(2.);
		
	fStiffness(0,0) = fStiffness(1,1) = fStiffness(0,1) = fStiffness(1,0) = new_stiffness(0,0);
	fStiffness(2,0) = fStiffness(2,1) = new_stiffness(1,0);
	fStiffness(0,2) = fStiffness(1,2) = new_stiffness(0,1);
	fStiffness(2,2) = new_stiffness(1,1);

	/* scale by chain rule */
	fStiffness(0,0) *= f0*f0;
	fStiffness(1,1) *= f1*f1;
	fStiffness(0,1) *= f0*f1;
	fStiffness(1,0) *= f0*f1;
	fStiffness(2,0) *= f0;
	fStiffness(0,2) *= f0;
	fStiffness(2,1) *= f1;
	fStiffness(1,2) *= f1;
	
	if (fabs(new_jump_u[0]) > kSmall)
	{	
		/* Add traction-dependent terms */
		double T2D_t;
		dArrayT new_state(state.Length(), state.Pointer());
		const dArrayT& tract = f2DModel->Traction(new_jump_u, new_state, sigma, false);
		T2D_t =  tract[0]/new_jump_u[0];
		
		fStiffness(0,0) += (1. - f0*f0)*T2D_t;
		fStiffness(1,1) += (1. - f1*f1)*T2D_t;
		T2D_t *= f0*f1;
		fStiffness(0,1) -= T2D_t;
		fStiffness(1,0) -= T2D_t;
	}
	
	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT From2Dto3DT::Status(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(state)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	dArrayT new_jump_u(2);
	new_jump_u[0] = sqrt(jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1]);
	new_jump_u[1] = jump_u[2];
	
	return f2DModel->Status(new_jump_u, state);
}

int From2Dto3DT::NumOutputVariables(void) const { return f2DModel->NumOutputVariables(); }
void From2Dto3DT::OutputLabels(ArrayT<StringT>& labels) const
{
	f2DModel->OutputLabels(labels);
}
		
void From2Dto3DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
	dArrayT new_jump_u(2);
	new_jump_u[0] = sqrt(jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1]);
	new_jump_u[1] = jump_u[2];
	
	f2DModel->ComputeOutput(new_jump_u, state, output);
}
void From2Dto3DT::PrintName(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	out << " 3D up-conversion of ";
	f2DModel->PrintName(out);
#else
#pragma unused(out)
#endif
}

/* print parameters to the output stream */
void From2Dto3DT::Print(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	f2DModel->Print(out);
#else
#pragma unused(out)
#endif
}
