/* $Id: TiedPotentialT.cpp,v 1.13 2003-03-26 20:00:08 cjkimme Exp $  */
/* created: cjkimme (10/23/2001) */

#include "TiedPotentialT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"
#include "fstreamT.h"
#include "StringT.h"
#include "iArrayT.h"

/* class parameters */

using namespace Tahoe;

const int    knumDOF = 2;
const double kExpMax = 100;

/* constructor */
TiedPotentialT::TiedPotentialT(ifstreamT& in): 
	SurfacePotentialT(knumDOF), TiedPotentialBaseT()
{
    in >> fnvec1; /* read in direction to sample stress state at */
    in >> fnvec2;
 
    /*make it a unit vector */
    double mag = sqrt(fnvec1*fnvec1+fnvec2*fnvec2);
    
    fnvec1 /= mag;
    fnvec2 /= mag;
 
    int nBulkGroups;
    in >> nBulkGroups; if (nBulkGroups < 1) throw ExceptionT::kBadInputValue;
    iBulkGroups.Dimension(nBulkGroups);
    for (int i = 0; i < nBulkGroups; i++)
    {
    	in >> iBulkGroups[i]; 
    	if (iBulkGroups[i] < 0) throw ExceptionT::kBadInputValue;
    	iBulkGroups[i]--;
    }
    
	in >> qTv; /* 0 for Xu-Needleman. 1 for TvergHutch */
	
	if (qTv)
	{
		in >> fsigma; if (fsigma < kSmall) throw ExceptionT::kBadInputValue;
		in >> d_n;  if (d_n < kSmall) throw ExceptionT::kBadInputValue;
		in >> d_t;  if (d_t < kSmall) throw ExceptionT::kBadInputValue;
		in >> fL_1; if (fL_1 < kSmall || fL_1 > 1.) throw ExceptionT::kBadInputValue;
		in >> fL_2; if (fL_2 > 1. || fL_2 < fL_1) throw ExceptionT::kBadInputValue;
		
		double fL_fail; 
		in >> fL_fail;
		
		in >> fL_0; if (fL_0 < 0. || fL_0 > 1.) throw ExceptionT::kBadInputValue;
		r_fail = 1.;
		if (fL_0 < fL_1)
			fsigma_critical = fsigma/fL_1*fL_0;
		else
			if (fL_0 < fL_2)
				fsigma_critical = fsigma;
			else
				fsigma_critical = fsigma*(1-fL_0)/(1-fL_2);
	}
	else
	{
	  	double q, r;
	  	in >> q;	// phi_t/phi_n
	  	in >> r; //delta_n* /d_n
	  	if (q < 0.0 || r < 0.0) throw ExceptionT::kBadInputValue;
	
		in >> d_n; if (d_n <= kSmall) throw ExceptionT::kBadInputValue;
		in >> d_t; if (d_t <= kSmall) throw ExceptionT::kBadInputValue;
		in >> phi_n; if (phi_n <= kSmall) throw ExceptionT::kBadInputValue;
		in >> r_fail; if (r_fail < 1.0) throw ExceptionT::kBadInputValue;
		fsigma_critical = phi_n/d_n/exp(1.); 
	}
	
	fsigma_critical *= fsigma_critical;
}

/* return the number of state variables needed by the model */
int TiedPotentialT::NumStateVariables(void) const { return 1; }

/* surface potential */ 
double TiedPotentialT::FractureEnergy(const ArrayT<double>& state) 
{
#pragma unused(state)

   	if (!qTv)
   		return phi_n;
   	else
   		if (fL_0 > fL_2)
   			return .5*fsigma*d_n*(1. - fL_0)*(1. - fL_0)/(1. - fL_2);
   		else
   			if (fL_0 > fL_1)
   				return .5*fsigma*d_n*(1. + fL_2 - 2.*fL_0);
   			else
   				return .5*fsigma*d_n*(1. + fL_2 - fL_1 - fL_0*fL_0/fL_1);
}

double TiedPotentialT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(jump_u)
#pragma unused(state)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	return 0.;
}
	
/* traction vector given displacement jump vector */	
const dArrayT& TiedPotentialT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, const bool& qIntegrate)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
	if (fTimeStep <= 0.0) {
#ifndef _SIERRA_TEST_
		cout << "\n TiedPotentialT::Traction: expecting positive time increment: "
		     << fTimeStep << endl;
#endif
		throw ExceptionT::kBadInputValue;
	}
#endif


	if (state[0] != 1. && state[0] != -10.)
	{
			fTraction = 0.;
	}
	else
	{
		if (!qTv)
		{
			double du_n = jump_u[1]+d_n;
			double fexpf = exp(-du_n/d_n)* exp(-jump_u[0]*jump_u[0]/d_t/d_t);
		
			fTraction[0] = 2.*phi_n*jump_u[0]/d_t*fexpf/d_t;
			fTraction[1] = phi_n/d_n*fexpf*du_n/d_n;
		}
		else
		{
		
			double effn = (jump_u[1]+fL_0*d_n);
			double fL = sqrt((jump_u[0]*jump_u[0]+effn*effn)/d_n/d_n);
			
			if (fL < fL_1)
			{	
				fTraction[0] = fsigma/fL_1*jump_u[0]/d_t*d_n/d_t;
				fTraction[1] = fsigma/fL_1*effn/d_n;
			}
			else 
				if (fL < fL_2) 
				{
					fTraction[0] = fsigma/fL*jump_u[0]/d_t*d_n/d_t;
					fTraction[1] = fsigma/fL*effn/d_n;
				}
				else
					if (fL < 1.)
					{
						fTraction[0] = fsigma*(1.-fL)/(1.-fL_2)*jump_u[0]/d_t*d_n/d_t;
						fTraction[1] = fsigma*(1.-fL)/(1.-fL_2)*effn/d_n;
					}
					else
						fTraction = 0.;
		}
		
		if (state[0] == -10. && qIntegrate)
			state[0] = 1.;
	}

	return fTraction;
	
}

/* potential stiffness */
const dMatrixT& TiedPotentialT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif 
	
	if (state[0] != -10. && state[0] != 1.)
	{
		fStiffness = 0.;
	}
	else
	{
		fStiffness = 0.;
		if (!qTv)
		{
			double du_n = jump_u[1]+d_n;
			double fexpf = exp(-du_n/d_n)* exp(-jump_u[0]*jump_u[0]/d_t/d_t);
		
			fStiffness[0] = 2.*phi_n/d_t/d_t*fexpf*(1.-2.*jump_u[0]*jump_u[0]/d_t/d_t);
			fStiffness[1] = -2.*phi_n*jump_u[0]/d_n*fexpf/d_t/d_t;
	    	fStiffness[2] = -2.*phi_n*du_n/d_n*fexpf*jump_u[0]/d_t/d_t/d_n;
			fStiffness[3] = phi_n/d_n*fexpf*(1.-du_n/d_n)/d_n;
		}
		else
		{				
			double u_t = jump_u[0];
			double u_n = jump_u[1]+fL_0*d_n;

			double dtm2 = 1./d_t/d_t;
			double dnm2 = 1./d_n/d_n;
			double L = sqrt(u_t*u_t*dtm2 + u_n*u_n*dnm2);
			
			if (L < fL_1)
			{
				fStiffness[0] = d_n/d_t*fsigma/fL_1/d_t;
				fStiffness[1] = 0.0;
				fStiffness[2] = 0.0;
				fStiffness[3] = fsigma/fL_1/d_n;
			}	
			else 
			{
				double lt_0 = jump_u[0]*dtm2;
				double lt_2 = jump_u[1]*dnm2;
				u_n *= dnm2;
		
				if (L < fL_2) // K2
				{
					double dijTerm = fsigma/L*d_n;
					
					fStiffness[0] = dijTerm*dtm2;
					fStiffness[3] = dijTerm*dnm2;
					dijTerm /= -L*L;
					fStiffness[0] += dijTerm*lt_0*lt_0;
					fStiffness[2] = fStiffness[1] = dijTerm*lt_0*lt_2;
					fStiffness[3] += dijTerm*u_n*u_n;
				}
				else 
					if (L < 1.) // K3
					{
						double dijTerm = fsigma*(1./L-1.)/(1.-fL_2)*d_n;
					
						fStiffness[0] = dijTerm*dtm2;
						fStiffness[3] = dijTerm*dnm2;
						dijTerm = -fsigma/(1.-fL_2)*d_n/L/L/L;
						fStiffness[0] += dijTerm*lt_0*lt_0;
						fStiffness[3] += dijTerm*u_n*u_n;
						fStiffness[1] = fStiffness[2] = dijTerm*lt_0*u_n;
					}
					else
					{
						fStiffness[0] = 0.0;
						fStiffness[1] = 0.0;
						fStiffness[2] = 0.0;
						fStiffness[3] = 0.0;	
					}
			}
		}
	}
	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT TiedPotentialT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(jump_u)
#pragma unused(state)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif
    
	double u_t1 = jump_u[0];
	double u_t  = sqrt(u_t1*u_t1);
	double u_n ;
	if (!qTv)
		u_n = jump_u[1]+d_n;
	else
		u_n = jump_u[1] + fL_0*d_n;
	
	/* square box for now */
	if (u_n > r_fail*d_n)
		return Failed;
	else 
		if (u_n > kSmall)
			return Critical;
		else
			return Precritical;
}

void TiedPotentialT::PrintName(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << "    TiedPotentialT (modified Xu-Needleman) 2D \n";
#endif
}

/* print parameters to the output stream */
void TiedPotentialT::Print(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << " Characteristic normal opening to failure. . . . = " << d_n     << '\n';
	out << " Characteristic tangential opening to failure. . = " << d_t     << '\n';
	out << " Mode I work to fracture (phi_n) . . . . . . . . = " << phi_n   << '\n';
	out << " Failure ratio (d_n/delta_n or d_t/delta_t). . . = " << r_fail   << '\n';
#endif
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int TiedPotentialT::NumOutputVariables(void) const { return 1; }

void TiedPotentialT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(1);
	labels[0] = "state[1]";
}

void TiedPotentialT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#pragma unused(state) 
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif	

	output[0] = 0.;//state[0];
}

bool TiedPotentialT::NeedsNodalInfo(void) { return true; }

int TiedPotentialT::NodalQuantityNeeded(void) 
{ 
        return kAverageCode; /*get stress tensor from bulk */ 
}

bool TiedPotentialT::InitiationQ(const double* sigma) 
{
	double t1 = sigma[0]*fnvec1+sigma[2]*fnvec2;
	double t2 = sigma[2]*fnvec1+sigma[1]*fnvec2;
	
	return t1*t1 + t2*t2 >= fsigma_critical;
}

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

