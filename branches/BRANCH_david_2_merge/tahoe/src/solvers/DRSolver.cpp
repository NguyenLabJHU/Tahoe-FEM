/* $Id: DRSolver.cpp,v 1.9 2004-01-05 07:07:19 paklein Exp $ */
/* created: PAK/CBH (10/03/1996) */
#include "DRSolver.h"

#include <iostream.h>
#include <math.h>

#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "FEManagerT.h"
#include "CCSMatrixT.h"

using namespace Tahoe;

/* constructor */
DRSolver::DRSolver(FEManagerT& fe_manager, int group): 
	NLSolver(fe_manager, group)
{
#ifdef __NO_RTTI__
	fCCSLHS = (CCSMatrixT*) fLHS;
#else
	fCCSLHS = TB_DYNAMIC_CAST(CCSMatrixT*, fLHS);
	if (!fCCSLHS) throw ExceptionT::kGeneralFail;
#endif
}

/* configure the global equation system */
void DRSolver::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	NLSolver::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* memory */
	fMass.Dimension(loc_num_eq);
	fVel.Dimension(loc_num_eq);
	fDisp.Dimension(loc_num_eq);
	fKd.Dimension(loc_num_eq);
	
	fNumEquations = loc_num_eq;
}

/* generate the solution for the current time sequence */
SolverT::SolutionStatusT DRSolver::Solve(int num_iterations)
{
	try
	{	  	

	/* form the residual force vector */
	fRHS = 0.0;
	fFEManager.FormRHS(Group());

	SolutionStatusT status = ExitIteration(fRHS.Magnitude(), fNumIteration);
	while (status != kConverged &&
		(num_iterations == -1 || IterationNumber() < num_iterations))
	{
		/* form the stiffness matrix */
		fLHS->Clear();				
		fFEManager.FormLHS(Group(), GlobalT::kDiagonal);
	
		/* compute mass for stability */
		ComputeMass();

		/* calculate velocity */
		ComputeVelocity();

		/* displacement update */
		fVel *= fFEManager.TimeStep();

		/* update displacements */
		fFEManager.Update(Group(), fVel);
				
		/* calculate critical damping */
		ComputeDamping();

		/* form the residual force vector */
		fRHS = 0.0;
		fFEManager.FormRHS(Group());
		
		/* check status */
		status = ExitIteration(fRHS.Magnitude(), fNumIteration);
	}  

	/* normal */
	return status;
	
	} /* end try */
	
	catch (ExceptionT::CodeT code) { return kFailed; }
}

/*************************************************************************
* Private
*************************************************************************/

/* compute the pseudo-mass */
void DRSolver::ComputeMass(void)
{
	double dt_sqr4 = 0.25*pow(1.1*fFEManager.TimeStep(), 2);

	for (int i = 0; i < fNumEquations; i++)
		fMass[i] = dt_sqr4*fCCSLHS->AbsRowSum(i);
}

void DRSolver::ComputeVelocity(void)
{
	if (fNumIteration == 0)
	{
		fVel = fRHS;
		
		double dt2 = 0.5*fFEManager.TimeStep();
		
		for (int i = 0; i < fNumEquations; i++)
			fVel[i] *= dt2/fMass[i];
	}
	else
	{
		double dt = fFEManager.TimeStep();
		double k1 = 2.0 - fDamp*dt;
		double k2 = 1.0/(2.0 + fDamp*dt);
		
		fVel *= k1;
		
		for (int i = 0; i < fNumEquations; i++)
			fVel[i] = (fVel[i] + 2.0*dt*fRHS[i]/fMass[i])*k2;
	}
}

void DRSolver::ComputeDamping(void)
{
	/* numerator */
	fFEManager.GetUnknowns(fGroup, 0, fDisp);
	fCCSLHS->Multx(fDisp, fKd);
	double numer = dArrayT::Dot(fDisp,fKd);
	
	/* denominator */
	double deno = 0.0;
	for (int i = 0; i < fNumEquations; i++)
	{	
		double di = fDisp[i];
		deno += di*di*fMass[i];
	}
	
	fDamp = 2.0*sqrt(numer/deno);
}

