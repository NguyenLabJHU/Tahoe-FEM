/* $Id: NLK0Solver.cpp,v 1.12 2004-06-17 07:42:05 paklein Exp $ */
/* created: paklein (10/01/1996) */
#include "NLK0Solver.h"
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "FEManagerT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

using namespace Tahoe;

/* line search parameters */
const double ks_max_factor = 5.0; //normally 1.0

/* constructor */
NLK0Solver::NLK0Solver(FEManagerT& fe_manager, int group):
	NLSolver(fe_manager, group),
	fFormTangent(1),
	fLastTangent(fe_manager.Output(), fLHS->CheckCode())
{
#ifdef __NO_RTTI__
	pCCSLHS = NULL; //not allowed for now
#else
	pCCSLHS = dynamic_cast<CCSMatrixT*>(fLHS);
#endif
	if (!pCCSLHS)
		ExceptionT::GeneralFail("NLK0Solver::NLK0Solver", 
			"solver requires matrix type: %d (symmetric, profile solver)", kProfileSolver);
}

/***********************************************************************
* Protected
***********************************************************************/

/* form and solve - returns the magnitude of the residual */
double NLK0Solver::SolveAndForm(int& iteration)
{		
	const char caller[] = "NLK0Solver::SolveAndForm";

	/* form the stiffness matrix */
	if (fFormTangent)
	{
		fLHS_lock = kOpen;
		if (fLHS_update) pCCSLHS->Clear();
		fFEManager.FormLHS(Group(), GlobalT::kSymmetric);
		fLHS_lock = kLocked;

		/* solve equation system */
		fUpdate = fRHS;
		if(!pCCSLHS->Solve(fRHS)) ExceptionT::BadJacobianDet(caller);
	
		/* check for positive definiteness */
		if ( pCCSLHS->HasNegativePivot() )
		{
			fFormTangent = 0; //no switching back, yet
		
			/* restore RHS */
			fRHS = fUpdate;
		}
		else /* store tangent */
			fLastTangent = (*pCCSLHS);
	}
	
	/* use last positive definite tangent */
	if (!fFormTangent)
		if (!fLastTangent.Solve(fRHS)) 
			ExceptionT::BadJacobianDet(caller);
			 		
	/* update system */
	fFEManager.Update(Group(), fRHS);
								
	/* compute new residual */
	iteration++;
	fRHS = 0.0;
	fFEManager.FormRHS(Group());

	/* combine residual magnitude with update magnitude */
	/* e = a1 |R| + a2 |delta_d|                        */
	//not implemented!
			
	return Residual(fRHS);
}
