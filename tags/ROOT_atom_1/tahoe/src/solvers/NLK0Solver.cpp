/* $Id: NLK0Solver.cpp,v 1.8 2002-11-28 17:30:31 paklein Exp $ */
/* created: paklein (10/01/1996) */

#include "NLK0Solver.h"
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "FEManagerT.h"
#include "fstreamT.h"

/* line search parameters */

using namespace Tahoe;

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
	if (!pCCSLHS) {
		cout << "\n NLK0Solver::NLK0Solver: solver requires matrix type: " << kProfileSolver 
		     << " (symmetric)" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/***********************************************************************
* Protected
***********************************************************************/

/* form and solve - returns the magnitude of the residual */
double NLK0Solver::SolveAndForm(bool junk)
{		
#pragma unused(junk)

	/* form the stiffness matrix */
	if (fFormTangent)
	{
		pCCSLHS->Clear();
		fFEManager.FormLHS(Group(), GlobalT::kNonSymmetric);

		/* solve equation system */
		fUpdate = fRHS;
		if(!pCCSLHS->Solve(fRHS)) throw ExceptionT::kBadJacobianDet;
	
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
		if (!fLastTangent.Solve(fRHS)) throw ExceptionT::kBadJacobianDet;
			 		
	/* update system */
	fFEManager.Update(Group(), fRHS);
								
	/* compute new residual */
	fRHS = 0.0;
	fFEManager.FormRHS(Group());

	/* combine residual magnitude with update magnitude */
	/* e = a1 |R| + a2 |delta_d|                        */
	//not implemented!
			
	return( fRHS.Magnitude() );
}
