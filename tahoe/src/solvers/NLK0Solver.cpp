/* $Id: NLK0Solver.cpp,v 1.2 2002-03-22 02:27:26 paklein Exp $ */
/* created: paklein (10/01/1996) */

#include "NLK0Solver.h"
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include "Constants.h"
#include "ExceptionCodes.h"
#include "FEManagerT.h"
#include "fstreamT.h"

/* line search parameters */
const double ks_max_factor = 5.0; //normally 1.0

/* constructor */
NLK0Solver::NLK0Solver(FEManagerT& fe_manager):
	NLSolver(fe_manager),
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
		throw eGeneralFail;
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
		fFEManager.FormLHS();

		/* solve equation system */
		fUpdate = fRHS;
		pCCSLHS->Solve(fRHS);
	
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
		fLastTangent.Solve(fRHS);
			 		
	/* update system */
	fFEManager.Update(fRHS);
								
	/* compute new residual */
	fRHS = 0.0;
	fFEManager.FormRHS();

	/* combine residual magnitude with update magnitude */
	/* e = a1 |R| + a2 |delta_d|                        */
	//not implemented!
			
	return( fRHS.Magnitude() );
}
