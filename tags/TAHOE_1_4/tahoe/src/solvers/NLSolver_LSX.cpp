/* $Id: NLSolver_LSX.cpp,v 1.4 2004-06-17 07:42:05 paklein Exp $ */
#include "NLSolver_LSX.h"
#include "FEManagerT.h"

#include <iostream.h>
#include "ifstreamT.h"
#include "ofstreamT.h"

using namespace Tahoe;

NLSolver_LSX::NLSolver_LSX(FEManagerT& fe_manager, int group):
	NLSolver_LS(fe_manager, group)
{
	ifstreamT& in = fFEManager.Input();
	
	/* read parameters */
	in >> fPuntTol;

	/* print parameters */
	ostream& out = fFEManager.Output();
	out << " Continuation tolerance. . . . . . . . . . . . . = " << fPuntTol << '\n';
}

/* start solution step */
void NLSolver_LSX::InitStep(void)
{
	/* inherited */
	NLSolver_LS::InitStep();
	
	/* reset minimum error */
	fMinStepRelError = 1.0;
}

/* allows continuation of NLSolver_LSX::fMinStepRelError is less
 * that NLSolver_LSX::fPuntTol */
NLSolver::SolutionStatusT NLSolver_LSX::ExitIteration(double error, int iteration)
{
	/* inherited */
	SolutionStatusT status = NLSolver_LS::ExitIteration(error, iteration);

	/* track smallest error */
	double rel_error = error/fError0;
	fMinStepRelError = (rel_error < fMinStepRelError) ? rel_error : fMinStepRelError;

	/* push through */
	if (status == kFailed &&         /* not converged */
		rel_error < fDivTolerance &&  /* not diverging */
		fMinStepRelError < fPuntTol) /* min error achieved */
	{
		cout << "\n NLSolver_LSX::ExitIteration: " << fMinStepRelError 
		     << " < " << fPuntTol << " (continuing tol)" << endl;
		status = kConverged;
	}
	
	return status;
}
