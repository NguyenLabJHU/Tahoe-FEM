/* $Id: NLSolverX.cpp,v 1.15 2004-09-14 18:16:10 paklein Exp $ */
/* created: paklein (08/25/1996) */
#include "NLSolverX.h"

#include <iostream.h>
#include <math.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "FEManagerT.h"
#include "CCSMatrixT.h"

/* control parameters */

using namespace Tahoe;

const int		kQuickConvTol = 6;
const int		kBiggerIncTol = 3;
const double	kRelDivTol    = 10.0;	

/* iteration status flags */
const int kContinue  = 0;
const int kConverged = 1;
const int kFailed    = 2;

/* constructor */
NLSolverX::NLSolverX(FEManagerT& fe_manager, int group):
	NLSolver(fe_manager, group)
{
ExceptionT::GeneralFail("NLSolverX::NLSolverX", "out of date");
#if 0
#ifdef __NO_RTTI__
	cout << "\n NLSolverX::Initialize: RTTI required" << endl;
	throw ExceptionT::kGeneralFail;
#else
	pCCS = TB_DYNAMIC_CAST(CCSMatrixT*, fLHS);
	if (!pCCS)
	{
		cout << "\n NLSolverX::Initialize: expecting CCS global matrix" << endl;
		throw ExceptionT::kBadInputValue;
	}
#endif		

	ifstreamT&  in = fFEManager.Input();
	ofstreamT& out = fFEManager.Output();

	in >> fMaxNewTangents;  if (fMaxNewTangents  < 1) throw ExceptionT::kBadInputValue;
	in >> fMaxTangentReuse; if (fMaxTangentReuse < 1) throw ExceptionT::kBadInputValue;
	in >> fMinFreshTangents;if (fMinFreshTangents != -1 && fMinFreshTangents < 1)
		throw ExceptionT::kBadInputValue;
	in >> fCheckNegPivots; if (fCheckNegPivots != 0 && fCheckNegPivots != 1)
		throw ExceptionT::kBadInputValue;
	int trust_exp;
	in >> trust_exp;
	if (trust_exp > 0)
	    cout << "\n NLSolverX::NLSolverX: WARNING: trust tolerance > 1" << endl;
	fTrustTol = pow(10.0,trust_exp);

	in >> fMinResRatio;     if (fMinResRatio > 1.0)   throw ExceptionT::kBadInputValue;

	out << " Maximum number new tangent matrices . . . . . . = " << fMaxNewTangents   << '\n';
	out << " Maximum number iterations with same tangent . . = " << fMaxTangentReuse  << '\n';
	out << " Number of new tangents at start (-1 for auto) . = " << fMinFreshTangents << '\n';
	out << " Automatic new tangent for negative pivots . . . = " << fCheckNegPivots   << '\n';
	out << " Trust tolerance for expanded iterations . . . . = " << fTrustTol         << '\n';
	out << " Mininum quasi-Newton reduction ratio. . . . . . = " << fMinResRatio      << '\n';
#endif
}

/* generate the solution for the current time sequence */
SolverT::SolutionStatusT NLSolverX::Solve(int num_iterations)
{
	try { 	

	/* form the first residual force vector */
	fRHS = 0.0;
	fFEManager.FormRHS(Group());	
	double      error = fRHS.Magnitude();
	double last_error = error;
			
	/* loop on error */
	fNumNewTangent     = 0;
	int into_trust     = 0;
	int negative_pivot = 0;
	int reuse_count    = 0;
	SolutionStatusT solutionflag = ExitIteration(error, fNumIteration);
	while(solutionflag == kContinue &&
		(num_iterations == -1 || IterationNumber() < num_iterations))
	{
		/* force new tangent */
		if (fMinFreshTangents > 0 &&
			fNumIteration < fMinFreshTangents)
		{
			fFormNewTangent = 1;
			fLHS_update = true;
		}
			
		if (fFormNewTangent)
		{
			fNumNewTangent++;
			reuse_count = 0;
			cout << " tangent:" << setw(3) << fNumNewTangent << ": ";				
		}
		else
		{
			reuse_count++;
			cout << " re-use :" << setw(3) << reuse_count <<  ": ";
		}

		/* reform the LHS matrix */
		if (fFormNewTangent) {
		
			/* recalculate */
			fLHS_lock = kOpen;
			fFEManager.FormLHS(Group(), fLHS->MatrixType());
			fLHS_lock = kLocked;
		
			/* compare with approxumate LHS */
			if (fLHS->CheckCode() == GlobalMatrixT::kCheckLHS) {
				const GlobalMatrixT* approx_LHS = ApproximateLHS(*fLHS);
				CompareLHS(*fLHS, *approx_LHS);
				delete approx_LHS;
			}
		}

		/* update solution */
		Iterate();
		fNumIteration++;

		/* tangent reformed next iteration? */
		if (fFormNewTangent)
			fLHS_update = true;

		/* recalculate residual */
		if (fLHS_update) {
			fLHS->Clear();
			fLHS_lock = kOpen; /* LHS open for assembly, too! */
		}
		else
			fLHS_lock = kIgnore; /* ignore assembled values */
		fRHS = 0.0;
		fFEManager.FormRHS(Group());	
		fLHS_lock = kLocked;
		fRHS_lock = kLocked;
		
		/* new error */
		error = Residual(fRHS);
		solutionflag = ExitIteration(error, fNumIteration);

		/* check for negative pivots */
		if (fCheckNegPivots && fFormNewTangent)
			negative_pivot = pCCS->HasNegativePivot();

		/* no fail on quasi-Newton step */
		if (solutionflag == kFailed && !fFormNewTangent)
		{
			cout << "\n NLSolverX: overriding kFailed on quasi-Newton step" << endl;
			solutionflag = kContinue;
			fFormNewTangent = 1;
			fLHS_update = true;
		}
		/* tangent re-use limit and fail recovery */
		else if (solutionflag == kContinue &&
				(error > fMinResRatio*last_error ||
				reuse_count == fMaxTangentReuse))
		{
			/* remove last quasi-Newton update */
			if (0 && error > fMinResRatio*last_error && !fFormNewTangent)
			{
				fLastUpdate *= -1.0;
				fFEManager.Update(Group(), fLastUpdate);
			}

			fFormNewTangent = 1;
			fLHS_update = true;				
		}
		else
			fFormNewTangent = 0;
			fLHS_update = false;
				
		/* don't re-use if negative pivot */
		if (negative_pivot)
		{
			fFormNewTangent = 1;
			fLHS_update = false;
			cout << "\n NLSolverX: no re-use for negative pivot" << endl;
		}
				
		/* too many new tangents */
		if (solutionflag == kContinue && fNumNewTangent == fMaxNewTangents)
		{
			if (error/fError0 < fTrustTol && !into_trust)
			{
				fNumNewTangent -= 3;
				into_trust      = 1;
				cout << "\n NLSolverX:: allowing 3 more tangents" << endl;
			}	
			else
			{
				cout << "\n NLSolverX:: too many new tangents" << endl;
				solutionflag = kFailed;
			}
		}
				
		/* store last error */
		last_error = error;
	}

	/* ensure good start for Relax() */
	if (into_trust) fFormNewTangent = 1;

	/* found solution - check relaxation */
	if (solutionflag == kConverged) solutionflag = DoConverged();

	/* normal */
	return solutionflag;

	} /* end try */
	
	/* exceptions */
	catch (ExceptionT::CodeT code) { return kFailed; }
}

/*************************************************************************
 * Protected
 *************************************************************************/

#if 0
/* relax system */
NLSolver::SolutionStatusT NLSolverX::Relax(int newtancount)
{	
#pragma unused (newtancount)

	cout <<   " Relaxation:" << '\n';

	/* reset counts */
	int iteration = -1;
		
	/* form the first residual force vector */
	fRHS = 0.0;
	fFEManager.FormRHS(Group());	
	double      error = fRHS.Magnitude();
	double last_error = error;
		
	/* loop on error */
	int negative_pivot = 0;
	int tangentcount   = 0;
	int reuse_count    = 0;
	SolutionStatusT solutionflag = ExitIteration(error, iteration);
	while (solutionflag == kContinue)
	{
		if (fFormNewTangent)
		{
			tangentcount++;
			reuse_count = 0;

			cout << " tangent:" << setw(3) << tangentcount << ": ";				
		}
		else
		{
			reuse_count++;
			cout << " re-use :" << setw(3) << reuse_count <<  ": ";
		}
			
		error = SolveAndForm(iteration);
		solutionflag = ExitIteration(error, iteration);

		/* check for negative pivots */
		if (fCheckNegPivots && fFormNewTangent)
			negative_pivot = pCCS->HasNegativePivot();

		/* tangent re-use limit */
		if (solutionflag == kContinue &&
			(error > fMinResRatio*last_error ||
			 reuse_count == fMaxTangentReuse))
		    fFormNewTangent = 1;
		else
		    fFormNewTangent = 0;
		    	
		/* don't re-use if negative pivot */
		if (negative_pivot) fFormNewTangent = 1;
			
		/* too many new tangents */
		if (solutionflag == kContinue && tangentcount == fMaxNewTangents)
		{
			cout << "\n NLSolverX:: too many new tangents\n" << endl;
			solutionflag = kFailed;
		}
			
		/* store last error */
		last_error = error;
	}
	return solutionflag;
}
#endif

/* handlers */
NLSolver::SolutionStatusT NLSolverX::DoConverged(void)
{
	/* increase time step ? (for multi-step sequences) */
	if (fQuickSolveTol > 1 && fNumNewTangent < fQuickSolveTol)
	{
		fQuickConvCount++;
		
//TEMP
cout << "\n NLSolverX::DoConverged: quick converged count: ";
cout << fQuickConvCount << "/" << fQuickSeriesTol << endl;
				
		if (fQuickConvCount >= fQuickSeriesTol)
		{
			fQuickConvCount = 0;
			fFEManager.IncreaseLoadStep();
		}
	}
	else
	{
		/* restart count if convergence is slow */
		fQuickConvCount = 0;	

//TEMP
cout << "\n NLSolverX::DoConverged: reset quick converged: ";
cout << fQuickConvCount << "/" << fQuickSeriesTol << endl;
	}

	/* success */
	return kConverged;
}

void NLSolverX::ResetStep(void)
{
	/* inherited */
	NLSolver::ResetStep();

	fFormNewTangent = 1;
}
