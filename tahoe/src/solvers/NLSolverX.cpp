/* $Id: NLSolverX.cpp,v 1.8 2002-12-13 02:42:55 paklein Exp $ */
/* created: paklein (08/25/1996) */
#include "NLSolverX.h"

#include <iostream.h>
#include <math.h>
#include "fstreamT.h"
#include "toolboxConstants.h"
#include "ExceptionT.h"
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
#ifdef __NO_RTTI__
	cout << "\n NLSolverX::Initialize: RTTI required" << endl;
	throw ExceptionT::kGeneralFail;
#else
	pCCS = dynamic_cast<CCSMatrixT*>(fLHS);
	if (!pCCS)
	{
		cout << "\n NLSolverX::Initialize: expecting CCS global matrix" << endl;
		throw ExceptionT::kBadInputValue;
	}
#endif		

	ifstreamT&  in = fFEManager.Input();
	ostream& out = fFEManager.Output();

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
	SolutionStatusT solutionflag = ExitIteration(error);
	while(solutionflag == kContinue &&
		(num_iterations == -1 || IterationNumber() < num_iterations))
	{
		/* force new tangent */
		if (fMinFreshTangents > 0 &&
			fNumIteration < fMinFreshTangents)
			fFormNewTangent = 1;
			
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
				
		error = SolveAndForm(fFormNewTangent, true);
		solutionflag = ExitIteration(error);

		/* check for negative pivots */
		if (fCheckNegPivots && fFormNewTangent)
			negative_pivot = pCCS->HasNegativePivot();

		/* no fail on quasi-Newton step */
		if (solutionflag == kFailed && !fFormNewTangent)
		{
			cout << "\n NLSolverX: overriding kFailed on quasi-Newton step" << endl;
			solutionflag = kContinue;
			fFormNewTangent = 1;
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
		}
		else
			fFormNewTangent = 0;
				
		/* don't re-use if negative pivot */
		if (negative_pivot)
		{
			fFormNewTangent = 1;
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


/* form and solve the equation system */
double NLSolverX::SolveAndForm(bool newtangent, bool clear_LHS)
{		
	/* form the stiffness matrix */
	if (newtangent)
	{
		if (clear_LHS) fLHS->Clear();
		fFEManager.FormLHS(Group(), GlobalT::kNonSymmetric);
	}
		 		
	/* solve equation system */
	if (!fLHS->Solve(fRHS)) throw ExceptionT::kBadJacobianDet;
	//double updatemag = fRHS.Magnitude(); //for alternate error measures

	/* update system */
	fFEManager.Update(Group(), fRHS);
	
	/* keep last update for reset */
	fLastUpdate = fRHS;
							
	/* compute new residual */
	fRHS = 0.0;
	fFEManager.FormRHS(Group());

	/* combine residual magnitude with update magnitude */
	/* e = a1 |R| + a2 |delta_d|                        */
	//not implemented!
			
	return( fRHS.Magnitude() );
}

/*************************************************************************
* Protected
*************************************************************************/

/* relax system */
NLSolver::SolutionStatusT NLSolverX::Relax(int newtancount)
{	
#pragma unused (newtancount)

	cout <<   " Relaxation:" << '\n';

	/* reset counts */
	fNumIteration    = -1;
		
	/* form the first residual force vector */
	fRHS = 0.0;
	fFEManager.FormRHS(Group());	
	double      error = fRHS.Magnitude();
	double last_error = error;
		
	/* loop on error */
	int negative_pivot = 0;
	int tangentcount   = 0;
	int reuse_count    = 0;
	SolutionStatusT solutionflag = ExitIteration(error);
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
			
		error = SolveAndForm(fFormNewTangent, true);
		solutionflag = ExitIteration(error);

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

	/* allow for multiple relaxation */
	GlobalT::RelaxCodeT relaxcode = fFEManager.RelaxSystem(Group());
	while (relaxcode != GlobalT::kNoRelax)
	{	
		/* reset global equations */
		if (relaxcode == GlobalT::kReEQ ||
		    relaxcode == GlobalT::kReEQRelax)
			fFEManager.SetEquationSystem(Group());
						
		/* new equilibrium */
		if (relaxcode == GlobalT::kRelax ||
		    relaxcode == GlobalT::kReEQRelax)
		{
			if (Relax() == kFailed)
				return kFailed;
	   		else
				relaxcode = fFEManager.RelaxSystem(Group());
		}
		else
			relaxcode = GlobalT::kNoRelax;
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
