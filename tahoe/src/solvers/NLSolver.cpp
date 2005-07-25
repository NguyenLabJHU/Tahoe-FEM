/* $Id: NLSolver.cpp,v 1.37.20.1 2005-07-25 02:37:25 paklein Exp $ */
/* created: paklein (07/09/1996) */
#include "NLSolver.h"

#include <iostream.h>
#include <math.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "FEManagerT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* constructor */
NLSolver::NLSolver(FEManagerT& fe_manager, int group):
	SolverT(fe_manager, group),
	fMaxIterations(-1),
	fMinIterations(0),
	fReformTangentIterations(1),
	fZeroTolerance(0.0),
	fTolerance(0.0),
	fDivTolerance(-1.0),
	fQuickConvCount(0),
	fIterationOutputCount(0),
	fVerbose(1),
	fQuickSolveTol(6),
	fQuickSeriesTol(3),
	fIterationOutputIncrement(0),
	fRestartIteration(0)
{
	SetName("nonlinear_solver");

	/* console variables */
	iAddVariable("max_iterations", fMaxIterations);
	iAddVariable("min_iterations", fMinIterations);
	iAddVariable("abs_tolerance", fZeroTolerance);
	iAddVariable("rel_tolerance", fTolerance);
	iAddVariable("div_tolerance", fDivTolerance);
	iAddVariable("iteration_output_inc", fIterationOutputIncrement);
}

/* start solution step */
void NLSolver::InitStep(void)
{
	/* inherited */
	SolverT::InitStep();

	/* open iteration output */
	InitIterationOutput();
	
	/* reset marker */
	fRestartIteration = IterationNumber();
}

/* generate the solution for the current time sequence */
SolverT::SolutionStatusT NLSolver::Solve(int max_iterations)
{
	/* write some header information */
	if (fLHS->CheckCode() != GlobalMatrixT::kNoCheck) {
		ofstreamT& out = fFEManager.Output();
		out << " NLSolver::Solve:\n"
		    << "      group = " << fGroup+1 << '\n'
		    << " iterations = " << max_iterations << '\n';	
	}

	try
	{ 	

	/* counters */
	int num_iterations = 0;
	int tan_iterations = 0;

	/* form the first residual force vector */
	fRHS_lock = kOpen;
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

	/* initial error */
	double error = Residual(fRHS);
	SolutionStatusT solutionflag = ExitIteration(error, fNumIteration);

#if 0
	/* check for relaxation */
	if (solutionflag == kConverged) {
		GlobalT::RelaxCodeT relaxcode = GetRelaxCode();

		/* reset global equations */
		if (relaxcode == GlobalT::kReEQ || relaxcode == GlobalT::kReEQRelax)
			ContinueStep(); // go back and renumber and go again... or do we just want to not renumber in here...
			
		/* recompute force and continue iterating */
		if (relaxcode == GlobalT::kRelax || relaxcode == GlobalT::kReEQRelax) {

			/* set marker */
			fRestartIteration = IterationNumber();

			/* tangent reformed next iteration? */
			if (tan_iterations + 1 == fReformTangentIterations)
				fLHS_update = true;
			else
				fLHS_update = false;

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

			/* test for convergence */
			solutionflag = ExitIteration(error, fNumIteration);
		}
	}
#endif
			
	/* loop on error */
	while (solutionflag == kContinue &&
		(max_iterations == -1 || num_iterations < max_iterations))
	{
		/* increment counters */
		num_iterations++;
		tan_iterations++;
	
		/* recompute LHS? */
		if (num_iterations == 1 || tan_iterations >= fReformTangentIterations) {
			fLHS_update = true;
			tan_iterations = 0;
		}
		else fLHS_update = false;

		/* reform the LHS matrix */
		if (fLHS_update) {
		
			/* recalculate */
			fLHS_lock = kOpen;
			fFEManager.FormLHS(Group(), fLHS->MatrixType());
			fLHS_lock = kLocked;
		
			/* compare with approximate LHS */
			if (fLHS->CheckCode() == GlobalMatrixT::kCheckLHS) {
				const GlobalMatrixT* approx_LHS = ApproximateLHS(*fLHS);
				CompareLHS(*fLHS, *approx_LHS);

				/* use the approximate matrix */
				GlobalMatrixT* tmp = fLHS;
				fLHS = (GlobalMatrixT*) approx_LHS;
				approx_LHS = tmp;

				delete approx_LHS;
			}
		}

		/* update solution */
		Iterate();
		fNumIteration++;

		/* tangent reformed next iteration? */
		if (tan_iterations + 1 == fReformTangentIterations)
			fLHS_update = true;
		else
			fLHS_update = false;

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

		/* test for convergence */
		solutionflag = ExitIteration(error, fNumIteration);

#if 0
		/* check for relaxation */
		if (solutionflag == kConverged) {
			GlobalT::RelaxCodeT relaxcode = GetRelaxCode();
		
			/* reset global equations */
			if (relaxcode == GlobalT::kReEQ || relaxcode == GlobalT::kReEQRelax)
				ContinueStep(); // go back and renumber and go again... or do we just want to not renumber in here...
			
			/* recompute force and continue iterating */
			if (relaxcode == GlobalT::kRelax || relaxcode == GlobalT::kReEQRelax) {

				/* set marker */
				fRestartIteration = IterationNumber();

				/* tangent reformed next iteration? */
				if (tan_iterations + 1 == fReformTangentIterations)
					fLHS_update = true;
				else
					fLHS_update = false;

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

				/* test for convergence */
				solutionflag = ExitIteration(error, fNumIteration);
			}
		}
#endif
	}

	/* found solution - check relaxation */
//	if (solutionflag == kConverged)
//		solutionflag = DoConverged();
#pragma message("delete me")
			
	return solutionflag;
	}
		
	/* abnormal ending */
	catch (ExceptionT::CodeT code)
	{
		cout << "\n NLSolver::Solve: exception at step number "
             << fFEManager.StepNumber() << " with step "
             << fFEManager.TimeStep() 
             << "\n     " << code << ": " << ExceptionT::ToString(code) << endl;

		/* error occurred here -> trip checksum */
		if (code != ExceptionT::kBadHeartBeat) fFEManager.Communicator().Sum(code);
		
		return kFailed;
	}
}

/* end solution step */
void NLSolver::CloseStep(void)
{
	/* inherited */
	SolverT::CloseStep();

	/* close iteration output */	
	CloseIterationOutput();

	/* increase time step ? (for multi-step sequences) */
	if (fQuickSolveTol > 1 && fNumIteration < fQuickSolveTol)
	{
		fQuickConvCount++;		
		cout << "\n NLSolver::CloseStep: quick converged count: ";
		cout << fQuickConvCount << "/" << fQuickSeriesTol << endl;

		if (fQuickConvCount >= fQuickSeriesTol)
			if (fFEManager.IncreaseLoadStep() == 1)
				fQuickConvCount = 0;
	}
	else
	{
		/* restart count if convergence is slow */
		fQuickConvCount = 0;	
		cout << "\n NLSolver::CloseStep: reset quick converged: ";
		cout << fQuickConvCount << "/" << fQuickSeriesTol << endl;	
	}
}

/* error handler */
void NLSolver::ResetStep(void)
{
	/* inherited */
	SolverT::ResetStep();

	/* reset count to increase load step */
	fQuickConvCount = 0;

	/* reset flag */
	fLHS_update = true;

	/* restore output */
	CloseIterationOutput();
}

/* renumbering handler */
void NLSolver::ContinueStep(void)
{
	// Calls InitStep to handle reconfiguration/renumbering
	InitStep();
}


/* (re-)set the reference error */
void NLSolver::SetReferenceError(double error)
{
	cout <<   "\n Group : " << fGroup+1 << '\n' << " Absolute error = " << error << endl;
	fError0 = error;
}

#if 0
/* handlers */
NLSolver::SolutionStatusT NLSolver::DoConverged(void)
{
	/* increase time step ? (for multi-step sequences) */
	if (fQuickSolveTol > 1 && fNumIteration < fQuickSolveTol)
	{
		fQuickConvCount++;		
		cout << "\n NLSolver::DoConverged: quick converged count: ";
		cout << fQuickConvCount << "/" << fQuickSeriesTol << endl;

		if (fQuickConvCount >= fQuickSeriesTol)
			if (fFEManager.IncreaseLoadStep() == 1)
				fQuickConvCount = 0;
	}
	else
	{
		/* restart count if convergence is slow */
		fQuickConvCount = 0;	
		cout << "\n NLSolver::DoConverged: reset quick converged: ";
		cout << fQuickConvCount << "/" << fQuickSeriesTol << endl;	
	}

	/* success */
	return kConverged;						
}
#endif

/* divert output for iterations */
void NLSolver::InitIterationOutput(void)
{
	if (fIterationOutputIncrement > 0)
	{
		/* root of output files */
		StringT root;
		root.Root(fFEManager.InputFile());
		
		/* remove processor designation */ 
		if (fFEManager.Size() > 1) root.Root();
		
		/* solver group */
		root.Append(".gp", Group());
		
		/* increment */
		root.Append(".", fFEManager.StepNumber());
		root.Append("of", fFEManager.NumberOfSteps());

		/* tack on processor designation */
		if (fFEManager.Size() > 1) root.Append(".p", fFEManager.Rank());

		/* set temporary output */
		fFEManager.DivertOutput(root);
		
		/* reset count */
		fIterationOutputCount = 0;
	}
}

void NLSolver::CloseIterationOutput(void)
{
	if (fIterationOutputIncrement > 0)
		fFEManager.RestoreOutput();
}

/* describe the parameters needed by the interface */
void NLSolver::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SolverT::DefineParameters(list);

	/* additional parameters */
	list.AddParameter(fMaxIterations, "max_iterations");

	ParameterT min_iterations(fMinIterations, "min_iterations");
	min_iterations.SetDefault(fMinIterations);
	list.AddParameter(min_iterations);

	ParameterT reform_tangent_its(fReformTangentIterations, "reform_tangent_iterations");
	reform_tangent_its.SetDefault(fReformTangentIterations);
	reform_tangent_its.AddLimit(fReformTangentIterations, LimitT::LowerInclusive);
	list.AddParameter(reform_tangent_its);

	list.AddParameter(fZeroTolerance, "abs_tolerance");
	list.AddParameter(fTolerance, "rel_tolerance");
	list.AddParameter(fDivTolerance, "divergence_tolerance");

	ParameterT quick_solve_iter(fQuickSolveTol, "quick_solve_iter");
	quick_solve_iter.SetDefault(fQuickSolveTol);
	list.AddParameter(quick_solve_iter);

	ParameterT quick_solve_count(fQuickSeriesTol, "quick_solve_count");
	quick_solve_count.SetDefault(fQuickSeriesTol);
	list.AddParameter(quick_solve_count);

	ParameterT output_inc(fIterationOutputIncrement, "output_inc");
	output_inc.SetDefault(fIterationOutputIncrement);
	output_inc.AddLimit(0, LimitT::LowerInclusive);
	list.AddParameter(output_inc);
}

/* accept parameter list */
void NLSolver::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SolverT::TakeParameterList(list);

	/* extract parameters */
	fMaxIterations = list.GetParameter("max_iterations");
	fMinIterations = list.GetParameter("min_iterations");
	fReformTangentIterations = list.GetParameter("reform_tangent_iterations");
	fZeroTolerance = list.GetParameter("abs_tolerance");
	fTolerance     = list.GetParameter("rel_tolerance");
	fDivTolerance  = list.GetParameter("divergence_tolerance");
	fQuickSolveTol = list.GetParameter("quick_solve_iter");
	fQuickSeriesTol= list.GetParameter("quick_solve_count");
	fIterationOutputIncrement = list.GetParameter("output_inc");
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* apply system update (socket for line searching) */
void NLSolver::Update(const dArrayT& update, const dArrayT* residual)
{
#pragma unused(residual)

	/* full Newton update */
	fFEManager.Update(Group(), update);
}

#if 0
/* relax system */
NLSolver::SolutionStatusT NLSolver::Relax(int newtancount)
{
#pragma unused(newtancount)

	cout <<   "\n Relaxation:" << '\n';

	/* keep current iteration count */
	int iteration = fNumIteration;
	fNumIteration = -1;

	/* re-solve equilibrium */
	SolutionStatusT status = Solve(-1);

	/* restore the iteration count */
	fNumIteration = iteration;

	return status;
}
#endif

/* returns 1 if the iteration loop should be left, otherwise
* returns 0.  The iteration loop can be exited for the
* following reasons:
*
*	(1) error meets the convergence criteria
*	(2) the iteration limit has been exceeded
*	(3) the error is diverging
*
* For (2) and (3), the load increment will be cut and the
* iteration re-entered with the next Step() call */
NLSolver::SolutionStatusT NLSolver::ExitIteration(double error, int iteration)
{
	int d_width = cout.precision() + kDoubleExtra;

	/* write convergence output */
	if (fIterationOutputIncrement > 0 && ++fIterationOutputCount >= fIterationOutputIncrement)
	{
		fFEManager.WriteOutput(double(iteration));
		fIterationOutputCount = 0;
	}
	
	/* return value */
	SolutionStatusT status = kContinue;
	
	/* first pass */
	if (iteration == -1)
	{
		cout <<   "\n Group : " << fGroup+1 << '\n';
		cout <<   " Absolute error = " << error << '\n';
		cout <<   " Relative error :\n\n";

		fError0 = error;
		
		/* hit on first try */
		if (fError0 < fZeroTolerance)
		{
			cout << setw(kIntWidth) << 0 << '\t' << error << endl;
			status = kConverged;
		}
		else
		{
			cout.flush();
			status = kContinue;
		}
	}
	/* interpret error */
	else
	{
		double relerror = error/fError0;
		if (fVerbose) 
			cout << setw(kIntWidth) << iteration  << ": Relative error = "
			     << setw(d_width) << relerror << endl;

		/* diverging solution */	
		if (relerror > fDivTolerance)
		{
			cout << "\n NLSolver::ExitIteration: diverging solution detected" << endl;			
			status = kFailed;
		}
		/* required number of iterations */
		else if (iteration < fMinIterations-1)
		{
			status = kContinue;
		}
		/* converged */
		else if (relerror < fTolerance || error < fZeroTolerance)
		{
			if (!fVerbose)
				cout << setw(kIntWidth) << iteration  << ": Relative error = " 
				     << setw(d_width) << relerror << " (converged)\n";
	
			fFEManager.Output() << "\n Converged at time = " << fFEManager.Time() << endl;
			status = kConverged;
		}
		/* iteration limit hit */
		else if (iteration >= fMaxIterations)
		{
			cout << setw(kIntWidth) << fNumIteration 
			     << ": Relative error = " << setw(d_width) << error/fError0 << '\n';	
			cout << "\n NLSolver::ExitIteration: max iterations hit" << endl;
			status = kFailed;
		}
		/* continue iterations */
		else
			status = kContinue;
	}

	return status;
}

/* do one iteration of the solution procedure */
void NLSolver::Iterate(void)
{
	/* solve equation system */
	if (!fLHS->Solve(fRHS)) ExceptionT::BadJacobianDet("NLSolver::Iterate");

	/* apply update to system */
	Update(fRHS, NULL);
}
