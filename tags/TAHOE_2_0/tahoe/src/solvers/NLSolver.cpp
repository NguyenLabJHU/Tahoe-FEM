/* $Id: NLSolver.cpp,v 1.31 2004-07-15 08:31:50 paklein Exp $ */
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
	fMinIterations(-1),
	fZeroTolerance(0.0),
	fTolerance(0.0),
	fDivTolerance(-1.0),
	fQuickConvCount(0),
	fIterationOutputCount(0),
	fVerbose(1),
	fQuickSolveTol(-1),
	fQuickSeriesTol(-1),
	fIterationOutputIncrement(0)
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
		
	/* loop on error */
	double error = Residual(fRHS);
	SolutionStatusT solutionflag = ExitIteration(error, fNumIteration);
	int num_iterations = 0;
	while (solutionflag == kContinue &&
		(max_iterations == -1 || num_iterations++ < max_iterations))
	{
		error = SolveAndForm(fNumIteration);
		solutionflag = ExitIteration(error, fNumIteration);
	}

	/* found solution - check relaxation */
	if (solutionflag == kConverged)
		solutionflag = DoConverged();
			
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
}

/* error handler */
void NLSolver::ResetStep(void)
{
	/* inherited */
	SolverT::ResetStep();

	/* reset count to increase load step */
	fQuickConvCount = 0;

	/* restore output */
	CloseIterationOutput();
}

/* (re-)set the reference error */
void NLSolver::SetReferenceError(double error)
{
	cout <<   "\n Group : " << fGroup+1 << '\n' << " Absolute error = " << error << endl;
	fError0 = error;
}

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
	min_iterations.SetDefault(0);
	list.AddParameter(min_iterations);

	list.AddParameter(fZeroTolerance, "abs_tolerance");
	list.AddParameter(fTolerance, "rel_tolerance");
	list.AddParameter(fDivTolerance, "divergence_tolerance");

	ParameterT quick_solve_iter(fQuickSolveTol, "quick_solve_iter");
	quick_solve_iter.SetDefault(6);
	list.AddParameter(quick_solve_iter);

	ParameterT quick_solve_count(fQuickSeriesTol, "quick_solve_count");
	quick_solve_count.SetDefault(3);
	list.AddParameter(quick_solve_count);

	ParameterT output_inc(fIterationOutputIncrement, "output_inc");
	output_inc.SetDefault(0);
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

/* relax system */
NLSolver::SolutionStatusT NLSolver::Relax(int newtancount)
{	
	cout <<   "\n Relaxation:" << '\n';

	/* iteration count */
	int iteration = -1;
	int count = newtancount - 1;

	/* form the first residual force vector */
	fLHS_lock = kOpen; /* LHS open for assembly, too! */
	fRHS_lock = kOpen;
	fRHS = 0.0;
	fFEManager.FormRHS(Group());	
	fLHS_lock = kLocked;
	fRHS_lock = kLocked;

	double error = Residual(fRHS);
		
	/* loop on error */
	SolutionStatusT solutionflag = ExitIteration(error, iteration);
	while (solutionflag == kContinue)
	{
		if (++count == newtancount) {	
			fLHS_update = true;
			count      = 0;
		}
		else fLHS_update = false;
			
		error = SolveAndForm(iteration);
		solutionflag = ExitIteration(error, iteration);
	}

	return solutionflag;
}

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
		else if (iteration > fMaxIterations)
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

/* form and solve the equation system */
double NLSolver::SolveAndForm(int& iteration)
{		
	/* form the stiffness matrix (must be cleared previously) */
	if (fLHS_update) {
		fLHS_lock = kOpen;
		fFEManager.FormLHS(Group(), GlobalT::kNonSymmetric);
		fLHS_lock = kLocked;
		
		/* compare with approxumate LHS */
		if (fLHS->CheckCode() == GlobalMatrixT::kCheckLHS) {
		
			/* compute approximate LHS */
			const GlobalMatrixT* approx_LHS = ApproximateLHS(*fLHS);
			
			/* compare */
			CompareLHS(*fLHS, *approx_LHS);
			
			/* clean-up */
			delete approx_LHS;
		}
	}

	/* solve equation system */
	if (!fLHS->Solve(fRHS)) ExceptionT::BadJacobianDet("NLSolver::SolveAndForm");

	/* apply update to system */
	Update(fRHS, NULL);

	/* recalculate residual */
	iteration++;
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

	/* could combine residual magnitude with update magnitude
	 * e = a1 |R| + a2 |delta_d|  --> not implemented */	
	return Residual(fRHS);
}
