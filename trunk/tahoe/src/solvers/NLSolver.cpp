/* $Id: NLSolver.cpp,v 1.17 2002-11-28 17:30:31 paklein Exp $ */
/* created: paklein (07/09/1996) */

#include "NLSolver.h"

#include <iostream.h>
#include <math.h>

#include "fstreamT.h"
#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "FEManagerT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* constructor */
NLSolver::NLSolver(FEManagerT& fe_manager, int group):
	SolverT(fe_manager, group),
	fMaxIterations(-1),
	fZeroTolerance(0.0),
	fTolerance(0.0),
	fDivTolerance(-1.0),
	fQuickConvCount(0),
	fIterationOutputCount(0),
	fVerbose(1)
{
	ifstreamT& in = fFEManager.Input();
	
	in >> fMaxIterations;	
	in >> fZeroTolerance;
	in >> fTolerance;
	in >> fDivTolerance;
	in >> fQuickSolveTol;
	in >> fQuickSeriesTol;
	in >> fIterationOutputIncrement;
	
	/* step increase disabled */
	if (fQuickSolveTol == -1) fQuickSeriesTol = -1;

	/* print parameters */
	ostream& out = fFEManager.Output();
	
	out << "\n O p t i m i z a t i o n   P a r a m e t e r s :\n\n";
	out << " Maximum number of iterations. . . . . . . . . . = " << fMaxIterations  << '\n';
	out << " Absolute convergence tolerance. . . . . . . . . = " << fZeroTolerance  << '\n';	
	out << " Relative convergence tolerance. . . . . . . . . = " << fTolerance      << '\n';	
	out << " Divergence tolerance. . . . . . . . . . . . . . = " << fDivTolerance   << '\n';	
	out << " Quick solution iteration count. (-1 to disable) = " << fQuickSolveTol  << '\n';	
	out << " Number of quick solutions before step increase. = " << fQuickSeriesTol << '\n';	
	out << " Iteration output print increment. . . . . . . . = " << fIterationOutputIncrement << endl;	
	
	/* checks */
	if (fMaxIterations < 0) throw ExceptionT::kBadInputValue;
	if (fZeroTolerance < 0.0 || fZeroTolerance > 1.0)
	{
		cout << "\n NLSolver::NLSolver: absolute convergence tolerance is out of\n"
		     <<   "    range: 0 <= tol <= 1: " << fZeroTolerance << endl;
		throw ExceptionT::kBadInputValue;
	}
	if (fTolerance < 0.0 || fTolerance > 1.0)
	{
		cout << "\n NLSolver::NLSolver: relative convergence tolerance is out of\n"
		     <<   "    range: 0 <= tol <= 1: " << fTolerance << endl;
		throw ExceptionT::kBadInputValue;
	}
	if (fDivTolerance < 0)  throw ExceptionT::kBadInputValue;
	if (fQuickSolveTol  != -1 && fQuickSolveTol  < 1) throw ExceptionT::kBadInputValue;
	if (fQuickSeriesTol != -1 && fQuickSeriesTol < 1) throw ExceptionT::kBadInputValue;
	if (fIterationOutputIncrement < 0)
	{
		cout << "\n NLSolver::NLSolver: expecting iteration output increment < 0: "
		     << fIterationOutputIncrement << endl;
		throw ExceptionT::kBadInputValue;
	}
	
	/* console variables */
	iAddVariable("max_iterations", fMaxIterations);
	iAddVariable("abs_tolerance", fZeroTolerance);
	iAddVariable("rel_tolerance", fTolerance);
	iAddVariable("div_tolerance", fDivTolerance);
	iAddVariable("iteration_output_inc", fIterationOutputIncrement);
}

/* generate the solution for the current time sequence */
SolverT::SolutionStatusT NLSolver::Solve(int num_iterations)
{
	try
	{ 	
	/* reset iteration count */
	fNumIteration = -1;

	/* open iteration output */
	InitIterationOutput();

	/* form the first residual force vector */
	fRHS = 0.0;
	fFEManager.FormRHS(Group());	
	double error = Residual(fRHS);
			
	/* loop on error */
	SolutionStatusT solutionflag = ExitIteration(error);
	while (solutionflag == kContinue &&
		(num_iterations == -1 || fNumIteration < num_iterations))
	{
		error = SolveAndForm(true);
		solutionflag = ExitIteration(error);
	}

	/* found solution - check relaxation */
	if (solutionflag == kConverged)
		solutionflag = DoConverged();
				
	/* close iteration output */	
	CloseIterationOutput();
			
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

/* error handler */
void NLSolver::ResetStep(void)
{
	/* inherited */
	SolverT::ResetStep();

	/* reset count to increase load step */
	fQuickConvCount = 0;

	/* restore output */
	if (fIterationOutputIncrement > 0) fFEManager.RestoreOutput();
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
		root.Root(fFEManager.Input().filename());
		
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

	/* reset iteration count */
	fNumIteration = -1;
		
	int count = newtancount - 1;

	/* form the first residual force vector */
	fRHS = 0.0;
	fFEManager.FormRHS(Group());	
	double error = Residual(fRHS);
		
	/* loop on error */
	SolutionStatusT solutionflag = ExitIteration(error);
	while (solutionflag == kContinue)
	{
	    int newtangent = 0;
		if (++count == newtancount)
		{	
			newtangent = 1;
			count      = 0;
		}
			
		error = SolveAndForm(newtangent);
		solutionflag = ExitIteration(error);
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
NLSolver::SolutionStatusT NLSolver::ExitIteration(double error)
{
	/* iteration count */
	++fNumIteration;

	/* write convergence output */
	if (++fIterationOutputCount == fIterationOutputIncrement)
	{
		fFEManager.WriteOutput(double(fNumIteration));
		fIterationOutputCount = 0;
	}
	
	/* return value */
	SolutionStatusT status = kContinue;
	
	/* first pass */
	if (fNumIteration == 0)
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
	/* iteration limit hit */
	else if (fNumIteration > fMaxIterations)
	{
		cout << "\n NLSolver::ExitIteration: max iterations hit" << endl;
		status = kFailed;
	}
	/* interpret error */
	else
	{
		int d_width = cout.precision() + kDoubleExtra;
		double relerror = error/fError0;
		if (fVerbose) cout << setw(kIntWidth) << fNumIteration 
		                   << ": Relative error = "
		                   << setw(d_width) << relerror << endl;

		/* diverging solution */	
		if (relerror > fDivTolerance)
		{
			cout << "\n NLSolver::ExitIteration: diverging solution detected" << endl;			
			status = kFailed;
		}
		/* converged */
		else if (relerror < fTolerance || error < fZeroTolerance)
		{
			if (!fVerbose)
			{
				cout << setw(kIntWidth) << fNumIteration;
				cout << setw(d_width)   << relerror << '\n';
			}
	
			fFEManager.Output() << "\n Converged at time = " << fFEManager.Time() << endl;
			status = kConverged;
		}
		/* continue iterations */
		else
			status = kContinue;
	}

	return status;
}

/* form and solve the equation system */
double NLSolver::SolveAndForm(bool newtangent)
{		
	/* form the stiffness matrix */
	if (newtangent)
	{
		fLHS->Clear();
		fFEManager.FormLHS(Group(), GlobalT::kNonSymmetric);
	}
		 		
	/* solve equation system */
	if (!fLHS->Solve(fRHS)) throw ExceptionT::kBadJacobianDet;

	/* apply update to system */
	Update(fRHS, NULL);
								
	/* compute new residual */
	fRHS = 0.0;
	fFEManager.FormRHS(Group());

	/* combine residual magnitude with update magnitude */
	/* e = a1 |R| + a2 |delta_d|                        */
	//not implemented!
			
	return Residual(fRHS);
}
