/* $Id: NLSolver.cpp,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (07/09/1996)                                          */

#include "NLSolver.h"

#include <iostream.h>
#include <math.h>

#include "fstreamT.h"
#include "Constants.h"
#include "ExceptionCodes.h"
#include "FEManagerT.h"

/* constructor */
NLSolver::NLSolver(FEManagerT& fe_manager):
	SolverT(fe_manager),
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
	if (fMaxIterations < 0) throw eBadInputValue;
	if (fZeroTolerance < 0.0 || fZeroTolerance > 1.0)
	{
		cout << "\n NLSolver::NLSolver: absolute convergence tolerance is out of\n"
		     <<   "    range: 0 <Êtol < 1: " << fZeroTolerance << endl;
		throw eBadInputValue;
	}
	if (fTolerance < 0.0 || fTolerance > 1.0)
	{
		cout << "\n NLSolver::NLSolver: relative convergence tolerance is out of\n"
		     <<   "    range: 0 <Êtol < 1: " << fTolerance << endl;
		throw eBadInputValue;
	}
	if (fDivTolerance < 0)  throw eBadInputValue;
	if (fQuickSolveTol  != -1 && fQuickSolveTol  < 1) throw eBadInputValue;
	if (fQuickSeriesTol != -1 && fQuickSeriesTol < 1) throw eBadInputValue;
	if (fIterationOutputIncrement < 0)
	{
		cout << "\n NLSolver::NLSolver: expecting iteration output increment < 0: "
		     << fIterationOutputIncrement << endl;
		throw eBadInputValue;
	}
	
	/* console variables */
	iAddVariable("max_iterations", fMaxIterations);
	iAddVariable("abs_tolerance", fZeroTolerance);
	iAddVariable("rel_tolerance", fTolerance);
	iAddVariable("div_tolerance", fDivTolerance);
	iAddVariable("iteration_output_inc", fIterationOutputIncrement);
}

/* generate the solution for the current time sequence */
void NLSolver::Run(void)
{
	/* solve load sequence */
	while (Step())
	{			
		/* residual loop */
		try
		{ 	
			/* apply boundary conditions */
			fFEManager.InitStep();

			/* open iteration output */
			InitIterationOutput();

			/* form the first residual force vector */
			fRHS = 0.0;
			fFEManager.FormRHS();	
			double error = Residual(fRHS);
			
			/* loop on error */
			IterationStatusT solutionflag = ExitIteration(error);
			while (solutionflag == kContinue)
			{
				error = SolveAndForm(true);
				solutionflag = ExitIteration(error);
			}

			/* found solution - check relaxation */
			if (solutionflag == kConverged)
				solutionflag = DoConverged();
				
			/* close iteration output */	
			CloseIterationOutput();
			
			/* close step */
			if (solutionflag == kConverged)
				fFEManager.CloseStep();
			/* solution procedure failed */
			else
				DoNotConverged();
		}
		
		catch (int code)
		{
			cout << "\n NLSolver::Run: exception at step number "
			     << fFEManager.StepNumber() << " with step "
			     << fFEManager.TimeStep() << endl;
			fFEManager.HandleException(code);
		}
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
NLSolver::IterationStatusT NLSolver::DoConverged(void)
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
	GlobalT::RelaxCodeT relaxcode = fFEManager.RelaxSystem();
	while (relaxcode != GlobalT::kNoRelax)
	{	
		/* reset global equations */
		if (relaxcode == GlobalT::kReEQ ||
		    relaxcode == GlobalT::kReEQRelax)
			fFEManager.Reinitialize();
						
		/* new equilibrium */
		if (relaxcode == GlobalT::kRelax ||
		    relaxcode == GlobalT::kReEQRelax)
		{
			if (Relax() == kFailed)
				return kFailed;
	   		else
				relaxcode = fFEManager.RelaxSystem();
		}
		else
			relaxcode = GlobalT::kNoRelax;
	}

	/* success */
	return kConverged;						
}

void NLSolver::DoNotConverged(void)
{
	/* step back to last converged */
	fFEManager.ResetStep();
	
	/* cut load increment */
	fFEManager.DecreaseLoadStep();
}

/* divert output for iterations */
void NLSolver::InitIterationOutput(void)
{
	if (fIterationOutputIncrement > 0)
	{
		/* root of output files */
		StringT root;
		root.Root(fFEManager.Input().filename());
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
	fFEManager.Update(update);
}

/* relax system */
NLSolver::IterationStatusT NLSolver::Relax(int newtancount)
{	
	cout <<   "\n Relaxation:" << '\n';

	/* reset iteration count */
	fNumIteration = -1;
		
	int count = newtancount - 1;

	/* form the first residual force vector */
	fRHS = 0.0;
	fFEManager.FormRHS();	
	double error = Residual(fRHS);
		
	/* loop on error */
	IterationStatusT solutionflag = ExitIteration(error);
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

/* advance to next load step. Returns 0 if there are no more
* steps. Overload to add class dependent initializations */
int NLSolver::Step(void)
{
	/* reset iteration count */
	fNumIteration = -1;
		
	/* inherited */
	return SolverT::Step();
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
NLSolver::IterationStatusT NLSolver::ExitIteration(double error)
{
	/* iteration count */
	++fNumIteration;

	/* write convergence output */
	if (++fIterationOutputCount == fIterationOutputIncrement)
	{
		fFEManager.WriteOutput(IOBaseT::kAtInc);
		fIterationOutputCount = 0;
	}
	
	/* return value */
	IterationStatusT status = kContinue;
	
	/* first pass */
	if (fNumIteration == 0)
	{
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
	else if (fNumIteration >= fMaxIterations)
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
		                   << setw(d_width)   << relerror << endl;

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
		fFEManager.FormLHS();
	}
		 		
	/* solve equation system */
	fLHS->Solve(fRHS);

	/* apply update to system */
	Update(fRHS, NULL);
								
	/* compute new residual */
	fRHS = 0.0;
	fFEManager.FormRHS();

	/* combine residual magnitude with update magnitude */
	/* e = a1 |R| + a2 |delta_d|                        */
	//not implemented!
			
	return Residual(fRHS);
}
