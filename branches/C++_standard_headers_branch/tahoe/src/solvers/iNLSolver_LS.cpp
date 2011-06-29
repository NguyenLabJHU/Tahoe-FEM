/* $Id: iNLSolver_LS.cpp,v 1.17 2004-09-09 23:54:55 paklein Exp $ */
/* created: paklein (01/01/2001) */
#include "iNLSolver_LS.h"

#include <iostream.h>
#include <math.h>


#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "FEManagerT.h"
#include "iConsoleT.h"
#include "CommandSpecT.h"
#include "ArgSpecT.h"

/* matrix types */
#include "CCSMatrixT.h"
#include "CCNSMatrixT.h"

/* constructor */

using namespace Tahoe;

iNLSolver_LS::iNLSolver_LS(FEManagerT& fe_manager, int group):
	NLSolver_LS(fe_manager, group),
	fFormTangent(true),
	fLineSearch(true)
{
	/* add variables */
	iAddVariable("form_tangent", fFormTangent);
	iAddVariable("line_search", fLineSearch);

	/* add console commands */
	iAddCommand(CommandSpecT("ResetStep"));
	iAddCommand(CommandSpecT("FormResidual"));
	iAddCommand(CommandSpecT("PivotInfo"));
	iAddCommand(CommandSpecT("InitStep"));

	CommandSpecT iterate("Iterate");
	ArgSpecT num_its(ArgSpecT::int_);
	num_its.SetPrompt("number of iterations");
	num_its.SetDefault(1);
	iterate.AddArgument(num_its);
	iAddCommand(iterate);

	CommandSpecT step("Step");
	ArgSpecT num_steps(ArgSpecT::int_);
	num_steps.SetPrompt("number of steps");
	num_steps.SetDefault(1);
	step.AddArgument(num_steps);
	iAddCommand(step);
}

/* interactive */
SolverT::SolutionStatusT iNLSolver_LS::Solve(int)
{
//TEMP - revised solvers means this interactive part needs to change
cout << "\n iNLSolver_LS::Solve: not updated for multifield" << endl;
return kFailed;

#if 0
	/* initial state */
	fIterationStatus = kConverged;
	
	/* run console */
	StringT log_file;
	log_file.Root(fFEManager.Input().filename());
	log_file.Append(".console.log");
	cout << "\n#### console input being appended to \"" << log_file
	     << "\" ####\n" << endl;
	iConsoleT(log_file, *this);
	
	/* finish time sequence */
	const int step_number = fFEManager.StepNumber();
	const int number_of_steps = fFEManager.NumberOfSteps();
	
	/* get the command spec */
	CommandSpecT* step_command = iCommand("Step");
	if (!step_command) throw ExceptionT::kGeneralFail;
	step_command->Argument(0).SetValue(number_of_steps - step_number);
	
	/* execute */
	StringT line;
	iDoCommand(*step_command, line);
#endif
}

/* console commands */
bool iNLSolver_LS::iDoCommand(const CommandSpecT& command, StringT& line)
{
	try
	{
		if (command.Name() == "Step")
		{
			/* resolve number of steps */
			int num_steps;
			command.Argument(0).GetValue(num_steps);

			/* run steps */
			return DoStep(num_steps);
		}
		else if (command.Name() == "InitStep")
			return DoInitStep();
		else if (command.Name() == "Iterate")
		{
			/* resolve argument */
			int num_iterations;
			command.Argument(0).GetValue(num_iterations);
			
			/* message */
			if (DoIterate(num_iterations) == kFailed)
				cout << "warning: iterations ended with FAIL" << endl;
			return true;
		}
		else if (command.Name() == "FormResidual")
		{
			/* compute new residual */
			fRHS = 0.0;
			fFEManager.FormRHS(Group());
			cout << "residual norm = " << fRHS.Magnitude() << endl;
			return true;
		}
		else if (command.Name() == "ResetStep")
		{
			/* step back to last converged */
//			fFEManager.ResetStep();

			/* initialize step */
			return DoInitStep();
		}
		else if (command.Name() == "PivotInfo")
		{
#ifdef __NO_RTTI__
			cout << "command not available: requires RTTI" << endl;
			return false;
#else
			/* get matrix pointer */
			const CCSMatrixT*   CCS_mat = TB_DYNAMIC_CAST(const CCSMatrixT*, fLHS);
			const CCNSMatrixT* CCNS_mat = TB_DYNAMIC_CAST(const CCNSMatrixT*, fLHS);
			double min, max, abs_min, abs_max;
			if (CCNS_mat) CCNS_mat->FindMinMaxPivot(min, max, abs_min, abs_max);
			else if (CCS_mat) CCS_mat->FindMinMaxPivot(min, max, abs_min, abs_max);
			else
			{
				cout << "requires matrix type: " << kProfileSolver << endl;
				return false;
			}
			
			/* write results */
			int d_width = OutputWidth(cout, &min);
			cout << "\n Matrix pivots:\n" 
			     <<   "     min = " << setw(d_width) << min << '\n' 
			     <<   "     max = " << setw(d_width) << max << '\n' 
			     <<   "   |min| = " << setw(d_width) << abs_min << '\n' 
			     <<   "   |max| = " << setw(d_width) << abs_max << '\n';

			/* initialize step */
			return true;
#endif
		}
		else
			/* inherited */
			return SolverT::iDoCommand(command, line);
	}
	
	catch (ExceptionT::CodeT code)
	{
		cout << "\n iNLSolver_LS::iDoCommand: exception at step number "
		     << fFEManager.StepNumber() << " with step "
		     << fFEManager.TimeStep() << endl;
//		fFEManager.HandleException(code);
		return false;
	}
}

/*************************************************************************
* Private
*************************************************************************/

/* apply flags */
void iNLSolver_LS::Update(const dArrayT& update, const dArrayT* residual)
{
	if (fLineSearch)
		NLSolver_LS::Update(update, residual);
	else
		NLSolver::Update(update, residual);
}

/* commands */
bool iNLSolver_LS::DoStep(int max_steps)
{
	/* close out current step */
	if (fIterationStatus == kContinue)
	{
		max_steps--;
		DoIterate(fMaxIterations);
	}

	/* continue */
	if (fIterationStatus != kConverged)
		return false;
	else
	{
		int count = 0;
		while (count++ < max_steps && DoInitStep())
		{
			if (DoIterate(fMaxIterations) != kConverged)
				return false;
		}
		
		/* finished command */
		if (count == max_steps + 1)
			return true;
		else
			return false;
	}
}

bool iNLSolver_LS::DoInitStep(void)
{
	/* close any iteration output */	
	CloseIterationOutput();

	return false; //TEMP
#if 0
	if (Step())
	{
		/* apply boundary conditions */
		fFEManager.InitStep();
		fIterationStatus = kContinue;
		return true;
	}
	else
	{
		cout << "reached end of time sequence" << endl;
		return false;
	}
#endif
}

NLSolver::SolutionStatusT iNLSolver_LS::DoIterate(int max_count)
{
	/* no action */
	if (max_count < 1) return fIterationStatus;
	
	switch (fIterationStatus)
	{
		case kConverged:
			cout << "solution is converged" << endl;
			break;

		case kFailed:
			cout << "solution procedure failed, reset step" << endl;
			break;

		case kContinue:
		{
			/* first iteration */
			if (fNumIteration == -1)
			{
				/* open iteration output */
				InitIterationOutput();
	
				/* clear all */
				fRHS = 0.0;
				fLHS->Clear();
				
				/* form residual */
				fFEManager.FormRHS(Group());
	
				/* initial error */
				double error = Residual(fRHS);
				fIterationStatus = ExitIteration(error, fNumIteration);
			}
				
			/* loop on error */
			int count = 0;
			while (fIterationStatus == kContinue && count++ < max_count)
			{
				fLHS_update = (fNumIteration == 0) ? true : fFormTangent;

				/* recalculate */
				if (fLHS_update) {
					fLHS->Clear();
					fFEManager.FormLHS(Group(), fLHS->MatrixType());
				}

				/* update solution */
				Iterate();
				fNumIteration++;

				/* new error */
				double error = Residual(fRHS);		

				/* test for convergence */
				fIterationStatus = ExitIteration(error, fNumIteration);
			}
		
			/* found solution - check relaxation */
			if (fIterationStatus == kConverged)
			{
				fIterationStatus = DoConverged();	
			}
			break;
		}
		default:
			cout << "unrecognized iteration status: " << fIterationStatus << endl;
	}
	
	/* return final status */
	return fIterationStatus;
}
