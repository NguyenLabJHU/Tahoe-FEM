/* $Id: iNLSolver_LS.cpp,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (01/01/2001)                                          */

#include "iNLSolver_LS.h"

#include <iostream.h>
#include <math.h>

#include "fstreamT.h"
#include "Constants.h"
#include "ExceptionCodes.h"
#include "FEManagerT.h"
#include "iConsoleT.h"

/* constructor */
iNLSolver_LS::iNLSolver_LS(FEManagerT& fe_manager):
	NLSolver_LS(fe_manager),
	fFormTangent(true),
	fLineSearch(true)
{
	/* add variables */
	iAddVariable("form_tangent", fFormTangent);
	iAddVariable("line_search", fLineSearch);

	/* add console commands */
	iAddCommand("ResetStep");
	iAddCommand("Iterate");
	iAddCommand("FormResidual");
	iAddCommand("InitStep");
	iAddCommand("Step");
}

/* interactive */
void iNLSolver_LS::Run(void)
{
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
	StringT arg;
	arg.Append("(", number_of_steps - step_number);
	arg.Append(")");
	iDoCommand("Step", arg);
}

/* console commands */
bool iNLSolver_LS::iDoCommand(const StringT& command, StringT& line)
{
	try
	{
		if (command == "Step")
		{
			/* resolve number of steps */
			int num_steps = 1;
			if (line[0] == '(' &&
			    !ResolveArgument(line, num_steps, &num_steps))
			{
				cout << "could not resolve integer argument from: \""
				     << line << '\"' << endl;
				return false;
			}
			
			/* run steps */
			return DoStep(num_steps);
		}
		else if (command == "InitStep")
			return DoInitStep();
		else if (command == "Iterate")
		{
			/* resolve argument */
			int num_iterations = 1;
			if (line[0] == '(' &&
			   !ResolveArgument(line, num_iterations, &num_iterations))
			{
				cout << "could not resolve integer argument from: \""
				     << line << '\"' << endl;
				return false;
			}
			
			/* message */
			if (DoIterate(num_iterations) == kFail)
				cout << "warning: iterations ended with FAIL" << endl;
			return true;
		}
		else if (command == "FormResidual")
		{
			/* compute new residual */
			fRHS = 0.0;
			fFEManager.FormRHS();
			cout << "residual norm = " << fRHS.Magnitude() << endl;
			return true;
		}
		else if (command == "ResetStep")
		{
			/* step back to last converged */
			fFEManager.ResetStep();

			/* initialize step */
			return DoInitStep();
		}
		else
			/* inherited */
			return SolverT::iDoCommand(command, line);
	}
	
	catch (int code)
	{
		cout << "\n iNLSolver_LS::iDoCommand: exception at step number "
		     << fFEManager.StepNumber() << " with step "
		     << fFEManager.TimeStep() << endl;
		fFEManager.HandleException(code);
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
		if (count == max_steps)
			return true;
		else
			return false;
	}
}

bool iNLSolver_LS::DoInitStep(void)
{
	/* close any iteration output */	
	CloseIterationOutput();

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
}

NLSolver::IterationStatusT iNLSolver_LS::DoIterate(int max_count)
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
	
				fRHS = 0.0;
				fFEManager.FormRHS();
	
				/* initial error */
				double error = Residual(fRHS);
				fIterationStatus = ExitIteration(error);
			}
				
			/* loop on error */
			int count = 0;
			while (fIterationStatus == kContinue && count++ < max_count)
			{
				bool form_tangent = (fNumIteration == 0) ? true : fFormTangent;
				double error = SolveAndForm(form_tangent);
				fIterationStatus = ExitIteration(error);
			}
		
			/* found solution - check relaxation */
			if (fIterationStatus == kConverged)
			{
				fIterationStatus = DoConverged();	
				fFEManager.CloseStep();
			}
			break;
		}
		default:
			cout << "unrecognized iteration status: " << fIterationStatus << endl;
	}
	
	/* return final status */
	return fIterationStatus;
}
