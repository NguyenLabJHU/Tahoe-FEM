/* $Id: NOXSolverT.cpp,v 1.1 2002-03-28 16:40:35 paklein Exp $ */
#include "NOXSolverT.h"
#include "fstreamT.h"
#include "FEManagerT.h"
#include "NOX_Solver_Manager.H"
#include "NOX_Tahoe_Group.h"

using namespace NOX::Status;
using namespace NOX::Tahoe;
using namespace NOX::Solver;

/* constructor */
NOXSolverT::NOXSolverT(FEManagerT& fe_manager):
	SolverT(fe_manager),
	fNOX(NULL),
	fGroup(NULL),
	fIterationOutputIncrement(0),
	fIterationOutputCount(0)
{

}

/* destructor */
NOXSolverT::~NOXSolverT(void)
{
	delete fNOX;
	delete fGroup;
}

/* generate the solution for the current time sequence */
void NOXSolverT::Run(void)
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
			SolutionStatusT solutionflag = kContinue;
#if 0
			while (solutionflag == kContinue)
			{
				error = SolveAndForm(true);
				solutionflag = ExitIteration(error);
			}
#endif
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
			cout << "\n NOXSolverT::Run: exception at step number "
			     << fFEManager.StepNumber() << " with step "
			     << fFEManager.TimeStep() << endl;
			fFEManager.HandleException(code);
		}
	}
}

/* error handler */
void NOXSolverT::ResetStep(void)
{
	// not implemented
	throw;
}

/* compute RHS for the given solution vector x */
bool NOXSolverT::computeRHS(const dArrayT& x, dArrayT& rhs)
{
#pragma unused(x)
#pragma unused(rhs)
	// compute update
	// apply to the system
	// set rhs as destination
	// evaluate residal
	return true;
}
  
/* compute the Jacobian given the specified solution vector x */
bool NOXSolverT::computeJacobian(const dArrayT& x, GlobalMatrixT& jacobian)
{
#pragma unused(x)
#pragma unused(jacobian)
	// compute update
	// apply to the system
	// set jacobian as destination
	// evaluate Jacobian
	return true;
}

/*************************************************************************
* Private
*************************************************************************/

/* success */
NOXSolverT::SolutionStatusT NOXSolverT::DoConverged(void)
{
	return kConverged;
}

/* solution failed */
void NOXSolverT::DoNotConverged(void)
{
	// not implemented
	throw;
}

/* divert output for iterations */
void NOXSolverT::InitIterationOutput(void)
{
	if (fIterationOutputIncrement > 0)
	{
		/* root of output files */
		StringT root;
		root.Root(fFEManager.Input().filename());
		
		/* remove processor designation */ 
		if (fFEManager.Size() > 1) root.Root();
		
		/* increment */
		root.Append(".", fFEManager.StepNumber());
		root.Append("of", fFEManager.NumberOfSteps());

		/* set temporary output */
		fFEManager.DivertOutput(root);
		
		/* reset count */
		fIterationOutputCount = 0;
	}
}

void NOXSolverT::CloseIterationOutput(void)
{
	if (fIterationOutputIncrement > 0)
		fFEManager.RestoreOutput();
}
