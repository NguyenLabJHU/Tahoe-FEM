/* $Id: NOXSolverT.cpp,v 1.2 2002-04-02 23:30:55 paklein Exp $ */
#include "NOXSolverT.h"

/* optional */
#ifdef __NOX__

#include "fstreamT.h"
#include "FEManagerT.h"
#include "ControllerT.h"
#include "NOX_Solver_Manager.H"
#include "NOX_Tahoe_Group.h"
#include "dArrayT.h"

/* NOX headers */
#include "NOX_Parameter_List.H"
#include "NOX_Status_AbsResid.H"
#include "NOX_Status_RelResid.H"
#include "NOX_Status_MaxIters.H"
#include "NOX_Status_Combo.H"
#include "NOX_Status_Combo.H"

using namespace NOX::Status;
using namespace NOX::Tahoe;
using namespace NOX::Solver;

inline static int Max(int a, int b) { return (a > b) ? a : b; };
inline static double Min(double a, double b) { return (a < b) ? a : b; };
inline static double Max(double a, double b) { return (a > b) ? a : b; };

/* constructor */
NOXSolverT::NOXSolverT(FEManagerT& fe_manager):
	SolverT(fe_manager),
	fNOXParameters(NULL),
	fMaxIterations(-1),
	fAbsResidual(1.1),
	fRelResidual(1.1),
	fQuickSolveTol(-1),
	fQuickSeriesTol(-1),
	fIterationOutputIncrement(-1),
	fIterationOutputCount(0)
{
	/* parameter stream */
	ifstreamT& in = fe_manager.Input();
	
	/* parameter lists */
	ArrayT<StringT> nox_parameters;
	ArrayT<StringT> nox_search_parameters;
	ArrayT<StringT> nox_direction_parameters;

	/* read */
	in >> fMaxIterations;
	in >> fAbsResidual;
	in >> fRelResidual;
	in >> fQuickSolveTol;
	in >> fQuickSeriesTol;
	in >> fIterationOutputIncrement;

	int num_parameters = -1;
	in >> num_parameters; num_parameters = Max(0, num_parameters); /* solver */
	if (num_parameters == 0) {
		num_parameters = 1;
		nox_parameters.Dimension(1);
		nox_parameters[0] = "Method Newton";
	} else {
		nox_parameters.Dimension(num_parameters);
		for (int i = 0; i < nox_parameters.Length(); i++)
			nox_parameters[i].GetLineFromStream(in);
	}

	in >> num_parameters; num_parameters = Max(0, num_parameters); /* direction */
	if (num_parameters == 0) {
		num_parameters = 1;
		nox_direction_parameters.Dimension(1);
		nox_direction_parameters[0] = "Method Newton";
	} else {
		nox_direction_parameters.Dimension(num_parameters);
		for (int i = 0; i < nox_direction_parameters.Length(); i++)
			nox_direction_parameters[i].GetLineFromStream(in);
	}

	in >> num_parameters; num_parameters = Max(0, num_parameters); /* line search */
	if (num_parameters == 0) {
		num_parameters = 1;
		nox_search_parameters.Dimension(1);
		nox_search_parameters[0] = "Method \"Full Step\"";
	} else {
		nox_search_parameters.Dimension(num_parameters);
		for (int i = 0; i < nox_search_parameters.Length(); i++)
			nox_search_parameters[i].GetLineFromStream(in);
	}

	/* filter */
	fMaxIterations = Max(1, fMaxIterations);
	fAbsResidual = Max(0.0, fAbsResidual);
	fRelResidual = Min(1.0, fRelResidual);
	fQuickSolveTol = Max(0, fQuickSolveTol);
	fQuickSeriesTol = Max(0, fQuickSeriesTol);
	fIterationOutputIncrement = Max(0, fIterationOutputIncrement);

	/* echo */
	ofstreamT& out = fe_manager.Output();
	out << " Maximum number of iterations. . . . . . . . . . = " << fMaxIterations << '\n';	
	out << " Absolute convergence tolerance. . . . . . . . . = " << fAbsResidual << '\n';	
	out << " Relative convergence tolerance. . . . . . . . . = " << fRelResidual << '\n';	
	out << " Quick solution iteration count. . . . . . . . . = " << fQuickSolveTol  << '\n';	
	out << " Number of quick solutions before step increase. = " << fQuickSeriesTol << '\n';	
	out << " Iteration output increment. . . . . . . . . . . = " << fIterationOutputIncrement << '\n';

	/* echo NOX parameters */
	out << " Number of NOX Solver parameters . . . . . . . . = " << nox_parameters.Length() << '\n';
	fNOXParameters = new NOX::Parameter::List();
	for (int i = 0; i < nox_parameters.Length(); i++) {

		StringT& line = nox_parameters[i];

		/* resolve */
		int count;
		StringT param, value;
		param.FirstWord(line, count, false);
		line.Drop(count);
		value.FirstWord(line, count, false);

		/* echo */
		out << "\t\"" << param << "\" = \"" << value << '\"' << endl;
		
		/* store */
		fNOXParameters->setParameter(param.Pointer(), value);
	}

	/* echo NOX direction parameters */
	out << " Number of NOX Direction parameters. . . . . . . = " << nox_direction_parameters.Length() << '\n';
	for (int i = 0; i < nox_direction_parameters.Length(); i++) {

		StringT& line = nox_direction_parameters[i];

		/* resolve */
		int count;
		StringT param, value;
		param.FirstWord(line, count, false);
		line.Drop(count);
		value.FirstWord(line, count, false);

		/* echo */
		out << "\t\"" << param << "\" = \"" << value << '\"' << endl;
		
		/* store */
		fNOXParameters->sublist("Direction").setParameter(param.Pointer(), value);
	}

	/* echo NOX search parameters */
	out << " Number of NOX Direction parameters. . . . . . . = " << nox_search_parameters.Length() << '\n';
	for (int i = 0; i < nox_search_parameters.Length(); i++) {

		StringT& line = nox_search_parameters[i];

		/* resolve */
		int count;
		StringT param, value;
		param.FirstWord(line, count, false);
		line.Drop(count);
		value.FirstWord(line, count, false);

		/* echo */
		out << "\t\"" << param << "\" = \"" << value << '\"' << endl;
		
		/* store */
		fNOXParameters->sublist("Line Search").setParameter(param.Pointer(), value);
	}
}

/* destructor */
NOXSolverT::~NOXSolverT(void)
{
	delete fNOXParameters;
}

/* generate the solution for the current time sequence */
void NOXSolverT::Run(void)
{
	/* time integrator */
	const ControllerT* controller = fFEManager.Controller();
	if (!controller) throw eGeneralFail;

	/* check */
	if (!fLHS) {
		cout << "\n NOXSolverT::NOXSolverT: global matrix not set" << endl;
		throw eGeneralFail;
	}	

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

			/* set up group */
			dArrayT u(fRHS.Length());
			fFEManager.GetUnknowns(controller->OrderOfUnknown(), u);
			NOX::Tahoe::Group group(*this, u, *fLHS);
			fLastSolution = u;
			
			/* compute initial residual */
			if (!group.computeRHS()) {
				cout << "\n NOXSolverT::Run: unable to compute initial residual" << endl;
				throw eGeneralFail;
			}
			double error0 = group.getNormRHS();
			cout << " Error = " << error0 << endl;
			
			/* construct combination stopping criteria */
			NOX::Status::AbsResid abs_resid(fAbsResidual);
			NOX::Status::RelResid rel_resid(error0, fRelResidual);
			NOX::Status::MaxIters max_iters(fMaxIterations);
			NOX::Status::Combo combo(NOX::Status::Combo::OR);
			combo.addTest(abs_resid);
			combo.addTest(rel_resid);
			combo.addTest(max_iters);

			/* set solver */
			NOX::Solver::Manager nox(group, combo, *fNOXParameters);

			/* solve */
			NOX::Status::StatusType nox_status;
			try {
				nox_status = nox.iterate();
				while (nox_status == NOX::Status::Unconverged) {
					cout << '\t' << group.getNormRHS()/error0 << endl;
					nox_status = nox.iterate();
				}
				cout << '\t' << group.getNormRHS()/error0 << endl;
			}
			catch (int error) { /* Tahoe throws int's */
				cout << "\n NOXSolverT::Run: exception during solve: " 
				     << fFEManager.Exception(error) << endl;
				throw error;
			}
			catch (const char* error) { /* NOX throws strings */
				cout << "\n NOXSolverT::Run: exception during solve: " << error << endl;
				throw eGeneralFail;
			}

			/* what to do next */
			SolutionStatusT status = kFailed;
			switch (nox_status) {				
				case NOX::Status::Converged:
					status = DoConverged(); /* "relax" */
					break;
					
				case NOX::Status::Unconverged:
				case NOX::Status::Failed:
					status = kFailed;
					break;
			
				default:
					cout << "\n NOXSolverT::Run: unrecognized exit status: " << status << endl;
					throw eGeneralFail;
			}
				
			/* close iteration output */	
			CloseIterationOutput();
			
			/* close step */
			if (status == kConverged)
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

/* (re-)configure the global equation system */
void NOXSolverT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	SolverT::Initialize(tot_num_eq, loc_num_eq, start_eq);

//TEMP - need this function for anything?
}

/* compute RHS for the given update to the solution vector dx */
bool NOXSolverT::computeRHS(const dArrayT& x, dArrayT& rhs)
{
	/* compute the update vector and apply */
	dArrayT update(x.Length());
	update.DiffOf(x, fLastSolution);
	fFEManager.Update(update);
	update.Free();
	
	/* store new configuration */
	fLastSolution = x;

	/* set target */
	rhs.Swap(fRHS);
	
	/* calculate */
	try {
		fRHS = 0.0;
		fFEManager.FormRHS();	
	}
	catch (int exception) {
		cout << "\n NOXSolverT::computeRHS: exception: "
		     << fFEManager.Exception(exception) << endl;	
		rhs.Swap(fRHS); /* restore target */
		return false;
	}
	
	/* restore target */
	rhs.Swap(fRHS);
	return true;
}
  
/* compute the Jacobian */
bool NOXSolverT::computeJacobian(GlobalMatrixT& jacobian)
{
	/* set target */
	GlobalMatrixT* temp = fLHS;
	fLHS = &jacobian;
	
	/* calculate */
	try {
		fLHS->Clear();
		fFEManager.FormLHS();	
	}
	catch (int exception) {
		cout << "\n NOXSolverT::computeJacobian: exception: "
		     << fFEManager.Exception(exception) << endl;	
		fLHS = temp; /* restore target */
		return false;
	}
	
	/* restore target */
	fLHS = temp;
	return true;
}

/*************************************************************************
* Private
*************************************************************************/

/* success */
NOXSolverT::SolutionStatusT NOXSolverT::DoConverged(void)
{
//TEMP - nothing yet. Need to add "relax" like NLSolver.
	return kConverged;
}

/* solution failed */
void NOXSolverT::DoNotConverged(void)
{
	/* message */
	cout << "\n NOXSolverT::DoNotConverged: resetting step, cutting load set" << endl;

	/* step back to last converged */
	fFEManager.ResetStep();
	
	/* cut load increment */
	fFEManager.DecreaseLoadStep();
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

#endif /* __NOX__ */
