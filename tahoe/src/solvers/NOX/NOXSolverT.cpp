/* $Id: NOXSolverT.cpp,v 1.9 2004-06-17 15:41:29 paklein Exp $ */
#include "NOXSolverT.h"

/* optional */
#ifdef __NOX__

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "FEManagerT.h"
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


using namespace Tahoe;
using namespace NOX::Status;
using namespace NOX::Solver;

inline static int Max(int a, int b) { return (a > b) ? a : b; };
inline static double Min(double a, double b) { return (a < b) ? a : b; };
inline static double Max(double a, double b) { return (a > b) ? a : b; };

/* constructor */
NOXSolverT::NOXSolverT(FEManagerT& fe_manager, int group, int unknowns_order):
	SolverT(fe_manager, group),
	fUnknownsOrder(unknowns_order),
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
SolverT::SolutionStatusT NOXSolverT::Solve(int num_iterations)
{
	try {

	/* check */
	if (!fLHS) {
		cout << "\n NOXSolverT::NOXSolverT: global matrix not set" << endl;
		throw ExceptionT::kGeneralFail;
	}	

	/* open iteration output */
	InitIterationOutput();

	/* set up group */
	dArrayT u(fRHS.Length());
	fFEManager.GetUnknowns(fGroup, fUnknownsOrder, u);
	Tahoe::Group group(*this, u, *fLHS);
	fLastSolution = u;
			
	/* compute initial residual */
	if (!group.computeRHS()) {
		cout << "\n NOXSolverT::Run: unable to compute initial residual" << endl;
		throw ExceptionT::kGeneralFail;
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
	fNumIteration = -1;
	try {
		nox_status = nox.iterate(); 
		fNumIteration++;
		while (nox_status == NOX::Status::Unconverged &&
			(num_iterations == -1 || fNumIteration < num_iterations)) {
			cout << '\t' << group.getNormRHS()/error0 << endl;
			nox_status = nox.iterate();
			fNumIteration++;
		}
		cout << '\t' << group.getNormRHS()/error0 << endl;
	}
	catch (ExceptionT::CodeT error) { /* Tahoe throws int's */
		cout << "\n NOXSolverT::Run: exception during solve: " 
	         << ExceptionT::ToString(error) << endl;
		throw error;
	}
	catch (const char* error) { /* NOX throws strings */
		cout << "\n NOXSolverT::Run: exception during solve: " << error << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* what to do next */
	SolutionStatusT status = kFailed;
	switch (nox_status) {				
		case NOX::Status::Converged:
			status = kConverged; /* "relax"? */
		break;
					
		case NOX::Status::Unconverged:
		case NOX::Status::Failed:
			status = kFailed;
			break;
			
		default:
			cout << "\n NOXSolverT::Run: unrecognized exit status: " << status << endl;
			throw ExceptionT::kGeneralFail;
	}
				
	/* close iteration output */	
	CloseIterationOutput();
			
	/* close step */
	return status;
	} /* end try */
		
	/* exception */
	catch (ExceptionT::CodeT code)
	{
		cout << "\n NOXSolverT::Run: exception at step number "
			 << fFEManager.StepNumber() << " with step "
			 << fFEManager.TimeStep() << endl;
		return kFailed;
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
	fFEManager.Update(Group(), update);
	update.Free();
	
	/* store new configuration */
	fLastSolution = x;

	/* set target */
	rhs.Swap(fRHS);
	
	/* calculate */
	try {
		/* open residual to assembly */
		fRHS_lock = kOpen;

		/* compute residual */
		fRHS = 0.0;
		fFEManager.FormRHS(Group());	

		/* lock residual */
		fRHS_lock = kLocked;
	}
	catch (ExceptionT::CodeT exception) {
		cout << "\n NOXSolverT::computeRHS: exception: "
		     << ExceptionT::ToString(exception) << endl;	
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
		/* open Jacobian to assembly */
		fLHS_lock = kOpen;

		/* compute Jacobian */
		fLHS->Clear();
		fFEManager.FormLHS(Group(), GlobalT::kNonSymmetric);	

		/* lock Jacobian */
		fLHS_lock = kLocked;
	}
	catch (ExceptionT::CodeT exception) {
		cout << "\n NOXSolverT::computeJacobian: exception: "
		     << ExceptionT::ToString(exception) << endl;	
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
