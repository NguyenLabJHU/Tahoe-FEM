/* $Id: MultiManagerT.cpp,v 1.9.12.1 2004-04-08 07:33:48 paklein Exp $ */
#include "MultiManagerT.h"

#ifdef BRIDGING_ELEMENT

#include "SolverT.h"
#include "DiagonalMatrixT.h"
#include "FEManagerT_bridging.h"
#include "NodeManagerT.h"
#include "OutputSetT.h"
#include "TimeManagerT.h"
#include "FieldT.h"

using namespace Tahoe;

/* constructor */
MultiManagerT::MultiManagerT(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	FEManagerT_bridging* fine, FEManagerT_bridging* coarse):
	FEManagerT(input, output, comm),
	fFine(fine),
	fCoarse(coarse),
	fDivertOutput(false)
{
	/* borrow parameters from coarse scale solver */
	fAnalysisCode = fCoarse->Analysis();
	fTimeManager = fCoarse->TimeManager();
	fOutputFormat = fCoarse->OutputFormat();

	/* don't compute initial conditions */
	fFine->SetComputeInitialCondition(false);
	fCoarse->SetComputeInitialCondition(false);
}

/* destructor */
MultiManagerT::~MultiManagerT(void)
{
	fTimeManager = NULL;
}

/* initialize members */
void MultiManagerT::Initialize(InitCodeT)
{
	const char caller[] = "MultiManagerT::Initialize";

	/* state */
	fStatus = GlobalT::kInitialization;

	/* configure projection/interpolation */
	int group = 0;
	int order1 = 0;
	StringT bridging_field = "displacement";
	bool make_inactive = true;
	//dArrayT mdmass;
	fFine->InitGhostNodes(fCoarse->ProjectImagePoints());
	fCoarse->InitInterpolation(fFine->GhostNodes(), bridging_field, *(fFine->NodeManager()));
	//fFine->LumpedMass(fFine->NonGhostNodes(), mdmass);
	fCoarse->InitProjection(*(fFine->CommManager()), fFine->NonGhostNodes(), bridging_field, *(fFine->NodeManager()), make_inactive);

	/* send coarse/fine output through the fFine output */
	int ndof = fFine->NodeManager()->NumDOF(group);
	ArrayT<StringT> labels(2*ndof);
	const char* coarse_labels[] = {"UC_X", "UC_Y", "UC_Z"};
	const char* fine_labels[] = {"UF_X", "UF_Y", "UF_Z"};
	int dex = 0;
	for (int i = 0; i < ndof; i++) labels[dex++] = coarse_labels[i];
	for (int i = 0; i < ndof; i++) labels[dex++] = fine_labels[i];
	const iArrayT& non_ghost_nodes = fFine->NonGhostNodes();
	fAtomConnectivities.Alias(non_ghost_nodes.Length(), 1, non_ghost_nodes.Pointer());
	OutputSetT output_set(GeometryT::kPoint, fAtomConnectivities, labels, false);
	fOutputID = fFine->RegisterOutput(output_set);

	/* set joint solver */
	int n1 = fFine->NumGroups();
	int n2 = fCoarse->NumGroups();
	if (n1 != n2) ExceptionT::GeneralFail(caller, "number of groups must match: %d != %d", n1, n2);
	fSolvers.Dimension(n1);
	fSolvers = NULL;
	SetSolver();
}

/* (re-)set the equation number for the given group */
void MultiManagerT::SetEquationSystem(int group)
{
	fFine->SetEquationSystem(group);
	fCoarse->SetEquationSystem(group);

	int neq1 = fFine->GlobalNumEquations(group);
	int neq2 = fCoarse->GlobalNumEquations(group);
	fGlobalNumEquations[group] = neq1 + neq2;

	/* set total equation numbers */
	fEqnos1.Dimension(neq1);
	fEqnos2.Dimension(neq2);
	fEqnos1.SetValueToPosition();
	fEqnos2.SetValueToPosition();
	fEqnos1 += 1;
	fEqnos2 += (1+neq1);

	/* final step in solver configuration */
	fSolvers[group]->Initialize(
		fGlobalNumEquations[group],
		fGlobalNumEquations[group],
		1);
}

/* initialize the current time increment for all groups */
ExceptionT::CodeT MultiManagerT::InitStep(void)
{
	ExceptionT::CodeT error = ExceptionT::kNoError;
	if (error == ExceptionT::kNoError) 
		error = fFine->InitStep();
	if (error == ExceptionT::kNoError) 
		error = fCoarse->InitStep();

	/* loop over solvers only */
	for (fCurrentGroup = 0; fCurrentGroup < NumGroups(); fCurrentGroup++)
		fSolvers[fCurrentGroup]->InitStep();
	fCurrentGroup = -1;

	/* OK */
	return error;
}

/* close the current time increment for all groups */
ExceptionT::CodeT MultiManagerT::CloseStep(void)
{
	/* write output from sub's */
	ExceptionT::CodeT error = ExceptionT::kNoError;
	if (error == ExceptionT::kNoError) 
		error = fFine->CloseStep();
	if (error == ExceptionT::kNoError)
		error = fCoarse->CloseStep();

	/* write coarse/fine output */
	if (fTimeManager->WriteOutput())
		WriteOutput(fTimeManager->Time());

	/* loop over solvers only */
	for (fCurrentGroup = 0; fCurrentGroup < NumGroups(); fCurrentGroup++)
		fSolvers[fCurrentGroup]->CloseStep();
	fCurrentGroup = -1;

	/* OK */
	return error;
}

/* called for all groups if the solution procedure for any group fails */
ExceptionT::CodeT MultiManagerT::ResetStep(void)
{
	ExceptionT::GeneralFail("MultiManagerT::ResetStep", "not implemented");
	return ExceptionT::kNoError;
}

/* compute LHS-side matrix and assemble to solver */
void MultiManagerT::FormLHS(int group, GlobalT::SystemTypeT sys_type) const
{
	const char caller[] = "MultiManagerT::FormLHS";

	/* state */
	SetStatus(GlobalT::kFormLHS);

	/* fine scale */
	SolverT* fine_solver = fFine->Solver(group);
	const GlobalMatrixT& lhs_1 = fFine->LHS(group);
	DiagonalMatrixT* diag_1 = TB_DYNAMIC_CAST(DiagonalMatrixT*, (DiagonalMatrixT*) &lhs_1);
	if (!diag_1) ExceptionT::GeneralFail(caller);
	diag_1->Clear();
	fine_solver->UnlockLHS();
	fFine->FormLHS(group, sys_type);
	fine_solver->LockLHS();
	dArrayT& vec_1 = diag_1->TheMatrix();
	fSolvers[group]->AssembleLHS(vec_1, fEqnos1);

	/* coarse scale */
	SolverT* coarse_solver = fCoarse->Solver(group);
	const GlobalMatrixT& lhs_2 = fCoarse->LHS(group);
	DiagonalMatrixT* diag_2 = TB_DYNAMIC_CAST(DiagonalMatrixT*, (DiagonalMatrixT*) &lhs_2);
	if (!diag_2) ExceptionT::GeneralFail(caller);
	diag_2->Clear();
	coarse_solver->UnlockLHS();
	fCoarse->FormLHS(group, sys_type);
	coarse_solver->LockLHS();
	dArrayT& vec_2 = diag_2->TheMatrix();
	fSolvers[group]->AssembleLHS(vec_2, fEqnos2);
}

/* compute RHS-side */
void MultiManagerT::FormRHS(int group) const
{
	const char caller[] = "MultiManagerT::FormLHS";

	/* state */
	SetStatus(GlobalT::kFormRHS);

	/* fine scale */
	SolverT* fine_solver = fFine->Solver(group);
	dArrayT& fine_rhs = (dArrayT&) fFine->RHS(group);
	fine_rhs = 0.0;
	fine_solver->UnlockRHS();
	fFine->FormRHS(group);
	fine_solver->LockRHS();
	fSolvers[group]->AssembleRHS(fine_rhs, fEqnos1);

	/* coarse scale */
	SolverT* coarse_solver = fCoarse->Solver(group);
	dArrayT& coarse_rhs = (dArrayT&) fCoarse->RHS(group);
	coarse_rhs = 0.0;
	coarse_solver->UnlockRHS();
	fCoarse->FormRHS(group);
	coarse_solver->LockRHS();
	fSolvers[group]->AssembleRHS(coarse_rhs, fEqnos2);
}

/* send update of the solution to the NodeManagerT */
void MultiManagerT::Update(int group, const dArrayT& update)
{
	StringT bridging_field = "displacement";
	int order = 0;

	int neq1 = fEqnos1.Length();
	int neq2 = fEqnos2.Length();

	/* update subs */
	dArrayT update_tmp;
	update_tmp.Alias(neq1, update.Pointer());
	fFine->Update(group, update_tmp);
	update_tmp.Alias(neq2, update.Pointer(neq1));
	fCoarse->Update(group, update_tmp);

	/* project fine scale solution on to the coarse grid */
	fCoarse->ProjectField(bridging_field, *(fFine->NodeManager()), order);

	/* interpolate coarse scale solution to the fine */
	fCoarse->InterpolateField(bridging_field, order, fFieldAtGhosts);
	fFine->SetFieldValues(bridging_field, fFine->GhostNodes(), order, fFieldAtGhosts);
}

/* system relaxation */
GlobalT::RelaxCodeT MultiManagerT::RelaxSystem(int group) const
{
	/* just call relax */
	GlobalT::RelaxCodeT relax = GlobalT::kNoRelax;
	relax = GlobalT::MaxPrecedence(relax, fCoarse->RelaxSystem(group));
	relax = GlobalT::MaxPrecedence(relax, fFine->RelaxSystem(group));

	return relax;
}

/* (temporarily) direct output away from main out */
void MultiManagerT::DivertOutput(const StringT& outfile)
{
	/* divert output file fine scale */
	StringT fine_outfile(outfile);
	fine_outfile.Append(".fine");
	fFine->DivertOutput(fine_outfile);

	/* divert output file coarse scale */
	StringT coarse_outfile(outfile);
	coarse_outfile.Append(".coarse");
	fCoarse->DivertOutput(coarse_outfile);
	
	fDivertOutput = true;
}

/* restore outputs to their regular destinations */
void MultiManagerT::RestoreOutput(void)
{
	fFine->RestoreOutput();
	fCoarse->RestoreOutput();
	fDivertOutput = false;
}

/* initiate the process of writing output */
void MultiManagerT::WriteOutput(double time)
{
	const char caller[] = "MultiManagerT::WriteOutput";

	/* compute the coarse scale part of the field */
	const NodeManagerT& nodes = *(fFine->NodeManager());
	int order = 0;
	int group = 0;
	StringT field = "displacement";
	int ndof = nodes.NumDOF(group);
	const iArrayT& source_points = fFine->NonGhostNodes();
	dArray2DT n_values(source_points.Length(), 2*ndof), coarse;
	fCoarse->CoarseField(field, nodes, order, coarse);
	if (source_points.Length() != coarse.MajorDim())
		ExceptionT::GeneralFail(caller);

	/* get the total field */
	const FieldT* source_field = nodes.Field(field);
	if (!source_field) ExceptionT::GeneralFail(caller, "could not resolve source field \"%s\"", field.Pointer());
	const dArray2DT& total = (*source_field)[order];

	/* compute the fine scale part of the field */
	for (int i = 0; i < source_points.Length(); i++)
	{
		const double* p_total = total(source_points[i]);
		double* p_crse = coarse(i);
		double* p_crse_out = n_values(i);
		double* p_fine_out = p_crse_out + ndof;
		for (int j = 0; j < ndof; j++)
		{
			p_crse_out[j] = p_crse[j];
			p_fine_out[j] = p_total[j] - p_crse[j]; 
		}
	}
	
	/* send result through output of fine scale solver */
	dArray2DT e_values;
	fFine->WriteOutput(fOutputID, n_values, e_values);
	
	/* write iteration output for sub's */
	if (fDivertOutput) {
		fFine->WriteOutput(time);
		fCoarse->WriteOutput(time);
	}
	
}

#endif /* BRIDGING_ELEMENT */
