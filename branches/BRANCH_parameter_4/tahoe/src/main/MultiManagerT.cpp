/* $Id: MultiManagerT.cpp,v 1.11 2004-06-26 18:38:09 paklein Exp $ */
#include "MultiManagerT.h"

#ifdef BRIDGING_ELEMENT

#include "SolverT.h"
#include "DiagonalMatrixT.h"
#include "FEManagerT_bridging.h"
#include "NodeManagerT.h"
#include "OutputSetT.h"
#include "TimeManagerT.h"
#include "ParticlePairT.h"
#include "ifstreamT.h"

using namespace Tahoe;

/* constructor */
MultiManagerT::MultiManagerT(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	FEManagerT_bridging* fine, FEManagerT_bridging* coarse):
	FEManagerT(input, output, comm, fArgv),
	fFine(fine),
	fCoarse(coarse),
	fDivertOutput(false),
	fFineField(NULL),
	fCoarseField(NULL),
	fFineToCoarse(true),
	fCoarseToFine(true),
	fCorrectOverlap(true),
	fCBTikhonov(0.0),
	fK2(0.0)
{
	const char caller[] = "MultiManagerT::MultiManagerT";

	/* borrow parameters from coarse scale solver */
	fAnalysisCode = fCoarse->Analysis();
	fTimeManager = fCoarse->TimeManager();
	fOutputFormat = fCoarse->OutputFormat();

	/* don't compute initial conditions */
	fFine->SetComputeInitialCondition(false);
	fCoarse->SetComputeInitialCondition(false);

	StringT bridging_field = "displacement";
	fFineField = fFine->NodeManager()->Field(bridging_field);
	if (!fFineField) ExceptionT::GeneralFail(caller, "could not resolve fine scale \"%s\" field", bridging_field.Pointer());
	fCoarseField = fCoarse->NodeManager()->Field(bridging_field);
	if (!fFineField) ExceptionT::GeneralFail(caller, "could not resolve coarse scale \"%s\" field", bridging_field.Pointer());
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

	/* fine scale node manager */
	NodeManagerT& fine_node_manager = *(fFine->NodeManager());

	/* configure projection/interpolation */
	int group = 0;
	int order1 = 0;
	bool make_inactive = true;
	fFine->InitGhostNodes(fCoarse->ProjectImagePoints());
	fCoarse->InitInterpolation(fFine->GhostNodes(), fFineField->Name(), fine_node_manager);
	fCoarse->InitProjection(*(fFine->CommManager()), fFine->NonGhostNodes(), fFineField->Name(), fine_node_manager, make_inactive);

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

	/* read the cross term flags */
	ifstreamT& in = Input();
	in >> fFineToCoarse
	   >> fCoarseToFine
	   >> fCorrectOverlap;
	   
	/* enforce zero bond density in projected cells */
	if (fCoarseToFine)
		fCoarse->DeactivateFollowerCells();

	/* correct overlap */
	if (fCorrectOverlap) {
	
		fCBTikhonov = fK2 = -99;
		int nip = -99;
		double r = -99; /* augmented lagrangian regularization */
		in >> fCBTikhonov
		   >> fK2
		   >> nip
		   >> r;
		if (fCBTikhonov < 0.0 || fK2 < 0.0 || r < 0.0)
			ExceptionT::GeneralFail(caller, "regularization must be >= 0.0: %g, %g, %g", fCBTikhonov, fK2, r);
	
		const dArray2DT& fine_init_coords = fine_node_manager.InitialCoordinates();
		const ParticlePairT* particle_pair = fFine->ParticlePair();
		if (!particle_pair) ExceptionT::GeneralFail(caller, "could not resolve ParticlePairT");
		
		if (fCorrectOverlap == 1)
			fCoarse->CorrectOverlap_1(particle_pair->Neighbors(), fine_init_coords, fCBTikhonov, fK2);
		else if (fCorrectOverlap == 2)
			fCoarse->CorrectOverlap_2(particle_pair->Neighbors(), fine_init_coords, fCBTikhonov, fK2, r, nip);
		else if (fCorrectOverlap == 22)
			fCoarse->CorrectOverlap_22(particle_pair->Neighbors(), fine_init_coords, fCBTikhonov, fK2, r, nip);
		else if (fCorrectOverlap == 3)
			fCoarse->CorrectOverlap_3(particle_pair->Neighbors(), fine_init_coords, fCBTikhonov, fK2, nip);
		else if (fCorrectOverlap == 4)
			fCoarse->CorrectOverlap_4(particle_pair->Neighbors(), fine_init_coords, fCBTikhonov, fK2, r, nip);
		else
			ExceptionT::GeneralFail(caller, "unrecognized overlap correction method %d", fCorrectOverlap);
	}

//TEMP - debugging
#if 0
if (1) {
	ofstreamT& out = Output();
	int prec = out.precision();
	out.precision(12);
	iArrayT r, c;
	dArrayT v;
	
	/* projection data */
	const PointInCellDataT& projection_data = fCoarse->ProjectionData();
	const InterpolationDataT& interp = projection_data.PointToNode();
	out << "projection:\n" << '\n';
	interp.ToMatrix(r, c, v);
	for (int i = 0; i < r.Length(); i++)
		out << r[i]+1 << " " << c[i]+1 << " " << v[i] << '\n';

	/* interpolation data */
	const PointInCellDataT& interpolation_data = fCoarse->InterpolationData();
	interpolation_data.InterpolationDataToMatrix(r, c, v);
	out << "interpolation:\n" << '\n';
	for (int i = 0; i < r.Length(); i++)
		out << r[i]+1 << " " << c[i]+1 << " " << v[i] << '\n';

	out.flush();
	out.precision(prec);
}
#endif
}

/* (re-)set the equation number for the given group */
void MultiManagerT::SetEquationSystem(int group, int start_eq_shift)
{
	fFine->SetEquationSystem(group, start_eq_shift);
	int neq1 = fFine->GlobalNumEquations(group);
	
	fCoarse->SetEquationSystem(group, neq1 + start_eq_shift);
	int neq2 = fCoarse->GlobalNumEquations(group);

	fGlobalNumEquations[group] = neq1 + neq2;

	/* set total equation numbers */
	fEqnos1.Dimension(neq1);
	fEqnos2.Dimension(neq2);
	fEqnos1.SetValueToPosition();
	fEqnos2.SetValueToPosition();
	fEqnos1 += (1 + start_eq_shift);
	fEqnos2 += (1 + start_eq_shift + neq1);

	/* final step in solver configuration */
	fSolvers[group]->Initialize(
		fGlobalNumEquations[group],
		fGlobalNumEquations[group],
		1);
	
	/* equations for assembly of cross terms */	
	if (fFineToCoarse) {
		const iArray2DT& eq = fCoarseField->Equations();
		fR_U_eqnos.Dimension(eq.Length());
		for (int i = 0; i < fR_U_eqnos.Length(); i++) {
			fR_U_eqnos[i] = eq[i];
			if (fR_U_eqnos[i] > 0)
				fR_U_eqnos[i] += start_eq_shift + neq1;
		}
	}
	if (fCoarseToFine)
		fR_Q_eqnos.Alias(fFineField->Equations());
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
	ExceptionT::CodeT error = ExceptionT::kNoError;
	if (error == ExceptionT::kNoError) 
		error = fFine->ResetStep();
	if (error == ExceptionT::kNoError) 
		error = fCoarse->ResetStep();

	/* loop over solvers only */
	for (fCurrentGroup = 0; fCurrentGroup < NumGroups(); fCurrentGroup++)
		fSolvers[fCurrentGroup]->ResetStep();
	fCurrentGroup = -1;

	/* OK */
	return error;
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
	DiagonalMatrixT* diag_1 = TB_DYNAMIC_CAST(DiagonalMatrixT*, const_cast<GlobalMatrixT*>(&lhs_1));
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
	DiagonalMatrixT* diag_2 = TB_DYNAMIC_CAST(DiagonalMatrixT*, const_cast<GlobalMatrixT*>(&lhs_2));
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
	
	/* skip all cross terms */
	if (!fFineToCoarse && !fCoarseToFine) return;

	/* total internal force vectors */
	int atoms_group = 0;
	const dArray2DT& resid_fine = fFine->InternalForce(group);
	int continuum_group = 0;
	const dArray2DT& resid_coarse = fCoarse->InternalForce(continuum_group);

//TEMP -debugging
#if 0
if (1) {
	const dArray2DT& u_fine = (*fFineField)[0];
	const dArray2DT& u_coarse = (*fCoarseField)[0];

	ofstreamT& out = Output();
	int prec = out.precision();
	out.precision(12);
	int iteration = fSolvers[group]->IterationNumber();
	out << "iteration = " << iteration << '\n';
	out << "u_fine =\n" << u_fine << '\n';
	out << "f_fine =\n" << resid_fine << '\n';

	out << "u_coarse =\n" << u_coarse << '\n';
	out << "f_coarse =\n" << resid_coarse << '\n';
	out.precision(prec);
}
#endif

	/* fine scale contribution to the coarse scale residual */
	dArray2DT& R_U = const_cast<dArray2DT&>(fR_U);
	R_U.Dimension(resid_coarse.MajorDim(), resid_coarse.MinorDim());
	R_U = 0.0;
	if (fFineToCoarse) {
		const iArrayT& ghost_atoms = fFine->GhostNodes();
		const PointInCellDataT& interpolation_data = fCoarse->InterpolationData();
		fCoarse->MultNTf(interpolation_data, resid_fine, ghost_atoms, R_U);
		fSolvers[group]->AssembleRHS(R_U, fR_U_eqnos);
	}

	/* mixed contribution to the fine scale residual */
	if (fCoarseToFine) {
		dArray2DT& R_Q = const_cast<dArray2DT&>(fR_Q);
		R_Q.Dimension(resid_fine.MajorDim(), resid_fine.MinorDim());
		R_Q = 0.0;
		R_U += resid_coarse;
		const PointInCellDataT& projection_data = fCoarse->ProjectionData();
		fCoarse->MultNTf(projection_data.PointToNode(), R_U, fCoarse->ProjectedNodes(), R_Q);	
		fSolvers[group]->AssembleRHS(R_Q, fR_Q_eqnos);	
	}
	
//TEMP - debugging
#if 0
if (1) {
	const dArrayT& rhs = fSolvers[group]->RHS();

	ofstreamT& out = Output();
	int prec = out.precision();
	out.precision(12);
	out << "R =\n" << rhs << '\n';	
	out.precision(prec);
}
#endif
}

/* send update of the solution to the NodeManagerT */
void MultiManagerT::Update(int group, const dArrayT& update)
{
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
	fCoarse->ProjectField(fFineField->Name(), *(fFine->NodeManager()), order);

	/* interpolate coarse scale solution to the fine */
	fCoarse->InterpolateField(fFineField->Name(), order, fFieldAtGhosts);
	fFine->SetFieldValues(fFineField->Name(), fFine->GhostNodes(), order, fFieldAtGhosts);
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
	int ndof = nodes.NumDOF(group);
	const iArrayT& source_points = fFine->NonGhostNodes();
	dArray2DT n_values(source_points.Length(), 2*ndof), coarse;
	fCoarse->CoarseField(fFineField->Name(), nodes, order, coarse);
	if (source_points.Length() != coarse.MajorDim())
		ExceptionT::GeneralFail(caller);

	/* get the total field */
	const dArray2DT& total = (*fFineField)[order];

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

bool MultiManagerT::DecreaseLoadStep(void) {
	return fCoarse->DecreaseLoadStep() && fFine->DecreaseLoadStep();
}

bool MultiManagerT::IncreaseLoadStep(void) {
	return fCoarse->IncreaseLoadStep() && fFine->IncreaseLoadStep();
}


#endif /* BRIDGING_ELEMENT */
