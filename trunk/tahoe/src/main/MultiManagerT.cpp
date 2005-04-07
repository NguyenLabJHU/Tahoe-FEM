/* $Id: MultiManagerT.cpp,v 1.22 2005-04-07 19:18:05 d-farrell2 Exp $ */
#include "MultiManagerT.h"

#ifdef BRIDGING_ELEMENT

#include "ifstreamT.h"
#include "SolverT.h"
#include "DiagonalMatrixT.h"
#include "FEManagerT_bridging.h"
#include "NodeManagerT.h"
#include "OutputSetT.h"
#include "TimeManagerT.h"
#include "ParticlePairT.h"
#include "FieldT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* constructor */
MultiManagerT::MultiManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
	const ArrayT<StringT>& argv, TaskT task):
	FEManagerT(input_file, output, comm, argv, task),
	fFineComm(&fComm),
	fFine(NULL),
	fCoarseComm(&fComm),
	fCoarse(NULL),
	fDivertOutput(false),
	fFineField(NULL),
	fCoarseField(NULL),
	fFineToCoarse(true),
	fCoarseToFine(true),
	fImpExp(IntegratorT::kImplicit)
{
	SetName("tahoe_multi");
}

/* destructor */
MultiManagerT::~MultiManagerT(void)
{
	/* clear pointer because base class will free and this time manager is taken from
	 * the coarse scale solver */
	fTimeManager = NULL;

	/* free subs */
	delete fFine;
	delete fCoarse;
}

/** (re-)set system to initial conditions */
ExceptionT::CodeT MultiManagerT::InitialCondition(void)
{
	ExceptionT::CodeT error = ExceptionT::kNoError;
	if (error == ExceptionT::kNoError) 
		error = fFine->InitialCondition();
	if (fignore == false)	// if we are to use the continuum
	{
		if (error == ExceptionT::kNoError) 
			error = fCoarse->InitialCondition();
	}
	return error;
}

/* (re-)set the equation number for the given group */
void MultiManagerT::Solve(void)
{
	const char caller[] = "FEExecutionManagerT::Solve";
	    
	/* time managers */
	TimeManagerT* atom_time = fFine->TimeManager();
	TimeManagerT* continuum_time = fCoarse->TimeManager();
	if (atom_time->NumberOfSteps() - continuum_time->NumberOfSteps() != 0) 
		ExceptionT::GeneralFail(caller, "coarse/fine number of steps mismatch: %d != %d",
			atom_time->NumberOfSteps(), continuum_time->NumberOfSteps());

	/* set to initial condition */
	ExceptionT::CodeT error = InitialCondition();

	/* loop over time increments */
	while (atom_time->Step() && continuum_time->Step()) /* same clock */
	{
		/* consistency check */
		if (fabs(atom_time->TimeStep() - continuum_time->TimeStep()) > kSmall)
			ExceptionT::GeneralFail(caller, "coarse/fine time step mismatch: %g != %g", 
				atom_time->TimeStep(), continuum_time->TimeStep());

		/* initialize the current time step */
		if (error == ExceptionT::kNoError) 
			error = InitStep();

		/* solve the current time step */
		if (error == ExceptionT::kNoError) 
			error = SolveStep();
			
		/* close the current time step */
		if (error == ExceptionT::kNoError)
			error = CloseStep();
				
		/* handle errors */
		switch (error)
		{
			case ExceptionT::kNoError:
				/* nothing to do */
				break;
			case ExceptionT::kGeneralFail:					
			case ExceptionT::kBadJacobianDet:
			{
				cout << '\n' << caller << ": trying to recover from error: " << ExceptionT::ToString(error) << endl;
				
				/* reset system configuration */
				error = ResetStep();
					
				/* cut time step */
				if (error == ExceptionT::kNoError) {
					if (!DecreaseLoadStep())
						ExceptionT::GeneralFail(caller, "could not decrease load step");
				}
				else
					ExceptionT::GeneralFail(caller, "could not reset step");
				break;
			}
			default: 
				ExceptionT::GeneralFail(caller, "no recovery from \"%s\"", ExceptionT::ToString(error));
		}
	}
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

/* execute the solution procedure */
ExceptionT::CodeT MultiManagerT::SolveStep(void)
{
	/* monolithic solution scheme - inherited */
	if (fSolvers.Length() == 1)
		return FEManagerT::SolveStep();
	else if (fSolvers.Length() == 0) /* staggered solution scheme */
		return SolveStep_Staggered();
	else
		ExceptionT::GeneralFail("MultiManagerT::SolveStep", "number of solvers should be 0 or 1: %d",
			fSolvers.Length());

	return ExceptionT::kGeneralFail;
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
	fCoarse->ProjectField(fFineField->FieldName(), *(fFine->NodeManager()), order);

	/* interpolate coarse scale solution to the fine */
	fCoarse->InterpolateField(fFineField->FieldName(), order, fFieldAtGhosts);
	fFine->SetFieldValues(fFineField->FieldName(), fFine->GhostNodes(), order, fFieldAtGhosts);
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
	fCoarse->CoarseField(fFineField->FieldName(), nodes, order, coarse);
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

/* describe the parameters needed by the interface */
void MultiManagerT::DefineParameters(ParameterListT& list) const
{
	/* inherited - don't call direct base class method */
	ParameterInterfaceT::DefineParameters(list);

	/* paths to input files for sub-tahoe's */
	list.AddParameter(ParameterT::String, "atom_input");
	list.AddParameter(ParameterT::String, "continuum_input",ParameterListT::ZeroOrOnce);
	
	// set boolean to see if continuum input should be ignored
	ParameterT ignore_continuum(ParameterT::Boolean, "ignore_continuum_input");
	ignore_continuum.SetDefault(false);
	list.AddParameter(ignore_continuum);

	/* name of the bridging field */
	list.AddParameter(ParameterT::Word, "bridging_field");

	/* file name for the overlap information */
	ParameterT overlap_file(ParameterT::String, "overlap_file");
	list.AddParameter(overlap_file, ParameterListT::ZeroOrOnce);

	/* force balance options */
	ParameterT fine_to_coarse(fFineToCoarse, "fine_to_coarse");
	fine_to_coarse.SetDefault(fFineToCoarse);
	list.AddParameter(fine_to_coarse);

	ParameterT coarse_to_fine(fCoarseToFine, "coarse_to_fine");
	coarse_to_fine.SetDefault(fCoarseToFine);
	list.AddParameter(coarse_to_fine);
}

/* information about subordinate parameter lists */
void MultiManagerT::DefineSubs(SubListT& sub_list) const
{
	/* inherited - don't call direct base class method */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* ghost atoms */
	sub_list.AddSub("ghost_atom_ID_list", ParameterListT::ZeroOrOnce);

	/* overlap correction method */
	sub_list.AddSub("overlap_correction_choice", ParameterListT::Once, true);
	
	/* single solver for both tahoe's */
	sub_list.AddSub("multi_solver_choice", ParameterListT::ZeroOrOnce, true);
}

/* return the description of the given inline subordinate parameter list */
void MultiManagerT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "multi_solver_choice")
	{
		order = ParameterListT::Choice;

		/* only nonlinear PCG solver for now */
		sub_lists.AddSub("PCG_solver");
	}
	else /* inherited - skip FEManagerT definitions */
		ParameterInterfaceT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MultiManagerT::NewSub(const StringT& name) const
{
	/* try to construct solver */
	FEManagerT* non_const_this = (FEManagerT*) this;	
	SolverT* solver = SolverT::New(*non_const_this, name, -1);
	if (solver)
		return solver;
	else if (name == "overlap_correction_choice")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);

		/* prescribe density in overlap */
		ParameterContainerT prescribe_overlap("prescribe_overlap");
		prescribe_overlap.AddParameter(ParameterT::Double, "bond_density");
		choice->AddSub(prescribe_overlap);

		/* common parameters */
		ParameterT smoothing(ParameterT::Double, "smoothing");
		smoothing.SetDefault(0.0);
		smoothing.AddLimit(0.0, LimitT::LowerInclusive);

		ParameterT bind_1(ParameterT::Double, "bind_to_1.0");
		bind_1.SetDefault(0.0);
		bind_1.AddLimit(0.0, LimitT::LowerInclusive);

		ParameterT init_stiffness(ParameterT::Double, "init_stiffness_factor");
		init_stiffness.SetDefault(1.0e-04);
		init_stiffness.AddLimit(0.0, LimitT::LowerInclusive);

		ParameterT reg(ParameterT::Double, "constraint_regularization");
		reg.SetDefault(0.0);
		reg.AddLimit(0.0, LimitT::LowerInclusive);

		ParameterT density_nip(ParameterT::Enumeration, "density_nip");
		density_nip.AddEnumeration("1", 1);
		density_nip.AddEnumeration("default", -1);
		density_nip.SetDefault(-1);

		/* solve shell-at-a-time (no constraints) */
		ParameterContainerT by_bond("by-bond");
		by_bond.AddParameter(smoothing);
		by_bond.AddParameter(bind_1);
		by_bond.AddParameter(reg);
		by_bond.AddParameter(density_nip);
		choice->AddSub(by_bond);

		/* solve bond-at-a-time with Augmented Lagrangian constraints */
		ParameterContainerT by_bond_L("by-bond_multiplier");
		by_bond_L.AddParameter(smoothing);
		by_bond_L.AddParameter(bind_1);
		by_bond_L.AddParameter(reg);
		by_bond_L.AddParameter(density_nip);
		ParameterT init_bound(ParameterT::Double, "init_bound_width");
		init_bound.SetDefault(1.0);
		init_bound.AddLimit(1.0, LimitT::LowerInclusive);
		by_bond_L.AddParameter(init_bound);
		choice->AddSub(by_bond_L);

		/* solve bond-at-a-time with penalized constraints */
		ParameterContainerT by_bond_p("by-bond_penalty");
		by_bond_p.AddParameter(smoothing);
		by_bond_p.AddParameter(bind_1);
//		by_bond_p.AddParameter(init_stiffness);
		by_bond_p.AddParameter(density_nip);
		ParameterT bound_tolerance(ParameterT::Double, "bound_tolerance");
		bound_tolerance.SetDefault(0.02);
		bound_tolerance.AddLimit(0.0, LimitT::Lower);
		by_bond_p.AddParameter(bound_tolerance);
		ParameterT stiffness_jump(ParameterT::Double, "stiffness_jump");
		stiffness_jump.SetDefault(2.0);
		stiffness_jump.AddLimit(0.0, LimitT::Lower);
		by_bond_p.AddParameter(stiffness_jump);
		choice->AddSub(by_bond_p);
			
		return choice;
	}	
	else /* inherited - skip FEManagerT definitions */
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void MultiManagerT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MultiManagerT::TakeParameterList";

	/* inherited - don't call direct base class method */
	ParameterInterfaceT::TakeParameterList(list);

	/* path to parameters file */
	StringT path;
	path.FilePath(fInputFile);
	TaskT task = kRun;
#pragma message("Ignore continuum only implemented for Dynamic Bridging Scale calculation, ask dave")	
	// check to see if continuum should be ignored
	ParameterT ignore_continuum = list.GetParameter("ignore_continuum_input");
	fignore = false;
	if (ignore_continuum)
		fignore = ignore_continuum;
	
	if (fignore == false)	// if we are to use the continuum
	{
		/* parse/validate continuum input */
		StringT continuum_input = list.GetParameter("continuum_input");
		continuum_input.ToNativePathName();
		continuum_input.Prepend(path);
		ParameterListT continuum_params;
		ParseInput(continuum_input, continuum_params, true, true, true, fArgv);
				
		/* construct continuum solver */
		if (fCoarseComm->Size() != 1)
			ExceptionT::GeneralFail(caller, "parallel execution error");
		if (Size() > 1) /* change file name so output files are unique */
		{
			StringT suffix;
			suffix.Suffix(continuum_input);
			continuum_input.Root();
			continuum_input.Append(".p", Rank());
			continuum_input.Append(suffix);	
		}
		
		StringT continuum_output_file;
		continuum_output_file.Root(continuum_input);
		continuum_output_file.Append(".out");
		fCoarseOut.open(continuum_output_file);
		fCoarse = TB_DYNAMIC_CAST(FEManagerT_bridging*, FEManagerT::New(continuum_params.Name(), continuum_input, fCoarseOut, *fCoarseComm, fArgv, task));
		if (!fCoarse) ExceptionT::GeneralFail(caller, "could not construct continuum solver");
		fCoarse->TakeParameterList(continuum_params);		
	}

	/* parse/validate atomistic input */
	StringT atom_input = list.GetParameter("atom_input");
	atom_input.ToNativePathName();
	atom_input.Prepend(path);
	ParameterListT atom_params;
	ParseInput(atom_input, atom_params, true, true, true, fArgv);

	/* construct atomistic solver */
	if (Size() != fFineComm->Size())
		ExceptionT::GeneralFail(caller, "parallel execution error");
	StringT atom_output_file;
	atom_output_file.Root(atom_input);
	if (Size() > 1) atom_output_file.Append(".p", Rank());
	atom_output_file.Append(".out");
	fFineOut.open(atom_output_file);
	fFine = TB_DYNAMIC_CAST(FEManagerT_bridging*, FEManagerT::New(atom_params.Name(), atom_input, fFineOut, *fFineComm, fArgv, task));
	if (!fFine) ExceptionT::GeneralFail(caller, "could not construct atomistic solver");
	fFine->TakeParameterList(atom_params);
	
	if (fignore == false)	// if we are to use the continuum
	{
		/* check consistency between time managers */
		TimeManagerT* atom_time = fFine->TimeManager();
		TimeManagerT* continuum_time = fCoarse->TimeManager();

		/* use parameters from coarse scale solver */
		fTimeManager = fCoarse->TimeManager();
		fOutputFormat = fCoarse->OutputFormat();

		/* don't compute initial conditions */
		fFine->SetComputeInitialCondition(false);
		fCoarse->SetComputeInitialCondition(false);

		/* resolve bridging fields */
		const StringT& bridging_field = list.GetParameter("bridging_field");
		fFineField = fFine->NodeManager()->Field(bridging_field);
		if (!fFineField) ExceptionT::GeneralFail(caller, "could not resolve fine scale \"%s\" field", bridging_field.Pointer());
		fCoarseField = fCoarse->NodeManager()->Field(bridging_field);
		if (!fFineField) ExceptionT::GeneralFail(caller, "could not resolve coarse scale \"%s\" field", bridging_field.Pointer());
		
		/* resolve integrator types */
		if (fFineField->Integrator().ImplicitExplicit() != fCoarseField->Integrator().ImplicitExplicit())
			ExceptionT::GeneralFail(caller, "time integrator mismatch");
		fImpExp = fFineField->Integrator().ImplicitExplicit();

		/* collect the ghost atom ID list */
		ArrayT<StringT> ghost_atom_ID;
		const ParameterListT* ghosts = list.List("ghost_atom_ID_list");
		if (ghosts)	StringListT::Extract(*ghosts, ghost_atom_ID);

		/* configure projection/interpolation */
		NodeManagerT& fine_node_manager = *(fFine->NodeManager());	
		int group = 0;
		int order1 = 0;
		bool make_inactive = true;
		fFine->InitGhostNodes(fFineField->FieldName(), ghost_atom_ID, fCoarse->ProjectImagePoints());
		fCoarse->InitInterpolation(fFineField->FieldName(), fFine->GhostNodes(), fine_node_manager.InitialCoordinates());
		fCoarse->InitProjection(fFineField->FieldName(), *(fFine->CommManager()), fFine->NonGhostNodes(), fine_node_manager, make_inactive);
	
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

		/* construct solver */
		int n1 = fFine->NumGroups();
		int n2 = fCoarse->NumGroups();
		if (n1 != n2) ExceptionT::GeneralFail(caller, "number of groups must match: %d != %d", n1, n2);
		const ParameterListT* multi_solver = list.ListChoice(*this, "multi_solver_choice");
		if (multi_solver)
		{
			/* construct */
			SolverT* solver = SolverT::New(*this, multi_solver->Name(), 0);
			if (!solver) ExceptionT::GeneralFail(caller, "could not construct solver \"%s\"",
				multi_solver->Name().Pointer());
			solver->TakeParameterList(*multi_solver);
			
			/* store */
			fSolvers.Dimension(1);
			fSolvers[0] = solver;
		}
		SetSolver();

		/* default solver phases */
		fMaxSolverLoops = 1;
		fSolverPhases.Dimension(1, 3);
		fSolverPhases(0,0) = 0;
		fSolverPhases(0,1) =-1;
		fSolverPhases(0,2) =-1;
		fSolverPhasesStatus.Dimension(fSolverPhases.MajorDim(), kNumStatusFlags);
		fSolverPhasesStatus = 0;

		/* terms to include in the equilibrium equations */
		fFineToCoarse = list.GetParameter("fine_to_coarse");
		fCoarseToFine = list.GetParameter("coarse_to_fine");
		if (fCoarseToFine) /* enforce zero bond density in projected cells */
			fCoarse->DeactivateFollowerCells();
		
		/* needed to solve overlap */
		const dArray2DT& fine_init_coords = fine_node_manager.InitialCoordinates();
		const ParticlePairT* particle_pair = fFine->ParticlePair();
		if (!particle_pair) ExceptionT::GeneralFail(caller, "could not resolve ParticlePairT");
		
		/* overlap correction method */
		StringT overlap_path;
		const ParameterT* overlap_file = list.Parameter("overlap_file");
		if (overlap_file) {
			overlap_path = *overlap_file;
			overlap_path.ToNativePathName();
			
			StringT path;
			path.FilePath(fInputFile);
			overlap_path.Prepend(path);
		} /* default name */
		else {
			overlap_path.Root(fInputFile);
			overlap_path.Append(".overlap");
		}
		
		const ParameterListT& overlap = list.GetListChoice(*this, "overlap_correction_choice");
		if (overlap.Name() == "prescribe_overlap")
		{
			double density = overlap.GetParameter("bond_density");
			cout << "\n " << caller << ": \"prescribe_overlap\" not implemented. p = 1.0" << endl;
		}
		else if (overlap.Name() == "by-bond_multiplier")
		{
			/* extract parameters */
			double smoothing = overlap.GetParameter("smoothing");
			double bind_1 = overlap.GetParameter("bind_to_1.0");
			double reg = overlap.GetParameter("constraint_regularization");
			int nip = overlap.GetParameter("density_nip");
			double init_bound_width = overlap.GetParameter("init_bound_width");
			
			/* compute overlap correction */
			double bound_0 = init_bound_width/2.0;
			fCoarse->CorrectOverlap_2(particle_pair->Neighbors(), fine_init_coords, 
				overlap_path, smoothing, bind_1, reg, bound_0, nip);
		}
		else if (overlap.Name() == "by-bond_penalty")
		{
			/* extract parameters */
			double smoothing = overlap.GetParameter("smoothing");
			double bind_1 = overlap.GetParameter("bind_to_1.0");
			double bound_tol = overlap.GetParameter("bound_tolerance");
			double stiffness_jump = overlap.GetParameter("stiffness_jump");
			int nip = overlap.GetParameter("density_nip");
			
			/* compute overlap correction */
			fCoarse->CorrectOverlap_22(particle_pair->Neighbors(), fine_init_coords, 
				overlap_path, smoothing, bind_1, bound_tol, stiffness_jump, nip);
		}
		else
			ExceptionT::GeneralFail(caller, "unrecognized overlap correction method \"%s\"",
				overlap.Name().Pointer());
	}
	else if(fignore == true)	// if we are to ignore the continuum
	{
#pragma message("Ignore continuum only implemented for Dynamic Bridging Scale calculation, ask dave")
#pragma message("Can get rid of the redundancy here, not quite sure how to do it well")
		/* check consistency between time managers */
		TimeManagerT* atom_time = fFine->TimeManager();

		/* use parameters from fine scale solver */
		fTimeManager = fFine->TimeManager();
		fOutputFormat = fFine->OutputFormat();

		/* don't compute initial conditions */
		fFine->SetComputeInitialCondition(false);

		/* resolve bridging fields */
		const StringT& bridging_field = list.GetParameter("bridging_field");
		fFineField = fFine->NodeManager()->Field(bridging_field);
		if (!fFineField) ExceptionT::GeneralFail(caller, "could not resolve fine scale \"%s\" field", bridging_field.Pointer());

		/* resolve integrator types */
		fImpExp = fFineField->Integrator().ImplicitExplicit();

		/* collect the ghost atom ID list */
		ArrayT<StringT> ghost_atom_ID;
		const ParameterListT* ghosts = list.List("ghost_atom_ID_list");
		if (ghosts)	StringListT::Extract(*ghosts, ghost_atom_ID);

		/* configure projection/interpolation */
		NodeManagerT& fine_node_manager = *(fFine->NodeManager());	
		int group = 0;
		int order1 = 0;
		bool make_inactive = true;
		fFine->InitGhostNodes(fFineField->FieldName(), ghost_atom_ID, false);

		/* send fine output through the fFine output */
		int ndof = fFine->NodeManager()->NumDOF(group);
		ArrayT<StringT> labels(ndof);
		const char* fine_labels[] = {"UF_X", "UF_Y", "UF_Z"};
		int dex = 0;
		for (int i = 0; i < ndof; i++) labels[dex++] = fine_labels[i];
		const iArrayT& non_ghost_nodes = fFine->NonGhostNodes();
		fAtomConnectivities.Alias(non_ghost_nodes.Length(), 1, non_ghost_nodes.Pointer());
		OutputSetT output_set(GeometryT::kPoint, fAtomConnectivities, labels, false);
		fOutputID = fFine->RegisterOutput(output_set);

		/* construct solver */
		int n1 = fFine->NumGroups();
		const ParameterListT* multi_solver = list.ListChoice(*this, "multi_solver_choice");
		if (multi_solver)
		{
			/* construct */
			SolverT* solver = SolverT::New(*this, multi_solver->Name(), 0);
			if (!solver) ExceptionT::GeneralFail(caller, "could not construct solver \"%s\"",
				multi_solver->Name().Pointer());
			solver->TakeParameterList(*multi_solver);
			
			/* store */
			fSolvers.Dimension(1);
			fSolvers[0] = solver;
		}
		SetSolver();

		/* default solver phases */
		fMaxSolverLoops = 1;
		fSolverPhases.Dimension(1, 3);
		fSolverPhases(0,0) = 0;
		fSolverPhases(0,1) =-1;
		fSolverPhases(0,2) =-1;
		fSolverPhasesStatus.Dimension(fSolverPhases.MajorDim(), kNumStatusFlags);
		fSolverPhasesStatus = 0;
	}
	else	// something is wrong
	{
		ExceptionT::GeneralFail(caller, "Option to ignore continuum not set properly, check input file");
	}
}

/* driver for staggered solution with single clock for both systems */
ExceptionT::CodeT MultiManagerT::SolveStep_Staggered(void)
{
	const char caller[] = "MultiManagerT::SolveStep_Staggered";

	/* solver phase status */
	const iArray2DT& atom_phase_status = fFine->SolverPhasesStatus();
	const iArray2DT& continuum_phase_status = fCoarse->SolverPhasesStatus();

	/* running error */
	ExceptionT::CodeT error = ExceptionT::kNoError;

	/* loop until both solved */
	dArray2DT field_at_ghosts;
	int d_width = OutputWidth(cout, field_at_ghosts.Pointer());	
	int group_num = 0;
	int order1 = 0; 
	double atoms_res, continuum_res, combined_res_0 = 0.0;
	int count = 0;
	int atom_last_iter, atom_iter, continuum_last_iter, continuum_iter;
	atom_last_iter = atom_iter = continuum_last_iter = continuum_iter = 0;
	while (count == 0 || (atom_iter > 0 || continuum_iter > 0)) //TEMP - assume just one phase
//	while (1 || error == ExceptionT::kNoError &&
//		(atom_phase_status(0, FEManagerT::kIteration) > 0 ||
//		continuum_phase_status(0, FEManagerT::kIteration) > 0)) //TEMP - assume just one phase
	{
		count++;

		/* solve atoms */
		if (1 || error == ExceptionT::kNoError) {
			fFine->ResetCumulativeUpdate(group_num);
			error = fFine->SolveStep();
		}

		/* apply solution to continuum */
		fCoarse->ProjectField(fFineField->FieldName(), *fFine->NodeManager(), order1);
		fCoarse->FormRHS(group_num);
		continuum_res = fCoarse->RHS(group_num).Magnitude(); //serial
					
		/* solve continuum */
		if (1 || error == ExceptionT::kNoError) {
			fCoarse->ResetCumulativeUpdate(group_num);
			error = fCoarse->SolveStep();
		}
				
		/* apply solution to atoms */
		int order = 0;  // displacement only for static case
		fCoarse->InterpolateField(fFineField->FieldName(), order, field_at_ghosts);
		fFine->SetFieldValues(fFineField->FieldName(), fFine->GhostNodes(), order, field_at_ghosts);
		fFine->FormRHS(group_num);
		atoms_res = fFine->RHS(group_num).Magnitude(); //serial

		/* reset the reference errors */
		if (count == 1) {
			combined_res_0 = atoms_res + continuum_res;
			fFine->SetReferenceError(group_num, combined_res_0);
			fCoarse->SetReferenceError(group_num, combined_res_0);
		}
					
		/* log residual */
		double tot_rel_error = (fabs(combined_res_0) > kSmall) ? 
			(atoms_res + continuum_res)/combined_res_0 : 0.0;
		cout << setw(kIntWidth) << count << ": "
			 << setw(d_width) << atoms_res << " (A) | "
			 << setw(d_width) << continuum_res << " (C) | "
			 << setw(d_width) << tot_rel_error << endl;

		/* number of interations in last pass */
		int atom_total_iter = atom_phase_status(0, FEManagerT::kIteration);
		int continuum_total_iter = continuum_phase_status(0, FEManagerT::kIteration);
		atom_iter = atom_total_iter - atom_last_iter;
		continuum_iter = continuum_total_iter - continuum_last_iter;
		atom_last_iter = atom_total_iter;
		continuum_last_iter = continuum_total_iter;
	}

	/* do not return error codes */
	return ExceptionT::kNoError;
}

#endif /* BRIDGING_ELEMENT */
