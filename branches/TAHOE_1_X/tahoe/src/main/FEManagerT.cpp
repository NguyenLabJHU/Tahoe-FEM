/* $Id: FEManagerT.cpp,v 1.74 2004-06-26 18:38:08 paklein Exp $ */
/* created: paklein (05/22/1996) */
#include "FEManagerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <string.h>
#include <float.h>
#include <ctype.h>

#include "toolboxConstants.h"
#include "ExceptionT.h"

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "TimeManagerT.h"
#include "ModelManagerT.h"
#include "ElementBaseT.h"
#include "IOManager.h"
#include "OutputSetT.h"
#include "CommandSpecT.h"
#include "ArgSpecT.h"
#include "eIntegratorT.h"
#include "nIntegratorT.h"
#include "CommunicatorT.h"
#include "CommManagerT.h"

/* nodes */
#include "NodeManagerT.h"

/* solvers */
#include "LinearSolver.h"
#include "NLSolver.h"
#include "DRSolver.h"
#include "NLK0Solver.h"
#include "ExpCD_DRSolver.h"
#include "NLSolverX.h"
#include "NLSolver_LS.h"
#include "PCGSolver_LS.h"
#include "iNLSolver_LS.h"
#include "NOXSolverT.h"
#include "NLSolver_LSX.h" //TEMP

using namespace Tahoe;

/* File/Version Control */
const char kCurrentVersion[] = "v3.4.1";
const char kProgramName[] = "tahoe";

/* static methods */
const char* FEManagerT::Version(void) { return kCurrentVersion; }

/* constructor */
FEManagerT::FEManagerT(ifstreamT& input, ofstreamT& output, CommunicatorT& comm, const ArrayT<StringT>& argv):
	ParameterInterfaceT("tahoe"),
	fArgv(argv),
	fMainIn(input),
	fMainOut(output),
	fComm(comm),
	fPrintInput(false),
	fComputeInitialCondition(true),
	fStatus(GlobalT::kConstruction),
	fTimeManager(NULL),
	fNodeManager(NULL),
	fElementGroups(NULL),
	fIOManager(NULL),
	fModelManager(NULL),
	fCommManager(NULL),
	fMaxSolverLoops(0),
	fGlobalEquationStart(0),
	fActiveEquationStart(0),
	fGlobalNumEquations(0),
	fCurrentGroup(-1),
	fAnalysisCode(GlobalT::kNoAnalysis)
{
	/* console name */
	iSetName("FE_manager");

	/* add console variables */
	iAddVariable("title", *((const StringT*) &fTitle));
	iAddVariable("restart_inc", fWriteRestart);
	
	/* console commands */
	ArgSpecT rs_file(ArgSpecT::string_);
	rs_file.SetPrompt("restart file");

	CommandSpecT read_rs("ReadRestart");
	read_rs.AddArgument(rs_file);
	iAddCommand(read_rs);

	CommandSpecT write_rs("WriteRestart");
	write_rs.AddArgument(rs_file);
	iAddCommand(write_rs);
	
	iAddCommand(CommandSpecT("WriteOutput"));
}

/* destructor */
FEManagerT::~FEManagerT(void)
{
	fStatus = GlobalT::kDestruction;

	delete fTimeManager;
	delete fNodeManager;
	delete fElementGroups;
	
	for (int i = 0; i < fSolvers.Length(); i++)
		delete fSolvers[i];
	
	for (int i = 0; i < fIntegrators.Length(); i++)
		delete fIntegrators[i];

	delete fIOManager;
	delete fModelManager;
	delete fCommManager;
	
	fStatus = GlobalT::kNone;	
}

/* initialize members */
void FEManagerT::Initialize(InitCodeT init)
{
	bool verbose = false;
#ifdef __CPLANT__
	verbose = true;
#endif

	/* state */
	fStatus = GlobalT::kInitialization;

	/* verify file */
	CheckInputFile();
	
	/* get title */
	fMainIn.next_char();
	fTitle.GetLineFromStream(fMainIn);
	cout << "\n Title: " << fTitle << endl;
	
	/* set model manager */
	fModelManager = new ModelManagerT(fMainOut);
	if (!fModelManager) throw ExceptionT::kOutOfMemory;
	if (verbose) cout << "    FEManagerT::Initialize: input" << endl;

	/* main parameters */
	ReadParameters(init);
	if (init == kParametersOnly) return;
	WriteParameters();
	if (verbose) cout << "    FEManagerT::Initialize: execution parameters" << endl;

	/* construct IO manager - configure in SetOutput below */
	StringT file_name(fMainIn.filename());
	fIOManager = new IOManager(fMainOut, kProgramName, kCurrentVersion, fTitle,
		file_name, fOutputFormat);	
	if (!fIOManager) throw ExceptionT::kOutOfMemory;

	/* set communication manager */
	fCommManager = New_CommManager();
	if (!fCommManager) throw ExceptionT::kOutOfMemory;
	if (verbose) cout << "    FEManagerT::Initialize: comm manager" << endl;
	
	/* construct the managers */
	fTimeManager = new TimeManagerT(*this);
	fTimeManager->Initialize();
	if (!fTimeManager) throw ExceptionT::kOutOfMemory;
	iAddSub(*fTimeManager);
	if (verbose) cout << "    FEManagerT::Initialize: time" << endl;

	/* set time integration integrator */
	SetIntegrator();
	if (verbose) cout << "    FEManagerT::Initialize: integrator" << endl;

	/* initial configuration of communication manager */
	if (fModelManager->DatabaseFormat() != IOBaseT::kTahoe)
		fCommManager->Configure();

	/* set fields */
	SetNodeManager();
	fCommManager->SetNodeManager(fNodeManager);
	if (verbose) cout << "    FEManagerT::Initialize: nodal data" << endl;
	
	/* construct element groups */
	SetElementGroups();
	if (verbose) cout << "    FEManagerT::Initialize: element groups" << endl;

	/* initial configuration of communication manager */
	if (0 && fModelManager->DatabaseFormat() == IOBaseT::kTahoe)
	{
		//TEMP
		ExceptionT::Stop("FEManagerT::Initialize", "kTahoe format not supported");
	
		fCommManager->Configure();
		//call to tell everything to reconfigure
	}

	/* set output manager */
	SetOutput();
	if (verbose) cout << "    FEManagerT::Initialize: io" << endl;
	if (init == kAllButSolver) return;

	/* set solution driver */
	SetSolver();
	if (verbose) cout << "    FEManagerT::Initialize: solver" << endl;
}	

/* solve all the time sequences */
void FEManagerT::Solve(void)
{
	const char caller[] = "FEManagerT::Solve";

	fTimeManager->Top();
	while (fTimeManager->NextSequence())
	{	
		/* set to initial condition */
		ExceptionT::CodeT error = InitialCondition();

		/* loop over time increments */
		while (error == ExceptionT::kNoError && fTimeManager->Step())
		{
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
					if (error == ExceptionT::kNoError)
						if (!DecreaseLoadStep())
							error = ExceptionT::kGeneralFail;

					break;
				}
				default:
					cout << '\n' << caller <<  ": no recovery for error: " << ExceptionT::ToString(error) << endl;
			}
		}
	}
}

/* manager messaging */
const ScheduleT* FEManagerT::Schedule(int num) const
{
	return fTimeManager->Schedule(num);
}

GlobalT::AnalysisCodeT FEManagerT::Analysis(void) const { return fAnalysisCode; }
bool FEManagerT::PrintInput(void) const { return fPrintInput; }

void FEManagerT::WriteEquationNumbers(int group) const
{
	fNodeManager->WriteEquationNumbers(group, fMainOut);
	fMainOut.flush();
}

GlobalT::SystemTypeT FEManagerT::GlobalSystemType(int group) const
{
	GlobalT::SystemTypeT type = fNodeManager->TangentType(group);
	for (int i = 0 ; i < fElementGroups->Length(); i++)
	{
		GlobalT::SystemTypeT e_type = (*fElementGroups)[i]->TangentType();

		/* using precedence */
		type = (e_type > type) ? e_type : type;
	}
	
	return type;
}

/* load control functions */
bool FEManagerT::DecreaseLoadStep(void) { return fTimeManager->DecreaseLoadStep(); }
bool FEManagerT::IncreaseLoadStep(void) { return fTimeManager->IncreaseLoadStep(); }

ExceptionT::CodeT FEManagerT::ResetStep(void)
{
	ExceptionT::CodeT error = ExceptionT::kNoError;
	try{
		/* state */
		fStatus = GlobalT::kResetStep;	

		/* time */
		fTimeManager->ResetStep();
	
		/* check group flag */
		if (fCurrentGroup != -1) throw;

		/* do group-by-group */
		for (fCurrentGroup = 0; fCurrentGroup < NumGroups(); fCurrentGroup++)
		{
			/* relaxation code */
			GlobalT::RelaxCodeT relax = GlobalT::kNoRelax;
	
			/* node manager */
			relax = GlobalT::MaxPrecedence(relax, fNodeManager->ResetStep(fCurrentGroup));
	
			/* elements */
			for (int i = 0 ; i < fElementGroups->Length(); i++)
				if ((*fElementGroups)[i]->InGroup(fCurrentGroup))
					relax = GlobalT::MaxPrecedence(relax, (*fElementGroups)[i]->ResetStep());
	
			/* solver */
			fSolvers[fCurrentGroup]->ResetStep();
		
			/* check to see if the equation system needs to be reset */
			if (relax == GlobalT::kReEQ)
				SetEquationSystem(fCurrentGroup);
			else if (relax != GlobalT::kNoRelax)
				ExceptionT::GeneralFail("FEManagerT::ResetStep", "unsupported relaxation code %d", relax);
		}
	}
	
	catch (ExceptionT::CodeT exc) {
		cout << "\n FEManagerT::ResetStep: caught exception: " 
		     << ExceptionT::ToString(exc) << endl;
		return exc;
	}
	
	/* reset group flag */
	fCurrentGroup = -1;
	
	/* OK */
	return ExceptionT::kNoError;
}

const double& FEManagerT::Time(void) const { return fTimeManager->Time(); }
const double& FEManagerT::TimeStep(void) const { return fTimeManager->TimeStep(); }
const int& FEManagerT::StepNumber(void) const { return fTimeManager->StepNumber() ; }
const int& FEManagerT::NumberOfSteps(void) const { return fTimeManager->NumberOfSteps(); }
int FEManagerT::SequenceNumber(void) const { return fTimeManager->SequenceNumber(); }
int FEManagerT::NumSequences(void) const { return fTimeManager->NumSequences(); }
const int& FEManagerT::IterationNumber(int group) const 
{
#if __option(extended_errorcheck)
	/* range check */
	if (group < 0 || group >= fSolvers.Length()) 
		ExceptionT::OutOfRange("FEManagerT::IterationNumber", "group is out of range: %d", group);
#endif
	return fSolvers[group]->IterationNumber(); 
}

/* solution messaging */
void FEManagerT::FormLHS(int group, GlobalT::SystemTypeT sys_type) const
{
	/* state */
	SetStatus(GlobalT::kFormLHS);
	
	/* nodal contributions - from special BC's */
	fNodeManager->FormLHS(group, sys_type);

	/* element contributions */
	for (int i = 0 ; i < fElementGroups->Length(); i++)
		if ((*fElementGroups)[i]->InGroup(group))
			(*fElementGroups)[i]->FormLHS(sys_type);
}

void FEManagerT::FormRHS(int group) const
{
	/* state */
	SetStatus(GlobalT::kFormRHS);

	/* unlock assembly into RHS */
	fSolvers[group]->UnlockRHS();

	/* nodal force contribution - F(t) */
	fNodeManager->FormRHS(group);

	/* element contribution */
	for (int i = 0 ; i < fElementGroups->Length(); i++)
		if ((*fElementGroups)[i]->InGroup(group))
			(*fElementGroups)[i]->FormRHS();
		
	/* output system info (debugging) */
	if (fSolvers[group]->Check() == GlobalMatrixT::kPrintRHS)
		WriteSystemConfig(fMainOut, group);

	/* lock assembly into RHS */
	fSolvers[group]->LockRHS();
}

/* the residual for the given group */
const dArrayT& FEManagerT::RHS(int group) const 
{
	return fSolvers[group]->RHS(); 
}

/* the residual for the given group */
const GlobalMatrixT& FEManagerT::LHS(int group) const 
{
	return fSolvers[group]->LHS(); 
}

/* collect the internal force on the specified node */
void FEManagerT::InternalForceOnNode(const FieldT& field, int node, dArrayT& force) const
{
	/* initialize */
	force = 0.0;

	/* element contribution */
	for (int i = 0 ; i < fElementGroups->Length(); i++)
		(*fElementGroups)[i]->AddNodalForce(field, node, force);
}

ExceptionT::CodeT FEManagerT::InitStep(void)
{
	try {
	/* state */
	SetStatus(GlobalT::kInitStep);
	
	/* set the default value for the output time stamp */
	fIOManager->SetOutputTime(Time());	

	/* loop over solver groups */
	if (fCurrentGroup != -1) throw ExceptionT::kGeneralFail;
	for (fCurrentGroup = 0; fCurrentGroup < NumGroups(); fCurrentGroup++)
	{
		/* solver */
		fSolvers[fCurrentGroup]->InitStep();

		/* nodes */
		fNodeManager->InitStep(fCurrentGroup);
	}
	fCurrentGroup = -1;

	/* elements */
	for (int i = 0 ; i < fElementGroups->Length(); i++)
		(*fElementGroups)[i]->InitStep();
	}
	
	catch (ExceptionT::CodeT exc) {
		cout << "\n FEManagerT::InitStep: caught exception: " 
		     << ExceptionT::ToString(exc) << endl;
		return exc;
	}
	
	/* OK */
	return ExceptionT::kNoError;
}

ExceptionT::CodeT FEManagerT::SolveStep(void)
{
	ExceptionT::CodeT error = ExceptionT::kNoError;
	try {
	
		/* check group flag */
		if (fCurrentGroup != -1) throw ExceptionT::kGeneralFail;
		
		int loop_count = 0;
		bool all_pass = false;
		SolverT::SolutionStatusT status = SolverT::kContinue;

		while (!all_pass && 
			(fMaxSolverLoops == -1 || loop_count < fMaxSolverLoops) &&
			status != SolverT::kFailed)
		{
			/* clear status */
			fSolverPhasesStatus = 0;
		 
			/* one solver after the next */
			all_pass = true;
			for (int i = 0; status != SolverT::kFailed && i < fSolverPhases.MajorDim(); i++)
			{
				/* group parameters */
				fCurrentGroup = fSolverPhases(i,0);
				int iter = fSolverPhases(i,1);
				int pass = fSolverPhases(i,2);
			
				/* call solver */
				status = fSolvers[fCurrentGroup]->Solve(iter);
				
				/* check result */
				fSolverPhasesStatus(i, kGroup) = fCurrentGroup;
				fSolverPhasesStatus(i, kIteration) = fSolvers[fCurrentGroup]->IterationNumber();
				if (status == SolverT::kFailed) {
					all_pass = false;
					fSolverPhasesStatus(i, kPass) = -1;					
				}
				else if (status == SolverT::kConverged && 
					(pass == -1 || fSolverPhasesStatus(i, kPass) <= pass))
				{
					all_pass = all_pass && true; /* must all be true */
					fSolverPhasesStatus(i, kPass) = 1;
				}
				else
				{
					all_pass = false;
					fSolverPhasesStatus(i, kPass) = 0;
				}
			}

			/* count */
			loop_count++;
			
			/* write status */
			if (fSolverPhases.MajorDim() > 1)
			{
				cout << "\n Solver status: pass " << loop_count << '\n';
				cout << setw(kIntWidth) << "#"
				     << setw(kIntWidth) << "solver"
				     << setw(kIntWidth) << "its."
				     << setw(kIntWidth) << "pass" << '\n';
				fSolverPhasesStatus.WriteNumbered(cout);
				cout << endl;
			}			
		}

		/* solver looping not successful */
		if (!all_pass) {
			cout << "\n Solver loop not converged after " << loop_count << " iterations" << endl;
			error = ExceptionT::kGeneralFail;
		} else if (fSolverPhases.MajorDim() > 1) {
			cout << "\n Solver loop converged in " << loop_count << " iterations" << endl;
		}
	}	
	catch (ExceptionT::CodeT exc) {
		cout << "\n FEManagerT::SolveStep: caught exception: " 
		     << ExceptionT::ToString(exc) << endl;
		error = exc;
	}
	
	/* reset group flag */
	fCurrentGroup = -1;
	
	/* done */
	return error;
}

ExceptionT::CodeT FEManagerT::CloseStep(void)
{
	try {
	/* state */
	SetStatus(GlobalT::kCloseStep);

	/* solvers - need to be called first because they may be diverted the
	 * output to a temp file and need to switch it back */
	if (fCurrentGroup != -1) throw ExceptionT::kGeneralFail;
	for (fCurrentGroup = 0; fCurrentGroup < NumGroups(); fCurrentGroup++)
		fSolvers[fCurrentGroup]->CloseStep();
	fCurrentGroup = -1;

	/* write output BEFORE closing nodes and elements */
	if (fTimeManager->WriteOutput())
		WriteOutput(Time());

	/* nodes - loop over all groups */
	if (fCurrentGroup != -1) throw ExceptionT::kGeneralFail;
	for (fCurrentGroup = 0; fCurrentGroup < NumGroups(); fCurrentGroup++)
		fNodeManager->CloseStep(fCurrentGroup);
	fCurrentGroup = -1;

	/* elements */
	for (int i = 0 ; i < fElementGroups->Length(); i++)
		(*fElementGroups)[i]->CloseStep();

	/* write restart file */
	WriteRestart();
	}
	catch (ExceptionT::CodeT exc) {
		cout << "\n FEManagerT::CloseStep: caught exception: " << ExceptionT::ToString(exc) << endl;
		return exc;
	}
	
	/* OK */
	return ExceptionT::kNoError;
}

void FEManagerT::Update(int group, const dArrayT& update)
{
	fNodeManager->Update(group, update);
}

void FEManagerT::GetUnknowns(int group, int order, dArrayT& unknowns) const
{
	fNodeManager->GetUnknowns(group, order, unknowns);
}

GlobalT::RelaxCodeT FEManagerT::RelaxSystem(int group) const
{
	/* state */
	SetStatus(GlobalT::kRelaxSystem);
	
	/* check node manager */
	GlobalT::RelaxCodeT relax = GlobalT::kNoRelax;
	relax = GlobalT::MaxPrecedence(relax, fNodeManager->RelaxSystem(group));
		
	/* check element groups - must touch all of them to reset */
	for (int i = 0 ; i < fElementGroups->Length(); i++)
		if ((*fElementGroups)[i]->InGroup(group))
			relax = GlobalT::MaxPrecedence(relax, (*fElementGroups)[i]->RelaxSystem());

	return relax;
}

/* global equation functions */
void FEManagerT::AssembleLHS(int group, const ElementMatrixT& elMat,
	const nArrayT<int>& eqnos) const
{
	fSolvers[group]->AssembleLHS(elMat, eqnos);
}

void FEManagerT::AssembleLHS(int group, const ElementMatrixT& elMat,
	const nArrayT<int>& row_eqnos, const nArrayT<int>& col_eqnos) const
{
	fSolvers[group]->AssembleLHS(elMat, row_eqnos, col_eqnos);
}

void FEManagerT::AssembleLHS(int group, const nArrayT<double>& diagonal_elMat, 
	const nArrayT<int>& eqnos) const
{
	fSolvers[group]->AssembleLHS(diagonal_elMat, eqnos);
}

void FEManagerT::OverWriteLHS(int group, const ElementMatrixT& elMat,
	const nArrayT<int>& eqnos) const
{
	fSolvers[group]->OverWriteLHS(elMat, eqnos);
}

void FEManagerT::DisassembleLHS(int group, dMatrixT& elMat, const nArrayT<int>& eqnos) const
{
	fSolvers[group]->DisassembleLHS(elMat, eqnos);
}

void FEManagerT::DisassembleLHSDiagonal(int group, dArrayT& diagonals, const nArrayT<int>& eqnos) const
{
	fSolvers[group]->DisassembleLHSDiagonal(diagonals, eqnos);
}

void FEManagerT::AssembleRHS(int group, const nArrayT<double>& elRes,
	const nArrayT<int>& eqnos) const
{
	fSolvers[group]->AssembleRHS(elRes, eqnos);
}

void FEManagerT::OverWriteRHS(int group, const dArrayT& elRes, const nArrayT<int>& eqnos) const
{
	fSolvers[group]->OverWriteRHS(elRes, eqnos);
}

void FEManagerT::DisassembleRHS(int group, dArrayT& elRes, const nArrayT<int>& eqnos) const
{
	fSolvers[group]->DisassembleRHS(elRes, eqnos);
}

/* writing results */
void FEManagerT::WriteOutput(double time)
{
	try
	{
		/* state */
		SetStatus(GlobalT::kWriteOutput);

		/* set output time */
		fIOManager->SetOutputTime(time);

		/* write marker */
		ofstreamT& out = Output();
		out << "\n Time = " << time << '\n';
		out << " Step " << fTimeManager->StepNumber() << " of " << fTimeManager->NumberOfSteps() << '\n';

		/* nodes */
		fNodeManager->WriteOutput();

		/* elements */
		for (int i = 0; i < fElementGroups->Length(); i++)
			(*fElementGroups)[i]->WriteOutput();

		/* write info */
		if (fPrintInput) {
			const ArrayT<int>* pbc_nodes = fCommManager->NodesWithGhosts();
			if (pbc_nodes) {
				iArrayT tmp, uni;
				tmp.Alias(*pbc_nodes);
				uni.Union(tmp);
				uni++;
				out << " number of nodes with ghosts = " << pbc_nodes->Length() << '\n';
				out << uni.wrap(5) << endl;
			}
		}
		
		/* system output */
		for (int group = 0; group < fSO_OutputID.Length(); group++)
			if (fSO_DivertOutput[group] && fSO_OutputID[group] != -1) 
			{
				/* get the output set */
				int ndof = fNodeManager->NumDOF(group);
				int nnd  = fSO_Connects.MajorDim();
				dArray2DT n_values(nnd, 2*ndof);
				n_values = 0.0;
		
				/* collect field values */
				ArrayT<FieldT*> fields;
				fNodeManager->CollectFields(group, fields);
				int column = 0;
				for (int i = 0; i < fields.Length(); i++)
				{
					/* displacement */
					const dArray2DT& disp = (*fields[i])[0];

					/* collect field values for the output nodes */
					dArray2DT disp_out;
					if (disp.MajorDim() == nnd)
						disp_out.Alias(disp);
					else {
						disp_out.Dimension(nnd, disp.MinorDim());
						disp_out.RowCollect(fSO_Connects, disp);
					}
				
					for (int j = 0; j < disp.MinorDim(); j++)
						n_values.ColumnCopy(column++, disp_out, j);
				}
			
				/* collect forces */
				column = ndof;
				int shift = ActiveEquationStart(group);
				const dArrayT& RHS = fSolvers[group]->RHS();		
				int num_eq = RHS.Length();
				for (int i = 0; i < fields.Length(); i++)
				{
					/* equations numbers */
					const iArray2DT& eqnos = fields[i]->Equations();
					int ndof_i = eqnos.MinorDim();

					/* loop over nodes */
					for (int j = 0; j < nnd; j++)
					{
						double* f = n_values(j) + column;
						int node = fSO_Connects[j];
						const int* eqno = eqnos(node);
						for (int k = 0; k < ndof_i; k++)
						{
							int eq_k = eqno[k] - shift;
							if (eq_k > -1 && eq_k < num_eq) /* active */
								f[k] = RHS[eq_k];
						}
					}

					/* next field */
					column += ndof_i;
				}
	
				/* write output */
				dArray2DT e_values;
				WriteOutput(fSO_OutputID[group], n_values, e_values);
			}
	}
	
	catch (ExceptionT::CodeT error) { ExceptionT::Throw(error, "FEManagerT::WriteOutput"); }
}

void FEManagerT::WriteOutput(int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
	fIOManager->WriteOutput(ID, n_values, e_values);
}

/* write a snapshot */
void FEManagerT::WriteOutput(const StringT& file, const dArray2DT& coords, const iArrayT& node_map,
	const dArray2DT& values, const ArrayT<StringT>& labels) const
{
	fIOManager->WriteOutput(file, coords, node_map, values, labels);
}

int FEManagerT::RegisterOutput(const OutputSetT& output_set) const
{
	/* check */
	if (!fIOManager) 
		ExceptionT::GeneralFail("FEManagerT::RegisterOutput", "I/O manager not initialized");

	/* limit registering output to initialization stage */
	if (fStatus != GlobalT::kInitialization) 
		ExceptionT::GeneralFail("FEManagerT::RegisterOutput", "output sets can only be registered during initialization");

	int ID = fIOManager->AddElementSet(output_set);
	if (Size() > 1 && Rank() == 0)
	{
		/* file name */
		StringT io_file;
		io_file.Root(fMainIn.filename()); // drop ".in"
		io_file.Root();                   // drop ".p0"
		io_file.Append(".io.ID");

		/* open stream */
		ofstreamT io;
		if (fIOManager->ElementSets().Length() == 1)
		{
			io.open(io_file);
			
			/* write header information */
			io << "# element block ID's for each output ID\n";
			io << "# [output ID] [num blocks] [list of block ID's]\n";
		}
		else
			io.open_append(io_file);

		/* set mode */
		switch (output_set.Mode())
		{
			case  OutputSetT::kElementBlock: 
			{	
				/* write block ID information */
				const ArrayT<StringT>& block_ID = output_set.BlockID();
				io << ID << " " << block_ID.Length();
				for (int i = 0; i < block_ID.Length(); i++)
					io << " " << block_ID[i];
				io << '\n';
				
				break;
			}
			case OutputSetT::kFreeSet: 
			{
				/* no ID's for free sets */
				io << ID << " 0\n";
				
				break;
			}
			case OutputSetT::kElementFromSideSet:
			{
				/* write block ID information */
				const ArrayT<StringT>& block_ID = output_set.BlockID();
				io << ID << " " << block_ID.Length();
				for (int i = 0; i < block_ID.Length(); i++)
					io << " " << block_ID[i];
				io << '\n';
				
				break;
			}
			default:
			{
				ExceptionT::GeneralFail("FEManagerT::RegisterOutput", "unrecognized output set mode: %d", output_set.Mode());
			}
		}
	}
	
	return ID;
}

void FEManagerT::WriteGeometryFile(const StringT& file_name,
	IOBaseT::FileTypeT output_format) const
{
	fIOManager->WriteGeometryFile(file_name, output_format);
}

const OutputSetT& FEManagerT::OutputSet(int ID) const
{
	/* check */
	if (!fIOManager) ExceptionT::GeneralFail("FEManagerT::OutputSet", "I/O manager not initialized");
	return fIOManager->OutputSet(ID);
}

/* (temporarily) direct output away from main out */
void FEManagerT::DivertOutput(const StringT& outfile)
{
	/* check */
	if (!fIOManager) ExceptionT::GeneralFail("FEManagerT::DivertOutput", "I/O manager not initialized");
	fIOManager->DivertOutput(outfile);
	
	/* resolved group */
	if (fCurrentGroup != -1)
		fSO_DivertOutput[fCurrentGroup]	= true;
	else
		fSO_DivertOutput = true;
}

void FEManagerT::RestoreOutput(void)
{
	/* check */
	if (!fIOManager) ExceptionT::GeneralFail("FEManagerT::RestoreOutput", "I/O manager not initialized");
	fIOManager->RestoreOutput();

	/* resolved group */
	if (fCurrentGroup != -1)
		fSO_DivertOutput[fCurrentGroup]	= false;
	else
		fSO_DivertOutput = false;
}

/* cross-linking */
ElementBaseT* FEManagerT::ElementGroup(int groupnumber) const
{
	/* check range */
	if (groupnumber > -1 && groupnumber < fElementGroups->Length())
		return (*fElementGroups)[groupnumber];
	else
		return NULL;
}

int FEManagerT::ElementGroupNumber(const ElementBaseT* pgroup) const
{
	int groupnum = -1;
	for (int i = 0; i < fElementGroups->Length() && groupnum == -1; i++)
		if ((*fElementGroups)[i] == pgroup) groupnum = i;

	return groupnum;
}

int FEManagerT::Rank(void) const { return fComm.Rank(); }
int FEManagerT::Size(void) const { return fComm.Size(); }

const ArrayT<int>* FEManagerT::ProcessorMap(void) const
{
	return fCommManager->ProcessorMap();
}

const ArrayT<int>* FEManagerT::NodeMap(void) const
{
	return fCommManager->NodeMap();
}

const ArrayT<int>* FEManagerT::PartitionNodes(void) const
{
	return fCommManager->PartitionNodes();
}

void FEManagerT::Wait(void) { fComm.Barrier(); }

/* global number of first local equation */
GlobalT::EquationNumberScopeT FEManagerT::EquationNumberScope(int group) const
{
	return fSolvers[group]->EquationNumberScope();
}

int FEManagerT::GetGlobalEquationStart(int group, int start_eq_shift) const
{
#pragma unused(group)

	/* no other equations */
	return 1 + start_eq_shift;
}

int FEManagerT::GetGlobalNumEquations(int group) const
{
	/* no other equations */
	return fNodeManager->NumEquations(group);
}

/* access to integrators */
const eIntegratorT* FEManagerT::eIntegrator(int index) const
{
	/* cast to eIntegratorT */
  return &(fIntegrators[index]->eIntegrator());
}

const nIntegratorT* FEManagerT::nIntegrator(int index) const
{
	/* cast to nIntegratorT */
  return &(fIntegrators[index]->nIntegrator());
}

void FEManagerT::SetTimeStep(double dt) const
{
	/* update the time manager */
	fTimeManager->SetTimeStep(dt);

	/* update all integrators */
	for (int i = 0; i < fIntegrators.Length(); i++)
		fIntegrators[i]->SetTimeStep(dt);

	/* update all solvers */
	for (int i = 0; i < fSolvers.Length(); i++)
		fSolvers[i]->SetTimeStep(dt);
}

/* returns 1 of ALL element groups have interpolant DOF's */
int FEManagerT::InterpolantDOFs(void) const
{
	return fElementGroups->InterpolantDOFs();
}

/* debugging */
void FEManagerT::WriteSystemConfig(ostream& out, int group) const
{
	int old_precision = out.precision();
	out.precision(DBL_DIG);
	
	/* node map */
	const ArrayT<int>* node_map = NodeMap();

	/* nodal data */
	const dArray2DT& coords = fNodeManager->InitialCoordinates();
	ArrayT<FieldT*> fields;
	fNodeManager->CollectFields(group, fields);

	/* dimensions */
	int nnd = coords.MajorDim();
	int nsd = coords.MinorDim();
	int ndf = fNodeManager->NumDOF(group);

	/* force vector */
	const dArrayT& RHS = fSolvers[group]->RHS();

	/* header */
	out << "\n time = " << Time() << '\n';
	out <<   " nodes = " << nnd << '\n';
	int d_width = OutputWidth(out, coords.Pointer());
	out << setw(kIntWidth) << "node"
	    << setw(kIntWidth) << "mapped";
	for (int i0 = 0; i0 < ndf; i0++)
		out << setw(kIntWidth - 2) << "eq[" << i0 + 1 << "]";
	for (int i1 = 0; i1 < nsd; i1++)
		out << setw(d_width - 2) << "x[" << i1 + 1 << "]";

	/* loop over fields */
	for (int i = 0; i < fields.Length(); i++)
	{
		const FieldT& field = *(fields[i]);
	
		/* labels */
		const ArrayT<StringT>& labels = field.Labels();
		
		/* loop over time derivatives */
		StringT suffix = "_";
		for (int k = 0; k <= field.Order(); k++)
		{
			for (int j = 0; j < labels.Length(); j++)
			{
				StringT label = labels[j];
				if (k > 0) label.Append(suffix);
				out << setw(d_width) << label;
			}
			suffix.Append("t");
		}
	}

	for (int i3 = 0; i3 < ndf; i3++)
		out << setw(d_width - 2) << "f[" << i3 + 1 << "]";
	out << '\n';

	/* loop over nodes */
	int shift = ActiveEquationStart(group);
	int num_eq = RHS.Length();
	iArrayT eq(ndf);
	for (int i = 0; i < nnd; i++)
	{
		/* local node number */
		out << setw(kIntWidth) << i+1;
		
		/* mapped node number */
		out << setw(kIntWidth) << ((node_map != NULL) ? (*node_map)[i] : i) + 1;
		
		/* (local) equation numbers - loop over fields */
		int i0 = 0;
		for (int k = 0; k < fields.Length(); k++)
			for (int j = 0; j < fields[k]->NumDOF(); j++)
			{
				eq[i0] = fields[k]->EquationNumber(i,j);
				out << setw(kIntWidth) << eq[i0];
				i0++;
			}
			
		/* coordinates */
		for (int i1 = 0; i1 < nsd; i1++)
			out << setw(d_width) << coords(i, i1);

		/* displacements - loop over fields */
		for (int k = 0; k < fields.Length(); k++)
		{
			const FieldT& field = *(fields[k]);
			for (int l = 0; l <= field.Order(); l++)
			{
				const dArray2DT& u = field[l];
				for (int j = 0; j < u.MinorDim(); j++)
					out << setw(d_width) << u(i,j);
			}
		}

		/* force */
		for (int i3 = 0; i3 < ndf; i3++)
		{
			int eq_i = eq[i3] - shift;
			if (eq_i > -1 && eq_i < num_eq)
				out << setw(d_width) << RHS[eq_i];
			else
				out << setw(d_width) << 0.0;
		}
		
		out << '\n';
	}
	out.flush();

	/* restore */
	out.precision(old_precision);
}

void FEManagerT::RegisterSystemOutput(int group)
{
#pragma unused(group)
	const char caller[] = "FEManagerT::RegisterSystemOutput";

	if (fSO_OutputID.Length() == 0)
		fSO_OutputID.Dimension(NumGroups());
	else if (fSO_OutputID.Length() != NumGroups())
		ExceptionT::GeneralFail(caller);
}

/* interactive */
bool FEManagerT::iDoCommand(const CommandSpecT& command, StringT& line)
{
	try
	{
		if (command.Name() == "ReadRestart")
		{
			StringT file_name;
			command.Argument(0).GetValue(file_name);
			ReadRestart(&file_name);
		}
		else if (command.Name() == "WriteRestart")
		{
			StringT file_name;
			command.Argument(0).GetValue(file_name);
			WriteRestart(&file_name);
		}
		else if (command.Name() == "WriteOutput")
		{
			WriteOutput(Time());
		}
		else
			/* inherited */
			return iConsoleObjectT::iDoCommand(command, line);
	}
	
	catch (ExceptionT::CodeT error)
	{
		cout << "caught exception: " << ExceptionT::ToString(error) << endl;
		return false;
	}
	
	return true;
}

/* describe the parameters needed by the interface */
void FEManagerT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list); 

	/* geometry file format */
	ParameterT geometry_format(ParameterT::Enumeration, "geometry_format");
	geometry_format.AddEnumeration("TahoeII", IOBaseT::kTahoeII);
#ifdef __ACCESS__ 
	geometry_format.AddEnumeration("ExodusII", IOBaseT::kExodusII);
#endif
	geometry_format.SetDefault(IOBaseT::kTahoeII);
	list.AddParameter(geometry_format);

	/* geometry file */
	list.AddParameter(ParameterT(ParameterT::String, "geometry_file"));

	/* output format */
	ParameterT output_format(ParameterT::Enumeration, "output_format");
	output_format.AddEnumeration("Tahoe", IOBaseT::kTahoe);
	output_format.AddEnumeration("TecPlot", IOBaseT::kTecPlot);
	output_format.AddEnumeration("EnSight", IOBaseT::kEnSight);
#ifdef __ACCESS__ 
	output_format.AddEnumeration("ExodusII", IOBaseT::kExodusII);
#endif
	output_format.SetDefault(IOBaseT::kTahoe);
	list.AddParameter(output_format);

	/* restart file name */
	list.AddParameter(ParameterT(ParameterT::String, "restart_file"), ParameterListT::ZeroOrOnce);

	/* restart output increment */
	ParameterT restart_output_inc(ParameterT::Integer, "restart_output_inc");
	restart_output_inc.SetDefault(0);
	restart_output_inc.AddLimit(0, LimitT::LowerInclusive);
	list.AddParameter(restart_output_inc);

	/* verbose echo of input */
	ParameterT echo_input(ParameterT::Boolean, "echo_input");
	echo_input.SetDefault(false);
	list.AddParameter(echo_input);

	/* solve for initial conditions */
	ParameterT compute_IC(ParameterT::Boolean, "compute_IC");
	compute_IC.SetDefault(false);
	list.AddParameter(compute_IC);
}

/* accept parameter list */
void FEManagerT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);
	
	//not implemented
}

/* information about subordinate parameter lists */
void FEManagerT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* time manager */
	sub_list.AddSub("time");

	/* node manager */
	sub_list.AddSub("nodes");

	/* element list */
	sub_list.AddSub("element_list");

	/* node manager */
	sub_list.AddSub("solvers", ParameterListT::OnePlus, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FEManagerT::NewSub(const StringT& list_name) const
{
	FEManagerT* non_const_this = (FEManagerT*) this;

	if (list_name == "time")
		return new TimeManagerT(*non_const_this);
	else if (list_name == "nodes")
		return new NodeManagerT(*non_const_this, *fCommManager);
	else if (list_name == "element_list")
		return new ElementListT(*non_const_this);
	else if (list_name == "linear_solver")
		return non_const_this->New_Solver(GlobalT::kLinearSolver);
	else if (list_name == "nonlinear_solver")
		return non_const_this->New_Solver(GlobalT::kNewtonSolver);	
	else if (list_name == "PCG_solver")
		return non_const_this->New_Solver(GlobalT::kPCGSolver_LS);	
	else /* inherited */
		return ParameterInterfaceT::NewSub(list_name);
}

/* return the description of the given inline subordinate parameter list */
void FEManagerT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	if (sub == "solvers")
	{
		order = ParameterListT::Choice;

		/* linear solver */
		sub_sub_list.AddSub("linear_solver");

		/* nonlinear solver */
		sub_sub_list.AddSub("nonlinear_solver");

		/* nonlinear PCG solver */
		sub_sub_list.AddSub("PCG_solver");
	}
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(sub, order, sub_sub_list);
}

/* returns true if the option was passed on the command line */
bool FEManagerT::CommandLineOption(const char* str, int& index) const
{
	for (int i = 0; i < fArgv.Length(); i++)
		if (fArgv[i] == str) {
			index = i;
			return true;
		}

	/* dummy */
	index = -1;
	return false;
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* "const" function that sets the status flag */
void FEManagerT::SetStatus(GlobalT::StateT status) const
{
	/* non-const this */
	FEManagerT* non_const_this = (FEManagerT*) this;
	non_const_this->fStatus = status;
}

/* look for input file key and check file version */
void FEManagerT::CheckInputFile(void)
{
	/* version check */
	StringT version;
	fMainIn >> version;		
	if (strcmp(version, kCurrentVersion) != 0)
	{
//TEMP - until input is parsed, we will only have spotty support for backward
//       compatibility of input files
		cout << "\n !!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		cout << " FEManagerT::CheckInputFile: input file version is not current. See\n"
		     <<   "     VERSION_NOTES for description of changes:\n";
		cout << "     file version: " << version << '\n';
		cout << "  current version: " << kCurrentVersion << '\n';
		cout << " WARNING: backward compatibility is not completely supported\n";
		cout << " !!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	}	
	else
		cout    << "\n Input file version: " << version << '\n';

	fMainOut << "\n Input file version: " << version << '\n';
}

void FEManagerT::WriteParameters(void) const
{
	fMainOut << "\n T i t l e : " << fTitle << '\n';

	fMainOut << "\n E x e c u t i o n   C o n t r o l   I n f o r m a t i o n :\n\n";
	fMainOut << " Analysis code . . . . . . . . . . . . . . . . . = " << fAnalysisCode << '\n';
	fMainOut << "    eq. " << GlobalT::kLinStatic       << ", linear elastostatic\n";
	fMainOut << "    eq. " << GlobalT::kLinDynamic      << ", linear elastodynamic\n";
	fMainOut << "    eq. " << GlobalT::kNLStatic        << ", nonlinear elastostatic\n";
	fMainOut << "    eq. " << GlobalT::kNLDynamic       << ", nonlinear elastodynamic\n";   	
	fMainOut << "    eq. " << GlobalT::kDR              << ", dynamic relaxation\n";   	
	fMainOut << "    eq. " << GlobalT::kLinExpDynamic   << ", linear explicit dynamic\n";   	
	fMainOut << "    eq. " << GlobalT::kNLExpDynamic    << ", nonlinear explicit dynamic\n";   	
	fMainOut << "    eq. " << GlobalT::kPML             << ", perfectly matched layer (PML)\n";   	
	fMainOut << "    eq. " << GlobalT::kMultiField      << ", general multiple field problem\n";   	

	fModelManager->EchoData (fMainOut);
	IOBaseT temp (fMainOut);
	fMainOut << " Output format . . . . . . . . . . . . . . . . . = " << fOutputFormat  << '\n';
	temp.OutputFormats (fMainOut);

	fMainOut << " Read restart file code  . . . . . . . . . . . . = " << fReadRestart << '\n';
	fMainOut << "    eq. 0, do not read restart file\n";
	fMainOut << "    eq. 1, read restart file\n";
	if (fReadRestart)
		fMainOut << " Restart file. . . . . . . . . . . . . . . . . . = " << fRestartFile << '\n';
	fMainOut << " Restart file increment (at current step size) . = " << fWriteRestart << '\n';
	fMainOut << " Input data print code . . . . . . . . . . . . . = " << fPrintInput << '\n';
	fMainOut << "    eq. 0, non verbose echo of input data\n";
	fMainOut << "    eq. 1, echo all input data\n";
	fMainOut << " Number of solver groups . . . . . . . . . . . . = " << fSolvers.Length()  << '\n';
	fMainOut << endl;
}

/* set the correct fNodeManager type */
void FEManagerT::SetNodeManager(void)
{
	const char caller[] = "FEManagerT::SetNodeManager";

	/* construct */
	try {
		if (!fCommManager) ExceptionT::GeneralFail(caller);
		fNodeManager = new NodeManagerT(*this, *fCommManager);
		fNodeManager->Initialize();			
	
		/* add to console */
		iAddSub(*fNodeManager);
	}
	catch (ExceptionT::CodeT code) {
		ExceptionT::Throw(code, caller, "exception constructing node manager");
	}
}

/* construct element groups */
void FEManagerT::SetElementGroups(void)
{
	/* construct element list */
	fElementGroups = new ElementListT(*this);

	/* echo element data */
	int num_groups;
	fMainIn >> num_groups;
	if (num_groups < 0) ExceptionT::BadInputValue("FEManagerT::SetElementGroups");
	fElementGroups->Dimension(num_groups);
	fElementGroups->EchoElementData(fMainIn, fMainOut);
		
	/* set console */
	for (int i = 0; i < fElementGroups->Length(); i++)
		iAddSub(*((*fElementGroups)[i]));
}

/* set the correct fSolutionDriver type */
void FEManagerT::SetSolver(void)
{
	const char caller[] = "FEManagerT::SetSolver";

	/* equation info */
	int num_groups = fSolvers.Length();
	fGlobalEquationStart.Dimension(num_groups);
	fActiveEquationStart.Dimension(num_groups);
	fGlobalNumEquations.Dimension(num_groups);

	/* no predefined solvers */ 
	if (fAnalysisCode == GlobalT::kMultiField)
	{
		if (fCurrentGroup != -1) throw;
		for (fCurrentGroup = 0; fCurrentGroup < fSolvers.Length(); fCurrentGroup++)
		{
			int index = -1;
			int type = -1;
			fMainIn >> index >> type;
			index--;
			if (fSolvers[index] != NULL)
				ExceptionT::BadInputValue(caller, "solver at index %d is already set", index+1);
	
			/* construct solver */
			fSolvers[index] = New_Solver(type, fCurrentGroup);
		}
		fCurrentGroup = -1;
	}
	else /* support for legacy analysis codes */
	{
		/* should have just one solver */
		if (fSolvers.Length() != 1) ExceptionT::SizeMismatch();
	
		/* solver set by analysis code */
		switch (fAnalysisCode)
		{
			case GlobalT::kLinStatic:
			case GlobalT::kLinDynamic:
			case GlobalT::kLinExpDynamic:
			case GlobalT::kNLExpDynamic:
			case GlobalT::kNLExpDynKfield:
			case GlobalT::kLinStaticHeat:
			case GlobalT::kLinTransHeat:
			case GlobalT::kPML:
				fSolvers[0] = New_Solver(GlobalT::kLinearSolver, 0);
				break;

			case GlobalT::kDR:
				fSolvers[0] = New_Solver(GlobalT::kDRSolver, 0);
				break;

			case GlobalT::kNLStatic:
			case GlobalT::kNLDynamic:
			case GlobalT::kNLStaticKfield:
			case GlobalT::kVarNodeNLStatic:
			case GlobalT::kNLStaticHeat:
			case GlobalT::kNLTransHeat:
			{
				int NL_solver_code;
				fMainIn >> NL_solver_code;
				fSolvers[0] = New_Solver(NL_solver_code, 0);
				break;
			}
			default:
				ExceptionT::BadInputValue(caller, "unknown analysis type: %d", fAnalysisCode);
		}
	}
	
	/* read solver phases */
	AutoArrayT<int> solver_list;
	if (fAnalysisCode == GlobalT::kMultiField) {
	
		int num_phases = -99;
		fMainIn >> num_phases;
		if (num_phases < fSolvers.Length())
			ExceptionT::BadInputValue(caller, "expecting at least %d solver phases: %d",
				fSolvers.Length(), num_phases);

		fSolverPhases.Dimension(num_phases, 3);
		fSolverPhases = -99;

		for (int i = 0; i < fSolverPhases.MajorDim(); i++) {
			int solver = -99; 
			int iters = -99;
			int pass_iters = -99;
			fMainIn >> solver >> iters >> pass_iters;
			solver--;

			/* checks */
			if (solver < 0 || solver >= fSolverPhases.MajorDim()) 
				ExceptionT::BadInputValue(caller, "solver number is out of range: %d", solver+1); 

			if (iters < 1 && iters != -1)
				ExceptionT::BadInputValue(caller, "solver iterations for solver %d must be -1 or > 0: %d", solver+1, iters);

			if (pass_iters < 0 && pass_iters != -1) 
				ExceptionT::BadInputValue(caller, "iterations to pass solver %d must be -1 or >= 0: %d", solver+1, iters);
			
			/* store values */
			solver_list.AppendUnique(solver);
			fSolverPhases(i,0) = solver;
			fSolverPhases(i,1) = iters;
			fSolverPhases(i,2) = pass_iters;
		}
		
		/* number of loops through the solvers */
		fMaxSolverLoops = -99;
		fMainIn >> fMaxSolverLoops;
		if (fMaxSolverLoops < 1 && fMaxSolverLoops != -1)
			ExceptionT::BadInputValue(caller, "expecting -1 or > 0: %d", fMaxSolverLoops);
	}
	else /* set default values for single solver */
	{
		fSolverPhases.Dimension(1, 3);	
		fSolverPhases(0,0) = 0;
		fSolverPhases(0,1) =-1;
		fSolverPhases(0,2) =-1;
		fMaxSolverLoops  = 1;
		solver_list.Append(0);
	}
	
	/* dimension solver phase status array */
	fSolverPhasesStatus.Dimension(fSolverPhases.MajorDim(), kNumStatusFlags);
	fSolverPhasesStatus = 0;
		
	/* echo solver information */
	fMainOut << "\n Multi-solver parameters: " << fSolverPhases.MajorDim() << '\n';
	fMainOut << setw(kIntWidth) << "rank"
	         << setw(kIntWidth) << "group"
	         << setw(2*kIntWidth) << "loop its."
	         << setw(2*kIntWidth) << "pass its." << endl;
	for (int i = 0; i < fSolverPhases.MajorDim(); i++)
		fMainOut << setw(kIntWidth) << i+1
	             << setw(kIntWidth) << fSolverPhases(i,0)
	             << setw(2*kIntWidth) << fSolverPhases(i,1)
	             << setw(2*kIntWidth) << fSolverPhases(i,2) << endl;
	fMainOut << " Maximum number of solver loops. . . . . . . . . = " << fMaxSolverLoops << '\n';
	
	/* check that all solvers hit at least once */
	if (solver_list.Length() != fSolvers.Length())
		ExceptionT::BadInputValue(caller, "must have at least one phase per solver");
	
	/* initialize and register system output */
	fSO_DivertOutput.Dimension(fSolvers.Length());
	fSO_DivertOutput = false;
	fSO_OutputID.Dimension(fSolvers.Length());
	fSO_OutputID = -1;
	if (fCurrentGroup != -1) throw;
	for (fCurrentGroup = 0; fCurrentGroup < fSolvers.Length(); fCurrentGroup++)
	{
		/* reset equation structure */
		SetEquationSystem(fCurrentGroup);

		/* console hierarchy */
		iAddSub(*(fSolvers[fCurrentGroup]));	

		/* writing output */
		if (fSolvers[fCurrentGroup]->Check() != 0)
		{
			/* point connectivities for all nodes (owned by this processor) */
			const ArrayT<int>* partition_nodes = fCommManager->PartitionNodes();
			int nnd = (partition_nodes != NULL) ? partition_nodes->Length() : fNodeManager->NumNodes();
			if (fSO_Connects.MajorDim() != nnd) {
				if (partition_nodes)
					fSO_Connects.Alias(nnd, 1, partition_nodes->Pointer());
				else {
					fSO_Connects.Dimension(fNodeManager->NumNodes(), 1);
					fSO_Connects.SetValueToPosition();
				}
			}
		
			/* collect the fields */
			ArrayT<FieldT*> fields;
			fNodeManager->CollectFields(fCurrentGroup, fields);
		
			/* loop over field and collect labels */
			int ndof = fNodeManager->NumDOF(fCurrentGroup);
			ArrayT<StringT> n_labels(2*ndof);
			int u_dex = 0, f_dex = u_dex + ndof;
			for (int i = 0; i < fields.Length(); i++)
			{
				const FieldT* field = fields[i];
				const ArrayT<StringT>& labels = field->Labels();
				for (int j = 0; j < labels.Length(); j++)
				{
					StringT f_label = "F_";
					f_label.Append(labels[j]);
					n_labels[u_dex++] = labels[j];
					n_labels[f_dex++] = f_label;
				}
			}

			/* register output set */
			OutputSetT output_set(GeometryT::kPoint, fSO_Connects, n_labels, false);
			fSO_OutputID[fCurrentGroup] = RegisterOutput(output_set);
		}
	}
	fCurrentGroup = -1;
}

void FEManagerT::ReadParameters(InitCodeT init)
{
	/* path to parameters file */
	StringT path;
	path.FilePath(fMainIn.filename());

	IOBaseT::FileTypeT format;
	StringT database;

	/* read */
	fMainIn >> fAnalysisCode;
	fMainIn >> format;
 	if (format != IOBaseT::kTahoe)
   	{
		fMainIn >> database;

		/* prepend full path name to database file */
		database.ToNativePathName();      
		database.Prepend(path);
    }
	fMainIn >> fOutputFormat;
	fMainIn >> fReadRestart;
	if (fReadRestart == 1)
	{
		fMainIn >> fRestartFile;
	    
	    /* prepend path */
		fRestartFile.ToNativePathName();
	    fRestartFile.Prepend(path);
	}
	fMainIn >> fWriteRestart;
	fMainIn >> fPrintInput;

	/* check */
	if (fReadRestart  != 0 && fReadRestart  != 1) throw ExceptionT::kBadInputValue;
	if (fWriteRestart < 0) throw ExceptionT::kBadInputValue;
	if (fPrintInput   != 0 && fPrintInput   != 1) throw ExceptionT::kBadInputValue;

	/* initialize the model manager */
	if (!fModelManager->Initialize(format, database, 
		init == kFull || init == kAllButSolver)) /* conditions under which to scan model */
	{
		cout << "\n FEManagerT::ReadParameters: error initializing model manager" << endl;
		throw ExceptionT::kBadInputValue;
	}

	/* read number of equation groups */	
	int num_groups = -1;
	if (fAnalysisCode == GlobalT::kMultiField)
		fMainIn >> num_groups;
	/* support for legacy analysis */
	else
		num_groups = 1;

	/* allocate to that NumGroups is correct */
	fSolvers.Dimension(num_groups);
	fSolvers = NULL;
}

/* set the execution integrator and send to nodes and elements.
* This function constructs the proper drived class integrator
* and passes it to the nodes and elements.  The integrator is
* then cast to a integrator to extract only the FEManagerT
* interface. */
void FEManagerT::SetIntegrator(void)
{
	fMainOut << "\n T i m e   I n t e g r a t o r s:\n";
	
	/* no predefined integrators */
	if (fAnalysisCode == GlobalT::kMultiField)
	{
		/* construct from stream */
		ifstreamT& in = Input();

		int n_int = -1;
		in >> n_int;
		
		fIntegrators.Dimension(n_int);
		fIntegrators = NULL;
		for (int i = 0; i < fIntegrators.Length(); i++)
		{
			int dex = -1;
			TimeManagerT::CodeT code;
			in >> dex >> code;
			
			IntegratorT* integrator = fTimeManager->New_Integrator(code);
			if (!integrator) {
				cout << "\n FEManagerT::SetIntegrator: exception constructing integrator " 
				     << i+1 << " of " << fIntegrators.Length() << endl;
				throw ExceptionT::kBadInputValue;
			}
			fIntegrators[i] = integrator;
		}
	}
	else /* legacy code - a single predefined integrator */
	{
		/* just one */
		fIntegrators.Dimension(1);
		fIntegrators = NULL;
		
		/* set by analysis type */
		IntegratorT* integrator = NULL;
		switch (fAnalysisCode)
		{
			case GlobalT::kLinStatic:
			{
				integrator = fTimeManager->New_Integrator(TimeManagerT::kLinearStatic);
				break;
			}
			case GlobalT::kNLStatic:
			case GlobalT::kNLStaticKfield:
			case GlobalT::kLinStaticHeat:
			case GlobalT::kNLStaticHeat:
			{
				integrator = fTimeManager->New_Integrator(TimeManagerT::kStatic);
				break;
			}
			case GlobalT::kLinTransHeat:
			case GlobalT::kNLTransHeat:
			{
				integrator = fTimeManager->New_Integrator(TimeManagerT::kTrapezoid);
				break;
			}
			case GlobalT::kLinDynamic:
			{
				integrator = fTimeManager->New_Integrator(TimeManagerT::kLinearHHT);
				break;
			}
			case GlobalT::kNLDynamic:
			{
				integrator = fTimeManager->New_Integrator(TimeManagerT::kNonlinearHHT);
				break;
			}
			case GlobalT::kLinExpDynamic:
			case GlobalT::kNLExpDynamic:
			case GlobalT::kNLExpDynKfield:
			case GlobalT::kPML:
			{
				integrator = fTimeManager->New_Integrator(TimeManagerT::kExplicitCD);
				break;
			}			
			default:
				cout << "\nFEManagerT::SetIntegrator: unknown integrator type\n" << endl;
				throw ExceptionT::kBadInputValue;
		}
		
		if (!integrator) throw ExceptionT::kGeneralFail;
		fIntegrators[0] = integrator;
	}
}

/* construct output */
void FEManagerT::SetOutput(void)
{	
	/* set global coordinates */
	fIOManager->SetCoordinates(fNodeManager->InitialCoordinates(), NULL);
	
	/* element groups register output data */
	for (int i = 0; i < fElementGroups->Length(); i++)
		(*fElementGroups)[i]->RegisterOutput();
		
	/* register output from nodes */		
	fNodeManager->RegisterOutput();
}

/* (re-)set system to initial conditions */
ExceptionT::CodeT FEManagerT::InitialCondition(void)
{
	const char caller[] = "FEManagerT::InitialCondition";

	/* state */
	fStatus = GlobalT::kInitialCondition;	

	/* set I/O */		
	fIOManager->NextTimeSequence(SequenceNumber());

	/* set system to initial state */
	fNodeManager->InitialCondition();
	for (int i = 0 ; i < fElementGroups->Length(); i++)
		(*fElementGroups)[i]->InitialCondition();
		
	/* initialize state: solve (t = 0) unless restarted */
	ExceptionT::CodeT error = ExceptionT::kNoError;
	if (!ReadRestart() && fComputeInitialCondition)
	{
		cout << "\n " << caller << ": computing initial conditions" << endl;
	
		/* current step size */
		double dt = TimeStep();
		
		/* override time step */
		SetTimeStep(0.0);
	
		/* initialize the current time step */
		if (error == ExceptionT::kNoError) error = InitStep();
			
		/* solve the current time step */
		if (error == ExceptionT::kNoError) error = SolveStep();
			
		/* close the current time step */
		if (error == ExceptionT::kNoError) error = CloseStep();

		/* restore time step */
		SetTimeStep(dt);

		if (error != ExceptionT::kNoError)
			cout << "\n " << caller << ": error encountered" << endl;

		cout << "\n " << caller << ": computing initial conditions: DONE" << endl;
	}

	return error;
}

/* restart file functions */
bool FEManagerT::ReadRestart(const StringT* file_name)
{
	/* state */
	fStatus = GlobalT::kReadRestart;	

	/* do the read */
	if (fReadRestart || file_name != NULL)
	{
		const StringT& rs_file = (file_name != NULL) ?
			*file_name : fRestartFile;
	
		cout << "\n Restart for sequence: ";
		cout << fTimeManager->SequenceNumber() + 1 << '\n';
		cout <<   "         Restart file: " << rs_file << endl;

		ifstreamT restart(rs_file);			
		if (restart.is_open())
		{
			StringT title;
			title.GetLineFromStream(restart);
			cout << '\n' << title << '\n';
			
			fTimeManager->ReadRestart(restart);		  	                  		  	
			fNodeManager->ReadRestart(restart);
			for (int i = 0 ; i < fElementGroups->Length(); i++)
			{
				/* open new stream for each group */
				StringT file_name = restart.filename();
				file_name.Append(".elem", i);
				ifstreamT restart_elem(file_name);
				restart_elem.precision(DBL_DIG - 1);

				
				/* write element restart data */
				(*fElementGroups)[i]->ReadRestart(restart_elem);
			}

			cout <<   "         Restart file: " << rs_file
			     << ": DONE"<< endl;
		}
		else
			ExceptionT::BadInputValue("FEManagerT::ReadRestart", "could not open file: \"%s\"",
				rs_file.Pointer());
		
		/* relax system with new configuration */
		for (fCurrentGroup = 0; fCurrentGroup < NumGroups(); fCurrentGroup++)
		{
			/* check group */
			GlobalT::RelaxCodeT relax_code = RelaxSystem(fCurrentGroup);
			
			/* reset the equation system */
			if (relax_code == GlobalT::kReEQ || relax_code == GlobalT::kReEQRelax)
		    	SetEquationSystem(fCurrentGroup);

			/* will not resolve the group */
			if (relax_code == GlobalT::kRelax || relax_code == GlobalT::kReEQRelax)
				cout << "\n FEManagerT::ReadRestart: will not resolve group " << fCurrentGroup+1 << endl;
		}
		fCurrentGroup = -1;
		return true;
	}
	else /* no file read */
		return false;
}

bool FEManagerT::WriteRestart(const StringT* file_name) const
{
	/* state */
	SetStatus(GlobalT::kWriteRestart);
	
	/* regular write */
	if (file_name == NULL)
	{
		/* resolve restart flag */
		bool write = false;
		if (fabs(TimeStep()) > kSmall && /* no output if clock is not running */
			fWriteRestart > 0 && StepNumber()%fWriteRestart == 0) write = true;

		/* write file */
		if (write)
		{
			StringT rs_file;
			rs_file.Root(fMainIn.filename());
			rs_file.Append(".rs", StepNumber());
			rs_file.Append("of", NumberOfSteps());
			ofstreamT restart(rs_file);

			/* skip on fail */
			if (restart.is_open())
			{
				/* format the stream */										
				restart.precision(DBL_DIG - 1); //full precision
				restart << fTitle << '\n';
			
				fTimeManager->WriteRestart(restart);
				fNodeManager->WriteRestart(restart);			
				for (int i = 0 ; i < fElementGroups->Length(); i++)
				{
					/* open new stream for each group */
					StringT file_name = restart.filename();
					file_name.Append(".elem", i);
					ofstreamT restart_elem(file_name);
					restart_elem.precision(DBL_DIG - 1);
				
					/* write element restart data */
					(*fElementGroups)[i]->WriteRestart(restart_elem);
				}
					
				return true;
			}
			else
			{
				cout <<  "\n FEManagerT::WriteRestart: could not open restart file: "
				     << rs_file << endl;	     
				return false;
			}
		}
		else /* no write */
			return false;
	}
	else
	{
		ofstreamT restart(*file_name);

		/* skip on fail */
		if (restart.is_open())
		{
			/* format the stream */										
			restart.precision(DBL_DIG - 1); //full precision
			restart << fTitle << '\n';
			
			fTimeManager->WriteRestart(restart);
			fNodeManager->WriteRestart(restart);			
			for (int i = 0 ; i < fElementGroups->Length(); i++)
			{
				/* open new stream for each group */
				StringT file_name = restart.filename();
				file_name.Append(".elem", i);
				ofstreamT restart_elem(file_name);
				restart_elem.precision(DBL_DIG - 1);
				
				/* write element restart data */
				(*fElementGroups)[i]->WriteRestart(restart_elem);
			}
				
			return true;
		}
		else
		{
			cout <<  "\n FEManagerT::WriteRestart: could not open restart file: "
			     << *file_name << endl;	
			return false;
		}
	}	
}

/* final steps in solver configuration
* (1) signal nodes to assign equation numbers
* (2) renumber if needed
* (3) set numbering scope
* (4) collect equations and send to solver
* (5) signal solver for final configuration */
void FEManagerT::SetEquationSystem(int group, int start_eq_shift)
{
	/* equation number scope */
	GlobalT::EquationNumberScopeT equation_scope = 
		fSolvers[group]->EquationNumberScope();

	/* assign (local) equation numbers */
	fNodeManager->SetEquationNumbers(group);
	fGlobalEquationStart[group] = GetGlobalEquationStart(group, start_eq_shift);
	fActiveEquationStart[group] = (equation_scope == GlobalT::kGlobal) ? 
		fGlobalEquationStart[group] : 1;
	fGlobalNumEquations[group]  = GetGlobalNumEquations(group);

	/* renumber locally */
	if (fSolvers[group]->RenumberEquations())
	{
		int num_fields = fNodeManager->NumFields(group);
		if (num_fields == 1)
		{
			/* lists of connectivities */
			AutoArrayT<const iArray2DT*> connects_1;
			AutoArrayT<const RaggedArray2DT<int>*> connects_2;
			AutoArrayT<const iArray2DT*> equivalent_nodes;

			/* collect nodally generated DOF's */
			fNodeManager->ConnectsU(group, connects_1, connects_2, equivalent_nodes);

			/* collect element groups */
			for (int i = 0 ; i < fElementGroups->Length(); i++)
				if ((*fElementGroups)[i]->InGroup(group))
					(*fElementGroups)[i]->ConnectsU(connects_1, connects_2);		

			/* renumber equations */
			try { fNodeManager->RenumberEquations(group, connects_1, connects_2); }
			catch (ExceptionT::CodeT exception) {
				cout << "\n FEManagerT::SetEquationSystem: could not renumber equations: exception: " 
				     << exception << endl;
			}
		}
		else /* renumbering does not support multiple fields in the same group
		      * because each row in the equations arrays is assumed to correspond
		      * to a unique tag */
		{
			cout << "\n FEManagerT::SetEquationSystem: equations could not be renumbered\n"
			     <<   "     because group " << group+1 << " contains " << num_fields << " fields." << endl;
		}
	}

	/* set equation number scope */
	fNodeManager->SetEquationNumberScope(group, equation_scope);
	
	/* collect interaction equations and send to solver */
	SendEqnsToSolver(group);
	
	/* final step in solver configuration */
	fSolvers[group]->Initialize(
		fGlobalNumEquations[group], 
		fNodeManager->NumEquations(group),
		fActiveEquationStart[group]);
}

void FEManagerT::SendEqnsToSolver(int group) const
{
	/* dynamic arrays */
	AutoArrayT<const iArray2DT*> eq_1;
	AutoArrayT<const RaggedArray2DT<int>*> eq_2;
	
	/* collect equation sets */
	fNodeManager->Equations(group, eq_1, eq_2);
	for (int i = 0 ; i < fElementGroups->Length(); i++)
		if ((*fElementGroups)[i]->InGroup(group))
			(*fElementGroups)[i]->Equations(eq_1, eq_2);

	/* send lists to solver */
	for (int j = 0; j < eq_1.Length(); j++)
		fSolvers[group]->ReceiveEqns(*(eq_1[j]));

	for (int k = 0; k < eq_2.Length(); k++)
		fSolvers[group]->ReceiveEqns(*(eq_2[k]));
}

/* construct a new CommManagerT */
CommManagerT* FEManagerT::New_CommManager(void) const
{
	if (!fModelManager) 
		ExceptionT::GeneralFail("FEManagerT::New_CommManager", "need ModelManagerT");

	CommManagerT* comm_man = new CommManagerT(fComm, *fModelManager);
	return comm_man;
}

/*************************************************************************
 * Private
 *************************************************************************/

SolverT* FEManagerT::New_Solver(int code, int group)
{
	const char caller[] = "FEManagerT::New_Solver";

	/* construct solver */
	SolverT* solver = NULL;
	switch (code)
	{
		case GlobalT::kLinearSolver:
			solver = new LinearSolver(*this, group);
			break;

		case GlobalT::kDRSolver:
			solver = new DRSolver(*this, group);
			break;
	
		case GlobalT::kNewtonSolver:
			solver = new NLSolver(*this, group);
			break;

		case GlobalT::kK0_NewtonSolver:
			solver = new NLK0Solver(*this, group);
			break;

		case GlobalT::kModNewtonSolver:
			solver = new NLSolverX(*this, group);
			break;

		case GlobalT::kExpCD_DRSolver:
			solver = new ExpCD_DRSolver(*this, group);
			break;

		case GlobalT::kNewtonSolver_LS:				
			solver = new NLSolver_LS(*this, group);
			break;

		//TEMP
		case GlobalT::kNewtonSolver_LSX:				
			solver = new NLSolver_LSX(*this, group);
			break;

		case GlobalT::kPCGSolver_LS:				
			solver = new PCGSolver_LS(*this, group);
			break;

		case GlobalT::kiNewtonSolver_LS:				
			solver = new iNLSolver_LS(*this, group);
			break;

		case GlobalT::kNOXSolver:
		{
#ifdef __NOX__
			//TEMP - need to figure out which set of DOF's to send to the solver.
			//       This gets complicated since each group could have more than
			//       one time integrator, which sets the order of the DOF. For now,
			//       grab the integrator from the first field in this group and
			//       let its integrator decide.
			
			/* all fields in the group */
			ArrayT<FieldT*> fields;
			fNodeManager->CollectFields(group, fields);
			if (fields.Length() < 1) ExceptionT::GeneralFail(caller);
			
			const nIntegratorT& integrator = fields[0]->nIntegrator();
			solver = new NOXSolverT(*this, group, integrator.OrderOfUnknown());
			break;
#else
			ExceptionT::GeneralFail(caller, "NOX not installed: %d", GlobalT::kNOXSolver);
#endif	
		}		
		default:
			ExceptionT::BadInputValue(caller, "unknown nonlinear solver code: %d", code);
	}

	/* fail */
	if (!solver) ExceptionT::GeneralFail(caller, "failed");

	return solver;
}

SolverT* FEManagerT::New_Solver(GlobalT::SolverTypeT solver_type)
{
	const char caller[] = "FEManagerT::New_Solver";

	/* construct solver */
	SolverT* solver = NULL;
	switch (solver_type)
	{
		case GlobalT::kLinearSolver:
			solver = new LinearSolver(*this);
			break;
			
		case GlobalT::kNewtonSolver:
			solver = new NLSolver(*this);
			break;

		case GlobalT::kPCGSolver_LS:				
			solver = new PCGSolver_LS(*this);
			break;

		default:
			ExceptionT::BadInputValue(caller, "unknown nonlinear solver code: %d", solver_type);
	}
	return solver;
}
