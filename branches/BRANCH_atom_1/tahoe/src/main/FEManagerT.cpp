/* $Id: FEManagerT.cpp,v 1.50.2.3 2002-12-16 09:23:51 paklein Exp $ */
/* created: paklein (05/22/1996) */
#include "FEManagerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <string.h>
#include <float.h>
#include <ctype.h>

#include "toolboxConstants.h"
#include "ExceptionT.h"

#include "fstreamT.h"
#include "TimeManagerT.h"
#include "ModelManagerT.h"
#include "ElementBaseT.h"
#include "IOManager.h"
#include "OutputSetT.h"
#include "CommandSpecT.h"
#include "ArgSpecT.h"
#include "eControllerT.h"
#include "nControllerT.h"
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

using namespace Tahoe;

/* File/Version Control */
const char kCurrentVersion[] = "v3.4.1";
const char kProgramName[]    = "tahoe";

/* static methods */
const char* FEManagerT::Version(void) { return kCurrentVersion; }

/* constructor */
FEManagerT::FEManagerT(ifstreamT& input, ofstreamT& output, CommunicatorT& comm):
	fMainIn(input),
	fMainOut(output),
	fComm(comm),
	fStatus(GlobalT::kConstruction),
	fTimeManager(NULL),
	fNodeManager(NULL),
	fIOManager(NULL),
	fCommManager(NULL),
	fMaxSolverLoops(0),
	fRestartCount(0),
	fGlobalEquationStart(0),
	fActiveEquationStart(0),
	fGlobalNumEquations(0),
	fCurrentGroup(-1)
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
	
	for (int i = 0; i < fSolvers.Length(); i++)
		delete fSolvers[i];
	
	for (int i = 0; i < fControllers.Length(); i++)
		delete fControllers[i];

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

	/* set communication manager */
	fCommManager = New_CommManager();
	if (!fCommManager) throw ExceptionT::kOutOfMemory;
	if (verbose) cout << "    FEManagerT::Initialize: comm manager" << endl;
	
	/* construct the managers */
	fTimeManager = new TimeManagerT(*this);
	if (!fTimeManager) throw ExceptionT::kOutOfMemory;
	iAddSub(*fTimeManager);
	if (verbose) cout << "    FEManagerT::Initialize: time" << endl;

	/* set time integration controller */
	SetController();
	if (verbose) cout << "    FEManagerT::Initialize: controller" << endl;

	/* initial configuration of communication manager */
	if (fModelManager->DatabaseFormat() != IOBaseT::kTahoe)
	{
		fComm.SetLog(cout);
		fCommManager->Configure();
	}

//TEMP
else ExceptionT::Stop("FEManagerT::Initialize", "kTahoe format not supported");

	/* resolve node type here */
	SetNodeManager();
	if (verbose) cout << "    FEManagerT::Initialize: nodal data" << endl;
	
	/* construct element groups */
	SetElementGroups();
	if (verbose) cout << "    FEManagerT::Initialize: element groups" << endl;

	/* initial configuration of communication manager */
	if (fModelManager->DatabaseFormat() == IOBaseT::kTahoe)
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
	fTimeManager->Top();
	while (fTimeManager->NextSequence())
	{	
		/* set to initial condition */
		InitialCondition();

		/* read kinematic data from restart file */
		ReadRestart();

		/* loop over time increments */
		bool seq_OK = true;
		while (seq_OK && fTimeManager->Step())
		{
			/* running status flag */
			ExceptionT::CodeT error = ExceptionT::kNoError;		

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
					cout << "\n FEManagerT::Solve: trying to recover: error: " << ExceptionT::ToString(error) << endl;
				
					/* reset system configuration */
					error = ResetStep();
					
					if (error == ExceptionT::kNoError)
						/* cut time increment */
						seq_OK = DecreaseLoadStep();
					else
						seq_OK = false;
					break;
				}
				default:
					cout << "\n FEManagerT::Solve: no recovery for error: " << ExceptionT::ToString(error) << endl;
					seq_OK = false;
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
	for (int i = 0 ; i < fElementGroups.Length(); i++)
	{
		GlobalT::SystemTypeT e_type = fElementGroups[i]->TangentType();

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

	/* nodes - ALL groups */
	for (fCurrentGroup = 0; fCurrentGroup < NumGroups(); fCurrentGroup++)
		fNodeManager->ResetStep(fCurrentGroup);
	fCurrentGroup = -1;
	
	/* elements */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->ResetStep();
		
	/* solver - ALL groups */
	for (fCurrentGroup = 0; fCurrentGroup < NumGroups(); fCurrentGroup++)
		fSolvers[fCurrentGroup]->ResetStep();
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
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		if (fElementGroups[i]->InGroup(group))
			fElementGroups[i]->FormLHS(sys_type);
}

void FEManagerT::FormRHS(int group) const
{
	/* state */
	SetStatus(GlobalT::kFormRHS);

	/* nodal force contribution - F(t) */
	fNodeManager->FormRHS(group);

	/* element contribution */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		if (fElementGroups[i]->InGroup(group))
			fElementGroups[i]->FormRHS();
		
	/* output system info (debugging) */
	if (fSolvers[group]->Check() == GlobalMatrixT::kPrintRHS)
		WriteSystemConfig(fMainOut, group);
}

/* collect the internal force on the specified node */
void FEManagerT::InternalForceOnNode(const FieldT& field, int node, dArrayT& force) const
{
	/* initialize */
	force = 0.0;

	/* element contribution */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->AddNodalForce(field, node, force);
}

ExceptionT::CodeT FEManagerT::InitStep(void)
{
	try {
	/* state */
	SetStatus(GlobalT::kInitStep);
	
	/* set the default value for the output time stamp */
	fIOManager->SetOutputTime(Time());

	/* check group flag */
	if (fCurrentGroup != -1) throw ExceptionT::kGeneralFail;

	/* nodes - ALL groups*/
	for (fCurrentGroup = 0; fCurrentGroup < NumGroups(); fCurrentGroup++)
		fNodeManager->InitStep(fCurrentGroup);
	fCurrentGroup = -1;

	/* elements */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->InitStep();
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

		/* status */
		iArray2DT solve_status(fSolverPhases.MajorDim(), 3);

		while (!all_pass && 
			loop_count < fMaxSolverLoops &&
			status != SolverT::kFailed)
		{
			/* clear status */
			solve_status = 0;
		 
			/* one solver after the next */
			all_pass = true;
			for (int i = 0; status != SolverT::kFailed && i < fSolverPhases.MajorDim(); i++)
			{
				/* group parameters */
				fCurrentGroup = fSolverPhases(i,0);
				int iter  = fSolverPhases(i,1);
				int pass  = fSolverPhases(i,2);
			
				/* call solver */
				status = fSolvers[fCurrentGroup]->Solve(iter);
				
				/* check result */
				solve_status(i,0) = fCurrentGroup;
				solve_status(i,1) = fSolvers[fCurrentGroup]->IterationNumber();
				if (status == SolverT::kFailed) {
					all_pass = false;
					solve_status(i,2) = -1;					
				}
				else if (status == SolverT::kConverged && 
					(pass == -1 || solve_status(i,1) <= pass))
				{
					all_pass = all_pass && true; /* must all be true */
					solve_status(i,2) = 1;
				}
				else
				{
					all_pass = false;
					solve_status(i,2) = 0;
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
				solve_status.WriteNumbered(cout);
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

	/* write output BEFORE closing nodes and elements */
	fTimeManager->CloseStep();

	/* check group flag */
	if (fCurrentGroup != -1) throw ExceptionT::kGeneralFail;

	/* nodes - ALL groups */
	for (fCurrentGroup = 0; fCurrentGroup < NumGroups(); fCurrentGroup++)
		fNodeManager->CloseStep(fCurrentGroup);
	fCurrentGroup = -1;

	/* elements */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->CloseStep();
		
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
	GlobalT::RelaxCodeT relax = GlobalT::kNoRelax;
	
	/* check node manager */
	relax = GlobalT::MaxPrecedence(relax, fNodeManager->RelaxSystem(group));
		
	/* check element groups - must touch all of them to reset */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		if (fElementGroups[i]->InGroup(group))
			relax = GlobalT::MaxPrecedence(relax, fElementGroups[i]->RelaxSystem());

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

		/* nodes */
		fNodeManager->WriteOutput();

		/* elements */
		for (int i = 0; i < fElementGroups.Length(); i++)
			fElementGroups[i]->WriteOutput();
	}
	
	catch (ExceptionT::CodeT error) { ExceptionT::Throw(error, "FEManagerT::WriteOutput"); }
}

void FEManagerT::WriteOutput(int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
	fIOManager->WriteOutput(ID, n_values, e_values);
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
		if (output_set.Mode() == OutputSetT::kElementBlock) {	
			/* write block ID information */
			const ArrayT<StringT>& block_ID = output_set.BlockID();
			io << ID << " " << block_ID.Length();
			for (int i = 0; i < block_ID.Length(); i++)
				io << " " << block_ID[i];
			io << '\n';
		}
		else if (output_set.Mode() == OutputSetT::kFreeSet) {
			/* no ID's for free sets */
			io << ID << " 0\n";
		}
		else ExceptionT::GeneralFail("FEManagerT::RegisterOutput", "unrecognized output set mode: %d", output_set.Mode());
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
}

void FEManagerT::RestoreOutput(void)
{
	/* check */
	if (!fIOManager) ExceptionT::GeneralFail("FEManagerT::RestoreOutput", "I/O manager not initialized");
	fIOManager->RestoreOutput();
}

/* cross-linking */
ElementBaseT* FEManagerT::ElementGroup(int groupnumber) const
{
	/* check range */
	if (groupnumber > -1 && groupnumber < fElementGroups.Length())
		return fElementGroups[groupnumber];
	else
		return NULL;
}

int FEManagerT::ElementGroupNumber(const ElementBaseT* pgroup) const
{
	int groupnum = -1;
	for (int i = 0; i < fElementGroups.Length() && groupnum == -1; i++)
		if (fElementGroups[i] == pgroup) groupnum = i;

	return groupnum;
}

#if 0
int FEManagerT::GlobalEquationNumber(int nodenum, int dofnum) const
{
	return fNodeManager->GlobalEquationNumber(nodenum, dofnum);
}
#endif

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

void FEManagerT::IncomingNodes(iArrayT& nodes_in ) const {  nodes_in.Free(); }
void FEManagerT::OutgoingNodes(iArrayT& nodes_out) const { nodes_out.Free(); }

void FEManagerT::RecvExternalData(dArray2DT& external_data)
{
#pragma unused(external_data)
	ExceptionT::GeneralFail("FEManagerT::RecvExternalData", "invalid request for external data");
}

void FEManagerT::SendExternalData(const dArray2DT& all_out_data)
{
#pragma unused(all_out_data)
	ExceptionT::GeneralFail("FEManagerT::RecvExternalData", "invalid send of external data");
}

void FEManagerT::SendRecvExternalData(const iArray2DT& all_out_data, iArray2DT& external_data)
{
#pragma unused(all_out_data)
#pragma unused(external_data)
	ExceptionT::GeneralFail("FEManagerT::RecvExternalData", "invalid exchange of external data");
}

void FEManagerT::Wait(void)
{
// do nothing
}

/* global number of first local equation */
GlobalT::EquationNumberScopeT FEManagerT::EquationNumberScope(int group) const
{
	return fSolvers[group]->EquationNumberScope();
}

int FEManagerT::GetGlobalEquationStart(int group) const
{
#pragma unused(group)

	/* no other equations */
	return 1;
}

int FEManagerT::GetGlobalNumEquations(int group) const
{
	/* no other equations */
	return fNodeManager->NumEquations(group);
}

/* access to controllers */
eControllerT* FEManagerT::eController(int index) const
{
	/* cast to eControllerT */
#ifdef __NO_RTTI__
	eControllerT* e_controller = (eControllerT*) fControllers[index];
		//NOTE: cast should be safe for all cases
#else
	eControllerT* e_controller = dynamic_cast<eControllerT*>(fControllers[index]);
	if (!e_controller) ExceptionT::GeneralFail();
#endif

	return e_controller;
}

nControllerT* FEManagerT::nController(int index) const
{
	/* cast to eControllerT */
#ifdef __NO_RTTI__
	nControllerT* n_controller = (nControllerT*) fControllers[index];
		//NOTE: cast should be safe for all cases
#else
	nControllerT* n_controller = dynamic_cast<nControllerT*>(fControllers[index]);
	if (!n_controller) ExceptionT::GeneralFail();
#endif

	return n_controller;
}

void FEManagerT::SetTimeStep(double dt) const
{
	//TEMP - for ALL controllers
	for (int i = 0; i < fControllers.Length(); i++)
		fControllers[i]->SetTimeStep(dt);
}

/* returns 1 of ALL element groups have interpolant DOF's */
int FEManagerT::InterpolantDOFs(void) const
{
	return fElementGroups.InterpolantDOFs();
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
	int ndf = 0;
	for (int i = 0; i < fields.Length(); i++)
		ndf += fields[i]->NumDOF();

	/* force vector */
	const dArrayT& RHS = fSolvers[group]->RHS();

	/* header */
	out << "\n time = " << Time() << '\n';
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
	/* construct */
	fNodeManager = new NodeManagerT(*this);
	if (!fNodeManager) ExceptionT::OutOfMemory();	
	fNodeManager->Initialize();			

	/* add to console */
	iAddSub(*fNodeManager);	
}

	/* construct element groups */
void FEManagerT::SetElementGroups(void)
{
	/* echo element data */
	int num_groups;
	fMainIn >> num_groups;
	if (num_groups < 1) ExceptionT::BadInputValue();
	fElementGroups.Dimension(num_groups);
	fElementGroups.EchoElementData(fMainIn, fMainOut, *this);
		
	/* set console */
	for (int i = 0; i < fElementGroups.Length(); i++)
		iAddSub(*(fElementGroups[i]));
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
				fSolvers[0] = New_Solver(SolverT::kLinear, 0);
				break;

			case GlobalT::kDR:
				fSolvers[0] = New_Solver(SolverT::kDR, 0);
				break;

			case GlobalT::kNLStatic:
			case GlobalT::kNLDynamic:
			case GlobalT::kNLStaticKfield:
			case GlobalT::kVarNodeNLStatic:
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
	if (fSolvers.Length() > 1) {
	
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
	
	/* initialize */
	if (fCurrentGroup != -1) throw;
	for (fCurrentGroup = 0; fCurrentGroup < fSolvers.Length(); fCurrentGroup++)
	{
		/* reset equation structure */
		SetEquationSystem(fCurrentGroup);

		/* console hierarchy */
		iAddSub(*(fSolvers[fCurrentGroup]));	
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

/* set the execution controller and send to nodes and elements.
* This function constructs the proper drived class controller
* and passes it to the nodes and elements.  The controller is
* then cast to a controller to extract only the FEManagerT
* interface. */
void FEManagerT::SetController(void)
{
	fMainOut << "\n T i m e   I n t e g r a t o r s:\n";
	
	/* no predefined integrators */
	if (fAnalysisCode == GlobalT::kMultiField)
	{
		/* construct from stream */
		ifstreamT& in = Input();

		int n_int = -1;
		in >> n_int;
		
		fControllers.Dimension(n_int);
		fControllers = NULL;
		for (int i = 0; i < fControllers.Length(); i++)
		{
			int dex = -1;
			TimeManagerT::CodeT code;
			in >> dex >> code;
			
			ControllerT* controller = fTimeManager->New_Controller(code);
			if (!controller) {
				cout << "\n FEManagerT::SetController: exception constructing controller " 
				     << i+1 << " of " << fControllers.Length() << endl;
				throw ExceptionT::kBadInputValue;
			}
			fControllers[i] = controller;
		}
	}
	else /* legacy code - a single predefined integrator */
	{
		/* just one */
		fControllers.Dimension(1);
		fControllers = NULL;
		
		/* set by analysis type */
		ControllerT* controller = NULL;
		switch (fAnalysisCode)
		{
			case GlobalT::kLinStatic:
			{
				controller = fTimeManager->New_Controller(TimeManagerT::kLinearStatic);
				break;
			}
			case GlobalT::kNLStatic:
			case GlobalT::kNLStaticKfield:
			case GlobalT::kLinStaticHeat:
			{
				controller = fTimeManager->New_Controller(TimeManagerT::kStatic);
				break;
			}
			case GlobalT::kLinTransHeat:
			{
				controller = fTimeManager->New_Controller(TimeManagerT::kTrapezoid);
				break;
			}
			case GlobalT::kLinDynamic:
			{
				controller = fTimeManager->New_Controller(TimeManagerT::kLinearHHT);
				break;
			}
			case GlobalT::kNLDynamic:
			{
				controller = fTimeManager->New_Controller(TimeManagerT::kNonlinearHHT);
				break;
			}
			case GlobalT::kLinExpDynamic:
			case GlobalT::kNLExpDynamic:
			case GlobalT::kNLExpDynKfield:
			case GlobalT::kPML:
			{
				controller = fTimeManager->New_Controller(TimeManagerT::kExplicitCD);
				break;
			}			
			default:
				cout << "\nFEManagerT::SetController: unknown controller type\n" << endl;
				throw ExceptionT::kBadInputValue;
		}
		
		if (!controller) throw ExceptionT::kGeneralFail;
		fControllers[0] = controller;
	}
}

/* construct output */
void FEManagerT::SetOutput(void)
{
	StringT file_name(fMainIn.filename());
	fIOManager = new IOManager(fMainOut, kProgramName, kCurrentVersion, fTitle,
		file_name, fOutputFormat);	
	if (!fIOManager) throw ExceptionT::kOutOfMemory;
	
	/* set global coordinates */
	fIOManager->SetCoordinates(fNodeManager->InitialCoordinates(), NULL);
	
	/* element groups register output data */
	for (int i = 0; i < fElementGroups.Length(); i++)
		fElementGroups[i]->RegisterOutput();
		
	/* register output from nodes */		
	fNodeManager->RegisterOutput();
}

/* (re-)set system to initial conditions */
void FEManagerT::InitialCondition(void)
{
	/* state */
	fStatus = GlobalT::kInitialCondition;	

	/* set I/O */		
	fIOManager->NextTimeSequence(SequenceNumber());

	/* set system to initial state */
	fNodeManager->InitialCondition();
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->InitialCondition();
}

/* restart file functions */
void FEManagerT::ReadRestart(const StringT* file_name)
{
	/* state */
	fStatus = GlobalT::kReadRestart;	

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
			for (int i = 0 ; i < fElementGroups.Length(); i++)
				fElementGroups[i]->ReadRestart(restart);

			cout <<   "         Restart file: " << rs_file
			     << ": DONE"<< endl;
		}
		else
		{
			cout << "\n FEManagerT::ReadRestart: could not open file: "
			     << rs_file << endl;
			throw ExceptionT::kBadInputValue;
		}
	}
}

void FEManagerT::WriteRestart(const StringT* file_name) const
{
	/* state */
	SetStatus(GlobalT::kWriteRestart);
	
	/* regular write */
	if (file_name == NULL)
	{
		/* non-const counter */
		FEManagerT* non_const = (FEManagerT*) this;
		int& counter = non_const->fRestartCount;

		/* resolve restart flag */
		bool write = false;
		if (fWriteRestart > 0 && ++counter == fWriteRestart) write = true;

		/* write file */
		if (write)
		{
			/* reset */
			counter = 0;

			StringT rs_file;
			rs_file.Root(fMainIn.filename());
			rs_file.Append(".rs", StepNumber());
			ofstreamT restart(rs_file);

			/* skip on fail */
			if (restart.is_open())
			{
				/* format the stream */										
				restart.precision(DBL_DIG - 1); //full precision
				restart << fTitle << '\n';
			
				fTimeManager->WriteRestart(restart);
				fNodeManager->WriteRestart(restart);			
				for (int i = 0 ; i < fElementGroups.Length(); i++)
					fElementGroups[i]->WriteRestart(restart);
			}
			else
				cout <<  "\n FEManagerT::WriteRestart: could not open restart file: "
				     << rs_file << endl;
		}
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
			for (int i = 0 ; i < fElementGroups.Length(); i++)
				fElementGroups[i]->WriteRestart(restart);
		}
		else
			cout <<  "\n FEManagerT::WriteRestart: could not open restart file: "
			     << *file_name << endl;	
	}	
}

/* final steps in solver configuration
* (1) signal nodes to assign equation numbers
* (2) renumber if needed
* (3) set numbering scope
* (4) collect equations and send to solver
* (5) signal solver for final configuration */
void FEManagerT::SetEquationSystem(int group)
{
	/* equation number scope */
	GlobalT::EquationNumberScopeT equation_scope = 
		fSolvers[group]->EquationNumberScope();

	/* assign (local) equation numbers */
	fNodeManager->SetEquationNumbers(group);
	fGlobalEquationStart[group] = GetGlobalEquationStart(group);
	fActiveEquationStart[group] = (equation_scope == GlobalT::kGlobal) ? 
		fGlobalEquationStart[group] : 1;
	fGlobalNumEquations[group]  = GetGlobalNumEquations(group);

	/* renumber locally */
	if (fSolvers[group]->RenumberEquations())
	{
		/* lists of connectivities */
		AutoArrayT<const iArray2DT*> connects_1;
		AutoArrayT<const RaggedArray2DT<int>*> connects_2;
	
		/* collect nodally generated DOF's */
		fNodeManager->ConnectsU(group, connects_1, connects_2);
	
		/* collect element groups */
		for (int i = 0 ; i < fElementGroups.Length(); i++)
			if (fElementGroups[i]->InGroup(group))
				fElementGroups[i]->ConnectsU(connects_1, connects_2);		
	
		/* renumber equations */
		try { fNodeManager->RenumberEquations(group, connects_1, connects_2); }
		catch (ExceptionT::CodeT exception) {
			cout << "\n FEManagerT::SetEquationSystem: could not renumber equations: exception: " 
			     << exception << endl;
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
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		if (fElementGroups[i]->InGroup(group))
			fElementGroups[i]->Equations(eq_1, eq_2);

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
		case SolverT::kLinear:
			solver = new LinearSolver(*this, group);
			break;

		case SolverT::kDR:
			solver = new DRSolver(*this, group);
			break;
	
		case SolverT::kNewtonSolver:
			solver = new NLSolver(*this, group);
			break;

		case SolverT::kK0_NewtonSolver:
			solver = new NLK0Solver(*this, group);
			break;

		case SolverT::kModNewtonSolver:
			solver = new NLSolverX(*this, group);
			break;

		case SolverT::kExpCD_DRSolver:
			solver = new ExpCD_DRSolver(*this, group);
			break;

		case SolverT::kNewtonSolver_LS:				
			solver = new NLSolver_LS(*this, group);
			break;

		case SolverT::kPCGSolver_LS:				
			solver = new PCGSolver_LS(*this, group);
			break;

		case SolverT::kiNewtonSolver_LS:				
			solver = new iNLSolver_LS(*this, group);
			break;

		case SolverT::kNOX:
		{
#ifdef __NOX__
			//TEMP - need to figure out which set of DOF's to send to the solver.
			//       This gets complicated since each group could have more than
			//       one time integrator, which sets the order of the DOF. For now,
			//       grab the controller from the first field in this group and
			//       let its integrator decide.
			
			/* all fields in the group */
			ArrayT<FieldT*> fields;
			fNodeManager->CollectFields(group, fields);
			if (fields.Length() < 1) throw ExceptionT::kGeneralFail;
			
			const nControllerT& controller = fields[0]->nController();
			solver = new NOXSolverT(*this, group, controller.OrderOfUnknown());
			break;
#else
			ExceptionT::GeneralFail(caller, "NOX not installed: %d", SolverT::kNOX);
#endif	
		}		
		default:
			ExceptionT::BadInputValue(caller, "unknown nonlinear solver code: %d", code);
	}

	/* fail */
	if (!solver) ExceptionT::GeneralFail(caller, "failed");

	return solver;
}
