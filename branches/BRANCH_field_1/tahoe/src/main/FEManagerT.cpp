/* $Id: FEManagerT.cpp,v 1.32.2.11 2002-05-11 20:59:30 paklein Exp $ */
/* created: paklein (05/22/1996) */
#include "FEManagerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <string.h>
#include <float.h>
#include <ctype.h>

#include "Constants.h"
#include "ExceptionCodes.h"

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

/* File/Version Control */
const char* kCurrentVersion = "v3.4.1";
const char* kProgramName    = "tahoe";

/* exception strings */
const char* eExceptionStrings[] = {
/* 0 */ "no error",
/* 1 */ "general fail",
/* 2 */ "stop",
/* 3 */ "out of memory",
/* 4 */ "index out of range",
/* 5 */ "dimension mismatch",
/* 6 */ "invalid value read from input",
/* 7 */ "zero or negative jacobian",
/* 8 */ "MPI message passing error",
/* 9 */ "database read failure",
/*10 */ "unknown"};	

/* constructor */
FEManagerT::FEManagerT(ifstreamT& input, ofstreamT& output):
	fMainIn(input),
	fMainOut(output),
	fStatus(GlobalT::kConstruction),
	fTimeManager(NULL),
	fNodeManager(NULL),
	fIOManager(NULL),
	fRestartCount(0),
	fGlobalEquationStart(0),
	fActiveEquationStart(0),
	fGlobalNumEquations(0)
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
	if (!fModelManager) throw eOutOfMemory;
	if (verbose) cout << "    FEManagerT::Initialize: input" << endl;

	/* main parameters */
	ReadParameters(init);
	if (init == kParametersOnly) return;
	WriteParameters();
	if (verbose) cout << "    FEManagerT::Initialize: execution parameters" << endl;
	
	/* construct the managers */
	fTimeManager = new TimeManagerT(*this);
	if (!fTimeManager) throw eOutOfMemory;
	iAddSub(*fTimeManager);
	if (verbose) cout << "    FEManagerT::Initialize: time" << endl;

	/* set time integration controller */
	SetController();
	if (verbose) cout << "    FEManagerT::Initialize: controller" << endl;

	/* resolve node type here */
	SetNodeManager();
	if (verbose) cout << "    FEManagerT::Initialize: nodal data" << endl;
	
	/* construct element groups */
	SetElementGroups();
	if (verbose) cout << "    FEManagerT::Initialize: element groups" << endl;

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
			int error = eNoError;		

			/* initialize the current time step */
			if (error == eNoError) 
				error = InitStep();
			
			/* solve the current time step */
			if (error == eNoError) 
				error = SolveStep();
			
			/* close the current time step */
			if (error == eNoError) 
				error = CloseStep();
			
			/* handle errors */
			switch (error)
			{
				case eNoError:
					/* nothing to do */
					break;

				case eGeneralFail:
				case eBadJacobianDet:
				{
					/* reset system configuration */
					error = ResetStep();
					
					if (error == eNoError)
						/* cut time increment */
						seq_OK = DecreaseLoadStep();
					else
						seq_OK = false;
					break;
				}
				default:
					cout << "FEManagerT::Solve: no recovery for error: " << Exception(error) << endl;
					seq_OK = false;
			}
		}
	}
}

/* signal that references to external data are stale - i.e. equ numbers */
void FEManagerT::Reinitialize(int group)
{
	/* reset equation structure */
	SetEquationSystem(group);		
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

#if 0
/* exception handling */
void FEManagerT::HandleException(int exception)
{
	/* state */
	fStatus = GlobalT::kException;

	switch (exception)
	{
		case eBadJacobianDet:
		{
			cout << "\n FEManagerT::HandleException: detected bad jacobian determinant" << endl;
		
			if (fAnalysisCode == GlobalT::kLinExpDynamic   ||
			    fAnalysisCode == GlobalT::kNLExpDynamic    ||
			    fAnalysisCode == GlobalT::kPML)
			{
				cout << " FEManagerT::HandleException: no adaptive step for analysis code "
				     << fAnalysisCode << endl;
				throw eGeneralFail;
			}
			else
			{
				ResetStep();
				DecreaseLoadStep();
			}
			break;
		}
		default:
			cout << "\n FEManagerT::HandleException: unrecoverable exception: "
			     << Exception(exception) << '\n';
			cout <<   "      time: " << Time() << '\n';
			cout <<   "      step: " << StepNumber() << endl;
			throw eGeneralFail;
	}
}
#endif

void FEManagerT::WriteExceptionCodes(ostream& out) const
{
	out << "\nE x c e p t i o n   c o d e s :\n\n";
	
	for (int i = 0; i < eNumExceptions; i++)
	{
		out << setw(kIntWidth) << i << " : ";
		out << eExceptionStrings[i] << '\n';
		out << endl;
	}	
}

const char* FEManagerT::Exception(int code) const
{
	if (code >= 0 && code <= 8)
		return eExceptionStrings[code];
	else
		return eExceptionStrings[9];
}

/* load control functions */
bool FEManagerT::DecreaseLoadStep(void) { return fTimeManager->DecreaseLoadStep(); }
bool FEManagerT::IncreaseLoadStep(void) { return fTimeManager->IncreaseLoadStep(); }

int FEManagerT::ResetStep(void)
{
	int error = eNoError;
	try{
	/* state */
	fStatus = GlobalT::kResetStep;	

	/* time */
	fTimeManager->ResetStep();

	/* nodes - ALL groups */
	for (int i = 0; i < NumGroups(); i++)
		fNodeManager->ResetStep(i);
	
	/* elements */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->ResetStep();
		
	/* solver - ALL groups */
	for (int i = 0; i < NumGroups(); i++)
		fSolvers[i]->ResetStep();
	}
	
	catch (int exc) {
		cout << "\n FEManagerT::ResetStep: caught exception: " 
		     << Exception(exc) << endl;
		return exc;
	}
	
	/* OK */
	return eNoError;
}

const double& FEManagerT::Time(void) const { return fTimeManager->Time(); }
const double& FEManagerT::TimeStep(void) const { return fTimeManager->TimeStep(); }
const int& FEManagerT::StepNumber(void) const { return fTimeManager->StepNumber() ; }
const int& FEManagerT::NumberOfSteps(void) const { return fTimeManager->NumberOfSteps(); }
int FEManagerT::SequenceNumber(void) const { return fTimeManager->SequenceNumber(); }
int FEManagerT::NumSequences(void) const { return fTimeManager->NumSequences(); }
const int& FEManagerT::IterationNumber(int group) const 
{ 
	return fSolvers[group]->IterationNumber(); 
}

/* solution messaging */
void FEManagerT::FormLHS(int group) const
{
	/* state */
	SetStatus(GlobalT::kFormLHS);
	
	/* nodal contributions - from F(x) BC's */
	fNodeManager->FormLHS(group);

	/* element contributions */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		if (fElementGroups[i]->Group() == group)
			fElementGroups[i]->FormLHS();
}

void FEManagerT::FormRHS(int group) const
{
	/* state */
	SetStatus(GlobalT::kFormRHS);

	/* nodal force contribution - F(t) */
	//fController->FormNodalForce(fNodeManager);
	fNodeManager->FormRHS(group);

	/* element contribution */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		if (fElementGroups[i]->Group() == group)
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

int FEManagerT::InitStep(void)
{
	try {
	/* state */
	SetStatus(GlobalT::kInitStep);
	
	/* set the default value for the output time stamp */
	fIOManager->SetOutputTime(Time());

	/* nodes - ALL groups*/
	for (int i = 0; i < NumGroups(); i++)
		fNodeManager->InitStep(i);

	/* elements */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->InitStep();
	}
	
	catch (int exc) {
		cout << "\n FEManagerT::InitStep: caught exception: " 
		     << Exception(exc) << endl;
		return exc;
	}
	
	/* OK */
	return eNoError;
}

int FEManagerT::SolveStep(void)
{
	int error = eNoError;
	try {

		/* one solver after the next */
		for (int i = 0; error == eNoError && i < fSolvers.Length(); i++)
		{
			error = fSolvers[i]->Solve(); 
		}
	}

	catch (int exc) {
		cout << "\n FEManagerT::SolveStep: caught exception: " 
		     << Exception(exc) << endl;
	}
	
	/* done */
	return error;
}

int FEManagerT::CloseStep(void)
{
	try {
	/* state */
	SetStatus(GlobalT::kCloseStep);

	/* write output BEFORE closing nodes and elements */
	fTimeManager->CloseStep();

	/* nodes - ALL groups */
	for (int i = 0; i < NumGroups(); i++)
		fNodeManager->CloseStep(i);

	/* elements */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->CloseStep();
		
	/* write restart file */
	WriteRestart();
	}
	catch (int exc) {
		cout << "\n FEManagerT::CloseStep: caught exception: " << Exception(exc) << endl;
		return exc;
	}
	
	/* OK */
	return eNoError;
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
		if (fElementGroups[i]->Group() == group)
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

void FEManagerT::AssembleRHS(int group, const dArrayT& elRes,
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
void FEManagerT::WriteOutput(double time, IOBaseT::OutputModeT mode)
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
			fElementGroups[i]->WriteOutput(mode);
	}
	
	catch (int error) { 
		cout << "\n FEManagerT::WriteOutput: caught exception: " << Exception(error) << endl;
		throw error; 
	}
}

void FEManagerT::WriteOutput(int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
	fIOManager->WriteOutput(ID, n_values, e_values);
}

int FEManagerT::RegisterOutput(const OutputSetT& output_set) const
{
	/* check */
	if (!fIOManager) {
		cout << "\n FEManagerT::RegisterOutput: I/O manager not yet initialized" << endl;
		throw eGeneralFail;
	}

	/* limit registering output to initialization stage */
	if (fStatus != GlobalT::kInitialization) {
		cout << "\n FEManagerT::RegisterOutput: output sets can only be registered\n"
		     <<   "     during initialization" << endl;
		throw eGeneralFail;
	}

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
		else {
			cout << "\n FEManagerT::RegisterOutput: unrecognized output set mode: "
			     << output_set.Mode() << endl;
			throw eGeneralFail;
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
	if (!fIOManager) {
		cout << "\n FEManagerT::OutputSet: I/O manager not yet initialized" << endl;
		throw eGeneralFail;
	}

	return fIOManager->OutputSet(ID);
}

/* (temporarily) direct output away from main out */
void FEManagerT::DivertOutput(const StringT& outfile)
{
	/* check */
	if (!fIOManager) {
		cout << "\n FEManagerT::DivertOutput: I/O manager not yet initialized" << endl;
		throw eGeneralFail;
	}

	fIOManager->DivertOutput(outfile);
}

void FEManagerT::RestoreOutput(void)
{
	/* check */
	if (!fIOManager) {
		cout << "\n FEManagerT::RestoreOutput: I/O manager not yet initialized" << endl;
		throw eGeneralFail;
	}

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

void FEManagerT::IncomingNodes(iArrayT& nodes_in ) const {  nodes_in.Free(); }
void FEManagerT::OutgoingNodes(iArrayT& nodes_out) const { nodes_out.Free(); }

void FEManagerT::RecvExternalData(dArray2DT& external_data)
{
#pragma unused(external_data)
	cout << "\n FEManagerT::RecvExternalData: invalid request for external data" << endl;
	throw eGeneralFail;
}

void FEManagerT::SendExternalData(const dArray2DT& all_out_data)
{
#pragma unused(all_out_data)
	cout << "\n FEManagerT::RecvExternalData: invalid send of external data" << endl;
	throw eGeneralFail;
}

void FEManagerT::SendRecvExternalData(const iArray2DT& all_out_data, iArray2DT& external_data)
{
#pragma unused(all_out_data)
#pragma unused(external_data)
	cout << "\n FEManagerT::RecvExternalData: invalid exchange of external data" << endl;
	throw eGeneralFail;
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
	if (!e_controller) throw eGeneralFail;
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
	if (!n_controller) throw eGeneralFail;
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
	const iArrayT* node_map = NodeMap();

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
			WriteOutput(Time(), IOBaseT::kAtInc);			
		}
		else
			/* inherited */
			return iConsoleObjectT::iDoCommand(command, line);
	}
	
	catch (int error)
	{
		cout << "caught exception: " << Exception(error) << endl;
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
	fMainIn >> fVersion;		
	if (strcmp(fVersion, kCurrentVersion) != 0)
	{
//TEMP - until input is parsed, we will only have spotty support for backward
//       compatibility of input files
		cout << "\n !!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		cout << " FEManagerT::CheckInputFile: input file version is not current. See\n"
		     <<   "     VERSION_NOTES for description of changes:\n";
		cout << "     file version: " << fVersion << '\n';
		cout << "  current version: " << kCurrentVersion << '\n';
		cout << " WARNING: backward compatibility is not completely supported\n";
		cout << " !!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	}	
	else
		cout    << "\n Input file version: " << fVersion << '\n';

	fMainOut << "\n Input file version: " << fVersion << '\n';
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
	if (!fNodeManager) throw eOutOfMemory;	
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
	if (num_groups < 1) throw eBadInputValue;
	fElementGroups.Allocate(num_groups);
	fElementGroups.EchoElementData(fMainIn, fMainOut, *this);
		
	/* set console */
	for (int i = 0; i < fElementGroups.Length(); i++)
		iAddSub(*(fElementGroups[i]));
}

/* set the correct fSolutionDriver type */
void FEManagerT::SetSolver(void)
{
	/* equation info */
	int num_groups = fSolvers.Length();
	fGlobalEquationStart.Dimension(num_groups);
	fActiveEquationStart.Dimension(num_groups);
	fGlobalNumEquations.Dimension(num_groups);

	/* no predefined solvers */ 
	if (fAnalysisCode == GlobalT::kMultiField)
	{
		for (int i = 0; i < fSolvers.Length(); i++)
		{
			int index = -1;
			int type = -1;
			fMainIn >> index >> type;
			index--;
			if (fSolvers[index] != NULL) {
				cout << "\n FEManagerT::SetSolver: solver at index "
				     << index+1 << " is already set" << endl;
				throw eBadInputValue;
			}
	
			/* construct solver */
			fSolvers[index] = New_Solver(type, i);
		}
	}
	else /* support for legacy analysis codes */
	{
		/* should have just one solver */
		if (fSolvers.Length() != 1) throw eSizeMismatch;
	
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
				cout << "\n FEManagerT::SetSolver: unknown analysis type: " << fAnalysisCode << endl;
				throw eBadInputValue;
		}
	}
	
	/* initialize */
	for (int i = 0; i < fSolvers.Length(); i++)
	{
		/* reset equation structure */
		SetEquationSystem(i);

		/* console hierarchy */
		iAddSub(*(fSolvers[i]));	
	}
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
	if (fReadRestart  != 0 && fReadRestart  != 1) throw eBadInputValue;
	if (fWriteRestart < 0) throw eBadInputValue;
	if (fPrintInput   != 0 && fPrintInput   != 1) throw eBadInputValue;

	/* initialize the model manager */
	if (!fModelManager->Initialize(format, database, 
		init == kFull || init == kAllButSolver)) /* conditions under which to scan model */
	{
		cout << "\n FEManagerT::ReadParameters: error initializing model manager" << endl;
		throw eBadInputValue;
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
				throw eBadInputValue;
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
				throw eBadInputValue;
		}
		
		if (!controller) throw eGeneralFail;
		fControllers[0] = controller;
	}
}

/* construct output */
void FEManagerT::SetOutput(void)
{
	StringT file_name(fMainIn.filename());
	fIOManager = new IOManager(fMainOut, kProgramName, kCurrentVersion, fTitle,
		file_name, fOutputFormat);	
	if (!fIOManager) throw eOutOfMemory;
	
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
			throw eBadInputValue;
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
			if (fElementGroups[i]->Group() == group)
				fElementGroups[i]->ConnectsU(connects_1, connects_2);		
	
		/* renumber equations */
		fNodeManager->RenumberEquations(group, connects_1, connects_2);
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
		if (fElementGroups[i]->Group() == group)
			fElementGroups[i]->Equations(eq_1, eq_2);

	/* send lists to solver */
	for (int j = 0; j < eq_1.Length(); j++)
		fSolvers[group]->ReceiveEqns(*(eq_1[j]));

	for (int k = 0; k < eq_2.Length(); k++)
		fSolvers[group]->ReceiveEqns(*(eq_2[k]));
}

SolverT* FEManagerT::New_Solver(int code, int group)
{
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
			if (fields.Length() < 1) throw eGeneralFail;
			
			const nControllerT& controller = fields[0]->nController();
			solver = new NOXSolverT(*this, group, controller.OrderOfUnknown());
			break;
#else
			cout << "\n FEManagerT::New_Solver: NOX not installed: " << SolverT::kNOX << endl;
			throw eGeneralFail;
#endif	
		}		
		default:			
			cout << "\n FEManagerT::New_Solver: unknown nonlinear solver code: ";
			cout << code << endl;
			throw eBadInputValue;
	}

	/* fail */
	if (!solver) {
		cout << "\n FEManagerT::New_Solver: failed" << endl;
		throw eGeneralFail;	
	}

	return solver;
}