/* $Id: FEManagerT.cpp,v 1.1.1.1 2001-01-29 08:20:21 paklein Exp $ */
/* created: paklein (05/22/1996)                                          */

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
#include "ElementBaseT.h"
#include "IOManager.h"

/* nodes */
#include "NodeManagerT.h"
#include "LinStaticNodes.h"
#include "FDNodeManager.h"
#include "DynNodeManager.h"
#include "FDDynNodeManagerT.h"
#include "VariFDNodesT.h"
#include "VariFDDynNodesT.h"
#include "XDOF_FDNodesT.h"
#include "DuNodeManager.h"

/* controllers */
#include "StaticController.h"
#include "Trapezoid.h"
#include "LinearHHTalpha.h"
#include "NLHHTalpha.h"
#include "ExplicitCDController.h"

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

/* File/Version Control */
const char* kCurrentVersion = "v3.01"; //version marker
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
/* 9 */ "unknown"};	

/* constructor */
FEManagerT::FEManagerT(ifstreamT& input, ofstreamT& output):
	fMainIn(input),
	fMainOut(output),
	fStatus(GlobalT::kConstruction),
	fTimeManager(NULL),
	fNodeManager(NULL),
	fElementGroups(*this),
	fSolutionDriver(NULL),
	fController(NULL),
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
	iAddVariable("model_file", *((const StringT*) &fModelFile));
	iAddVariable("restart_inc", fWriteRestart);
	
	/* console commands */
	iAddCommand("ReadRestart");
	iAddCommand("WriteRestart");
	iAddCommand("WriteOutput");
}

/* destructor */
FEManagerT::~FEManagerT(void)
{
	fStatus = GlobalT::kDestruction;
	delete fTimeManager;
	delete fNodeManager;
	delete fSolutionDriver;
	delete fController;
	delete fIOManager;
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
	
	/* main parameters */
	ReadParameters();
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

	/* resolve node and element types here */
	SetNodeManager();
	if (verbose) cout << "    FEManagerT::Initialize: nodal data" << endl;
	
	/* construct element groups */
	SetElementGroups();
	if (verbose) cout << "    FEManagerT::Initialize: element groups" << endl;

	/* set I/O manager */
	SetIO();
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

		/* step through sequence */
		fSolutionDriver->Run();
		
		/* write restart file */
		//WriteRestart();
	}
}

/* signal that references to external data are stale - i.e. equ numbers */
void FEManagerT::Reinitialize(void)
{
	/* node */
	fNodeManager->Reinitialize();
	
	/* elements */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->Reinitialize();

	/* reset equation structure */
	SetEquationSystem();		
}

ifstreamT& FEManagerT::OpenExternal(ifstreamT& in,  ifstreamT& in2, ostream& out,
	bool verbose, const char* fail) const
{
	/* check for external file */
	char nextchar = in.next_char();
	if (isdigit(nextchar))
		return in;
	else
	{
		/* open external file */
		StringT file;
		in >> file;
		if (verbose) out << " external file: " << file << '\n';
			
		/* open stream */
		file.ToNativePathName();
		in2.open(file);
		if (!in2.is_open())
		{
			if (verbose && fail) cout << "\n " << fail << ": " << file << endl;
			throw eBadInputValue;
		}

		/* set comments */
		if (in.skip_comments()) in2.set_marker(in.comment_marker());

		return in2;
	}
}

/* manager messaging */
LoadTime* FEManagerT::GetLTfPtr(int num) const
{
	return fTimeManager->GetLTf(num);
}

double FEManagerT::LoadFactor(int nLTf) const
{
	return fTimeManager->LoadFactor(nLTf);
}

int FEManagerT::NumberOfLTf(void) const
{
	return fTimeManager->NumberOfLTf();
}

GlobalT::AnalysisCodeT FEManagerT::Analysis(void) const { return fAnalysisCode; }
bool FEManagerT::PrintInput(void) const { return fPrintInput; }

void FEManagerT::WriteEquationNumbers(void) const
{
	fNodeManager->WriteEquationNumbers(fMainOut);
	fMainOut.flush();
}

GlobalT::SystemTypeT FEManagerT::GlobalSystemType(void) const
{
	GlobalT::SystemTypeT type = fNodeManager->TangentType();
	for (int i = 0 ; i < fElementGroups.Length(); i++)
	{
		GlobalT::SystemTypeT e_type = fElementGroups[i]->TangentType();

		/* using precedence */
		type = (e_type > type) ? e_type : type;
	}
	return type;
}

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
			    fAnalysisCode == GlobalT::kVarNodeNLExpDyn ||
			    fAnalysisCode == GlobalT::kNLExpDynKfield)
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

/* time sequence messaging */
bool FEManagerT::Step(void)
{
	/* more steps */
	return fTimeManager->Step();
}

void FEManagerT::ResetStep(void)
{
	/* state */
	fStatus = GlobalT::kResetStep;	

	/* time */
	fTimeManager->ResetStep();

	/* nodes */
	fNodeManager->ResetStep();
	
	/* elements */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->ResetStep();
		
	/* solver */
	fSolutionDriver->ResetStep();
}

const double& FEManagerT::Time(void) const { return fTimeManager->Time(); }
const double& FEManagerT::TimeStep(void) const { return fTimeManager->TimeStep(); }
const int& FEManagerT::StepNumber(void) const { return fTimeManager->StepNumber() ; }
const int& FEManagerT::NumberOfSteps(void) const { return fTimeManager->NumberOfSteps(); }
int FEManagerT::SequenceNumber(void) const { return fTimeManager->SequenceNumber(); }
int FEManagerT::NumSequences(void) const { return fTimeManager->NumSequences(); }
const int& FEManagerT::IterationNumber(void) const { return fSolutionDriver->IterationNumber(); }

void FEManagerT::SetLocalEqnos(const iArray2DT& nodes,
	iArray2DT& eqnos) const
{
	fNodeManager->SetLocalEqnos(nodes, eqnos);
}

void FEManagerT::RegisterLocal(LocalArrayT& array) const
{
	fNodeManager->RegisterLocal(array);
}

/* solution messaging */
void FEManagerT::FormLHS(void) const
{
	/* state */
	SetStatus(GlobalT::kFormLHS);
	
	/* nodal contributions - from F(x) BC's */
	fNodeManager->FormLHS();

	/* element contributions */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->FormLHS();
}

void FEManagerT::FormRHS(void) const
{
	/* state */
	SetStatus(GlobalT::kFormRHS);

	/* nodal force contribution - F(t) */
	fController->FormNodalForce(fNodeManager);

	/* element contribution */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->FormRHS();
		
	/* output system info (debugging) */
	if (fSolutionDriver->Check() == GlobalMatrixT::kPrintRHS)
		WriteSystemConfig(fMainOut);
}

/* collect the internal force on the specified node */
void FEManagerT::InternalForceOnNode(int node, dArrayT& force) const
{
	/* initialize */
	force = 0.0;

	/* element contribution */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->AddNodalForce(node, force);
}

void FEManagerT::InitStep(void) const
{
	/* state */
	SetStatus(GlobalT::kInitStep);

	/* nodes */
	fNodeManager->InitStep();

	/* elements */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->InitStep();
}

void FEManagerT::CloseStep(void) const
{
	/* state */
	SetStatus(GlobalT::kCloseStep);

	/* nodes */
	fNodeManager->CloseStep();

	/* elements */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->CloseStep();
		
	/* write output */
	fTimeManager->CloseStep(); //TEMP will move control here later
	
	/* write restart file */
	WriteRestart();
}

void FEManagerT::Update(const dArrayT& update)
{
	fNodeManager->Update(update);
}

void FEManagerT::ActiveDisplacements(dArrayT& activedisp) const
{
	fNodeManager->ActiveDisplacements(activedisp);
}

GlobalT::RelaxCodeT FEManagerT::RelaxSystem(void) const
{
	GlobalT::RelaxCodeT relax = GlobalT::kNoRelax;
	
	/* check node manager */
	relax = GlobalT::MaxPrecedence(relax, fNodeManager->RelaxSystem());
		
	/* check element groups - must touch all of them to reset */
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		relax = GlobalT::MaxPrecedence(relax, fElementGroups[i]->RelaxSystem());

	return relax;
}

/* global equation functions */
void FEManagerT::AssembleLHS(const ElementMatrixT& elMat,
	const iArrayT& eqnos) const
{
	fSolutionDriver->AssembleLHS(elMat, eqnos);
}

void FEManagerT::OverWriteLHS(const ElementMatrixT& elMat,
	const iArrayT& eqnos) const
{
	fSolutionDriver->OverWriteLHS(elMat, eqnos);
}

void FEManagerT::DisassembleLHS(dMatrixT& elMat, const iArrayT& eqnos) const
{
	fSolutionDriver->DisassembleLHS(elMat, eqnos);
}

void FEManagerT::AssembleRHS(const dArrayT& elRes,
	const iArrayT& eqnos) const
{
	fSolutionDriver->AssembleRHS(elRes, eqnos);
}

void FEManagerT::OverWriteRHS(const dArrayT& elRes, const iArrayT& eqnos) const
{
	fSolutionDriver->OverWriteRHS(elRes, eqnos);
}

void FEManagerT::DisassembleRHS(dArrayT& elRes, const iArrayT& eqnos) const
{
	fSolutionDriver->DisassembleRHS(elRes, eqnos);
}

/* writing results */
void FEManagerT::WriteOutput(IOBaseT::OutputModeT mode)
{
	try
	{
		/* state */
		SetStatus(GlobalT::kWriteOutput);

		/* nodes */
		fNodeManager->WriteOutput(mode);
		
		/* elements */
		for (int i = 0; i < fElementGroups.Length(); i++)
			fElementGroups[i]->WriteOutput(mode);

//TEMP - skip
#if 0	
		cout << "\n FEManagerT::WriteOutput: computing total linear momentum " << endl;
		dArrayT momentum((fNodeManager->Displacements()).MinorDim());
		momentum = 0.0;
		for (int j = 0 ; j < fElementGroups.Length(); j++)
		{
			/* get element group */
			ElementBaseT* pbase = fElementGroups[j];
			ElasticT* pelastic = dynamic_cast<ElasticT*>(pbase);
			if (!pelastic) throw eGeneralFail;
	
			/* accumulate from elements */
			pelastic->AddLinearMomentum(momentum);	
		}
		
		fMainOut << "\n momentum = ";
		momentum.Write(fMainOut);
		fMainOut << '\n';
#endif
	}
	
	catch (int error) { HandleException(error); }
}

void FEManagerT::WriteOutput(int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
	fIOManager->WriteOutput(Time(), ID, n_values, e_values);
}

int FEManagerT::RegisterOutput(const OutputSetT& output_set)
{
	return fIOManager->AddElementSet(output_set);
}

void FEManagerT::WriteGeometryFile(const StringT& file_name,
	IOBaseT::FileTypeT output_format) const
{
	fIOManager->WriteGeometryFile(file_name, output_format);
}

const OutputSetT& FEManagerT::OutputSet(int ID) const
{
	return fIOManager->OutputSet(ID);
}

/* (temporarily) direct output away from main out */
void FEManagerT::DivertOutput(const StringT& outfile)
{
	fIOManager->DivertOutput(outfile);
}

void FEManagerT::RestoreOutput(void)
{
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

int FEManagerT::GlobalEquationNumber(int nodenum, int dofnum) const
{
	return fNodeManager->GlobalEquationNumber(nodenum, dofnum);
}

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
GlobalT::EquationNumberScopeT FEManagerT::EquationNumberScope(void) const
{
	return fSolutionDriver->EquationNumberScope();
}

int FEManagerT::GetGlobalEquationStart(void) const
{
	/* no other equations */
	return 1;
}

int FEManagerT::GetGlobalNumEquations(void) const
{
	/* no other equations */
	return fNodeManager->NumEquations();
}

/* access to controllers */
eControllerT* FEManagerT::eController(void) const
{
	/* cast to eControllerT */
#ifdef __NO_RTTI__
	eControllerT* e_controller = (eControllerT*) fController;
		//NOTE: cast should be safe for all cases
#else
	eControllerT* e_controller = dynamic_cast<eControllerT*>(fController);
	if (!e_controller) throw eGeneralFail;
#endif

	return e_controller;
}

nControllerT* FEManagerT::nController(void) const
{
	/* cast to eControllerT */
#ifdef __NO_RTTI__
	nControllerT* n_controller = (nControllerT*) fController;
		//NOTE: cast should be safe for all cases
#else
	nControllerT* n_controller = dynamic_cast<nControllerT*>(fController);
	if (!n_controller) throw eGeneralFail;
#endif

	return n_controller;
}

void FEManagerT::SetTimeStep(double dt) const
{
	if (!fController) throw eGeneralFail;

	fController->SetTimeStep(dt);
}

// TEMP - adding nodes to the Vari_ nodemanagers
// returns the global node number for the new node
int FEManagerT::AddNode(const ArrayT<LocalArrayT::TypeT>& types,
	const dArray2DT& values)
{
	//could use RTTI here
	switch (fAnalysisCode)
	{
		case GlobalT::kVarNodeNLStatic:
		{
			VariFDNodesT* pnodes = (VariFDNodesT*) fNodeManager;
			return  pnodes->AddNode(types,values);
			break;
		}
		case GlobalT::kVarNodeNLExpDyn:
		{
			VariFDDynNodesT* pnodes = (VariFDDynNodesT*) fNodeManager;
			return pnodes->AddNode(types,values);
			break;
		}
		default:
		
			cout << "\n FEManagerT::AddNode: incorrect analysis code: ";
			cout << fAnalysisCode << endl;
			throw eGeneralFail;
	}

	/* catch all */
	return 0;
}

/* returns 1 of ALL element groups have interpolant DOF's */
int FEManagerT::InterpolantDOFs(void) const
{
	return fElementGroups.InterpolantDOFs();
}

/* debugging */
void FEManagerT::WriteSystemConfig(ostream& out) const
{
	int old_precision = out.precision();
	out.precision(DBL_DIG);
	
	/* node map */
	const iArrayT* node_map = NodeMap();

	/* nodal data */
	const dArray2DT& coords = fNodeManager->InitialCoordinates();
	const dArray2DT&   disp = fNodeManager->Displacements();

	/* dimensions */
	int nnd = coords.MajorDim();
	int nsd = coords.MinorDim();
	int ndf = disp.MinorDim();

	/* force vector */
	const dArrayT& RHS = fSolutionDriver->RHS();

	/* header */
	out << "\n time = " << Time() << '\n';
	int d_width = OutputWidth(out, coords.Pointer());
	out << setw(kIntWidth) << "node"
	    << setw(kIntWidth) << "mapped";
	for (int i0 = 0; i0 < ndf; i0++)
		out << setw(kIntWidth - 2) << "eq[" << i0 + 1 << "]";
	for (int i1 = 0; i1 < nsd; i1++)
		out << setw(d_width - 2) << "x[" << i1 + 1 << "]";
	for (int i2 = 0; i2 < ndf; i2++)
		out << setw(d_width - 2) << "d[" << i2 + 1 << "]";
	for (int i3 = 0; i3 < ndf; i3++)
		out << setw(d_width - 2) << "f[" << i3 + 1 << "]";
	out << '\n';

	/* loop over nodes */
	int shift = ActiveEquationStart();
	int num_eq = RHS.Length();
	iArrayT eq(ndf);
	for (int i = 0; i < nnd; i++)
	{
		/* local node number */
		out << setw(kIntWidth) << i+1;
		
		/* mapped node number */
		out << setw(kIntWidth) << ((node_map != NULL) ? (*node_map)[i] : i) + 1;
		
		/* (local) equation numbers */
		for (int i0 = 0; i0 < ndf; i0++)
		{
			eq[i0] = GlobalEquationNumber(i, i0);
			out << setw(kIntWidth) << eq[i0];
		}
			
		/* coordinates */
		for (int i1 = 0; i1 < nsd; i1++)
			out << setw(d_width) << coords(i, i1);

		/* displacement */
		for (int i2 = 0; i2 < ndf; i2++)
			out << setw(d_width) << disp(i, i2);

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
bool FEManagerT::iDoCommand(const StringT& command, StringT& line)
{
	try
	{
		if (command == "ReadRestart")
		{
			StringT file_name;
			if (ResolveArgument(line, file_name, NULL))
				ReadRestart(&file_name);
			else
			{
				cout << "could not resolve file name from \"" << line
				     << '\"'<< endl;
			}
		}
		else if (command == "WriteRestart")
		{
			StringT file_name;
			if (ResolveArgument(line, file_name, NULL))
				WriteRestart(&file_name);
			else
			{
				cout << "could not resolve file name from \"" << line
				     << '\"'<< endl;
			}
		}
		else if (command == "WriteOutput")
		{
			WriteOutput(IOBaseT::kAtInc);			
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
	StringT version;
	fMainIn >> version;
		
	if (strcmp(version, kCurrentVersion) != 0)
	{
		cout << "\n FEManagerT::CheckInputFile: input file version is not current. See\n"
		     <<   "     VERSION_NOTES for description of changes:\n";
		cout << "     file version: " << version << '\n';
		cout << "  current version: " << kCurrentVersion << endl;
		throw eBadInputValue;
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
	fMainOut << "    eq. " << GlobalT::kAugLagStatic    << ", nonlinear static with element level DOF's\n";   	
	fMainOut << " Input format. . . . . . . . . . . . . . . . . . = " << fInputFormat  << '\n';
	fMainOut << "    eq. " << IOBaseT::kTahoe         << ", standard ASCII\n";
	fMainOut << "    eq. " << IOBaseT::kTahoeII       << ", random access ASCII\n";
	fMainOut << "    eq. " << IOBaseT::kExodusII      << ", ExodusII\n";
	if (fInputFormat == IOBaseT::kTahoeII || fInputFormat == IOBaseT::kExodusII)
		fMainOut << " Geometry file . . . . . . . . . . . . . . . . . = " << fModelFile  << '\n';
	fMainOut << " Output format . . . . . . . . . . . . . . . . . = " << fOutputFormat << '\n';
	fMainOut << "    eq. " << IOBaseT::kTahoe         << ", standard ASCII\n";
	fMainOut << "    eq. " << IOBaseT::kTecPlot       << ", TecPlot\n";
	fMainOut << "    eq. " << IOBaseT::kEnSight       << ", Ensight 6 ASCII\n";
	fMainOut << "    eq. " << IOBaseT::kEnSightBinary << ", Ensight 6 binary\n";
	fMainOut << "    eq. " << IOBaseT::kExodusII      << ", ExodusII\n";
	fMainOut << " Read restart file code  . . . . . . . . . . . . = " << fReadRestart << '\n';
	fMainOut << "    eq. 0, do not read restart file\n";
	fMainOut << "    eq. 1, read restart file\n";
	if (fReadRestart)
		fMainOut << " Restart file. . . . . . . . . . . . . . . . . . = " << fRestartFile << '\n';
	fMainOut << " Restart file increment (at current step size) . = " << fWriteRestart << '\n';
	fMainOut << " Input data print code . . . . . . . . . . . . . = " << fPrintInput << '\n';
	fMainOut << "    eq. 0, non verbose echo of input data\n";
	fMainOut << "    eq. 1, echo all input data\n";
	fMainOut << endl;
}

/* set the correct fNodeManager type */
void FEManagerT::SetNodeManager(void)
{
	switch (fAnalysisCode)
	{
		case GlobalT::kLinStaticHeat:
		case GlobalT::kLinStatic:

			fNodeManager = new LinStaticNodes(*this);
			break;

		case GlobalT::kLinTransHeat:
		
			fNodeManager = new DuNodeManager(*this);
			break;

		case GlobalT::kNLStatic:
		case GlobalT::kDR:

			fNodeManager = new FDNodeManager(*this);
			break;

		case GlobalT::kLinDynamic:
		case GlobalT::kLinExpDynamic:

			fNodeManager = new DynNodeManager(*this);
			break;

		case GlobalT::kNLDynamic:
		case GlobalT::kNLExpDynamic:

			fNodeManager = new FDDynNodeManagerT(*this);
			break;

		case GlobalT::kVarNodeNLStatic:
		
			fNodeManager = new VariFDNodesT(*this);
			break;
	
		case GlobalT::kVarNodeNLExpDyn:
		
			fNodeManager = new VariFDDynNodesT(*this);
			break;

		case GlobalT::kAugLagStatic:
		
			fNodeManager = new XDOF_FDNodesT(*this);
			break;
	
		default:

			cout << "FEManagerT::SetNodeManager: unknown analysis type." << endl;
			throw eBadInputValue;
	}
	
	if (!fNodeManager) throw eOutOfMemory;
	
	fNodeManager->Initialize();			

	/* cast to eControllerT */
#ifdef __NO_RTTI__
	nControllerT* n_controller = (nControllerT*) fController;
		//NOTE: cast should be safe for all cases
#else
	nControllerT* n_controller = dynamic_cast<nControllerT*>(fController);
	if (!n_controller) throw eGeneralFail;
#endif
	
	/* set controller */
	fNodeManager->SetController(n_controller);
	
	/* add to console */
	iAddSub(*fNodeManager);	
}

	/* construct element groups */
void FEManagerT::SetElementGroups(void)
{
	/* cast to eControllerT */
#ifdef __NO_RTTI__
	eControllerT* e_controller = (eControllerT*) fController;
		//NOTE: cast should be safe for all cases
#else
	eControllerT* e_controller = dynamic_cast<eControllerT*>(fController);
	if (!e_controller) throw eGeneralFail;
#endif
	
	/* echo element data */
	int num_groups;
	fMainIn >> num_groups;
	if (num_groups < 1) throw eBadInputValue;
	fElementGroups.Allocate(num_groups);
	fElementGroups.EchoElementData(fMainIn, fMainOut,
		e_controller);
		
	/* set console */
	for (int i = 0; i < fElementGroups.Length(); i++)
		iAddSub(*(fElementGroups[i]));
}

/* set the correct fSolutionDriver type */
void FEManagerT::SetSolver(void)
{	
	switch (fAnalysisCode)
	{
		case GlobalT::kLinStatic:
		case GlobalT::kLinDynamic:
		case GlobalT::kLinExpDynamic:
		case GlobalT::kNLExpDynamic:
		case GlobalT::kVarNodeNLExpDyn:
		case GlobalT::kNLExpDynKfield:
		case GlobalT::kLinStaticHeat:
		case GlobalT::kLinTransHeat:

			fSolutionDriver = new LinearSolver(*this);
			break;

		case GlobalT::kDR:

			fSolutionDriver = new DRSolver(*this);
			break;

		case GlobalT::kNLStatic:
		case GlobalT::kNLDynamic:
		case GlobalT::kNLStaticKfield:
		case GlobalT::kVarNodeNLStatic:
		case GlobalT::kAugLagStatic:
		{
			int NL_solver_code;
			fMainIn >> NL_solver_code;
			
			/* construct nonlinear solver */
			switch (NL_solver_code)
			{
				case SolverT::kNewtonSolver:			
					fSolutionDriver = new NLSolver(*this);	
					break;

				case SolverT::kK0_NewtonSolver:				
					fSolutionDriver = new NLK0Solver(*this);
					break;

				case SolverT::kModNewtonSolver:				
					fSolutionDriver = new NLSolverX(*this);
					break;

				case SolverT::kExpCD_DRSolver:				
					fSolutionDriver = new ExpCD_DRSolver(*this);
					break;

				case SolverT::kNewtonSolver_LS:				
					fSolutionDriver = new NLSolver_LS(*this);
					break;

				case SolverT::kPCGSolver_LS:				
					fSolutionDriver = new PCGSolver_LS(*this);
					break;

				case SolverT::kiNewtonSolver_LS:				
					fSolutionDriver = new iNLSolver_LS(*this);
					break;
			
				default:			
					cout << "\n FEManagerT::SetSolver: unknown nonlinear solver code: ";
					cout << NL_solver_code << endl;
			}
			break;
		}

		default:

			cout << "\n FEManagerT::SetSolver: unknown analysis type: " << fAnalysisCode << endl;
			throw eBadInputValue;
	}

	if (!fSolutionDriver) throw eOutOfMemory;

	/* reset equation structure */
	SetEquationSystem();
	
	/* console */
	iAddSub(*fSolutionDriver);
}

void FEManagerT::ReadParameters(void)
{
	/* read */
	fMainIn >> fAnalysisCode;
	fMainIn >> fInputFormat;
	if (fInputFormat == IOBaseT::kTahoeII ||
	    fInputFormat == IOBaseT::kExodusII) fMainIn >> fModelFile;
	fMainIn >> fOutputFormat;
	fMainIn >> fReadRestart;
	if (fReadRestart == 1) fMainIn >> fRestartFile;
	fMainIn >> fWriteRestart;
	fMainIn >> fPrintInput;

	/* check */
	if (fInputFormat  != IOBaseT::kTahoe   &&
	    fInputFormat  != IOBaseT::kTahoeII &&
	    fInputFormat  != IOBaseT::kExodusII) throw eBadInputValue;
	if (fOutputFormat != IOBaseT::kTahoe  &&
	    fOutputFormat != IOBaseT::kTecPlot   &&
	    fOutputFormat != IOBaseT::kEnSight   &&
	    fOutputFormat != IOBaseT::kEnSightBinary &&
	    fOutputFormat != IOBaseT::kExodusII) throw eBadInputValue;
	if (fReadRestart  != 0 && fReadRestart  != 1) throw eBadInputValue;
	if (fWriteRestart < 0) throw eBadInputValue;
	if (fPrintInput   != 0 && fPrintInput   != 1) throw eBadInputValue;
}

/* set the execution controller and send to nodes and elements.
* This function constructs the proper drived class controller
* and passes it to the nodes and elements.  The controller is
* then cast to a controller to extract only the FEManagerT
* interface. */
void FEManagerT::SetController(void)
{
	fMainOut << "\n T i m e   I n t e g r a t i o n   C o n t r o l l e r:\n";

	switch (fAnalysisCode)
	{
		case GlobalT::kLinStatic:
		case GlobalT::kNLStatic:
		case GlobalT::kNLStaticKfield:
		case GlobalT::kVarNodeNLStatic:
		case GlobalT::kAugLagStatic:
		case GlobalT::kLinStaticHeat:
		{
			fController = new StaticController(fMainOut);
			break;
		}
		case GlobalT::kLinTransHeat:
		{
			fController = new Trapezoid(fMainOut);
			break;
		}
		case GlobalT::kLinDynamic:
		{
			fController = new LinearHHTalpha(*fTimeManager, fMainIn, fMainOut,
				kHHTalphaAuto_O2);
			break;
		}
		case GlobalT::kNLDynamic:
		{
			fController = new NLHHTalpha(*fTimeManager, fMainIn, fMainOut,
				kHHTalphaAuto_O2);
			break;
		}
		case GlobalT::kLinExpDynamic:
		case GlobalT::kNLExpDynamic:
		case GlobalT::kVarNodeNLExpDyn:
		case GlobalT::kNLExpDynKfield:
		{
			fController = new ExplicitCDController(fMainOut);
			break;
		}			
		default:
		
			cout << "\nFEManagerT::SetController: unknown controller type\n" << endl;
			throw eBadInputValue;
	}
	
	if (!fController) throw eOutOfMemory;
}

/* construct I/O */
void FEManagerT::SetIO(void)
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
		if (fWriteRestart == 1 && StepNumber() == fTimeManager->NumberOfSteps())
			write = true;
		else if (fWriteRestart> 1 && ++counter == fWriteRestart)
			write = true;

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
void FEManagerT::SetEquationSystem(void)
{
	/* check */
	if (!fSolutionDriver)
	{
		cout << "\n FEManagerT::SetEquationSystem: invalid solution manager" << endl;
		throw eGeneralFail;
	}

	/* equation number scope */
	GlobalT::EquationNumberScopeT equation_scope = fSolutionDriver->EquationNumberScope();

	/* assign (local) equation numbers */
	fNodeManager->SetEquationNumbers();
	fGlobalEquationStart = GetGlobalEquationStart();
	fActiveEquationStart = (equation_scope == GlobalT::kGlobal) ? fGlobalEquationStart : 1;
	fGlobalNumEquations  = GetGlobalNumEquations();

	/* renumber locally */
	if (fSolutionDriver->RenumberEquations())
	{
		/* lists of connectivities */
		AutoArrayT<const iArray2DT*> connects_1;
		AutoArrayT<const RaggedArray2DT<int>*> connects_2;
	
		/* collect nodally generated DOF's */
		fNodeManager->ConnectsU(connects_1, connects_2);
	
		/* collect element groups */
		for (int i = 0 ; i < fElementGroups.Length(); i++)
			fElementGroups[i]->ConnectsU(connects_1, connects_2);		
	
		/* renumber equations */
		fNodeManager->RenumberEquations(connects_1, connects_2);
	}

	/* set equation number scope */
	fNodeManager->SetEquationNumberScope(equation_scope);
	
	/* collect interaction equations and send to solver */
	SendEqnsToSolver();
	
	/* final step in solver configuration */
	fSolutionDriver->Initialize(fGlobalNumEquations, fNodeManager->NumEquations(),
		fActiveEquationStart);
}

void FEManagerT::SendEqnsToSolver(void) const
{
	/* dynamic arrays */
	AutoArrayT<const iArray2DT*> eq_1;
	AutoArrayT<const RaggedArray2DT<int>*> eq_2;
	
	/* collect equation sets */
	fNodeManager->Equations(eq_1, eq_2);
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->Equations(eq_1, eq_2);

	/* send lists to solver */
	for (int j = 0; j < eq_1.Length(); j++)
		fSolutionDriver->ReceiveEqns(*(eq_1[j]));

	for (int k = 0; k < eq_2.Length(); k++)
		fSolutionDriver->ReceiveEqns(*(eq_2[k]));
}
