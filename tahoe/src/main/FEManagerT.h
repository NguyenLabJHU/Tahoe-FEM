/* $Id: FEManagerT.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (05/22/1996)                                          */

#ifndef _FE_MANAGER_H_
#define _FE_MANAGER_H_

/* program parameters */
#include "GlobalT.h"

/* base class */
#include "iConsoleObjectT.h"

/* direct members */
#include "StringT.h"
#include "ElementListT.h"
#include "LocalArrayT.h"
#include "IOBaseT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class ofstreamT;
class TimeManagerT;
class NodeManagerT;
class ControllerT;
class nControllerT;
class eControllerT;
class LoadTime;
class SolverT;
class dMatrixT;
class LocalArrayT;
class iArrayT;
class dArrayT;
class iArray2DT;
class ElementMatrixT;
class dArray2DT;
class iAutoArrayT;
class IOManager;
class OutputSetT;

class FEManagerT: public iConsoleObjectT
{
public:

	/* degree of initialization */
	enum InitCodeT {kFull = 0,
	      kParametersOnly = 1,
	        kAllButSolver = 2};

	/* constructor */
	FEManagerT(ifstreamT& input, ofstreamT& output);

	/* destructor */
	virtual ~FEManagerT(void);

	/* initialize members */
	void Initialize(InitCodeT init = kFull);
	
	/* solve all the time sequences */
	void Solve(void);

	/* signal that references to external data are stale - ie eq numbers */
	void Reinitialize(void);
	
	/* accessors */
	ifstreamT& Input(void) const;
	ofstreamT& Output(void) const;
	GlobalT::AnalysisCodeT Analysis(void) const;
	bool PrintInput(void) const;
	GlobalT::SystemTypeT GlobalSystemType(void) const;
	const GlobalT::StateT& RunState(void) const;

	/* stream utils */
	ifstreamT& OpenExternal(ifstreamT& in,  ifstreamT& in2, ostream& out,
		bool echo_name, const char* fail) const;

	/* functions of time */
	LoadTime* GetLTfPtr(int num) const;
	double LoadFactor(int nLTf) const;
	int NumberOfLTf(void) const;

	/* equation system */
	void WriteEquationNumbers(void) const;
	GlobalT::EquationNumberScopeT EquationNumberScope(void) const;
	int GlobalEquationNumber(int nodenum, int dofnum) const;
	int GlobalEquationStart(void) const;
	int ActiveEquationStart(void) const;
	int GlobalNumEquations(void) const;

	/* exception handling */
	virtual void HandleException(int exception);
	void WriteExceptionCodes(ostream& out) const;
	const char* Exception(int code) const;

	/* load control functions (returns true if successful) */
	bool DecreaseLoadStep(void);
	bool IncreaseLoadStep(void);
	
	/* time sequence messaging */
	virtual bool Step(void);
	void ResetStep(void);
	
	/* solution accessors */
	const double& Time(void) const;
	const double& TimeStep(void) const;
	const int& StepNumber(void) const;
	const int& NumberOfSteps(void) const;
	void SetTimeStep(double dt) const;
	int SequenceNumber(void) const;
	int NumSequences(void) const;
	const int& IterationNumber(void) const;

	/* I/O info */
	IOBaseT::FileTypeT InputFormat(void) const;
	IOBaseT::FileTypeT OutputFormat(void) const;
	const StringT& ModelFile(void) const;

	/* local reordering */
	void SetLocalEqnos(const iArray2DT& nodes, iArray2DT& eqnos) const;
	void RegisterLocal(LocalArrayT& array) const;
	
	/* solution messaging */
	void FormLHS(void) const;
	void FormRHS(void) const;
	
	/* collect the internal force on the specified node */
	void InternalForceOnNode(int node, dArrayT& force) const;

	/* first/last functions called during a time increment */
	virtual void InitStep(void) const;
	virtual void CloseStep(void) const;

	/* solution update */
	virtual void Update(const dArrayT& update);

	/* intermediate updates at the current time step */	
	void ActiveDisplacements(dArrayT& activedisp) const;

	/* system relaxation */
	virtual GlobalT::RelaxCodeT RelaxSystem(void) const;

	/* assembling the global equation system */
	void AssembleLHS(const ElementMatrixT& elMat, const iArrayT& eqnos) const;
	void OverWriteLHS(const ElementMatrixT& elMat, const iArrayT& eqnos) const;
	void DisassembleLHS(dMatrixT& elMat, const iArrayT& eqnos) const;

	void AssembleRHS(const dArrayT& elRes, const iArrayT& eqnos) const;
	void OverWriteRHS(const dArrayT& elRes, const iArrayT& eqnos) const;
	void DisassembleRHS(dArrayT& elRes, const iArrayT& eqnos) const;
	
	/* writing results */
	void WriteOutput(IOBaseT::OutputModeT mode);
	virtual int RegisterOutput(const OutputSetT& output_set);
	virtual void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values);
	IOManager* OutputManager(void) const;
	void WriteGeometryFile(const StringT& file_name, IOBaseT::FileTypeT output_format) const;
	const OutputSetT& OutputSet(int ID) const;

	/* (temporarily) direct output away from main out */
	void DivertOutput(const StringT& outfile);
	void RestoreOutput(void);
	
	/* cross-linking - create your own trouble */
	NodeManagerT* NodeManager(void) const;
	ElementBaseT* ElementGroup(int groupnumber) const;
		// 1 <= groupnumber <= fNumElementGroups
		// returns NULL if out of range
	int ElementGroupNumber(const ElementBaseT* pgroup) const;
		// returns the element group number (0...) for pgroup,
		// or -1 if not found, result not valid until after
		// fElementGroups is fully constructed

	/* external nodes functions (parallel execution) */
	virtual void IncomingNodes(iArrayT& nodes_in) const;
	virtual void OutgoingNodes(iArrayT& nodes_out) const;
	virtual void SendExternalData(const dArray2DT& all_out_data);
	virtual void RecvExternalData(dArray2DT& external_data);
	virtual void SendRecvExternalData(const iArray2DT& all_out_data, iArray2DT& external_data);
	virtual void Wait(void);

	/* access to controllers */
	eControllerT* eController(void) const;
	nControllerT* nController(void) const;
	
	// TEMP - adding nodes to the Vari_ nodemanagers
	// returns the global node number for the new node
	int AddNode(const ArrayT<LocalArrayT::TypeT>& types, const dArray2DT& values);
		//NOTE: error to add node for non-Vari analysis codes

	/* returns 1 of ALL element groups have interpolant DOF's */
	int InterpolantDOFs(void) const;

	/* debugging */
	virtual void WriteSystemConfig(ostream& out) const;
	virtual const iArrayT* NodeMap(void) const { return NULL; }
	virtual const iArrayT* ElementMap(int blockID) const;

	/* basic MP support */
	virtual int Rank(void) const { return 0; }
	virtual int Size(void) const { return 1; }

	/* interactive */
	virtual bool iDoCommand(const StringT& command, StringT& line);

protected:

	/* "const" function that sets the status flag */
	void SetStatus(GlobalT::StateT status) const;

	/* look for input file key and check file version */
	void CheckInputFile(void);

	/* initialization functions */
	virtual void ReadParameters(void);
	void WriteParameters(void) const;
	void SetController(void);
	virtual void SetNodeManager(void);
	virtual void SetElementGroups(void);
	void SetSolver(void);
	virtual void SetIO(void);

	/* (re-)set system to initial conditions */
	virtual void InitialCondition(void);
	  	
	/* initialize/restart functions
	 *
	 * Initialize functions reset all kinematic data to the
	 * default initial system state.  The restart functions
	 * should read/write any data that overrides the default
	 * values */
	void ReadRestart(const StringT* file_name = NULL);
	void WriteRestart(const StringT* file_name = NULL) const;

	/* final steps in solver configuration
	 * (1) signal nodes to assign equation numbers
	 * (2) renumber if needed
	 * (3) set numbering scope
	 * (4) collect equations and send to solver
	 * (5) signal solver for final configuration */
	virtual void SetEquationSystem(void);
	virtual int GetGlobalEquationStart(void) const;
	virtual int GetGlobalNumEquations(void) const;
	
	/* collect element equations and send to solver */
	void SendEqnsToSolver(void) const;

private:

	/* no copies/assignment */
	FEManagerT(FEManagerT&);
	FEManagerT& operator=(FEManagerT&) const;

protected:

	/* I/O streams */
	ifstreamT&  fMainIn;
	StringT     fModelFile;
	ofstreamT& 	fMainOut;
	StringT		fTitle;
	StringT     fRestartFile;

	/* state */
	GlobalT::StateT fStatus;
	
	/* execution parameters */
	GlobalT::AnalysisCodeT fAnalysisCode;
	IOBaseT::FileTypeT  fInputFormat;
	IOBaseT::FileTypeT  fOutputFormat;
	bool fReadRestart;
	int  fWriteRestart;
	bool fPrintInput;
	
	/* the managers */
	TimeManagerT* fTimeManager;
	NodeManagerT* fNodeManager;
	ElementListT  fElementGroups;
	SolverT*      fSolutionDriver;
	ControllerT*  fController;
	IOManager*    fIOManager;
	
	/* restart file counter */
	int fRestartCount;
	int fGlobalEquationStart;
	int fActiveEquationStart;
	int fGlobalNumEquations;
};

/* inlines */
inline ifstreamT& FEManagerT::Input(void) const { return fMainIn;  }
inline ofstreamT& FEManagerT::Output(void) const { return fMainOut; }
inline const GlobalT::StateT& FEManagerT::RunState(void) const { return fStatus; }
inline IOBaseT::FileTypeT FEManagerT::InputFormat(void) const { return fInputFormat; }
inline IOBaseT::FileTypeT FEManagerT::OutputFormat(void) const { return fOutputFormat; }
inline const StringT& FEManagerT::ModelFile(void) const { return fModelFile; }

inline NodeManagerT* FEManagerT::NodeManager(void) const { return fNodeManager; }
inline IOManager* FEManagerT::OutputManager(void) const { return fIOManager; }
inline const iArrayT* FEManagerT::ElementMap(int blockID) const
{
#pragma unused(blockID)
	return NULL;
}

inline int FEManagerT::GlobalEquationStart(void) const { return fGlobalEquationStart; }
inline int FEManagerT::ActiveEquationStart(void) const { return fActiveEquationStart; }
inline int FEManagerT::GlobalNumEquations(void) const { return fGlobalNumEquations; }

#endif /* _FE_MANAGER_H_ */
