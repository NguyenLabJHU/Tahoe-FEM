/* $Id: FEManagerT.h,v 1.13.2.4 2002-04-26 02:24:21 paklein Exp $ */
/* created: paklein (05/22/1996) */

#ifndef _FE_MANAGER_H_
#define _FE_MANAGER_H_

/* program parameters */
#include "GlobalT.h"

/* base class */
#include "iConsoleObjectT.h"

/* direct members */
#include "StringT.h"
#include "ElementListT.h"
#include "IOBaseT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class ofstreamT;
class ModelManagerT;
class TimeManagerT;
class NodeManagerPrimitive; //TEMP - rename
class ControllerT;
class nControllerT;
class eControllerT;
class ScheduleT;
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
class FieldT;

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
	
	/** \name accessors */
	/*@{*/
	ifstreamT& Input(void) const;
	ofstreamT& Output(void) const;
	GlobalT::AnalysisCodeT Analysis(void) const;
	GlobalT::SystemTypeT GlobalSystemType(int group) const;
	const GlobalT::StateT& RunState(void) const;

	/** get schedule function */
	const ScheduleT* Schedule(int num) const;

	/** return the number of equation groups */
	int NumGroups(void) const { return fSolvers.Length(); };
	/*@}*/

	bool PrintInput(void) const;

//NOTE - NEED THESE?
//	int NumberOfLTf(void) const;
//	double LoadFactor(int nLTf) const;

	/** \name equation system */
	/*@{*/
	/** write the field equations to for the given group to the stream */
	void WriteEquationNumbers(int group) const;

	/** determine the numbering scope of the equations for the given group */
	GlobalT::EquationNumberScopeT EquationNumberScope(int group) const;

	/** return the global equation number */
	int EquationNumber(int field, int node, int dof) const;

	/** the global number of the first equation on this processor, regardless of
	 * the FEManagerT::EquationNumberScope for that group. */
	int GlobalEquationStart(int group) const;

	/** the first equation number owned by this processor */
	int ActiveEquationStart(int group) const;

	/** total number of equations in the specified group */
	int GlobalNumEquations(int group) const;
	/*@}*/

	/** \name exception handling */
	/*@{*/
	virtual void HandleException(int exception);
	void WriteExceptionCodes(ostream& out) const;
	const char* Exception(int code) const;
	/*@}*/

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
	const int& IterationNumber(int group) const;

	/* I/O info */
	const StringT& Version(void) const;
	IOBaseT::FileTypeT OutputFormat(void) const;
	ModelManagerT* ModelManager (void) const;
	const StringT& Title(void) const { return fTitle; };

	/* local reordering */
	void SetLocalEqnos(const iArray2DT& nodes, iArray2DT& eqnos) const;
	void RegisterLocal(LocalArrayT& array) const;
	
	/* solution messaging */
	void FormLHS(void) const;
	void FormRHS(void) const;
	
	/** collect the internal force on the specified node */
	void InternalForceOnNode(const FieldT& field, int node, dArrayT& force) const;

	/* first/last functions called during a time increment */
	virtual void InitStep(void) const;
	virtual void CloseStep(void) const;

	/** send update of the solution to the NodeManagerT */
	virtual void Update(const dArrayT& update);

	/** return the current values of the unknowns 
	 * \param order time derivative of the unknowns to collect. Must be
	 *        in range
	 * \param unknowns destination for the current values field values
	 *        for unprescribed degrees of freedom */
	virtual void GetUnknowns(int order, dArrayT& unknowns) const;

	/* system relaxation */
	virtual GlobalT::RelaxCodeT RelaxSystem(void) const;

	/** \name assembly methods 
	 * methods for assembling contributions to the global equations systems */
	/*@{*/
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& eqnos) const;
	
	
	
	
	
	
	
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& row_eqnos,
		const nArrayT<int>& col_eqnos) const;
	void OverWriteLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& eqnos) const;
	void DisassembleLHS(int group, dMatrixT& elMat, const nArrayT<int>& eqnos) const;
	void DisassembleLHSDiagonal(int group, dArrayT& diagonals, const nArrayT<int>& eqnos) const;

	void AssembleRHS(int group, const dArrayT& elRes, const nArrayT<int>& eqnos) const;
	void OverWriteRHS(int group, const dArrayT& elRes, const nArrayT<int>& eqnos) const;
	void DisassembleRHS(int group, dArrayT& elRes, const nArrayT<int>& eqnos) const;
	/*@}*/

	/** pointer to the I/O manager */
	IOManager* OutputManager(void) const;
	
	/** register an output set to write output data. See OutputSetT for more information.
	 * \return the ID for the output set. This value is needed to send data to the
	 *         correct destination with a subsequent call to FEManagerT::WriteOutput */
	virtual int RegisterOutput(const OutputSetT& output_set) const;

	/** return a reference to the output set with the given output ID
	 * \param ID ID of the output set that was returned when the set was
	 *        registered with FEManagerT::RegisterOutput */
	const OutputSetT& OutputSet(int ID) const;

	/** initiate the process of writing output from all output sets 
	 * \param time time label associated with the output data */
	virtual void WriteOutput(double time, IOBaseT::OutputModeT mode);
	
	/** write results for a single output set
	 * \param ID output set ID for the given data
	 * \param n_values nodal output values
	 * \param e_values element output values */
	virtual void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values);

	/** write a geometry file for the current model */
	void WriteGeometryFile(const StringT& file_name, IOBaseT::FileTypeT output_format) const;

	/* (temporarily) direct output away from main out */
	virtual void DivertOutput(const StringT& outfile);
	virtual void RestoreOutput(void);
	
	/* cross-linking - create your own trouble */
	NodeManagerPrimitive* NodeManager(void) const;
	ElementBaseT* ElementGroup(int groupnumber) const;
		// 1 <= groupnumber <= fNumElementGroups
		// returns NULL if out of range
	int ElementGroupNumber(const ElementBaseT* pgroup) const;
		// returns the element group number (0...) for pgroup,
		// or -1 if not found, result not valid until after
		// fElementGroups is fully constructed

	/** \name access to controllers */
	/*@{*/
	ControllerT* Controller(int index) { return fControllers[index]; };
	const ControllerT* Controller(int index) const { return fControllers[index]; };
	eControllerT* eController(int index) const;
	nControllerT* nController(int index) const;
	/*@}*/

	/* returns 1 of ALL element groups have interpolant DOF's */
	int InterpolantDOFs(void) const;

	/** debugging */
	virtual void WriteSystemConfig(ostream& out) const;
	virtual const iArrayT* NodeMap(void) const { return NULL; }
	virtual const iArrayT* ElementMap(const StringT& block_ID) const;

	/** \name basic MP info */
	/*@{*/
	virtual int Rank(void) const { return 0; }
	virtual int Size(void) const { return 1; }

	/* external nodes functions (parallel execution) */
	virtual void IncomingNodes(iArrayT& nodes_in) const;
	virtual void OutgoingNodes(iArrayT& nodes_out) const;
	virtual void SendExternalData(const dArray2DT& all_out_data);
	virtual void RecvExternalData(dArray2DT& external_data);
	virtual void SendRecvExternalData(const iArray2DT& all_out_data, iArray2DT& external_data);
	virtual void Wait(void);
	/*@}*/

	/** interactive */
	virtual bool iDoCommand(const CommandSpecT& command, StringT& line);

protected:

	/** "const" function that sets the status flag */
	void SetStatus(GlobalT::StateT status) const;

	/** look for input file key and check file version */
	void CheckInputFile(void);

	/** \name phases of FEManagerT::Initialize. */
	/*@{*/
	virtual void ReadParameters(InitCodeT init);
	void WriteParameters(void) const;
	void SetController(void);
	virtual void SetNodeManager(void);
	virtual void SetElementGroups(void);
	void SetSolver(void);
	virtual void SetOutput(void);
	/*@}*/

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

	/** \name disallowed */
	/*@{*/
	FEManagerT(FEManagerT&);
	FEManagerT& operator=(FEManagerT&) const;
	/*@}*/

	/** construct a solver of the specified type. */
	SolverT* New_Solver(int code) const;
	
protected:

	/** \name I/O streams */
	/*@{*/
	ifstreamT& fMainIn;
	ofstreamT& fMainOut;
	/*@}*/

	/** \name info strings */
	/*@{*/
	StringT fVersion;
	StringT fTitle;
	StringT fRestartFile;
	/*@}*/
	
	/** \name execution parameters */
	/*@{*/
	GlobalT::AnalysisCodeT fAnalysisCode;
	IOBaseT::FileTypeT  fOutputFormat;
	bool fReadRestart;
	int  fWriteRestart;
	bool fPrintInput;
	/*@}*/
	
	/** the managers */
	/*@{*/
	TimeManagerT* fTimeManager;
	NodeManagerPrimitive* fNodeManager;
	ElementListT fElementGroups;
	ArrayT<SolverT*> fSolvers;
	ArrayT<ControllerT*> fControllers;
	IOManager*    fIOManager;
	ModelManagerT* fModelManager;
	/*@}*/	

	/** \name run time information */
	/*@{*/
	GlobalT::StateT fStatus; /**< state */
	int fRestartCount; 	     /**< restart output counter */
	/*@}*/
	
	/** \name equation system
	 * information by group is determined during the call to 
	 * FEManagerT::SetEquationSystem */
	/*@{*/
	iArrayT fGlobalEquationStart;
	iArrayT fActiveEquationStart;
	iArrayT fGlobalNumEquations;
	/*@}*/
};

/* inlines */
inline const StringT& FEManagerT::Version(void) const { return fVersion; }
inline ifstreamT& FEManagerT::Input(void) const { return fMainIn;  }
inline ofstreamT& FEManagerT::Output(void) const { return fMainOut; }
inline const GlobalT::StateT& FEManagerT::RunState(void) const { return fStatus; }
inline IOBaseT::FileTypeT FEManagerT::OutputFormat(void) const { return fOutputFormat; }
inline ModelManagerT* FEManagerT::ModelManager (void) const { return fModelManager; }
inline NodeManagerPrimitive* FEManagerT::NodeManager(void) const { return fNodeManager; }
inline IOManager* FEManagerT::OutputManager(void) const { return fIOManager; }
inline const iArrayT* FEManagerT::ElementMap(const StringT& block_ID) const
{
#pragma unused(block_ID)
	return NULL;
}

inline int FEManagerT::GlobalEquationStart(int group) const { return fGlobalEquationStart[group]; }
inline int FEManagerT::ActiveEquationStart(int group) const { return fActiveEquationStart[group]; };
inline int FEManagerT::GlobalNumEquations(int group) const { return fGlobalNumEquations[group]; }

#endif /* _FE_MANAGER_H_ */
