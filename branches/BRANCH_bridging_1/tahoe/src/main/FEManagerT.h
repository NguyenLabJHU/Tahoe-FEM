/* $Id: FEManagerT.h,v 1.32.2.3 2003-02-14 02:49:04 paklein Exp $ */
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
#include "iArray2DT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ofstreamT;
class ModelManagerT;
class TimeManagerT;
class NodeManagerT;
class IntegratorT;
class nIntegratorT;
class eIntegratorT;
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
class CommunicatorT;
class CommManagerT;

class FEManagerT: public iConsoleObjectT
{
public:

	/** degree of initialization. Passed to FEManagerT::Initialize. */
	enum InitCodeT {kFull = 0, /**< initialize to solve */
	      kParametersOnly = 1, /**< read top-level parameters only */
	        kAllButSolver = 2  /**< do everything except initialize the equation system and solvers */
	        };

	/** constructor */
	FEManagerT(ifstreamT& input, ofstreamT& output, CommunicatorT& comm);

	/** destructor */
	virtual ~FEManagerT(void);

	/** initialize members */
	virtual void Initialize(InitCodeT init = kFull);
	
	/** solve all the time sequences */
	void Solve(void);
	
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

	/** the current group. When performing operations over a group, this 
	 * returns the group being operated on. At other times it will return -1. */
	const int& CurrentGroup(void) const { return fCurrentGroup; }

	/** returns true for verbose echo of input */
	bool PrintInput(void) const;

	/** version number */
	static const char* Version(void);
	IOBaseT::FileTypeT OutputFormat(void) const;
	const StringT& Title(void) const { return fTitle; };

	/** returns 1 of ALL element groups have interpolant DOF's */
	int InterpolantDOFs(void) const;

	/** pointer to the I/O manager */
	IOManager* OutputManager(void) const;

	/** the model database manager */
	ModelManagerT* ModelManager (void) const;

	/** the time manager */
	TimeManagerT* TimeManager(void) const;

	/** the node manager */
	NodeManagerT* NodeManager(void) const;

	/** the communication manager */
	CommManagerT* CommManager(void) const;

	/** pointer to an element group */
	ElementBaseT* ElementGroup(int groupnumber) const;

	/** resolve the index of the given element group */
	int ElementGroupNumber(const ElementBaseT* pgroup) const;

	/** pointer to an element group */
	SolverT* Solver(int group) { return fSolvers[group]; };
	
	/** the MP communicator */
	CommunicatorT& Communicator(void) const { return fComm; };
	/*@}*/

	/** \name equation system */
	/*@{*/
	/** (re-)set the equation number for the given group */
	virtual void SetEquationSystem(int group);

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
	
	/** \name time */
	/*@{*/
	/* load control functions (returns true if successful) */
	bool DecreaseLoadStep(void);
	bool IncreaseLoadStep(void);
	
	/* solution accessors */
	const double& Time(void) const;
	const double& TimeStep(void) const;
	const int& StepNumber(void) const;
	const int& NumberOfSteps(void) const;
	void SetTimeStep(double dt) const;
	int SequenceNumber(void) const;
	int NumSequences(void) const;
	/*@}*/

	/** \name solution messaging 
	 * Methods called by the solver during the solution process. Either can be called
	 * an arbitrary number of times per time increment. However, FEManagerT::FormRHS
	 * must be called before the corresponding call to FEManagerT::FormLHS during a
	 * given iteration. */
	/*@{*/
	/** iteration number for the solution of the given group over the
	 * current time increment */	
	const int& IterationNumber(int group) const;

	/** return the iteration number of the current group. Returns -1
	 * if no group is current. */
	int IterationNumber(void) const;

	/** compute LHS-side matrix and assemble to solver.
	 * \param group equation group to solve
	 * \param sys_type "maximum" LHS matrix type needed by the solver. The GlobalT::SystemTypeT
	 *        enum is ordered by generality. The solver should indicate the most general
	 *        system type that is actually needed. */
	void FormLHS(int group, GlobalT::SystemTypeT sys_type) const;

	/** compute RHS-side, residual force vector and assemble to solver
	 * \param group equation group to solve */
	void FormRHS(int group) const;

	/** send update of the solution to the NodeManagerT */
	virtual void Update(int group, const dArrayT& update);

	/** system relaxation */
	virtual GlobalT::RelaxCodeT RelaxSystem(int group) const;

	/** return the current values of the unknowns 
	 * \param group equation group 
	 * \param order time derivative of the unknowns to collect. Must be
	 *        in range
	 * \param unknowns destination for the current values field values
	 *        for unprescribed degrees of freedom */
	virtual void GetUnknowns(int group, int order, dArrayT& unknowns) const;
	/*@}*/

	/** \name assembly methods 
	 * methods for assembling contributions to the global equations systems */
	/*@{*/
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& eqnos) const;
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& row_eqnos,
		const nArrayT<int>& col_eqnos) const;
	void AssembleLHS(int group, const nArrayT<double>& diagonal_elMat, const nArrayT<int>& eqnos) const;
	void OverWriteLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& eqnos) const;
	void DisassembleLHS(int group, dMatrixT& elMat, const nArrayT<int>& eqnos) const;
	void DisassembleLHSDiagonal(int group, dArrayT& diagonals, const nArrayT<int>& eqnos) const;

	void AssembleRHS(int group, const nArrayT<double>& elRes, const nArrayT<int>& eqnos) const;
	void OverWriteRHS(int group, const dArrayT& elRes, const nArrayT<int>& eqnos) const;
	void DisassembleRHS(int group, dArrayT& elRes, const nArrayT<int>& eqnos) const;
	/*@}*/
	
	/** \name output */
	/*@{*/
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
	virtual void WriteOutput(double time);
	
	/** write results for a single output set
	 * \param ID output set ID for the given data
	 * \param n_values nodal output values
	 * \param e_values element output values */
	virtual void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values);

	/** write a geometry file for the current model */
	void WriteGeometryFile(const StringT& file_name, IOBaseT::FileTypeT output_format) const;

	/** (temporarily) direct output away from main out */
	virtual void DivertOutput(const StringT& outfile);

	/** restore outputs to their regular destinations */
	virtual void RestoreOutput(void);

	/** collect the internal force on the specified node */
	void InternalForceOnNode(const FieldT& field, int node, dArrayT& force) const;
	/*@}*/

	/** \name access to integrators */
	/*@{*/
	int NumIntegrators(void) const { return fIntegrators.Length(); };
	IntegratorT* Integrator(int index) { return fIntegrators[index]; };
	const IntegratorT* Integrator(int index) const { return fIntegrators[index]; };
	eIntegratorT* eIntegrator(int index) const;
	nIntegratorT* nIntegrator(int index) const;
	/*@}*/

	/** debugging */
	virtual void WriteSystemConfig(ostream& out, int group) const;

	/** \name basic MP info */
	/*@{*/
	int Rank(void) const;
	int Size(void) const;
	
	/** the local node to home processor map. Returns the home processor
	 * for each local node. Returns NULL if there is no map, indicating 
	 * that the native processor for all nodes is this one. */
	const ArrayT<int>* ProcessorMap(void) const;

	/** node numbering map. The global id of each local node. Returns
	 * NULL if there is no map, indicating the local and global node
	 * numbers are the same. */
	const ArrayT<int>* NodeMap(void) const;

	/** list of nodes owned by the partition. Returns NULL if there is no list,
	 * indicating \e all nodes are owned by this partition */
	const ArrayT<int>* PartitionNodes(void) const;

	virtual const iArrayT* ElementMap(const StringT& block_ID) const;

	void Wait(void);
	/*@}*/

	/** interactive */
	virtual bool iDoCommand(const CommandSpecT& command, StringT& line);

	/** \name solution steps 
	 * These methods are called during FEManagerT::Solve during the solution
	 * process. These would only be called if the solution process is going to
	 * be driven externally, i.e., without calling FEManagerT::Solve.
	 * All steps return ExceptionT::kNoError = 0 unless a problem occurs. */
	/*@{*/
	/** (re-)set system to initial conditions */
	virtual void InitialCondition(void);

	/** initialize the current time increment for all groups */
	virtual ExceptionT::CodeT InitStep(void);

	/** execute the solution procedure */
	virtual ExceptionT::CodeT SolveStep(void);

	/** close the current time increment for all groups */
	virtual ExceptionT::CodeT CloseStep(void);

	/** called for all groups if the solution procedure for any group fails */
	virtual ExceptionT::CodeT ResetStep(void);
	
	/** solver phase information. Results of the last call to FEManagerT::SolveStep */
	const iArray2DT& SolverPhasesStatus(void) const { return fSolverPhasesStatus; };

	/** enum defining the contents of the FEManagerT::SolverPhasesStatus array */
	enum PhaseStatusT {
          kGroup = 0, /**< solver group associated with the phase */
      kIteration = 1, /**< number of iterations executed by the phase */
           kPass = 2, /**< 1/0 if phase PASS/FAIL-ed test conditions */
 kNumStatusFlags = 3  /**< number of phase status flags */
	};
	/*@}*/

	/** \name initialize/restart functions */
	/*@{*/
	void ReadRestart(const StringT* file_name = NULL);
	void WriteRestart(const StringT* file_name = NULL) const;
	/*@}*/

protected:

	/** "const" function that sets the status flag */
	void SetStatus(GlobalT::StateT status) const;

	/** look for input file key and check file version */
	void CheckInputFile(void);

	/** \name phases of FEManagerT::Initialize. */
	/*@{*/
	virtual void ReadParameters(InitCodeT init);
	virtual void WriteParameters(void) const;
	void SetIntegrator(void);
	virtual void SetNodeManager(void);
	virtual void SetElementGroups(void);
	void SetSolver(void);
	virtual void SetOutput(void);
	/*@}*/

	/** (re-) set cached value of the first equation number for the given
	 * group on this processor. This value is cached because communication 
	 * is required. */
	virtual int GetGlobalEquationStart(int group) const;

	/** (re-) set cached value of the total number of equations for the given
	 * group. This value is cached because communication is required. */
	virtual int GetGlobalNumEquations(int group) const;
	
	/** collect element equations and send to solver */
	void SendEqnsToSolver(int group) const;

	/** construct a new CommManagerT. Should be called some time after the
	 * ModelManagerT has been constructed */
	virtual CommManagerT* New_CommManager(void) const;

private:

	/** \name disallowed */
	/*@{*/
	FEManagerT(FEManagerT&);
	FEManagerT& operator=(FEManagerT&) const;
	/*@}*/

	/** construct a solver of the specified type. This function cannot be
	 * const because a non-const reference to the FEManagerT is passed to
	 * the solvers. */
	SolverT* New_Solver(int code, int group);
		
protected:

	/** \name I/O streams */
	/*@{*/
	ifstreamT& fMainIn;
	ofstreamT& fMainOut;
	/*@}*/
	
	/** MP environment */
	CommunicatorT& fComm;

	/** \name info strings */
	/*@{*/
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
	
	/** \name the managers */
	/*@{*/
	TimeManagerT* fTimeManager;
	NodeManagerT* fNodeManager;
	ElementListT fElementGroups;
	ArrayT<SolverT*> fSolvers;
	ArrayT<IntegratorT*> fIntegrators;
	IOManager*    fIOManager;
	ModelManagerT* fModelManager;
	CommManagerT* fCommManager;
	/*@}*/
	
	/** \name multi-solver phases */
	/*@{*/
	/** multi-solver phases. For cases with more than one solver, this
	 * array contains information about how the multiple solvers should be
	 * handled to determine the final solution. When there is only one solver
	 * this array will be empty. */
	iArray2DT fSolverPhases;

	/** maximum number of loops through the solvers. This is either a number
	 * greater than zero or -1, for no limit. */
	int fMaxSolverLoops;
	
	/** status of solver phases. Updated during call to FEManagerT::SolveStep */
	iArray2DT fSolverPhasesStatus;
	/*@}*/

	/** \name run time information */
	/*@{*/
	GlobalT::StateT fStatus; /**< state */
	int fRestartCount; 	     /**< restart output counter */
	int fCurrentGroup;       /**< current group being operated on */
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
inline ifstreamT& FEManagerT::Input(void) const { return fMainIn;  }
inline ofstreamT& FEManagerT::Output(void) const { return fMainOut; }
inline const GlobalT::StateT& FEManagerT::RunState(void) const { return fStatus; }
inline IOBaseT::FileTypeT FEManagerT::OutputFormat(void) const { return fOutputFormat; }
inline ModelManagerT* FEManagerT::ModelManager (void) const { return fModelManager; }
inline TimeManagerT* FEManagerT::TimeManager(void) const { return fTimeManager; }
inline NodeManagerT* FEManagerT::NodeManager(void) const { return fNodeManager; }
inline CommManagerT* FEManagerT::CommManager(void) const { return fCommManager; }
inline IOManager* FEManagerT::OutputManager(void) const { return fIOManager; }
inline const iArrayT* FEManagerT::ElementMap(const StringT& block_ID) const
{
#pragma unused(block_ID)
	return NULL;
}

inline int FEManagerT::IterationNumber(void) const
{
	int curr_group = CurrentGroup();
	if (curr_group >= 0)
		return IterationNumber(curr_group);
	else
		return -1;
}

inline int FEManagerT::GlobalEquationStart(int group) const { return fGlobalEquationStart[group]; }
inline int FEManagerT::ActiveEquationStart(int group) const { return fActiveEquationStart[group]; };
inline int FEManagerT::GlobalNumEquations(int group) const { return fGlobalNumEquations[group]; }

} // namespace Tahoe 
#endif /* _FE_MANAGER_H_ */
