/* $Id: ElementSupportT.h,v 1.29 2004-06-26 18:29:48 paklein Exp $ */
#ifndef _ELEMENT_SUPPORT_T_H_
#define _ELEMENT_SUPPORT_T_H_

/* headers */
#include <iostream.h>
#include "ExceptionT.h"

/* direct members */
#include "IOBaseT.h"
#include "GlobalT.h"
#include "dArray2DT.h"
#ifndef _FRACTURE_INTERFACE_LIBRARY_
#include "FieldT.h"
#else
#include "StringT.h"
#include "GroupAverageT.h"
#endif

namespace Tahoe {

/* forward declarations */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
class FEManagerT;
class NodeManagerT;
class TimeManagerT;
class XDOF_ManagerT;
class FieldT;
class eIntegratorT;
#else
class dMatrixT;
#endif
class GroupAverageT;
class ElementMatrixT;
template <class TYPE> class nArrayT;
class dArrayT;
class ifstreamT;
class ofstreamT;
class ElementBaseT;
class ModelManagerT;
class iArrayT;
class ScheduleT;
class StringT;
class OutputSetT;
class LocalArrayT;
class CommManagerT;
class CommunicatorT;

/** support for the ElementBaseT class hierarchy. A limited interface to get 
 * information in and out of an ElementBaseT */
class ElementSupportT
{
public:

#ifdef _FRACTURE_INTERFACE_LIBRARY_
	/* Parameters normally read from input stream must be passed through ElementSupport */
	enum CodeT { kGeometryCode = 0, /**< Topology of surface element */
	    		    kNumIntPts = 1, /**< Number of integration points */
	             kCloseSurface = 2, /**< Initially close cohesive surfaces? */
				   kOutputArea = 3, /**< Output fracture area */
	             kMaterialCode = 4};/**< Which cohesive law to use */ 
#endif

	/** constructor */
	ElementSupportT(void);

#ifndef _FRACTURE_INTERFACE_LIBRARY_

	/** \name initialization 
	 * Cached values are reset when source are reset */
	/*@{*/
	/** (re-)set the FEManagerT */
	void SetFEManager(FEManagerT* fe);

	/** (re-)set the NodeManagerT */
	void SetNodes(NodeManagerT* nodes);
	/*@}*/

	/** \name accessors */
	/*@{*/
	/** Tahoe version string */
	const char* Version(void) const;

#endif // ndef _FRACTURE_INTERFACE_LIBRARY_

	/** verbose echo */
	bool PrintInput(void) const;
	
#ifdef _FRACTURE_INTERFACE_LIBRARY_

	/** set the number of nodes in the fracture interface */
	void SetNumNodes(int nn);	
	
	void SetInitialCoordinates(double *InitialCoords);

	/** set the memory used to hold the reference configuration */
	void SetInitialCoordinates(dArray2DT *InitialCoords);

	void SetCurrentCoordinates(double *CurrentCoords);
	
	/** set the memory used to hold the current configuration */
	void SetCurrentCoordinates(dArray2DT *CurrentCoords);

	/** use the displacements and reference configuration to update the current one */
	void UpdateCurrentCoordinates(double *displacements);

	/** set the time step in the fracture interface */
	void SetTimeStep(double dt);

	/** set the model manager in the fracture interface */
	void SetModelManager(ModelManagerT* modelManager);

	/** set the number of elements in the fracture interface */
	void SetNumElements(int nelem);

	/** accessor for the number of elements in the fracture interface */
	int NumElements(void) const;

	/** accessor for element floating point input when streams are not available */
	dArrayT *FloatInput(void) const;
	
	/** accessor for element integer input when streams are not available */
	iArrayT *IntInput(void) const;
	
	/** generate equation numbers based on connectivity information */
	void SetEqnos(int *conn, const int& nelem, const int& nElemNodes, const int&nNodes);

	dArrayT& Residual(void) const { return *fResidual; };
	
	dMatrixT& Stiffness(void) const { return *fStiffness; };
	
	void SetMaterialInput(double *inputFloats, int length);
	
	void SetElementInput(int *inputInts, int length);
		
	int ReturnInputInt(CodeT label);
	
	double *StateVariableArray(void);
	
	void SetStateVariableArray(double *incomingArray);
	
	void SetBlockID(StringT &Id);
	
	StringT& BlockID(void);
	
	void OutputSize(int& nNodeOutputVars, int& nElemOutputVars);
	
	void SetOutputCodes(iArrayT& fNodalOutputCodes, iArrayT& fElementOutputCodes);
	
	void SetOutputPointers(double *nodalOutput, double *elemOutput);
	
	void SetResidual(double *nodalForces);

#endif // def _FRACTURE_INTERFACE_LIBRARY_

	/** number of nodes */
	int NumNodes(void) const;
	
	/** number of spatial dimensions */
	int NumSD(void) const;

	/** initial coordinates */
	const dArray2DT& InitialCoordinates(void) const;
	
	/** current coordinates */
	const dArray2DT& CurrentCoordinates(void) const;

	/** register the local coordinate array with its source */
	void RegisterCoordinates(LocalArrayT& array) const;

	/** global analysis type */
	GlobalT::AnalysisCodeT Analysis(void) const { return fAnalysis; };

	/** return a const reference to the run state flag */
	const GlobalT::StateT& RunState(void) const;

	/** return a pointer to the specified schedule function. Returns
	 * NULL if the number is out of range. */
	const ScheduleT* Schedule(int num) const;

	/** solver iteration number for the specified group */
	const int& IterationNumber(int group) const;
	
	/** return the iteration number for the current solver group. Returns
	 * -1 of no solver group is current */
	int IterationNumber(void) const;

	/** the group number being solved or -1 if not defined */
	int CurrentGroup(void) const;
	
	/** exception string */
	const char* Exception(ExceptionT::CodeT exception) const;
	
	/** simulation time */
	const double& Time(void) const;
	
	/** simulation step number */
	const int& StepNumber(void) const;

	/** total number of simulations steps at the current step size */
	const int& NumberOfSteps(void) const;
	
	/** time increment */
	const double& TimeStep(void) const;
	
	/** index of the element group in the list of elements */
	int ElementGroupNumber(const ElementBaseT* element) const;

#ifndef _FRACTURE_INTERFACE_LIBRARY_

	/** the element group at the specified index in the element list */
	ElementBaseT& ElementGroup(int index) const;

	/** XDOF support */
	XDOF_ManagerT& XDOF_Manager(void) const;

	/** return a pointer to the field with the specified name. returns NULL
	 * if a field with the given name is not found. */
	const FieldT* Field(const char* name) const;

	/** return the element controller appropriate for the given field */
	const eIntegratorT* eIntegrator(const FieldT& field) const;
	/*@}*/

	/** \name basic MP support */
	/*@{*/
	/** total number of processes */
	int Size(void) const;

	/** rank of this process */
	int Rank(void) const;

	/** low-level global communicator */
	const CommunicatorT& Communicator(void) const;

	/** the nodes not native to this processor. Returns NULL if there is no 
	 * list, indicating \e all nodes are owned by this partition */
	const ArrayT<int>* ExternalNodes(void) const;

	/** the nodes native to this processor that appear on other processors.
	 * Returns NULL if there is no list, indicating \e all nodes are owned by 
	 * this partition */
	const ArrayT<int>* BorderNodes(void) const;
	/*@}*/	
#endif
	
	/** geometry information */
	ModelManagerT& Model(void) const;

	/** comm information */
	CommManagerT& CommManager(void) const;

	/** node number map. returns NULL if there is not map */
	const ArrayT<int>* NodeMap(void) const;
	
	/** element number map for the given block ID */
	const iArrayT* ElementMap(const StringT& block_ID) const;

	/** \name assembly functions */
	/*@{*/
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& eqnos) const;
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& row_eqnos,
		const nArrayT<int>& col_eqnos) const;
	void AssembleLHS(int group, const nArrayT<double>& diagonal_elMat, const nArrayT<int>& eqnos) const;
	void AssembleRHS(int group, const nArrayT<double>& elRes, const nArrayT<int>& eqnos) const;
	/*@}*/

	/** \name nodal averaging */
	/*@{*/
	/** initialize work space to the number of values to be averaged */
	void ResetAverage(int n_values) const;

	/** assemble values 
	 * \param nodes list of nodes for the values being assembled: [nnd] 
	 * \param vals values to be assembled: [nnd] x [nvals] */
	void AssembleAverage(const iArrayT& nodes, const dArray2DT& vals) const;

	/** average assembled values and return the array of averages 
	 * values: [nnd] x [nvals] */
	const dArray2DT& OutputAverage(void) const;

	/** return averaged values for the nodes with assembled values. Returned
	 * nodes are ordered by increasing node number */
	void OutputUsedAverage(dArray2DT& average_values) const;
	/*@}*/

	/** \name input/output */
	/*@{*/
	/** the parameters stream */
	ifstreamT& Input(void) const;

	/** the echo file */
	ofstreamT& Output(void) const;

	/** format of the output files */
	IOBaseT::FileTypeT OutputFormat(void) const;

	/** register the output set. returns the ID that should be used with
	 * ElementSupport::WriteOutput */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	int RegisterOutput(const OutputSetT& output_set) const;
#else
	/** for SIERRA interface, we're just holding the data till they want it */
	int RegisterOutput(ArrayT<StringT>& n_labels, ArrayT<StringT>& e_labels);
#endif

	/** write results for a single output set
	 * \param ID output set ID for the given data
	 * \param n_values nodal output values
	 * \param e_values element output values */
	void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values) const;
	
	/** return true if output is going to be written for the current time step */
	bool WriteOutput(void) const;

	/** write a snapshot */
	void WriteOutput(const StringT& file, const dArray2DT& coords, const iArrayT& node_map,
		const dArray2DT& values, const ArrayT<StringT>& labels) const;
	
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/** return a reference to the output set with the given ID */
	const OutputSetT& OutputSet(int ID) const;
#endif
	/*@}*/

#ifndef _FRACTURE_INTERFACE_LIBRARY_
 	/** \name verified access 
	 * Use these if you don't want to keep checking that the pointers
	 * have been initialized. */
	/*@{*/
	/** the top-level manager */
	FEManagerT& FEManager(void) const;

	/** the nodes */
	NodeManagerT& Nodes(void) const;

	/** the time manager */
	TimeManagerT& TimeManager(void) const;
	/*@}*/

	/** \name command line flags */
	/*@{*/
	const ArrayT<StringT>& Argv(void) const;
	bool CommandLineOption(const char* str) const;

	/** returns the index of the requested option or -1 if not bound */
	bool CommandLineOption(const char* str, int& index) const;
	/*@}*/
#endif

private:

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/** \name managers */
	/*@{*/
	/** the boss */
	FEManagerT* fFEManager;
	
	/** the nodes */
	NodeManagerT* fNodes;

	/** the time manager */
	TimeManagerT* fTimeManager;

	const dArray2DT *fInitialCoordinates;
	const dArray2DT *fCurrentCoordinates;
#endif

	/** the model manager */
 	ModelManagerT* fModelManager;	

	/** the communication manager */
	CommManagerT* fCommManager;
	/*@}*/
	
	/** \name cached parameters
	 * Pre-set to allow fast access */
	/*@{*/
	GlobalT::AnalysisCodeT fAnalysis;
	const GlobalT::StateT* fRunState;
	/*@}*/

#ifdef _FRACTURE_INTERFACE_LIBRARY_	
 
	dArrayT *fResidual;
	dMatrixT *fStiffness;
	
	dArray2DT *fInitialCoordinates;
	dArray2DT *fCurrentCoordinates;

	int fItNum, fElem;
	int fNumSD, fNumNodes;
	
	GroupAverageT* fGroupAverage;
	
	double fTimeStep;
	
	ifstreamT *ifst;
	ofstreamT *ofst;
	
	dArrayT *fparams;
	iArrayT *iparams;
	
	iArrayT *ieqnos;
	
	double *fStateVars, *fNodalOutput, *fElemOutput;

	StringT sBlockID;
	
	ArrayT<StringT> fNodeOutputLabels;
	ArrayT<StringT> fElemOutputLabels;
#endif
	
};

#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* the top-level manager */
inline FEManagerT& ElementSupportT::FEManager(void) const
{
	if (!fFEManager)
		ExceptionT::GeneralFail("ElementSupportT::FEManager", "pointer not set");
	return *fFEManager;
}

/* the nodes */
inline NodeManagerT& ElementSupportT::Nodes(void) const
{
	if (!fNodes) 
		ExceptionT::GeneralFail("ElementSupportT::Nodes", "pointer not set");
	return *fNodes;
}

/* the time manager */
inline TimeManagerT& ElementSupportT::TimeManager(void) const
{
	if (!fTimeManager) 
		ExceptionT::GeneralFail("ElementSupportT::TimeManager", "pointer not set");
	return *fTimeManager;
}
#else
inline int ElementSupportT::NumElements(void) const { return fElem; }
inline dArrayT *ElementSupportT::FloatInput(void) const { return fparams; }
inline iArrayT *ElementSupportT::IntInput(void) const { return iparams; }
#endif

/* return a const reference to the run state flag */
inline const GlobalT::StateT& ElementSupportT::RunState(void) const
{
	if (!fRunState)
#ifndef _FRACTURE_INTERFACE_LIBRARY_
		ExceptionT::GeneralFail("ElementSupportT::RunState", "not set");
#else
	        throw ExceptionT::kGeneralFail;
#endif
	return *fRunState;
}

/* geometry information */
inline ModelManagerT& ElementSupportT::Model(void) const
{
	if (!fModelManager) 
		ExceptionT::GeneralFail("ElementSupportT::Model", "pointer not set");
	return *fModelManager;
}

/* comm information */
inline CommManagerT& ElementSupportT::CommManager(void) const
{
	if (!fCommManager) 
#ifndef _FRACTURE_INTERFACE_LIBRARY_
		ExceptionT::GeneralFail("ElementSupportT::CommManager", "pointer not set");
#else
	        throw ExceptionT::kGeneralFail;
#endif
	return *fCommManager;
}

/* number of nodes */
inline int ElementSupportT::NumNodes(void) const
{
#if __option(extended_errorcheck)
	if (!fInitialCoordinates) ExceptionT::GeneralFail("ElementSupportT::NumNodes", 
		"no initial coordinates");
#endif
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return fInitialCoordinates->MajorDim();
#else
	return fNumNodes;
#endif
}
	
/* number of spatial dimensions */
inline int ElementSupportT::NumSD(void) const
{
#if __option(extended_errorcheck)
	if (!fInitialCoordinates) ExceptionT::GeneralFail("ElementSupportT::NumSD", 
		"no initial coordinates");
#endif
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return fInitialCoordinates->MinorDim();
#else
	return fNumSD;
#endif
}

} // namespace Tahoe 
#endif /* _ELEMENT_SUPPORT_T_H_ */
