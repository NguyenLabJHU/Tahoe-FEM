/* $Id: BasicSupportT.h,v 1.1.4.1 2004-04-08 07:33:55 paklein Exp $ */
#ifndef _TAHOE_SUPPORT_T_H_
#define _TAHOE_SUPPORT_T_H_

/* direct members */
#include "GlobalT.h"

namespace Tahoe {

/* forward declarations */
class FEManagerT;
class NodeManagerT;
class XDOF_ManagerT;
class dArray2DT;
class LocalArrayT;
class ScheduleT;
class ifstreamT;
class ofstreamT;
class FieldT;
class CommunicatorT;
class ElementMatrixT;
class OutputSetT;
class ModelManagerT;
class CommManagerT;
class TimeManagerT;
class ElementBaseT;
template <class TYPE> class nArrayT;

/** Base class for support within tahoe. Provides a limited interface to get 
 * information in and out components. */
class BasicSupportT
{
public:

	/** constructor */
	BasicSupportT(void);

	/** \name initialization 
	 * Cached values are reset when source are reset */
	/*@{*/
	/** (re-)set the FEManagerT */
	void SetFEManager(const FEManagerT* fe);

	/** (re-)set the NodeManagerT */
	void SetNodeManager(NodeManagerT* nodes);
	/*@}*/

	/** \name accessors */
	/*@{*/
	/** Tahoe version string */
	const char* Version(void) const;

	/** verbose echo */
	bool PrintInput(void) const;

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
	
	/** simulation time */
	const double& Time(void) const;
	
	/** simulation step number */
	const int& StepNumber(void) const;

	/** total number of simulations steps at the current step size */
	const int& NumberOfSteps(void) const;
	
	/** time increment */
	const double& TimeStep(void) const;
	
	/** return a pointer to the field with the specified name. returns NULL
	 * if a field with the given name is not found. */
	const FieldT* Field(const char* name) const;
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

	/** node number map. returns NULL if there is not map */
	const ArrayT<int>* NodeMap(void) const;
	
	/** element number map for the given block ID */
	const iArrayT* ElementMap(const StringT& block_ID) const;
	/*@}*/

	/** \name assembly functions */
	/*@{*/
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& eqnos) const;
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& row_eqnos,
		const nArrayT<int>& col_eqnos) const;
	void AssembleLHS(int group, const nArrayT<double>& diagonal_elMat, const nArrayT<int>& eqnos) const;
	void AssembleRHS(int group, const nArrayT<double>& elRes, const nArrayT<int>& eqnos) const;
	/*@}*/

	/** \name input/output */
	/*@{*/
	/** the parameters stream */
	ifstreamT& Input(void) const;

	/** the echo file */
	ofstreamT& Output(void) const;

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/** register the output set. returns the ID that should be used with
	 * ElementSupport::WriteOutput */
	int RegisterOutput(const OutputSetT& output_set) const;
#endif

	/** write results for a single output set
	 * \param ID output set ID for the given data
	 * \param n_values nodal output values
	 * \param e_values element output values */
	void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values) const;
	
	/** return true if output is going to be written for the current time step */
	bool WriteOutput(void) const;
	
	/** return a reference to the output set with the given ID */
	const OutputSetT& OutputSet(int ID) const;
	/*@}*/

 	/** \name verified access 
	 * Use these if you don't want to keep checking that the pointers
	 * have been initialized. */
	/*@{*/
	/** the top-level manager */
	const FEManagerT& FEManager(void) const;

	/** the nodes */
	NodeManagerT& NodeManager(void) const;

	/** XDOF support */
	XDOF_ManagerT& XDOF_Manager(void) const;

	/** the time manager */
	TimeManagerT& TimeManager(void) const;

	/** geometry information */
	ModelManagerT& ModelManager(void) const;

	/** comm information */
	CommManagerT& CommManager(void) const;

	/** the element group at the specified index in the element list */
	ElementBaseT& ElementGroup(int index) const;
	/*@}*/

private:

	/** \name cached pointers to components */
	/*@{*/
	/** the top level */
	const FEManagerT* fFEManager;
	
	/** the nodes */
	NodeManagerT* fNodeManager;

	/** the time manager */
	TimeManagerT* fTimeManager;

	/** the model manager */
 	ModelManagerT* fModelManager;	

	/** the communication manager */
	CommManagerT* fCommManager;

	/** low-level communicator */
	const CommunicatorT* fCommunicator;
	/*@}*/
	
	/** \name cached data */
	/*@{*/
	GlobalT::AnalysisCodeT fAnalysis;
	const GlobalT::StateT* fRunState;

#ifdef _FRACTURE_INTERFACE_LIBRARY_
	int fNumNodes;
	int fNumSD;
	double fTimeStep;
	int fIterationNumber;
	ifstreamT *ifst;
	ofstreamT *ofst;	
	dArray2DT *fInitialCoordinates;
	dArray2DT *fCurrentCoordinates;
#endif
	/*@}*/
};

/* the top-level manager */
inline const FEManagerT& BasicSupportT::FEManager(void) const {
	if (!fFEManager) ExceptionT::GeneralFail("BasicSupportT::FEManager", "pointer not set");
	return *fFEManager;
}

/* the nodes */
inline NodeManagerT& BasicSupportT::NodeManager(void) const {
	if (!fNodeManager) ExceptionT::GeneralFail("BasicSupportT::NodeManager", "pointer not set");
	return *fNodeManager;
}

/* the time manager */
inline TimeManagerT& BasicSupportT::TimeManager(void) const {
	if (!fTimeManager)  ExceptionT::GeneralFail("BasicSupportT::TimeManager", "pointer not set");
	return *fTimeManager;
}

/* return a const reference to the run state flag */
inline const GlobalT::StateT& BasicSupportT::RunState(void) const {
	if (!fRunState) ExceptionT::GeneralFail("BasicSupportT::RunState", "not set");
	return *fRunState;
}

/* geometry information */
inline ModelManagerT& BasicSupportT::ModelManager(void) const {
	if (!fModelManager) ExceptionT::GeneralFail("BasicSupportT::Model", "pointer not set");
	return *fModelManager;
}

/* comm information */
inline CommManagerT& BasicSupportT::CommManager(void) const {
	if (!fCommManager) ExceptionT::GeneralFail("BasicSupportT::CommManager", "pointer not set");
	return *fCommManager;
}

inline int BasicSupportT::NumNodes(void) const {
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return fNumNodes;
#else
	return InitialCoordinates().MajorDim();
#endif
}
	
inline int BasicSupportT::NumSD(void) const {
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return fNumSD;
#else
	return InitialCoordinates().MinorDim();
#endif
}

/* low-level communicator */
inline const CommunicatorT& BasicSupportT::Communicator(void) const {
	if (!fCommunicator) ExceptionT::GeneralFail("BasicSupportT::Communicator", "pointer not set");
	return *fCommunicator;
}

} /* namespace Tahoe */

#endif /* _TAHOE_SUPPORT_T_H_ */
