/* $Id: ElementSupportT.h,v 1.6 2002-10-20 22:48:15 paklein Exp $ */
#ifndef _ELEMENT_SUPPORT_T_H_
#define _ELEMENT_SUPPORT_T_H_

/* headers */
#include <iostream.h>
#include "ExceptionT.h"

/* direct members */
#include "GlobalT.h"

namespace Tahoe {

/* forward declarations */
class FEManagerT;
class NodeManagerT;
class GroupAverageT;
class XDOF_ManagerT;
class ElementMatrixT;
template <class TYPE> class nArrayT;
class dArrayT;
class ifstreamT;
class ofstreamT;
class ScheduleT;
class ElementBaseT;
class ModelManagerT;
class iArrayT;
class StringT;
class OutputSetT;
class dArray2DT;
class LocalArrayT;
class FieldT;
class eControllerT;

/** support for the ElementBaseT class hierarchy. A limited interface to get 
 * information in and out of an ElementBaseT */
class ElementSupportT
{
public:

	/** constructor */
	ElementSupportT(void);

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
	const StringT& Version(void) const;
	
	/** verbose echo */
	bool PrintInput(void) const;
	
	/** number of nodes */
	int NumNodes(void) const { return fNumNodes; };
	
	/** number of spatial dimensions */
	int NumSD(void) const { return fNumSD; };

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

	/** solver iteration numbers */
	const int& IterationNumber(int group) const;
	
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

	/** the element group at the specified index in the element list */
	ElementBaseT& ElementGroup(int index) const;

	/** geometry information */
	ModelManagerT& Model(void) const;
	
	/** XDOF support */
	XDOF_ManagerT& XDOF_Manager(void) const;
	
	/** node number map. returns NULL if there is not map */
	const iArrayT* NodeMap(void) const;
	
	/** element number map for the given block ID */
	const iArrayT* ElementMap(const StringT& block_ID) const;

	/** return a pointer to the field with the specified name. returns NULL
	 * if a field with the given name is not found. */
	const FieldT* Field(const char* name) const;

	/** return the element controller appropriate for the given field */
	const eControllerT* eController(const FieldT& field) const;
	/*@}*/

	/** \name basic MP support */
	/*@{*/
	/** total number of processes */
	int Size(void) const;

	/** rank of this process */
	int Rank(void) const;

	/** the nodes not native to this processor */
	void IncomingNodes(iArrayT& nodes_in) const;

	/** the nodes native to this processor that appear on other processors */
	void OutgoingNodes(iArrayT& nodes_out) const;

	/** send data out.
	 * \param all_out_data outgoing data for \e every node on this processor: [nnd] x [nvals] */
	void SendExternalData(const dArray2DT& all_out_data) const;

	/** receive incoming data.
	 * \param external_data data from other processors for each node in 
	 *        the ElementSupport::IncomingNodes array */
	void RecvExternalData(dArray2DT& external_data) const;
	/*@}*/

	/** \name assembly functions */
	/*@{*/
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& eqnos) const;
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& row_eqnos,
		const nArrayT<int>& col_eqnos) const;
	void AssembleRHS(int group, const dArrayT& elRes, const nArrayT<int>& eqnos) const;
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

	/** register the output set. returns the ID that should be used with
	 * ElementSupport::WriteOutput */
	int RegisterOutput(const OutputSetT& output_set) const;

	/** write results for a single output set
	 * \param ID output set ID for the given data
	 * \param n_values nodal output values
	 * \param e_values element output values */
	void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values) const;
	/*@}*/

private:

	/** \name verified access 
	 * Use these if you don't want to keep checking that the pointers
	 * have been initialized. */
	/*@{*/
	/** the top-level manager */
	FEManagerT& FEManager(void) const;

	/** the nodes */
	NodeManagerT& Nodes(void) const;
	/*@}*/

private:

	/** the boss */
	FEManagerT* fFEManager;
	
	/** the nodes */
	NodeManagerT* fNodes;
	
	/** \name cached parameters
	 * Pre-set to allow fast access */
	/*@{*/
	GlobalT::AnalysisCodeT fAnalysis;
	const GlobalT::StateT* fRunState;
	int fNumSD;
	int fNumNodes;		
	/*@}*/
};

/* the top-level manager */
inline FEManagerT& ElementSupportT::FEManager(void) const
{
	if (!fFEManager) {
		cout << "\n ElementSupportT::FEManager: pointer not set" << endl;
		throw ExceptionT::kGeneralFail;
	}
	return *fFEManager;
}

/* the nodes */
inline NodeManagerT& ElementSupportT::Nodes(void) const
{
	if (!fNodes) {
		cout << "\n ElementSupportT::Nodes: pointer not set" << endl;
		throw ExceptionT::kGeneralFail;
	}
	return *fNodes;
}

/* return a const reference to the run state flag */
inline const GlobalT::StateT& ElementSupportT::RunState(void) const
{
	if (!fRunState) {
		cout << "\n ElementSupportT::RunState: not set" << endl;
		throw ExceptionT::kGeneralFail;
	}
	return *fRunState;
}

} // namespace Tahoe 
#endif /* _ELEMENT_SUPPORT_T_H_ */
