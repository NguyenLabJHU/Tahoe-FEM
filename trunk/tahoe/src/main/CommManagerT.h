/* $Id: CommManagerT.h,v 1.8 2004-11-10 17:10:41 paklein Exp $ */
#ifndef _COMM_MANAGER_T_H_
#define _COMM_MANAGER_T_H_

/* direct members */
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "AutoArrayT.h"
#include "InverseMapT.h"
#include "MessageT.h" /* message enum's */
#include "dArrayT.h"
#include "nVariArray2DT.h"
#include "Array2DT.h"

namespace Tahoe {

/* forward declarations */
class PartitionT;
class CommunicatorT;
class ModelManagerT;
class NodeManagerT;
class iArrayT;
class MessageT;

/** manage processor to processor transactions. Manages partition information.
 * Creates ghosts nodes. Manages communication lists. Manipulates the
 * ModelManagerT to create the environment for the local processor. */
class CommManagerT
{
public:

	/** enum */
	enum CodesT { kNULLMessageID = -99 };

	/** constructor */
	CommManagerT(CommunicatorT& comm, ModelManagerT& model_manager);

	/** destructor */
	~CommManagerT(void);

	/** number of processes */
	int Size(void) const { return fSize; };

	/** rank of this process */
	int Rank(void) const { return fRank; };

	/** width of the communication skin */
	double Skin(void) const { return fSkin; };

	/** set the width of the communication skin */
	void SetSkin(double skin) { fSkin = skin; };

	/** set or clear partition information. Needs to be set before calling 
	 * CommManagerT::Configure. */
	void SetPartition(PartitionT* partition);

	/** set or clear node manager information */
	void SetNodeManager(NodeManagerT* node_manager);

	/** the communicator */
	const CommunicatorT& Communicator(void) const { return fComm; };

	/** \name periodic boundaries */
	/*@{*/
	void SetPeriodicBoundaries(int i, double x_i_min, double x_i_max);
	void ClearPeriodicBoundaries(int i);
	
	/** accessor */
	const dArray2DT& PeriodicBoundaries(void) { return fPeriodicBoundaries; };
	
	/** enforce the periodic boundary conditions */
	void EnforcePeriodicBoundaries(void);

	/** return the number of real nodes */
	int NumRealNodes(void) const { return fNumRealNodes; };
	/*@}*/

	/** configure the current local coordinate list and register it with the
	 * model manager. The first time this method is called, it will call
	 * CommManagerT::FirstConfigure before performing the usual operations. 
	 * Partition information, if it exists should be (un-)set with
	 * CommManagerT::SetPartition before calling this. */
	void Configure(void);

	/** \name numbering maps */
	/*@{*/
	/** the local node to home processor map. Returns the home processor
	 * for each local node. Returns NULL if there is no map, indicating 
	 * that the native processor for all nodes is this one. */
	const ArrayT<int>* ProcessorMap(void) const;

	/** node numbering map. The global id of each local node. Returns
	 * NULL if there is no map, indicating the local and global node
	 * numbers are the same. */
	const ArrayT<int>* NodeMap(void) const;

	/** list of nodes owned by the partition. These nodes are numbered locally
	 * within the nodes appearing on this processor. Returns NULL if there is no 
	 * list, indicating \e all nodes are owned by this partition */
	const ArrayT<int>* PartitionNodes(void) const;

	/** return true if the list of partition nodes may be changing */
	bool PartitionNodesChanging(void) const;

	/** inverse of fPartitionNodes list. Gives the index in CommManagerT::fPartitionNodes
	 * of the nodes listed in that array:
	 *
	 *    fPartitionNodes_inv[fPartitionNodes[i]] = i
	 *
	 * Returns NULL if there is no inverse map, indicating all nodes on this partition
	 * are owned by this partition.
	 */
	const InverseMapT* PartitionNodes_inv(void) const;

	/** list of nodes \e not owned by the partition. These are the nodes for which information
	 * is received from other processors. These nodes are numbered locally
	 * within the nodes appearing on this processor. Returns NULL if there is no 
	 * list, indicating \e all nodes are owned by this partition. */
	const ArrayT<int>* ExternalNodes(void) const;

	/** list of nodes adjacent to nodes in CommManagerT::ExternalNodes. These are the nodes 
	 * for which information is sent to other processors. These nodes are numbered locally
	 * within the nodes appearing on this processor. Returns NULL if there is no 
	 * list, indicating \e all nodes are owned by this partition */
	const ArrayT<int>* BorderNodes(void) const;
	/*@}*/

	/** \name ghost nodes 
	 * Lists of {real, ghost} node pairs created on this processor. Return NULL if there
	 * are no ghost nodes on this processor. */
	/*@{*/
	/** nodes with ghosts */
	const ArrayT<int>* NodesWithGhosts(void) const;

	/** images of the nodes listed in CommManagerT::GhostedNodes */
	const ArrayT<int>* GhostNodes(void) const;
	/*@}*/

	/** \name configuring persistent communications */
	/*@{*/
	/** set up a persistent all gather communication. Distribute the values per
	 * node for all nodes owned by this process and collect values from nodes
	 * on other processes. If called with an nArray2DT argument, the type and minor
	 * dimension of the array must match the subsequent calls to CommManagerT::AllGather.
	 * \param t data type to be transmitted
	 * \param num_vals number of values per node 
	 * \return ID for this communication */
	int Init_AllGather(MessageT::TypeT t, int num_vals);
	int Init_AllGather(const nArray2DT<int>& values);
	int Init_AllGather(const nArray2DT<double>& values);

	/** clear the persistent communication
	 * \param ID ID for the communication to be cleared obtained during
	 *        CommManagerT::Init_AllGather. 
	 * \note  Currently, no mechanism is implemented for recovering cleared
	 *        ID's. The ID's returned by CommManagerT::Init_AllGather will
	 *        increase sequentially regardless of whether any previous
	 *        communications have been cleared. */
	void Clear_AllGather(int id);

	/** perform the all gather. The values from this partition must already by
	 * in the appropriate location in the destination array */
	void AllGather(int id, nArray2DT<double>& values);
	void AllGather(int id, nArray2DT<int>& values);
	/*@}*/

private:

	/** collect partition nodes */
	void CollectPartitionNodes(const ArrayT<int>& n2p_map, int part, 
		AutoArrayT<int>& part_nodes) const;

	/** perform actions needed the first time CommManagerT::Configure is called. */
	void FirstConfigure(void);

	/** \name methods for configuring computation with a spatial decomposition */
	/*@{*/
	/** distribute nodes on spatial grid */
	void Distribute(void);
	
	/** set border information */
	void SetExchange(void);
	/*@}*/

	/** determine the local coordinate bounds 
	 * \param coords coordinate list
	 * \param local rows of coords to use in determining bounds
	 * \param bounds returns with the bounds of the given points */
	void GetBounds(const dArray2DT& coords, const iArrayT& local, dArray2DT& bounds) const;

	/** \name not allowed */
	/*@{*/
	/** copy constructor */
	CommManagerT(CommManagerT&);

	/** assignment operator */
	const CommManagerT& operator=(const CommManagerT&);
	/*@}*/

private:

	/** communicator */
	CommunicatorT& fComm;

	/** the model manager */
	ModelManagerT& fModelManager;

	int fSize;
	int fRank;

	/** width of communication layer for non-graph based decompositions */
	double fSkin;

	/** \name periodic boundaries */
	/*@{*/
	/** flags to indicate if periodic boundary conditions are imposed */
	ArrayT<bool> fIsPeriodic;
	
	/** rows give the lower and upper periodic bounds for that coordinate */
	dArray2DT fPeriodicBoundaries;

	/** periodic distance along each coordinate */
	dArrayT fPeriodicLength;

	int fNumRealNodes;
	AutoArrayT<int> fPBCNodes;
	AutoArrayT<int> fPBCNodes_ghost;
	AutoArrayT<int> fPBCNodes_face;
	/*@}*/
	
	/** processor bounds */
	dArray2DT fBounds;

	/** partition information */
	PartitionT* fPartition;

	/** node manager */
	NodeManagerT* fNodeManager;
	
	/** true if CommManagerT::Configure has not been called yet */
	bool fFirstConfigure;

	/** \name maps */
	/*@{*/
	/** native processor per node */
	AutoArrayT<int> fProcessor;

	/** local to global node map */
	AutoArrayT<int> fNodeMap;

	/** list of nodes owned by this partition */
	AutoArrayT<int> fPartitionNodes;

	/** inverse of CommManagerT::fPartitionNodes list */
	InverseMapT fPartitionNodes_inv;

	/** list of nodes \e not owned by this partition */
	AutoArrayT<int> fExternalNodes;

	/** list of nodes adjacent to any external nodes */
	AutoArrayT<int> fBorderNodes;
	/*@}*/
	
	/** \name persistent communications */
	/*@{*/
	/** number of values per node for each message */
	AutoArrayT<int> fNumValues;
	
	/** communications for nodal values */
	AutoArrayT<MessageT*> fCommunications;

	/** communications for ghost nodes associated with the communications
	 * in CommManagerT::fCommunications */
	AutoArrayT<MessageT*> fGhostCommunications;
	/*@}*/

	/** \name communication buffers */
	/*@{*/
	dArray2DT fd_send_buffer;
	dArray2DT fd_recv_buffer;
	nVariArray2DT<double> fd_send_buffer_man;
	nVariArray2DT<double> fd_recv_buffer_man;

	iArray2DT fi_send_buffer;
	iArray2DT fi_recv_buffer;
	nVariArray2DT<int> fi_send_buffer_man;
	nVariArray2DT<int> fi_recv_buffer_man;
	/*@}*/

	/** \name spatial decomposition */
	/*@{*/
	/** ID of surrounding processors needed for spatial decomposition only */
	iArray2DT fAdjacentCommID; /* [nsd] x {low, high} */

	/** exchange nodes */
	Array2DT<AutoArrayT<int> > fExchange;
	/*@}*/
};

/* processor map */
inline const ArrayT<int>* CommManagerT::ProcessorMap(void) const
{
	if (fProcessor.Length() > 0)
		return &fProcessor;
	else
		return NULL;
}

/* node numbering map */
inline const ArrayT<int>* CommManagerT::NodeMap(void) const
{
	if (fNodeMap.Length() > 0)
		return &fNodeMap;
	else
		return NULL;
}

/* list of nodes owned by the partition */
inline const ArrayT<int>* CommManagerT::PartitionNodes(void) const
{
	if (fPartitionNodes.Length() > 0)
		return &fPartitionNodes;
	else
		return NULL;
}

inline const ArrayT<int>* CommManagerT::ExternalNodes(void) const
{
	if (fExternalNodes.Length() > 0)
		return &fExternalNodes;
	else
		return NULL;
}

inline const ArrayT<int>* CommManagerT::BorderNodes(void) const
{
	if (fBorderNodes.Length() > 0)
		return &fBorderNodes;
	else
		return NULL;
}

/* nodes with ghosts */
inline const ArrayT<int>* CommManagerT::NodesWithGhosts(void) const
{
	if (fPBCNodes.Length() > 0)
		return &fPBCNodes;
	else
		return NULL;
}

/* the ghosts */
inline const ArrayT<int>* CommManagerT::GhostNodes(void) const
{
	if (fPBCNodes_ghost.Length() > 0)
		return &fPBCNodes_ghost;
	else
		return NULL;
}

inline int CommManagerT::Init_AllGather(const nArray2DT<int>& values)
{
	return Init_AllGather(MessageT::Integer, values.MinorDim());
}

inline int CommManagerT::Init_AllGather(const nArray2DT<double>& values)
{
	return Init_AllGather(MessageT::Double, values.MinorDim());
}

} /* namespace Tahoe */

#endif /* _COMM_MANAGER_T_H_ */
