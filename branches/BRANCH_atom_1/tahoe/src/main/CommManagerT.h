/* $Id: CommManagerT.h,v 1.1.2.3 2002-12-16 09:22:16 paklein Exp $ */
#ifndef _COMM_MANAGER_T_H_
#define _COMM_MANAGER_T_H_

/* direct members */
#include "dArray2DT.h"
#include "AutoArrayT.h"
#include "InverseMapT.h"

namespace Tahoe {

/* forward declarations */
class PartitionT;
class CommunicatorT;
class ModelManagerT;
class iArrayT;

/** manage processor to processor transactions. Manages partition information.
 * Creates ghosts nodes. Manages communication lists. Manipulates the
 * ModelManagerT to create the environment for the local processor. */
class CommManagerT
{
public:

	/** constructor */
	CommManagerT(CommunicatorT& comm, ModelManagerT& model_manager);

	/** set or clear partition information */
	void SetPartition(PartitionT* partition);

	/** \name setting periodic boundaries */
	/*@{*/
	void SetPeriodicBoundaries(int i, double x_i_min, double x_i_max);
	void ClearPeriodicBoundaries(int i);
	/*@}*/

	/** configure the current local coordinate list and register it with the
	 * model manager. The first time this method is called, it will call
	 * CommManagerT::FirstConfigure before performing the usual operations. */
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

	/** list of nodes owned by the partition. Returns NULL if there is no list,
	 * indicating \e all nodes are owned by this partition */
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
	/*@}*/

private:

	/** collect partition nodes */
	void CollectPartitionNodes(const ArrayT<int>& n2p_map, int part, 
		AutoArrayT<int>& part_nodes) const;

	/** perform actions needed the first time CommManagerT::Configure is called. */
	void FirstConfigure(void);

	/** determine the local coordinate bounds 
	 * \param coords coordinate list
	 * \param local rows of coords to use in determining bounds
	 * \param bounds returns with the bounds of the given points */
	void GetBounds(const dArray2DT& coords, const iArrayT& local, dArray2DT& bounds) const;

private:

	/** communicator */
	CommunicatorT& fComm;

	/** the model manager */
	ModelManagerT& fModelManager;

	/** \name periodic boundaries */
	/*@{*/
	/** flags to indicate if periodic boundary conditions are imposed */
	ArrayT<bool> fIsPeriodic;
	
	/** rows give the lower and upper periodic bounds for that coordinate */
	dArray2DT fPeriodicBoundaries;
	/*@}*/
	
	/** processor bounds */
	dArray2DT fBounds;

	/** partition information */
	PartitionT* fPartition;
	
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

} /* namespace Tahoe */

#endif /* _COMM_MANAGER_T_H_ */
