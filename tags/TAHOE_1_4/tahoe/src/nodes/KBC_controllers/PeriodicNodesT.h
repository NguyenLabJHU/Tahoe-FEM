/* $Id: PeriodicNodesT.h,v 1.2 2002-10-20 22:41:38 paklein Exp $ */
#ifndef _PERIODIC_NODES_T_H_
#define _PERIODIC_NODES_T_H_

/* base class */
#include "TiedNodesT.h"

namespace Tahoe {

/** Nodes tied across periodic boundaries. The function is the same
 * as TiedNodesT in that the degrees of freedom (and their derivatives)
 * from the leader are prescribed for the follower node. This class
 * differs in that pairs are determined by comparing coordinates
 * across periodic strides along the coordinate axes. */
class PeriodicNodesT: public TiedNodesT
{
public:	

	/** constructor */
	PeriodicNodesT(NodeManagerT& node_manager, BasicFieldT& field);

	/** write class parameters */
	void WriteParameters(ostream& out) const;

protected:

	/** \name called by TiedNodesT::Initialize */
	/*@{*/
	/** read class parameters. Called after nodes information is read.
	 * Reads peridic boundary information. */
	virtual void ReadParameters(ifstreamT& in);

	/** set initial tied node pairs. Initializes the data in TiedNodesT::fLeaderIds,
	 * TiedNodesT::fFollowerIds, TiedNodesT::fNodePairs, and TiedNodesT::fPairStatus.
	 * PeriodicNodesT::InitTiedNodePairs pairs nodes that coincide across the periodic
	 * strides. The routine is based on TiedNodesT::InitTiedNodePairs. The follower
	 * array may be resized during the operation. */
	virtual void InitTiedNodePairs(const iArrayT& leader, iArrayT& follower);
	/*@}*/

protected:

	/** true if the coordinate direction is periodic */
	ArrayT<bool> fIsPeriodic;

	/** periodic stride length along each coordinate axis */
	dArrayT fPeriodicStride;
};

} // namespace Tahoe 
#endif /* _PERIODIC_NODES_T_H_ */
