/* $Id: PeriodicNodesT.h,v 1.4 2004-09-01 10:46:27 paklein Exp $ */
#ifndef _PERIODIC_NODES_T_H_
#define _PERIODIC_NODES_T_H_

/* base class */
#include "TiedNodesT.h"

/* direct members */
#include "dArrayT.h"

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
	PeriodicNodesT(const BasicSupportT& support, BasicFieldT& field);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** set initial tied node pairs. Initializes the data in TiedNodesT::fLeaderIds,
	 * TiedNodesT::fFollowerIds, TiedNodesT::fNodePairs, and TiedNodesT::fPairStatus.
	 * PeriodicNodesT::InitTiedNodePairs pairs nodes that coincide across the periodic
	 * strides. The routine is based on TiedNodesT::InitTiedNodePairs. The follower
	 * array may be resized during the operation. */
	virtual void InitTiedNodePairs(const iArrayT& leader, iArrayT& follower);

protected:

	/** true if the coordinate direction is periodic */
	ArrayT<bool> fIsPeriodic;

	/** periodic stride length along each coordinate axis */
	dArrayT fPeriodicStride;
};

} // namespace Tahoe 
#endif /* _PERIODIC_NODES_T_H_ */
