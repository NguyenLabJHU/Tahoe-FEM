/* $Id: CartesianShiftT.h,v 1.2 2005-06-04 17:17:15 paklein Exp $ */
#ifndef _CARTESIAN_SHIFT_T_H_
#define _CARTESIAN_SHIFT_T_H_

/* base class */
#include "MessageT.h"

/* direct members */
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "CommunicatorT.h"
#include "nVariArray2DT.h"
#include "AutoArrayT.h"

namespace Tahoe {

/* forward declarations */
class PartitionT;

/** record for communications with shifts on a cartesian grid */
class CartesianShiftT: public MessageT
{
public:

	/** constructor 
	 * \param comm communicator
	 * \param tag message tag
	 * \param send_nodes list of incoming nodes in each direction: [2*nsd]
	 * \param send_nodes list of outgoing nodes in each direction: [2*nsd] */
	CartesianShiftT(CommunicatorT& comm, int tag,
		const iArray2DT& adjacent_comm_ID, const iArray2DT& swap,
		const ArrayT<AutoArrayT<int> >& send_nodes,
		const ArrayT<AutoArrayT<int> >& recv_nodes);

	/** perform the exchange. When called with only one
	 * argument, the data from this partition must already be in the
	 * appropriate place within the destination array. */
	/*@{*/
	/** perform double gather.
	 * \param gather source/destination array */
	void AllGather(nArray2DT<double>& gather);

	/** perform int gather.
	 * \param gather source/destination array */
	void AllGather(nArray2DT<int>& gather);

	/** specialization for gathering a single value per node. Can
	 * only be used if num_values is 1. */
	void AllGather(nArrayT<double>& gather);

	/** specialization for gathering a single value per node. Can
	 * only be used if num_values is 1. */
	void AllGather(nArrayT<int>& gather);
	/*@}*/

private:

	/** \name spatial decomposition */
	/*@{*/
	/** ID of surrounding processors needed for spatial decomposition only */
	const iArray2DT& fAdjacentCommID; /* [nsd] x {low, high} */

	/** communication patterns for shift: {+x, +y, +z,..., -x, -y, -z,...}*/
	const iArray2DT& fSwap;
	/*@}*/

	/** \name exchanged nodes
	 * send nodes for each phase of shift */
	/*@{*/
	const ArrayT<AutoArrayT<int> >& fSendNodes;
	const ArrayT<AutoArrayT<int> >& fRecvNodes;
	/*@}*/

	/** \name buffers */
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
};

/* specialization for gathering a single value per node */
inline void CartesianShiftT::AllGather(nArrayT<double>& gather) {
	nArray2DT<double> gather_2D(gather.Length(), 1, gather.Pointer());
	AllGather(gather_2D);
}

inline void CartesianShiftT::AllGather(nArrayT<int>& gather) {
	nArray2DT<int> gather_2D(gather.Length(), 1, gather.Pointer());
	AllGather(gather_2D);
}

} /* namespace Tahoe */

#endif /* _CARTESIAN_SHIFT_T_H_ */
