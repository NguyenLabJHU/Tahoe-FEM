/* $Id: PointToPointT.h,v 1.2 2003-01-27 06:42:48 paklein Exp $ */
#ifndef _POINT_TO_POINT_T_H_
#define _POINT_TO_POINT_T_H_

/* base class */
#include "MessageT.h"

/* direct members */
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "CommunicatorT.h"

namespace Tahoe {

/* forward declarations */
class PartitionT;

/** record for performing point to point communications with either
 * blocking or non-blocking communications */
class PointToPointT: public MessageT
{
public:

	/** constructor */
	PointToPointT(CommunicatorT& comm, const PartitionT& partition);

	/** destructor */
	virtual ~PointToPointT(void);

	/** allocate buffers 
	 * \param t data type being communicated
	 * \param num_values number of values per node
	 * \param gather array which sets the data type and values per node
	 *        for the communication */
	/*@{*/
	void Initialize(MessageT::TypeT t, int num_values);
	void Initialize(nArray2DT<int>& gather);
	void Initialize(nArray2DT<double>& gather);
	/*@}*/

	/** current maximum tag number */
	static int MaxTag(void) { return sMaxTag; };

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

	/** partition information */
	const PartitionT& fPartition;

	/** \name buffers */
	/*@{*/
	int fMinorDim;
	ArrayT<dArray2DT> fdRecvBuffer;
	ArrayT<dArray2DT> fdSendBuffer;

	ArrayT<iArray2DT> fiRecvBuffer;
	ArrayT<iArray2DT> fiSendBuffer;
	/*@}*/

	/** \name data for non-blocking communications */
	/*@{*/
	ArrayT<MPI_Request> fRecvRequest;
	ArrayT<MPI_Request> fSendRequest;
	/*@}*/
	
	/** \name message tags 
	 * Message tags are generated automatically based on the number of
	 * instantiated PointToPointT objects. */
	/*@{*/
	int fTag;
	static int sMaxTag;
	static int sTagCount;
	/*@}*/
};

/* specialization for gathering a single value per node */
inline void PointToPointT::AllGather(nArrayT<double>& gather)
{
	if (fMinorDim != 1) ExceptionT::SizeMismatch("PointToPointT::AllGather");
	nArray2DT<double> gather_2D(gather.Length(), 1, gather.Pointer());
	AllGather(gather_2D);
}

inline void PointToPointT::AllGather(nArrayT<int>& gather)
{
	if (fMinorDim != 1) ExceptionT::SizeMismatch("PointToPointT::AllGather");
	nArray2DT<int> gather_2D(gather.Length(), 1, gather.Pointer());
	AllGather(gather_2D);
}

} /* namespace Tahoe */

#endif /* _POINT_TO_POINT_T_H_ */
