/* $Id: PointToPointT.h,v 1.1.2.1 2002-12-19 03:09:14 paklein Exp $ */
#ifndef _POINT_TO_POINT_T_H_
#define _POINT_TO_POINT_T_H_

/* base class */
#include "MessageT.h"

/* direct members */
#include "dArray2DT.h"
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
	 * \param partition partition information
	 * \param num_values number of values per node */
	void Initialize(int num_values);

	/** current maximum tag number */
	static int MaxTag(void) { return sMaxTag; };

	/** perform the exchange */
	void AllGather(nArrayT<double>& gather);

private:

	/** partition information */
	const PartitionT& fPartition;

	/** \name buffers */
	/*@{*/
	ArrayT<dArray2DT> fRecvBuffer;
	ArrayT<dArray2DT> fSendBuffer;
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

} /* namespace Tahoe */

#endif /* _POINT_TO_POINT_T_H_ */
