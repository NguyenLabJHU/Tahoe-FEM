/* $Id: PointToPointT.cpp,v 1.1.2.1 2002-12-19 03:09:14 paklein Exp $ */
#include "PointToPointT.h"
#include "CommunicatorT.h"
#include "PartitionT.h"

using namespace Tahoe;

namespace Tahoe {
int PointToPointT::sMaxTag = 0;
int PointToPointT::sTagCount = 0;
}

/* constructor */
PointToPointT::PointToPointT(CommunicatorT& comm, const PartitionT& partition):
	MessageT(comm),
	fPartition(partition),
	fTag(sMaxTag++)
{
	sTagCount++;
}

/* destructor */
PointToPointT::~PointToPointT(void)
{
	/* try to clear any incomplete communications */
	fComm.FreeRequests(fRecvRequest);
	fComm.FreeRequests(fSendRequest);
	
	/* reset max tag */
	if (--sTagCount == 0) sMaxTag = 0;
}

/* allocate buffers */
void PointToPointT::Initialize(int num_values)
{
	/* communication list */
	const iArrayT& commID = fPartition.CommID();

	/* check if allocation is already OK */
	bool exit = true;
	if (fRecvRequest.Length() != commID.Length() ||
	    fSendRequest.Length() != commID.Length()) exit = false;
	for (int j = 0; j < commID.Length() && exit; j++)
		if (fRecvBuffer[j].MinorDim() != num_values ||
		    fSendBuffer[j].MinorDim() != num_values) exit = false;
	if (exit) return;

	/* was possibly initialize before */
	if (fRecvRequest.Length() > 0) {
		fComm.FreeRequests(fRecvRequest);
		fComm.FreeRequests(fSendRequest);
	}

	/* allocate requests */
	fRecvRequest.Dimension(commID.Length());
	fSendRequest.Dimension(commID.Length());

	/* allocate buffers */
	fRecvBuffer.Dimension(commID.Length());
	fSendBuffer.Dimension(commID.Length());
	for (int i = 0; i < commID.Length(); i++)
	{
		const iArrayT& nodes_in = *(fPartition.NodesIn(commID[i]));
		fRecvBuffer[i].Dimension(nodes_in.Length(), num_values);
		
		const iArrayT& nodes_out = *(fPartition.NodesOut(commID[i]));
		fSendBuffer[i].Dimension(nodes_out.Length(), num_values);
	}
}

/* perform the exchange */
void PointToPointT::AllGather(nArrayT<double>& gather)
{
	/* communication list */
	const iArrayT& commID = fPartition.CommID();

	/* post receives */
	int rank = fComm.Rank();
	for (int i = 0; i < commID.Length(); i++)
		fComm.PostReceive(fRecvBuffer[i], rank, fTag, fRecvRequest[i]);
		
	/* post sends */
	for (int i = 0; i < commID.Length(); i++)
	{
		/* collect outgoing data */
		dArray2DT& send = fSendBuffer[i];
		//
		// -------> do it
		//
	
		/* post */
		fComm.PostSend(send, commID[i], fTag, fSendRequest[i]);
	}
	
	/* process receives */
	for (int i = 0; i < commID.Length(); i++)
	{
		/* index of message */
		int index = fComm.WaitReceive(fRecvRequest);
	
		/* process receive */
		const dArray2DT& recv = fRecvBuffer[index];
		//
		// -------> do it
		//
	}

	/* complete sends */
	fComm.WaitSends(fSendRequest);		
}
