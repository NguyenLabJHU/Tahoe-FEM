/* $Id: PointToPointT.cpp,v 1.3 2003-11-10 22:14:39 cjkimme Exp $ */
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
	fMinorDim(0),
	fTag(sMaxTag++)
{
	sTagCount++;
}

/* destructor */
PointToPointT::~PointToPointT(void)
{
	/* try to clear any incomplete communications */
  //fComm.FreeRequests(fRecvRequest);
  //fComm.FreeRequests(fSendRequest);
	
	/* reset max tag */
	if (--sTagCount == 0) sMaxTag = 0;
}

/* allocate buffers */
void PointToPointT::Initialize(MessageT::TypeT t, int num_values)
{
	/* set type */
	fType = t;

	/* communication list */
	const iArrayT& commID = fPartition.CommID();

	/* check if allocation is already OK */
	bool exit = true;
	if (fRecvRequest.Length() != commID.Length() ||
	    fSendRequest.Length() != commID.Length()) exit = false;

	if (fType == Double) {
		for (int j = 0; j < commID.Length() && exit; j++)
			if (fdRecvBuffer[j].MinorDim() != num_values ||
			    fdSendBuffer[j].MinorDim() != num_values) exit = false;
		if (exit) return;
	} else if (fType == Integer) {
			for (int j = 0; j < commID.Length() && exit; j++)
			if (fiRecvBuffer[j].MinorDim() != num_values ||
			    fiSendBuffer[j].MinorDim() != num_values) exit = false;
		if (exit) return;
	} else
		ExceptionT::GeneralFail("PointToPointT::Initialize", "unrecognized type %d", fType);


	/* was possibly initialize before */
	if (fRecvRequest.Length() > 0) {
		fComm.FreeRequests(fRecvRequest);
		fComm.FreeRequests(fSendRequest);
	}

	/* allocate requests */
	fRecvRequest.Dimension(commID.Length());
	fSendRequest.Dimension(commID.Length());

	/* allocate buffers */
	if (fType == Double) {
		fdRecvBuffer.Dimension(commID.Length());
		fdSendBuffer.Dimension(commID.Length());
		for (int i = 0; i < commID.Length(); i++)
		{
			const iArrayT& nodes_in = *(fPartition.NodesIn(commID[i]));
			fdRecvBuffer[i].Dimension(nodes_in.Length(), num_values);
			
			const iArrayT& nodes_out = *(fPartition.NodesOut(commID[i]));
			fdSendBuffer[i].Dimension(nodes_out.Length(), num_values);
		}
		
		/* free other buffers */
		fiRecvBuffer.Free();
		fiSendBuffer.Free();
	}
	else if (fType == Integer) {
		fiRecvBuffer.Dimension(commID.Length());
		fiSendBuffer.Dimension(commID.Length());
		for (int i = 0; i < commID.Length(); i++)
		{
			const iArrayT& nodes_in = *(fPartition.NodesIn(commID[i]));
			fiRecvBuffer[i].Dimension(nodes_in.Length(), num_values);
			
			const iArrayT& nodes_out = *(fPartition.NodesOut(commID[i]));
			fiSendBuffer[i].Dimension(nodes_out.Length(), num_values);
		}

		/* free other buffers */
		fdRecvBuffer.Free();
		fdSendBuffer.Free();
	}
	else
		ExceptionT::GeneralFail("PointToPointT::Initialize", "unrecognized type %d", fType);
	
	/* keep dim */
	fMinorDim = num_values;
}

/* perform the exchange */
void PointToPointT::AllGather(nArray2DT<double>& gather)
{
	/* log received data */
	if (fComm.LogLevel() == CommunicatorT::kLow)
		fComm.Log() << " gather in:\n" << gather << endl;

	/* communication list */
	const iArrayT& commID = fPartition.CommID();

	/* post receives */
//	int rank = fComm.Rank();
	for (int i = 0; i < commID.Length(); i++)
		fComm.PostReceive(fdRecvBuffer[i], commID[i], fTag, fRecvRequest[i]);
		
	/* post sends */
	for (int i = 0; i < commID.Length(); i++)
	{
		/* destination process */
		int proc = commID[i];
	
		/* outgoing nodes */
		const iArrayT& out_nodes = *(fPartition.NodesOut(proc));
	
		/* collect outgoing data */
		dArray2DT& send = fdSendBuffer[i];
		send.RowCollect(out_nodes, gather);
	
		/* post */
		fComm.PostSend(send, proc, fTag, fSendRequest[i]);
	}
	
	/* process receives */
	for (int i = 0; i < commID.Length(); i++)
	{
		/* index of message */
		int index, source;
		fComm.WaitReceive(fRecvRequest, index, source);
		int proc = commID[index];
	
		/* incoming nodes */
		const iArrayT& in_nodes = *(fPartition.NodesIn(proc));
	
		/* process receive */
		const dArray2DT& recv = fdRecvBuffer[index];
		gather.Assemble(in_nodes, recv);
		
		/* log received data */
		if (fComm.LogLevel() == CommunicatorT::kLow)
			fComm.Log() << "received:\n" << recv << endl;
	}

	/* complete sends */
	fComm.WaitSends(fSendRequest);		

	/* log received data */
	if (fComm.LogLevel() == CommunicatorT::kLow)
		fComm.Log() << " gather out:\n" << gather << endl;
}

void PointToPointT::AllGather(nArray2DT<int>& gather)
{
	/* log received data */
	if (fComm.LogLevel() == CommunicatorT::kLow)
		fComm.Log() << " gather in:\n" << gather << endl;

	/* communication list */
	const iArrayT& commID = fPartition.CommID();

	/* post receives */
//	int rank = fComm.Rank();
	for (int i = 0; i < commID.Length(); i++)
		fComm.PostReceive(fiRecvBuffer[i], commID[i], fTag, fRecvRequest[i]);
		
	/* post sends */
	for (int i = 0; i < commID.Length(); i++)
	{
		/* destination process */
		int proc = commID[i];
	
		/* outgoing nodes */
		const iArrayT& out_nodes = *(fPartition.NodesOut(proc));
	
		/* collect outgoing data */
		iArray2DT& send = fiSendBuffer[i];
		send.RowCollect(out_nodes, gather);
	
		/* post */
		fComm.PostSend(send, proc, fTag, fSendRequest[i]);
	}
	
	/* process receives */
	for (int i = 0; i < commID.Length(); i++)
	{
		/* index of message */
		int index, source;
		fComm.WaitReceive(fRecvRequest, index, source);
		int proc = commID[index];
	
		/* incoming nodes */
		const iArrayT& in_nodes = *(fPartition.NodesIn(proc));
	
		/* process receive */
		const iArray2DT& recv = fiRecvBuffer[index];
		gather.Assemble(in_nodes, recv);
		
		/* log received data */
		if (fComm.LogLevel() == CommunicatorT::kLow)
			fComm.Log() << "received:\n" << recv << endl;
	}

	/* complete sends */
	fComm.WaitSends(fSendRequest);		

	/* log received data */
	if (fComm.LogLevel() == CommunicatorT::kLow)
		fComm.Log() << " gather out:\n" << gather << endl;
}
