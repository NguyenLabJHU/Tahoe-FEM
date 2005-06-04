/* $Id: CartesianShiftT.cpp,v 1.1 2005-06-04 16:58:35 paklein Exp $ */
#include "CartesianShiftT.h"
#include "CommunicatorT.h"
#include "PartitionT.h"

using namespace Tahoe;

/* constructor */
CartesianShiftT::CartesianShiftT(CommunicatorT& comm, int tag,
	const iArray2DT& adjacent_comm_ID, const iArray2DT& swap,
	const ArrayT<AutoArrayT<int> >& send_nodes,
	const ArrayT<AutoArrayT<int> >& recv_nodes):
	MessageT(comm, tag),
	fAdjacentCommID(adjacent_comm_ID),
	fSwap(swap),
	fSendNodes(send_nodes),
	fRecvNodes(recv_nodes),
	fd_send_buffer_man(fd_send_buffer),
	fd_recv_buffer_man(fd_recv_buffer),
	fi_send_buffer_man(fi_send_buffer),
	fi_recv_buffer_man(fi_recv_buffer)	
{
	const char caller[] = "CartesianShiftT::CartesianShiftT";
	if (fRecvNodes.Length() < 2 ||
		fRecvNodes.Length() != fSendNodes.Length() || 
		fAdjacentCommID.Length() != fRecvNodes.Length() ||
		fSwap.Length() != 2*fAdjacentCommID.Length())
		ExceptionT::SizeMismatch(caller, "incoming/outgoing node lists");
}

/* perform the exchange */
void CartesianShiftT::AllGather(nArray2DT<double>& gather)
{
	/* log received data */
	if (fComm.LogLevel() == CommunicatorT::kLow)
		fComm.Log() << " gather in:\n" << gather << endl;

	/* exhange with swaps */
	MPI_Request request;
	for (int i = 0; i < fSwap.MajorDim(); i++)
	{
		/* swap info */
		int dir = fSwap(i,0);
		int sgn = fSwap(i,1);

		/* neighboring process ID's */
		int ID_s = fAdjacentCommID(dir,sgn);
		int ID_r = fAdjacentCommID(dir,1-sgn);

		/* message tags */
		int tag_s = fComm.Rank();
		int tag_r = ID_r;

		/* exchange lists */
		const AutoArrayT<int>& send_nodes = fSendNodes[i];
		const AutoArrayT<int>& recv_nodes = fRecvNodes[i];
		
		/* dimensions */
		int n_s = send_nodes.Length();
		int n_r = recv_nodes.Length();		
		int n_v = gather.MinorDim();

		/* has neighboring processes different from self and non-zero exchange */
		bool do_send = (ID_s != -1 && ID_s != tag_s && n_s > 0);
		bool do_recv = (ID_r != -1 && ID_r != tag_s && n_r > 0);

		/* post receive */
		if (do_recv)
		{
			/* reset incoming buffer size */
			fd_recv_buffer_man.Dimension(n_r, n_v);

			/* post receive */
			fComm.PostReceive(fd_recv_buffer, ID_r, tag_r, request);
		}

		/* collect/send data */
		if (do_send)
		{
			/* reset incoming buffer size */
			fd_send_buffer_man.Dimension(n_s, n_v);

			/* collect values */
			fd_send_buffer.RowCollect(send_nodes, gather);
		
			/* post send */
			fComm.Send(fd_send_buffer, ID_s, tag_s);
		}

		/* complete receive */
		if (do_recv) 
		{
			/* complete receive */
			fComm.Wait(request);

			/* process incoming values */
			gather.Assemble(recv_nodes, fd_recv_buffer);
		}
	}

	/* log received data */
	if (fComm.LogLevel() == CommunicatorT::kLow)
		fComm.Log() << " gather out:\n" << gather << endl;
}

/* perform the exchange */
void CartesianShiftT::AllGather(nArray2DT<int>& gather)
{
	/* log received data */
	if (fComm.LogLevel() == CommunicatorT::kLow)
		fComm.Log() << " gather in:\n" << gather << endl;

	/* exhange with swaps */
	MPI_Request request;
	for (int i = 0; i < fSwap.MajorDim(); i++)
	{
		/* swap info */
		int dir = fSwap(i,0);
		int sgn = fSwap(i,1);

		/* neighboring process ID's */
		int ID_s = fAdjacentCommID(dir,sgn);
		int ID_r = fAdjacentCommID(dir,1-sgn);

		/* message tags */
		int tag_s = fComm.Rank();
		int tag_r = ID_r;

		/* exchange lists */
		const AutoArrayT<int>& send_nodes = fSendNodes[i];
		const AutoArrayT<int>& recv_nodes = fRecvNodes[i];
		
		/* dimensions */
		int n_s = send_nodes.Length();
		int n_r = recv_nodes.Length();		
		int n_v = gather.MinorDim();

		/* has neighboring processes different from self and non-zero exchange */
		bool do_send = (ID_s != -1 && ID_s != tag_s && n_s > 0);
		bool do_recv = (ID_r != -1 && ID_r != tag_s && n_r > 0);

		/* post receive */
		if (do_recv)
		{
			/* reset incoming buffer size */
			fi_recv_buffer_man.Dimension(n_r, n_v);

			/* post receive */
			fComm.PostReceive(fi_recv_buffer, ID_r, tag_r, request);
		}

		/* collect/send data */
		if (do_send)
		{
			/* reset incoming buffer size */
			fi_send_buffer_man.Dimension(n_s, n_v);

			/* collect values */
			fi_send_buffer.RowCollect(send_nodes, gather);
		
			/* post send */
			fComm.Send(fi_send_buffer, ID_s, tag_s);
		}

		/* complete receive */
		if (do_recv) 
		{
			/* complete receive */
			fComm.Wait(request);

			/* process incoming values */
			gather.Assemble(recv_nodes, fi_recv_buffer);
		}
	}

	/* log received data */
	if (fComm.LogLevel() == CommunicatorT::kLow)
		fComm.Log() << " gather out:\n" << gather << endl;
}
