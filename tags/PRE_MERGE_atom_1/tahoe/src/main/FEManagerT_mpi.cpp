/* $Id: FEManagerT_mpi.cpp,v 1.27 2003-01-03 03:31:04 paklein Exp $ */
/* created: paklein (01/12/2000) */
#include "FEManagerT_mpi.h"
#include <time.h>

#include "ModelManagerT.h"
#include "AutoArrayT.h"
#include "RaggedArray2DT.h"
#include "fstreamT.h"

#include "NodeManagerT.h"
#include "ElementBaseT.h"
#include "IOManager_mpi.h"
#include "GraphT.h"
#include "IOBaseT.h"
#include "PartitionT.h"
#include "CommunicatorT.h"
#include "CommManagerT.h"

#include "ModelFileT.h"
#include "ExodusT.h"

using namespace Tahoe;

/* constructor */
FEManagerT_mpi::FEManagerT_mpi(ifstreamT& input, ofstreamT& output, 
	CommunicatorT& comm, PartitionT* partition, TaskT task):
	FEManagerT(input, output, comm),
	fPartition(partition),
	fTask(task),
	fExternIOManager(NULL)
{	
	if (fTask == kRun)
	{
		/* checks */
		if (!fPartition)
		{
			cout << "\n FEManagerT_mpi::FEManagerT_mpi: partition information required if task is "
			     << kRun << endl;
			throw ExceptionT::kBadInputValue;
		}
		else if (fPartition->ID() != Rank())
		{
			cout << "\n FEManagerT_mpi::FEManagerT_mpi: partition ID " << fPartition->ID()
			     << " does not match process rank " << Rank() << endl;
			throw ExceptionT::kMPIFail;
		}
		
		/* initial time */
		flast_time = clock();

		StringT log_file;
		log_file.Root(input.filename());
		log_file.Append(".p", Rank());
		log_file.Append(".log");
		flog.open(log_file);
		
		/* redirect log messages */
		fComm.SetLog(flog);
		
		//TEMP
		TimeStamp("FEManagerT_mpi::FEManagerT_mpi");
	}
}

/* destructor */
FEManagerT_mpi::~FEManagerT_mpi(void)
{
	//TEMP
	TimeStamp("FEManagerT_mpi::~FEManagerT_mpi");

	/* restore log messages */
	fComm.SetLog(cout);
}

ExceptionT::CodeT FEManagerT_mpi::InitStep(void)
{
	/* give heartbeat */
	fComm.Log("FEManagerT_mpi::InitStep", "init", true);

	/* set default output time stamp */
	if (fExternIOManager) fExternIOManager->SetOutputTime(Time());

	/* inherited */
	ExceptionT::CodeT error = FEManagerT::InitStep();
	if (error != ExceptionT::kNoError) {
		cout << "\n FEManagerT_mpi::InitStep: error: " << error << endl;
	}
	return error;
}

ExceptionT::CodeT FEManagerT_mpi::SolveStep(void)
{
	/* inherited */
	ExceptionT::CodeT error = FEManagerT::SolveStep();
	if (error != ExceptionT::kNoError) {
		cout << "\n FEManagerT_mpi::SolveStep: return: " << error << endl;
	}
	return error;
}

/* system relaxation */
GlobalT::RelaxCodeT FEManagerT_mpi::RelaxSystem(int group) const
{
	/* inherited */
	GlobalT::RelaxCodeT relax = FEManagerT::RelaxSystem(group);

	/* gather codes */
	ArrayT<int> all_relax(Size());
	fComm.AllGather(relax, all_relax);

	/* code precedence */
	for (int i = 0; i < all_relax.Length(); i++)
		relax = GlobalT::MaxPrecedence(relax, GlobalT::RelaxCodeT(all_relax[i]));
	
	/* report */
	if (relax != GlobalT::kNoRelax)
	{
		cout << "\n Relaxation code at time = " << Time() << '\n';
		cout << setw(kIntWidth) << "proc";	
		cout << setw(kIntWidth) << "code" << '\n';	
		for (int i = 0; i < all_relax.Length(); i++)
		{
			cout << setw(kIntWidth) << i;	
			cout << setw(kIntWidth) << all_relax[i];
			cout << '\n';	
		}
	}

	return relax;
}

/* update solution */
void FEManagerT_mpi::Update(int group, const dArrayT& update)
{
	//TEMP
	TimeStamp("FEManagerT_mpi::Update");
	
	/* check sum */
	if (fComm.Sum(ExceptionT::kNoError) != 0) 
		throw ExceptionT::kBadHeartBeat; /* must trigger try block in FEManagerT::SolveStep */

	/* inherited */
	FEManagerT::Update(group, update);
}

#if 0
/* writing results */
const dArray2DT& FEManagerT_mpi::Coordinates(void) const
{
	return fNodeManager->InitialCoordinates();
}
#endif

/* initiate the process of writing output from all output sets */
void FEManagerT_mpi::WriteOutput(double time)
{
	/* set output time for the external IO manager */
	if (fExternIOManager) fExternIOManager->SetOutputTime(time);

	/* inherited */
	FEManagerT::WriteOutput(time);
}
	
void FEManagerT_mpi::WriteOutput(int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
	/* output assembly mode */
	if (!fExternIOManager)
		/* do local IO */
		FEManagerT::WriteOutput(ID, n_values, e_values);
	else
		/* distribute/assemble/write */
		fExternIOManager->WriteOutput(ID, n_values, e_values);
}

/* (temporarily) direct output away from main out */
void FEManagerT_mpi::DivertOutput(const StringT& outfile)
{
	/* do local IO */
	if (!fExternIOManager)
	  {
		/* need processor designation for split output */
		StringT outfile_p = outfile;
		if (Size() > 1) outfile_p.Append(".p", Rank());

		FEManagerT::DivertOutput(outfile_p);
	  }
	else /* external I/O */
		fExternIOManager->DivertOutput(outfile);
}

void FEManagerT_mpi::RestoreOutput(void)
{
	/* do local IO */
	if (!fExternIOManager)
		FEManagerT::RestoreOutput();
	else /* external I/O */
		fExternIOManager->RestoreOutput();
}

/* return list of ID's of external nodes */
void FEManagerT_mpi::IncomingNodes(iArrayT& nodes_in) const
{
	if (fTask == kRun)
	{
		if (!fPartition) throw ExceptionT::kGeneralFail;
		nodes_in = fPartition->Nodes_External();
	}
}

void FEManagerT_mpi::OutgoingNodes(iArrayT& nodes_out) const
{
	if (fTask == kRun)
	{
		if (!fPartition) throw ExceptionT::kGeneralFail;
		nodes_out = fPartition->Nodes_Border();
	}
}

/* get external nodal values */
void FEManagerT_mpi::RecvExternalData(dArray2DT& external_data)
{
#ifdef __TAHOE_MPI__
	//TEMP
	//TimeStamp("FEManagerT_mpi::RecvExternalData");

	int shift = 0;
	const iArrayT& nodes_ex = fPartition->Nodes_External();
	if (nodes_ex.Length() > 0) shift = nodes_ex[0];
	//NOTE - mapping from local node number to incoming
	//       node sequence assumes sequential numbering of
	//       incoming node numbers starting at in_nodes[0]
	//NOTE - really shouldn't be communicating if the number of
	//       external nodes is zero, but right now the communicator
	//       may include processes that aren't involved
	
	/* loop until all receives completed */
	const iArrayT& commID = fPartition->CommID();
	for (int i = 0; i < commID.Length(); i++)
	{
		/* grab completed receive */
		int index;
		MPI_Status status;
		if (MPI_Waitany(fRecvRequest.Length(), fRecvRequest.Pointer(),
			&index, &status) != MPI_SUCCESS) throw ExceptionT::kMPIFail;
		
		/* process receive */
		if (status.MPI_ERROR == MPI_SUCCESS)
		{
			const iArrayT& in_nodes = *(fPartition->NodesIn(commID[index]));
			const dArray2DT&   recv = fRecvBuffer[index];

			//TEMP
			//flog << "\n incoming from: " << commID[index] << '\n';
			//flog << " number of values = " << recv.MajorDim() << '\n';
			//recv.WriteNumbered(flog);
			//flog.flush();

			/* incoming nodes are always highest numbered */
			int num_nodes = in_nodes.Length();
			for (int j = 0; j < num_nodes; j++)
				external_data.SetRow(in_nodes[j] - shift, recv(j));
		}
		else
			throw ExceptionT::kMPIFail;
	}
	
	//TEMP
	//flog << " all external data:\n";
	//external_data.WriteNumbered(flog);

	/* complete all sends */
	for (int ii = 0; ii < commID.Length(); ii++)
	{
		/* grab completed receive */
		int index;
		MPI_Status status;
		if (MPI_Waitany(fSendRequest.Length(), fSendRequest.Pointer(),
			&index, &status) != MPI_SUCCESS) throw ExceptionT::kMPIFail;
			
		if (0 && status.MPI_ERROR != MPI_SUCCESS)
		{
			flog << "\n FEManagerT_mpi::RecvExternalData: error completing send\n"
			     <<   "     from " << Rank() << " to " << commID[ii] << endl;
			throw ExceptionT::kMPIFail;
		}
	}
#else
if (external_data.Length() > 0)
{
	cout << "\n FEManagerT_mpi::RecvExternalData: invalid request for external data" << endl;
	throw ExceptionT::kGeneralFail;
}
#endif /* __TAHOE_MPI__ */
}

/* send external data */
void FEManagerT_mpi::SendExternalData(const dArray2DT& all_out_data)
{
#ifdef __TAHOE_MPI__
	//TEMP
	//TimeStamp("FEManagerT_mpi::SendExternalData");

	/* check */
	if (all_out_data.MajorDim() != fNodeManager->NumNodes())
	{
		cout << "\n FEManagerT_mpi::SendExternalData: expecting outgoing array with\n"
		     <<   "    major dimension (number of nodes) " << fNodeManager->NumNodes() << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* allocate communication buffers */
	AllocateBuffers(all_out_data.MinorDim(), fRecvBuffer, fSendBuffer);

	/* communication list */
	const iArrayT& commID = fPartition->CommID();

	/* post non-blocking receives */
	for (int j = 0; j < commID.Length(); j++)
		if (MPI_Irecv(fRecvBuffer[j].Pointer(), fRecvBuffer[j].Length(),
			MPI_DOUBLE, commID[j], MPI_ANY_TAG, fComm, &fRecvRequest[j])
			!= MPI_SUCCESS) throw ExceptionT::kMPIFail;

	/* post non-blocking sends */
	for (int i = 0; i < commID.Length(); i++)
	{
		/* outgoing nodes */
		const iArrayT& out_nodes = *(fPartition->NodesOut(commID[i]));
	
		/* outgoing data */
		fSendBuffer[i].RowCollect(out_nodes, all_out_data);

		//TEMP
		//flog << "\n sending to: " << commID[i] << '\n';
		//flog << " number of values = " << fSendBuffer[i].MajorDim() << '\n';
		//fSendBuffer[i].WriteNumbered(flog);
		//flog.flush();
					
		/* post send */
		if (MPI_Isend(fSendBuffer[i].Pointer(), fSendBuffer[i].Length(),
			MPI_DOUBLE, commID[i], Rank(), fComm, &fSendRequest[i])
			!= MPI_SUCCESS) throw ExceptionT::kMPIFail;		
	}
#else
	if (all_out_data.Length() > 0)
	{
		cout << "\n FEManagerT_mpi::SendExternalData: invalid send of external data" << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif /* __TAHOE_MPI__ */
}

void FEManagerT_mpi::SendRecvExternalData(const iArray2DT& all_out_data,
	iArray2DT& external_data)
{
#ifndef __TAHOE_MPI__
#pragma unused(all_out_data)
#pragma unused(external_data)
	cout << "\n FEManagerT_mpi::SendRecvExternalData: invalid exchange of external data" << endl;
	throw ExceptionT::kGeneralFail;
#else

	/* checks */
	if (!fPartition)
	{
		cout << "\n FEManagerT_mpi::SendRecvExternalData: invalid pointer to partition data" << endl;
		throw ExceptionT::kGeneralFail;
	}

	if (all_out_data.MajorDim() != fNodeManager->NumNodes())
	{
		cout << "\n FEManagerT_mpi::SendRecvExternalData: expecting outgoing array with\n"
		     <<   "    major dimension (number of nodes) " << fNodeManager->NumNodes() << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* communication list */
	const iArrayT& commID = fPartition->CommID();

	/* requests */
	ArrayT<MPI_Request> recv_request(commID.Length());
	ArrayT<MPI_Request> send_request(commID.Length());

	/* allocate communication buffers */
	ArrayT<iArray2DT> recv;
	ArrayT<iArray2DT> send;
	AllocateBuffers(all_out_data.MinorDim(), recv, send);

	/* post non-blocking receives */
	for (int j = 0; j < commID.Length(); j++)
		if (MPI_Irecv(recv[j].Pointer(), recv[j].Length(),
			MPI_INT, commID[j], MPI_ANY_TAG, fComm, &recv_request[j])
			!= MPI_SUCCESS) throw ExceptionT::kMPIFail;

	/* post non-blocking sends */
	for (int k = 0; k < commID.Length(); k++)
	{
		/* outgoing nodes */
		const iArrayT& out_nodes = *(fPartition->NodesOut(commID[k]));
	
		/* outgoing data */
		send[k].RowCollect(out_nodes, all_out_data);
					
		/* post send */
		if (MPI_Isend(send[k].Pointer(), send[k].Length(),
			MPI_INT, commID[k], Rank(), fComm, &send_request[k])
			!= MPI_SUCCESS) throw ExceptionT::kMPIFail;		
	}

	int shift = (fPartition->Nodes_External())[0];
	//NOTE - mapping from local node number to incoming
	//       node sequence assumes sequential numbering of
	//       incoming node numbers starting at in_nodes[0]
	
	/* loop until all receives completed */
	for (int i = 0; i < commID.Length(); i++)
	{
		/* grab completed receive */
		int index;
		MPI_Status status;
		if (MPI_Waitany(recv_request.Length(), recv_request.Pointer(),
			&index, &status) != MPI_SUCCESS) throw ExceptionT::kMPIFail;
		
		/* process receive */
		if (status.MPI_ERROR == MPI_SUCCESS)
		{
			const iArrayT& in_nodes = *(fPartition->NodesIn(commID[index]));
			const iArray2DT& incoming = recv[index];

			/* incoming nodes are always highest numbered */
			int num_nodes = in_nodes.Length();
			for (int j = 0; j < num_nodes; j++)
				external_data.SetRow(in_nodes[j] - shift, incoming(j));
		}
		else
			throw ExceptionT::kMPIFail;
	}
	
	/* complete all sends */
	for (int ii = 0; ii < commID.Length(); ii++)
	{
		/* grab completed receive */
		int index;
		MPI_Status status;
		if (MPI_Waitany(send_request.Length(), send_request.Pointer(),
			&index, &status) != MPI_SUCCESS) throw ExceptionT::kMPIFail;
			
		if (0 && status.MPI_ERROR != MPI_SUCCESS)
		{
			flog << "\n FEManagerT_mpi::RecvExternalData: error completing send\n"
			     <<   "     from " << Rank() << " to " << commID[ii] << endl;
			throw ExceptionT::kMPIFail;
		}
	}
#endif
}

/* return the local node to processor map */
void FEManagerT_mpi::NodeToProcessorMap(const iArrayT& node, iArrayT& processor) const
{
	/* initialize */
	processor.Dimension(node);
	processor = -1;
	
	/* empty list */
	if (node.Length() == 0) return;
	
	/* no parition data */
	if (!fPartition) {
		processor = Rank();
		return;
	}
	
	/* range of node numbers */
	int shift, max;
	node.MinMax(shift, max);
	int range = max - shift + 1;

	/* node to index-in-node-array map */
	iArrayT index(range);
	index = -1;
	for (int i = 0; i < range; i++)
		index[node[i] - shift] = i;
	
	/* mark external */
	const iArrayT& comm_ID = Partition().CommID();
	for (int i = 0; i < comm_ID.Length(); i++)
	{	
		int proc = comm_ID[i];
		const iArrayT* comm_nodes = Partition().NodesIn(proc);
		if (!comm_nodes) throw ExceptionT::kGeneralFail;
		for (int j = 0; j < comm_nodes->Length(); j++)
		{
			int nd = (*comm_nodes)[j] - shift;
			if (nd > -1 && nd < range) /* in the list range */
			{
				int dex = index[nd];
				if (dex != -1) /* in the list */
					processor[dex] = proc;
			}
		}
	}
	
	/* assume all others are from this proc */
	int rank = Rank();
	for (int i = 0; i < range; i++)
		if (processor[i] == -1) processor[i] = rank;
}

void FEManagerT_mpi::Wait(void)
{
	/* synchronize */
	fComm.Barrier();
}

/* domain decomposition */
void FEManagerT_mpi::Decompose(ArrayT<PartitionT>& partition, GraphT& graphU,
	bool verbose, int method)
{
	//TEMP
	//TimeStamp("FEManagerT_mpi::Decompose");

	/* check */
	if (partition.Length() == 1)
	{
		cout << "\n FEManagerT_mpi::Decompose: expecting more than 1 partition" << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* geometry file must be ascii external */
	if (fInputFormat != IOBaseT::kTahoeII && fInputFormat != IOBaseT::kExodusII)
	{
		cout << "\n FEManagerT_mpi::Decompose: requires input format with external\n";
		cout <<   "     geometry information. Use code " << IOBaseT::kTahoeII
		     << " or " << IOBaseT::kExodusII << endl;
		throw ExceptionT::kBadInputValue;
	}	

	/* decomposition method */
	bool use_new_methods = false; //TEMP
	bool dual_graph = (InterpolantDOFs() == 0);
	if (dual_graph && use_new_methods)
		DoDecompose_2(partition, graphU, verbose, method);
	else
		DoDecompose_1(partition, graphU, verbose, method);

	if (fInputFormat == IOBaseT::kTahoeII)
	{
		/* label element sets in partition data */
		for (int j = 0; j < partition.Length(); j++)
		{
			/* original model file */
			ModelFileT model_ALL;
			model_ALL.OpenRead(fModelFile);
			
			/* set number of element sets */
			iArrayT elementID;
			if (model_ALL.GetElementSetID(elementID) != ModelFileT::kOK) throw ExceptionT::kGeneralFail;
			ArrayT<StringT> IDlist(elementID.Length());
			for (int i = 0; i < IDlist.Length(); i++)
				IDlist[i].Append(elementID[i]);
			partition[j].InitElementBlocks(IDlist);	
			for (int i = 0; i < elementID.Length(); i++)
			{
				/* get element set */
				iArray2DT set;
				if (model_ALL.GetElementSet(elementID[i], set) != ModelFileT::kOK)
					throw ExceptionT::kGeneralFail;
					
				/* correct node numbering offset */
				set--;	
				
				/* set partition */
				partition[j].SetElements(IDlist[i], set);
			}
		}
	}
	else if (fInputFormat == IOBaseT::kExodusII)
	{
		/* label element sets in partition data */
		for (int j = 0; j < partition.Length(); j++)
		{
			/* original model file */
			ExodusT model_ALL(cout);
			model_ALL.OpenRead(fModelFile);
			
			/* set number of element sets */
			iArrayT elementID(model_ALL.NumElementBlocks());
			model_ALL.ElementBlockID(elementID);
			ArrayT<StringT> IDlist(elementID.Length());
			for (int i = 0; i < IDlist.Length(); i++)
				IDlist[i].Append(elementID[i]);
			partition[j].InitElementBlocks(IDlist);	
			for (int i = 0; i < elementID.Length(); i++)
			{
				/* get set dimensions */
				int num_elems;
				int num_elem_nodes;
				model_ALL.ReadElementBlockDims(elementID[i], num_elems, num_elem_nodes);

				/* get element set */
				iArray2DT set(num_elems, num_elem_nodes);
				GeometryT::CodeT geometry_code;
				model_ALL.ReadConnectivities(elementID[i], geometry_code, set);
					
				/* correct node numbering offset */
				set--;	
				
				/* set partition */
				partition[j].SetElements(IDlist[i], set);
			}
		}
	}
	else throw ExceptionT::kGeneralFail;
}

/*************************************************************************
* Protected
*************************************************************************/

void FEManagerT_mpi::ReadParameters(InitCodeT init)
{
	//TEMP
	TimeStamp("FEManagerT_mpi::ReadParameters");

	/* inherited */
	FEManagerT::ReadParameters(init);

	/* collect model file and input format from ModelManager */
	fInputFormat = fModelManager->DatabaseFormat();
	fModelFile = fModelManager->DatabaseName();
	
	/* set for parallel execution */
	if (fTask == kRun)
	{
		StringT name, suffix;
		
		/* input file name */
		name.Root(fMainIn.filename());
		suffix.Suffix(fMainIn.filename());
		name.Append(".p", Rank());
		name.Append(suffix);
		fMainIn.set_filename(name);

		/* model file name */
		suffix.Suffix(fModelFile);
		fModelFile.Root();
		fModelFile.Append(".n", Size());
		fModelFile.Append(".p", Rank());
		fModelFile.Append(suffix);
		
		/* (re-)set model manager to partial geometry file */
		if (!fModelManager->Initialize(fInputFormat, fModelFile, true)) {
			cout << "\n FEManagerT_mpi::ReadParameters: error initializing model manager" << endl;
			throw ExceptionT::kBadInputValue;
		}
		
		/* restart file name */
		if (fReadRestart)
		{
			suffix.Suffix(fRestartFile);
			fRestartFile.Root();
			fRestartFile.Append(".p", Rank());
			fRestartFile.Append(suffix);
		}
	}
}

void FEManagerT_mpi::SetNodeManager(void)
{
	//TEMP
	TimeStamp("FEManagerT_mpi::SetNodeManager");

	/* inherited */
	FEManagerT::SetNodeManager();

#ifdef __TAHOE_MPI__
	/* set for parallel execution */
	if (fTask == kRun)
	{
		/* check */
		if (fPartition->ID() < 0) throw ExceptionT::kGeneralFail;
	
		/* communication ID list */
		const iArrayT& commID = fPartition->CommID();
		fRecvRequest.Dimension(commID.Length());
		fSendRequest.Dimension(commID.Length());
	}
#endif /* __TAHOE_MPI__ */
}

void FEManagerT_mpi::SetElementGroups(void)
{
	/* inherited */
	FEManagerT::SetElementGroups();
	
//TEMP - contact not yet supported in parallel
	if (fElementGroups.HasContact())
		cout << "\n FEManagerT_mpi::SetElementGroups: WARNING: no contact between partitions"
			 << endl;
}

/* (re-)set system to initial conditions */
void FEManagerT_mpi::InitialCondition(void)
{
	/* inherited */
	FEManagerT::InitialCondition();
	
	/* set I/O */
	if (fExternIOManager) fExternIOManager->NextTimeSequence(SequenceNumber());
}

#if 0
/* reduce single value */
int FEManagerT_mpi::AllReduce(MPI_Op operation, int value)
{
#ifndef __TAHOE_MPI__
#pragma unused(operation)
#pragma unused(value)
	cout << "\n FEManagerT_mpi::AllReduce: illegal request to reduce value" << endl;
	throw ExceptionT::kGeneralFail;
#else
	int reduction = 0;
	if (MPI_Allreduce(&value, &reduction, 1, MPI_INT, operation, fComm)
		!= MPI_SUCCESS) throw ExceptionT::kMPIFail;
	return reduction;
#endif
}
#endif

/* global number of first local equation */
int FEManagerT_mpi::GetGlobalEquationStart(int group) const
{
	if (Size() == 1)
		return FEManagerT::GetGlobalEquationStart(group);
	else
	{
		/* number of local equations */
		int num_eq = fNodeManager->NumEquations(group);

		/* collect from all */
		int size = Size();
		iArrayT all_num_eq(size);
		fComm.AllGather(num_eq, all_num_eq);

		/* compute offset to local equations */
		int offset = 0;
		for (int i = 0; i < Rank(); i++)
			offset += all_num_eq[i];

		/* equation start */
		return offset + 1; //OFFSET
	}
}

int FEManagerT_mpi::GetGlobalNumEquations(int group) const
{
	if (Size() == 1)
		return FEManagerT::GetGlobalNumEquations(group);
	else
	{
		int loc_num_eq = fNodeManager->NumEquations(group);
		iArrayT all_num_eq(Size());
		fComm.AllGather(loc_num_eq, all_num_eq);
		return all_num_eq.Sum();
	}
}

/* construct a new CommManagerT */
CommManagerT* FEManagerT_mpi::New_CommManager(void) const
{
	/* inherited */
	CommManagerT* comm_man = FEManagerT::New_CommManager();
	if (!comm_man) ExceptionT::GeneralFail();

	/* set the partition data */
	comm_man->SetPartition(fPartition);
	return comm_man;
}

/*************************************************************************
 * Private
 *************************************************************************/

void FEManagerT_mpi::AllocateBuffers(int minor_dim, ArrayT<dArray2DT>& recv,
	ArrayT<dArray2DT>& send)
{
	/* check */
	if (!fPartition)
	{
		cout << "\n FEManagerT_mpi::AllocateBuffers: invalid pointer to partition" << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* communication list */
	const iArrayT& commID = fPartition->CommID();

	/* check */
	bool exit = true;
	if (recv.Length() != commID.Length() ||
	    send.Length() != commID.Length()) exit = false;
	for (int j = 0; j < commID.Length() && exit; j++)
		if (recv[j].MinorDim() != minor_dim ||
		    send[j].MinorDim() != minor_dim) exit = false;
	if (exit) return;

	/* free any existing */
	recv.Free();
	send.Free();

	/* allocate buffers */
	recv.Dimension(commID.Length());
	send.Dimension(commID.Length());
	for (int i = 0; i < commID.Length(); i++)
	{
		const iArrayT& nodes_in = *(fPartition->NodesIn(commID[i]));
		recv[i].Dimension(nodes_in.Length(), minor_dim);
		
		const iArrayT& nodes_out = *(fPartition->NodesOut(commID[i]));
		send[i].Dimension(nodes_out.Length(), minor_dim);
	}
}

void FEManagerT_mpi::AllocateBuffers(int minor_dim, ArrayT<iArray2DT>& recv,
	ArrayT<iArray2DT>& send)
{
	/* check */
	if (!fPartition)
	{
		cout << "\n FEManagerT_mpi::AllocateBuffers: invalid pointer to partition" << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* communication list */
	const iArrayT& commID = fPartition->CommID();

	/* check */
	bool exit = true;
	if (recv.Length() != commID.Length() ||
	    send.Length() != commID.Length()) exit = false;
	for (int j = 0; j < commID.Length() && exit; j++)
		if (recv[j].MinorDim() != minor_dim ||
		    send[j].MinorDim() != minor_dim) exit = false;
	if (exit) return;

	/* free any existing */
	recv.Free();
	send.Free();

	/* allocate buffers */
	recv.Dimension(commID.Length());
	send.Dimension(commID.Length());
	for (int i = 0; i < commID.Length(); i++)
	{
		const iArrayT& nodes_in = *(fPartition->NodesIn(commID[i]));
		recv[i].Dimension(nodes_in.Length(), minor_dim);
		
		const iArrayT& nodes_out = *(fPartition->NodesOut(commID[i]));
		send[i].Dimension(nodes_out.Length(), minor_dim);
	}
}

/* collect computation effort for each node */
void FEManagerT_mpi::WeightNodalCost(iArrayT& weight) const
{
	weight.Dimension(fNodeManager->NumNodes());
	weight = 1;
	fNodeManager->WeightNodalCost(weight);
	for (int i = 0 ; i < fElementGroups.Length(); i++)
		fElementGroups[i]->WeightNodalCost(weight);
}

/* write time stamp to log file */
void FEManagerT_mpi::TimeStamp(const char* message, bool flush_stream) const
{
	/* get elasped time */
	clock_t t = clock();
	double elasped_time = double(t - flast_time)/CLOCKS_PER_SEC;

	/* cast away const-ness */
	FEManagerT_mpi* tmp = (FEManagerT_mpi*) this;

	tmp->flog << " " << message << ": " << WallTime();
	tmp->flog << " elapsed time: " << elasped_time << " sec.\n\n";

#if __option(extended_errorcheck)
	if (flush_stream) tmp->flog.flush();
#else
#pragma unused(flush_stream)	
#endif

	/* store last */
	tmp->flast_time = t;
}

/* returns the time string */
const char* FEManagerT_mpi::WallTime(void) const
{
	time_t t;
	time(&t);
	return ctime(&t);
}

/* decomposition methods */
void FEManagerT_mpi::DoDecompose_2(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int
	method)
{
	/* connectivities for partititioning */
	AutoArrayT<const iArray2DT*> connects_1;
	AutoArrayT<const RaggedArray2DT<int>*> connects_2;

	/* collect element groups */
	for (int s = 0 ; s < fElementGroups.Length(); s++)
		fElementGroups[s]->ConnectsU(connects_1, connects_2);		
	
	/* dual graph partitioning graph */
	AutoArrayT<const iArray2DT*> connectsX_1;
	
	/* collect minimal connects */
	for (int s = 0 ; s < fElementGroups.Length(); s++)
		fElementGroups[s]->ConnectsX(connectsX_1);

	/* initialize graph */
	GraphT& graphX = graph;
	for (int r = 0; r < connectsX_1.Length(); r++)
		graphX.AddGroup(*(connectsX_1[r]));
		
	/* make graph */
	if (verbose) cout << " FEManagerT_mpi::DoDecompose_2: constructing dual graph" << endl;
	clock_t t0 = clock();		
	graphX.MakeGraph();
	clock_t t1 = clock();
	if (verbose)
		cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
	     << " sec: FEManagerT_mpi::DoDecompose_2: construct graph" << endl;
	
	/* generate partition */
	iArrayT config(1); //TEMP - will be scalar soon?
	config[0] = partition.Length();	
	iArrayT weight;
	WeightNodalCost(weight);
	graphX.Partition(config, weight, connects_1, connects_2, partition, true, method);
}

void FEManagerT_mpi::DoDecompose_1(ArrayT<PartitionT>& partition, GraphT& graph, 
	bool verbose, int method)
{
	/* connectivities for partititioning */
	AutoArrayT<const iArray2DT*> connects_1;
	AutoArrayT<const RaggedArray2DT<int>*> connects_2;

	/* collect connectivies from all solver groups */
	for (int i = 0; i < NumGroups(); i++)
		fNodeManager->ConnectsU(i,connects_1,connects_2);

	/* collect element groups */
	for (int s = 0 ; s < fElementGroups.Length(); s++)
		fElementGroups[s]->ConnectsU(connects_1, connects_2);

	/* initialize graph */
	GraphT& graphU = graph;
	for (int r = 0; r < connects_1.Length(); r++)
		graphU.AddGroup(*(connects_1[r]));
	for (int k = 0; k < connects_2.Length(); k++)
		graphU.AddGroup(*(connects_2[k]));
		
	/* make graph */
	clock_t t0 = clock();
	if (verbose) cout << " FEManagerT_mpi::DoDecompose_1: constructing graph" << endl;
	graphU.MakeGraph();
	clock_t t1 = clock();
	if (verbose)
		cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
		     << " sec: FEManagerT_mpi::DoDecompose_1: construct graph" << endl;
	
	/* dual graph partitioning graph */
	int dual_graph = (InterpolantDOFs() == 0) ? 1 : 0;
	AutoArrayT<const iArray2DT*> connectsX_1;
	GraphT graphX;
	if (dual_graph == 1)
	{
		if (verbose) cout << " FEManagerT_mpi::DoDecompose_1: constructing dual graph" << endl;
		
		/* collect element groups */
		for (int s = 0 ; s < fElementGroups.Length(); s++)
			fElementGroups[s]->ConnectsX(connectsX_1);

		/* initialize graph */
		for (int r = 0; r < connectsX_1.Length(); r++)
			graphX.AddGroup(*(connectsX_1[r]));
		
		/* make graph */
		clock_t t0 = clock();		
		graphX.MakeGraph();
		clock_t t1 = clock();
		if (verbose)
			cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
		     << " sec: FEManagerT_mpi::DoDecompose_1: construct X graph" << endl;
	}
	
	/* generate partition */
	iArrayT config(1); //TEMP - will be scalar soon?
	config[0] = partition.Length();	
	iArrayT weight;
	WeightNodalCost(weight);
	if (dual_graph == 1)
		graphX.Partition(config, weight, graphU, partition, true, method);
	else
		graphU.Partition(config, weight, partition, true, method);
}
