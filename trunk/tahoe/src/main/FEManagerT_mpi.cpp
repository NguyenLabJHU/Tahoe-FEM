/* $Id: FEManagerT_mpi.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (01/12/2000)                                          */

#include "FEManagerT_mpi.h"

#include <time.h>

#ifdef __MPI__
#include "mpi.h"
#endif

#include "AutoArrayT.h"
#include "RaggedArray2DT.h"
#include "fstreamT.h"

#include "NodeManagerT.h"
#include "ElementBaseT.h"
#include "IOManager_mpi.h"
#include "GraphT.h"
#include "IOBaseT.h"
#include "PartitionT.h"
#include "ModelFileT.h"
#include "ExodusT.h"

/* MPI information */
static int rank(void)
{
#ifdef __MPI__
	int rank;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) throw eMPIFail;
	return rank;
#else
	return 0;
#endif
};

static int size(void)
{
#ifdef __MPI__
	int size;
	if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) throw eMPIFail;
	return size;
#else
	return 1;
#endif
}

/* constructor */
FEManagerT_mpi::FEManagerT_mpi(ifstreamT& input, ofstreamT& output, const PartitionT* partition,
	TaskT task):
	FEManagerT(input, output),
	fRank(rank()),
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
			throw eBadInputValue;
		}
		else if (fPartition->ID() != fRank)
		{
			cout << "\n FEManagerT_mpi::FEManagerT_mpi: partition ID " << fPartition->ID()
			     << " does not match process rank " << Rank() << endl;
			throw eMPIFail;
		}
		
		/* initial time */
		flast_time = clock();

		StringT log_file;
		log_file.Root(input.filename());
		log_file.Append(".p", fRank);
		log_file.Append(".log");
		flog.open(log_file);
		
		//TEMP
		TimeStamp("FEManagerT_mpi::FEManagerT_mpi");
	}
}

/* destructor */
FEManagerT_mpi::~FEManagerT_mpi(void)
{
	//TEMP
	TimeStamp("FEManagerT_mpi::~FEManagerT_mpi");

#ifdef __MPI__
	if (fTask == kRun)
	{
		const iArrayT& commID = fPartition->CommID();

		/* free any uncompleted receive requests */
		for (int i = 0; i < fRecvRequest.Length(); i++)
			if (fRecvRequest[i] != MPI_REQUEST_NULL)
			{
				flog << " FEManagerT_mpi::~FEManagerT_mpi: cancelling Recv from: "
				     << commID[i] << endl;
				
				/* cancel request */
				MPI_Cancel(&fRecvRequest[i]);
				MPI_Status status;
				MPI_Wait(&fRecvRequest[i], &status);
				int flag;
				MPI_Test_cancelled(&status, &flag);
				if (flag == true)
					flog << " FEManagerT_mpi::~FEManagerT_mpi: cancelling Recv from: "
					     << commID[i] << ": DONE" <<endl;
				else	
					flog << " FEManagerT_mpi::~FEManagerT_mpi: cancelling Recv from: "
					     << commID[i] << ": FAIL" <<endl;		
			}

		/* free any uncompleted send requests */
		for (int ii = 0; ii < fSendRequest.Length(); ii++)
			if (fSendRequest[ii] != MPI_REQUEST_NULL)
			{
				flog << " FEManagerT_mpi::~FEManagerT_mpi: cancelling Send from: "
				     << commID[ii] << endl;
				
				/* cancel request */
				MPI_Cancel(&fSendRequest[ii]);
				MPI_Status status;
				MPI_Wait(&fSendRequest[ii], &status);
				int flag;
				MPI_Test_cancelled(&status, &flag);
				if (flag == true)
					flog << " FEManagerT_mpi::~FEManagerT_mpi: cancelling Send from: "
					     << commID[ii] << ": DONE" <<endl;
				else	
					flog << " FEManagerT_mpi::~FEManagerT_mpi: cancelling Send from: "
					     << commID[ii] << ": FAIL" <<endl;		
			}
	}
#endif /* __MPI__ */
}

/* exception handling */
void FEManagerT_mpi::HandleException(int exception)
{
	//TEMP
	TimeStamp("FEManagerT_mpi::HandleException");

#ifdef __MPI__
	/* broadcast to running processes */
	if (exception != eNoError) AllReduce(MPI_SUM, exception);

	/* gather exceptions */
	int size = Size();
	iArrayT codes(size);
	if (MPI_Allgather(&exception, 1, MPI_INT, codes.Pointer(), 1, MPI_INT, MPI_COMM_WORLD)
		!= MPI_SUCCESS) throw eMPIFail;
			
	/* exception handling */
	if (codes.HasValue(eBadJacobianDet))
	{
		cout << "\n FEManagerT_mpi::HandleException: caught bad jacobian exception" << endl;
		
		/* override local code */
		exception = eBadJacobianDet;
	}
	/* all others are unrecoverable */
	else
	{
		/* output */
		cout << "\n FEManagerT_mpi::HandleException: unrecoverable exception: "
		     << Time() << endl;
		
		/* end job */
		exception = eGeneralFail;
	}

	/* output */
	if (Rank() == 0)
	{
		cout << setw(kIntWidth) << "proc" << "  code" << '\n';	
		for (int i = 0; i < size; i++)
			cout << setw(kIntWidth) << i << ": " << Exception(codes[i]) << '\n';	
	}
#endif /* __MPI__ */

	/* inherited */
	FEManagerT::HandleException(exception);
}

/* time sequence messaging */
bool FEManagerT_mpi::Step(void)
{
	//TEMP
	//TimeStamp("FEManagerT_mpi::Step");

	/* inherited */
	bool result = FEManagerT::Step();

//NOTE: posting sends here means they will be uncompleted at destruction
//      if one of the processes throws an exception. Cancelling the
//      uncompleted requests has been unreliable. Receives posted in
//      SendExternalData. PAK (04/05/2000)
//	if (result)
//	{
//		/* post non-blocking receives */
//		const iArrayT& commID = fPartition->CommID();
//		for (int i = 0; i < commID.Length(); i++)
//			if (MPI_Irecv(fRecvBuffer[i].Pointer(), fRecvBuffer[i].Length(),
//				MPI_DOUBLE, commID[i], MPI_ANY_TAG, MPI_COMM_WORLD, &fRecvRequest[i])
//				!= MPI_SUCCESS) throw eMPIFail;
//	}

	return result;
}

/* solution update */
void FEManagerT_mpi::Update(const dArrayT& update)
{
	//TEMP
	TimeStamp("FEManagerT_mpi::Update");

#if 0
	/* check for errors */
	if (AllReduce(MPI_SUM, eNoError) != 0) throw eNoError;
#endif

	/* inherited */
	FEManagerT::Update(update);
}

/* system relaxation */
GlobalT::RelaxCodeT FEManagerT_mpi::RelaxSystem(void) const
{
	//TEMP
	//TimeStamp("FEManagerT_mpi::RelaxSystem");

	/* inherited */
	GlobalT::RelaxCodeT relax = FEManagerT::RelaxSystem();
	
#ifdef __MPI__
	/* gather all codes */
	int size = Size();
	ArrayT<GlobalT::RelaxCodeT> all_relax(size);
	if (MPI_Allgather(&relax, 1, MPI_INT, all_relax.Pointer(), 1, MPI_INT, MPI_COMM_WORLD)
		!= MPI_SUCCESS) throw eMPIFail;
		
	/* code precedence */
	for (int i = 0; i < size; i++)
		relax = GlobalT::MaxPrecedence(relax, all_relax[i]);
	
	//TEMP
	if (relax != GlobalT::kNoRelax)
	{
		cout << "\n Relaxation code at time = " << Time() << '\n';
		cout << setw(kIntWidth) << "proc";	
		cout << setw(kIntWidth) << "code" << '\n';	
		for (int i = 0; i < size; i++)
		{
			cout << setw(kIntWidth) << i;	
			cout << setw(kIntWidth) << all_relax[i];
			cout << '\n';	
		}
	}
#endif /* __MPI__ */

	return relax;
}

/* writing results */
const dArray2DT& FEManagerT_mpi::Coordinates(void) const
{
	return fNodeManager->InitialCoordinates();
}

void FEManagerT_mpi::WriteOutput(int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
//cout << "\n FEManagerT_mpi::WriteOutput" << endl;

#ifdef __MPI__
	/* check for errors */
	if (AllReduce(MPI_SUM, eNoError) != 0) throw eNoError;
#endif

	if (!fExternIOManager)
		/* do local IO */
		FEManagerT::WriteOutput(ID, n_values, e_values);
	else
		/* distribute/assemble/write */
		fExternIOManager->WriteOutput(Time(), ID, n_values, e_values);
}

/* return list of ID's of external nodes */
void FEManagerT_mpi::IncomingNodes(iArrayT& nodes_in) const
{
	if (fTask == kRun)
	{
		if (!fPartition) throw eGeneralFail;
		nodes_in = fPartition->Nodes_External();
	}
}

void FEManagerT_mpi::OutgoingNodes(iArrayT& nodes_out) const
{
	if (fTask == kRun)
	{
		if (!fPartition) throw eGeneralFail;
		nodes_out = fPartition->Nodes_Border();
	}
}

/* get external nodal values */
void FEManagerT_mpi::RecvExternalData(dArray2DT& external_data)
{
#ifdef __MPI__
	//TEMP
	//TimeStamp("FEManagerT_mpi::RecvExternalData");

	int shift = (fPartition->Nodes_External())[0];
	//NOTE - mapping from local node number to incoming
	//       node sequence assumes sequential numbering of
	//       incoming node numbers starting at in_nodes[0]
	
	/* loop until all receives completed */
	const iArrayT& commID = fPartition->CommID();
	for (int i = 0; i < commID.Length(); i++)
	{
		/* grab completed receive */
		int index;
		MPI_Status status;
		if (MPI_Waitany(fRecvRequest.Length(), fRecvRequest.Pointer(),
			&index, &status) != MPI_SUCCESS) throw eMPIFail;
		
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
			throw eMPIFail;
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
			&index, &status) != MPI_SUCCESS) throw eMPIFail;
			
		if (0 && status.MPI_ERROR != MPI_SUCCESS)
		{
			flog << "\n FEManagerT_mpi::RecvExternalData: error completing send\n"
			     <<   "     from " << fRank << " to " << commID[ii] << endl;
			throw eMPIFail;
		}
	}
#else
if (external_data.Length() > 0)
{
	cout << "\n FEManagerT_mpi::RecvExternalData: invalid request for external data" << endl;
	throw eGeneralFail;
}
#endif /* __MPI__ */
}

/* send external data */
void FEManagerT_mpi::SendExternalData(const dArray2DT& all_out_data)
{
#ifdef __MPI__
	//TEMP
	//TimeStamp("FEManagerT_mpi::SendExternalData");

	/* check for errors */
	if (AllReduce(MPI_SUM, eNoError) != 0) throw eNoError;

	/* check */
	if (all_out_data.MajorDim() != fNodeManager->NumNodes())
	{
		cout << "\n FEManagerT_mpi::SendExternalData: expecting outgoing array with\n"
		     <<   "    major dimension (number of nodes) " << fNodeManager->NumNodes() << endl;
		throw eGeneralFail;
	}

	/* allocate communication buffers */
	AllocateBuffers(all_out_data.MinorDim(), fRecvBuffer, fSendBuffer);

	/* communication list */
	const iArrayT& commID = fPartition->CommID();

	/* post non-blocking receives */
	for (int j = 0; j < commID.Length(); j++)
		if (MPI_Irecv(fRecvBuffer[j].Pointer(), fRecvBuffer[j].Length(),
			MPI_DOUBLE, commID[j], MPI_ANY_TAG, MPI_COMM_WORLD, &fRecvRequest[j])
			!= MPI_SUCCESS) throw eMPIFail;

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
			MPI_DOUBLE, commID[i], fRank, MPI_COMM_WORLD, &fSendRequest[i])
			!= MPI_SUCCESS) throw eMPIFail;		
	}
#else
	if (all_out_data.Length() > 0)
	{
		cout << "\n FEManagerT_mpi::SendExternalData: invalid send of external data" << endl;
		throw eGeneralFail;
	}
#endif /* __MPI__ */
}

void FEManagerT_mpi::SendRecvExternalData(const iArray2DT& all_out_data,
	iArray2DT& external_data)
{
#ifndef __MPI__
#pragma unused(all_out_data)
#pragma unused(external_data)
	cout << "\n FEManagerT_mpi::SendRecvExternalData: invalid exchange of external data" << endl;
	throw eGeneralFail;
#else

	/* checks */
	if (!fPartition)
	{
		cout << "\n FEManagerT_mpi::SendRecvExternalData: invalid pointer to partition data" << endl;
		throw eGeneralFail;
	}

	if (all_out_data.MajorDim() != fNodeManager->NumNodes())
	{
		cout << "\n FEManagerT_mpi::SendRecvExternalData: expecting outgoing array with\n"
		     <<   "    major dimension (number of nodes) " << fNodeManager->NumNodes() << endl;
		throw eGeneralFail;
	}

	/* check for errors */
	if (AllReduce(MPI_SUM, eNoError) != 0) throw eNoError;

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
			MPI_INT, commID[j], MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request[j])
			!= MPI_SUCCESS) throw eMPIFail;

	/* post non-blocking sends */
	for (int k = 0; k < commID.Length(); k++)
	{
		/* outgoing nodes */
		const iArrayT& out_nodes = *(fPartition->NodesOut(commID[k]));
	
		/* outgoing data */
		send[k].RowCollect(out_nodes, all_out_data);
					
		/* post send */
		if (MPI_Isend(send[k].Pointer(), send[k].Length(),
			MPI_INT, commID[k], fRank, MPI_COMM_WORLD, &send_request[k])
			!= MPI_SUCCESS) throw eMPIFail;		
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
			&index, &status) != MPI_SUCCESS) throw eMPIFail;
		
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
			throw eMPIFail;
	}
	
	/* complete all sends */
	for (int ii = 0; ii < commID.Length(); ii++)
	{
		/* grab completed receive */
		int index;
		MPI_Status status;
		if (MPI_Waitany(send_request.Length(), send_request.Pointer(),
			&index, &status) != MPI_SUCCESS) throw eMPIFail;
			
		if (0 && status.MPI_ERROR != MPI_SUCCESS)
		{
			flog << "\n FEManagerT_mpi::RecvExternalData: error completing send\n"
			     <<   "     from " << fRank << " to " << commID[ii] << endl;
			throw eMPIFail;
		}
	}
#endif
}

void FEManagerT_mpi::Wait(void)
{
#ifdef __MPI__
	if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) throw eMPIFail;
#endif
}

/* domain decomposition */
void FEManagerT_mpi::Decompose(ArrayT<PartitionT>& partition, GraphT& graphU,
	bool verbose)
{
	//TEMP
	//TimeStamp("FEManagerT_mpi::Decompose");

	/* check */
	if (partition.Length() == 1)
	{
		cout << "\n FEManagerT_mpi::Decompose: expecting more than 1 partition" << endl;
		throw eGeneralFail;
	}

	/* geometry file must be ascii external */
	if (fInputFormat != IOBaseT::kTahoeII && fInputFormat != IOBaseT::kExodusII)
	{
		cout << "\n FEManagerT_mpi::Decompose: requires input format with external\n";
		cout <<   "     geometry information. Use code " << IOBaseT::kTahoeII
		     << " or " << IOBaseT::kExodusII << endl;
		throw eBadInputValue;
	}	

	/* connectivities for partititioning */
	AutoArrayT<const iArray2DT*> connects_1;
	AutoArrayT<const RaggedArray2DT<int>*> connects_2;

	/* collect element groups */
	for (int s = 0 ; s < fElementGroups.Length(); s++)
		fElementGroups[s]->ConnectsU(connects_1, connects_2);		

	/* initialize graph */
	for (int r = 0; r < connects_1.Length(); r++)
		graphU.AddGroup(*(connects_1[r]));
	for (int k = 0; k < connects_2.Length(); k++)
		graphU.AddGroup(*(connects_2[k]));
		
	/* make graph */
	clock_t t0 = clock();		
	graphU.MakeGraph();
	clock_t t1 = clock();
	if (verbose)
		cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
		     << " sec: FEManagerT_mpi::Decompose: construct graph" << endl;
	
	/* dual graph partitioning graph */
	int dual_graph = (InterpolantDOFs() == 0) ? 1 : 0;
	AutoArrayT<const iArray2DT*> connectsX_1;
	GraphT graphX;
	if (dual_graph == 1)
	{
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
		     << " sec: FEManagerT_mpi::Decompose: construct X graph" << endl;
	}
	
	/* generate partition */
	iArrayT config(1); //TEMP - will be scalar soon?
	config[0] = partition.Length();	
	iArrayT weight;
	WeightNodalCost(weight);
	if (dual_graph == 1)
		graphX.Partition(config, weight, graphU, partition, true);
	else
		graphU.Partition(config, weight, partition, true);

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
			if (model_ALL.GetElementSetID(elementID) != ModelFileT::kOK) throw eGeneralFail;
			partition[j].InitElementBlocks(elementID);	
			for (int i = 0; i < elementID.Length(); i++)
			{
				/* get element set */
				iArray2DT set;
				if (model_ALL.GetElementSet(elementID[i], set) != ModelFileT::kOK)
					throw eGeneralFail;
					
				/* correct node numbering offset */
				set--;	
				
				/* set partition */
				partition[j].SetElements(elementID[i], set);
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
			partition[j].InitElementBlocks(elementID);	
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
				partition[j].SetElements(elementID[i], set);
			}
		}
	}
	else throw eGeneralFail;

//TEMP
//#ifdef __MACOS__
#if 0
	cout << "\n FEManagerT_mpi::Decompose: writing graph to output\n" << endl;

	fMainOut << "//######################### U con ##########################\n";
	for (int i = 0; i < connects_2.Length(); i++)
	{
		fMainOut << "\n group: " << i+1 << '\n';
		iArrayT tmp(connects_2[i]->Length(), connects_2[i]->Pointer());
		tmp++;
		connects_2[i]->WriteNumbered(fMainOut);
		tmp--;
	}
	fMainOut << "//######################### U con ##########################\n";

	fMainOut << "//######################### graph ##########################\n";
	graph.Write(fMainOut);
	fMainOut << "//######################### graph ##########################\n";
	fMainOut.flush();

	/* write partition data */
	for (int q = 0; q < partition.Length(); q++)
	{
		StringT file_name;
		file_name.Root(fModelFile);
		file_name.Append(".n", partition.Length());
		file_name.Append(".global");
		file_name.Append(".part", q);
		
		ofstream out_q(file_name);
		out_q << "# data for partition: " << q << '\n';
		out_q << partition[q] << '\n';
		out_q.close();
	}
#endif
	
	/* write partition data */
	for (int q = 0; q < partition.Length(); q++)
	{
		/* set to local scope */
		partition[q].SetScope(PartitionT::kLocal);

		StringT file_name;
		file_name.Root(fModelFile);
		file_name.Append(".n", partition.Length());
		file_name.Append(".part", q);
		
		ofstream out_q(file_name);
		out_q << "# data for partition: " << q << '\n';
		out_q << partition[q] << '\n';
		out_q.close();
	}
}

/* basic MP support */
int FEManagerT_mpi::Rank(void) const
{
#ifdef __MPI__
	return rank();
#else
	return FEManagerT::Rank();
#endif
}

int FEManagerT_mpi::Size(void) const
{
#ifdef __MPI__
	return size();
#else
	return FEManagerT::Size();
#endif
}

/*************************************************************************
* Protected
*************************************************************************/

void FEManagerT_mpi::ReadParameters(void)
{
	//TEMP
	TimeStamp("FEManagerT_mpi::ReadParameters");

	/* inherited */
	FEManagerT::ReadParameters();
	
	/* set for parallel execution */
	if (fTask == kRun)
	{
		StringT name, suffix;
		
		/* input file name */
		name.Root(fMainIn.filename());
		suffix.Suffix(fMainIn.filename());
		name.Append(".p", fRank);
		name.Append(suffix);
		fMainIn.set_filename(name);

		/* model file name */
		suffix.Suffix(fModelFile);
		fModelFile.Root();
		fModelFile.Append(".n", Size());
		fModelFile.Append(".p", fRank);
		fModelFile.Append(suffix);
		
		/* restart file name */
		if (fReadRestart)
		{
			suffix.Suffix(fRestartFile);
			fRestartFile.Root();
			fRestartFile.Append(".p", fRank);
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

#ifdef __MPI__
	/* set for parallel execution */
	if (fTask == kRun)
	{
		/* check */
		if (fPartition->ID() < 0) throw eGeneralFail;
	
		/* communication ID list */
		const iArrayT& commID = fPartition->CommID();
		fRecvRequest.Allocate(commID.Length());
		fSendRequest.Allocate(commID.Length());
	}
#endif /* __MPI__ */
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

/* reduce single value */
int FEManagerT_mpi::AllReduce(MPI_Op operation, int value)
{
#ifndef __MPI__
#pragma unused(operation)
#pragma unused(value)
	cout << "\n FEManagerT_mpi::AllReduce: illegal request to reduce value" << endl;
	throw eGeneralFail;
#else
	int reduction = 0;
	if (MPI_Allreduce(&value, &reduction, 1, MPI_INT, operation, MPI_COMM_WORLD)
		!= MPI_SUCCESS) throw eMPIFail;
	return reduction;
#endif
}

/* global number of first local equation */
int FEManagerT_mpi::GetGlobalEquationStart(void) const
{
#ifdef __MPI__
	/* number of local equations */
	int num_eq = fNodeManager->NumEquations();

	/* collect from all */
	int size = Size();
	iArrayT all_num_eq(size);
	if (MPI_Allgather(&num_eq, 1, MPI_INT, all_num_eq.Pointer(), 1,
		MPI_INT, MPI_COMM_WORLD) != MPI_SUCCESS) throw eMPIFail;

	/* compute offset to local equations */
	int offset = 0;
	for (int i = 0; i < fRank; i++)
		offset += all_num_eq[i];

	/* equation start */
	return offset + 1; //OFFSET

#else
	return FEManagerT::GetGlobalEquationStart();
#endif /* __MPI__ */
}

int FEManagerT_mpi::GetGlobalNumEquations(void) const
{
#ifdef __MPI__
	int loc_num_eq = fNodeManager->NumEquations();
	iArrayT all_num_eq(Size());
	if (MPI_Allgather(&loc_num_eq, 1, MPI_INT, all_num_eq.Pointer(), 1,
		MPI_INT, MPI_COMM_WORLD) != MPI_SUCCESS) throw eMPIFail;
	return all_num_eq.Sum();
#else
	/* inherited */
	return FEManagerT::GetGlobalNumEquations();
#endif /* __MPI__ */
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
		throw eGeneralFail;
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
	recv.Allocate(commID.Length());
	send.Allocate(commID.Length());
	for (int i = 0; i < commID.Length(); i++)
	{
		const iArrayT& nodes_in = *(fPartition->NodesIn(commID[i]));
		recv[i].Allocate(nodes_in.Length(), minor_dim);
		
		const iArrayT& nodes_out = *(fPartition->NodesOut(commID[i]));
		send[i].Allocate(nodes_out.Length(), minor_dim);
	}
}

void FEManagerT_mpi::AllocateBuffers(int minor_dim, ArrayT<iArray2DT>& recv,
	ArrayT<iArray2DT>& send)
{
	/* check */
	if (!fPartition)
	{
		cout << "\n FEManagerT_mpi::AllocateBuffers: invalid pointer to partition" << endl;
		throw eGeneralFail;
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
	recv.Allocate(commID.Length());
	send.Allocate(commID.Length());
	for (int i = 0; i < commID.Length(); i++)
	{
		const iArrayT& nodes_in = *(fPartition->NodesIn(commID[i]));
		recv[i].Allocate(nodes_in.Length(), minor_dim);
		
		const iArrayT& nodes_out = *(fPartition->NodesOut(commID[i]));
		send[i].Allocate(nodes_out.Length(), minor_dim);
	}
}

/* collect computation effort for each node */
void FEManagerT_mpi::WeightNodalCost(iArrayT& weight) const
{
	weight.Allocate(fNodeManager->NumNodes());
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
