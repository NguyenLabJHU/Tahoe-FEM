/* $Id: IOManager_mpi.cpp,v 1.25 2002-11-28 17:06:31 paklein Exp $ */
/* created: paklein (03/14/2000) */
#include "IOManager_mpi.h"

#include "fstreamT.h"
#include "OutputBaseT.h"
#include "OutputSetT.h"
#include "PartitionT.h"
#include "ModelManagerT.h"

using namespace Tahoe;

/* constructor */
IOManager_mpi::IOManager_mpi(ifstreamT& in, CommunicatorT& comm,
	const iArrayT& io_map,
	const IOManager& local_IO, const PartitionT& partition,
	const StringT& model_file, IOBaseT::FileTypeT format):
	IOManager(in, local_IO),
	fComm(comm),
	fIO_map(io_map),
	fPartition(partition),
	fOutputGeometry(NULL)
{
	if (io_map.Length() != local_IO.ElementSets().Length()) {
		cout << "\n IOManager_mpi::IOManager_mpi: length of the io_map (" 
		     << io_map.Length() << ") does not\n" 
		     <<   "     match the number of output sets (" << 
		     local_IO.ElementSets().Length() << ")" << endl;
		throw ExceptionT::kSizeMismatch;	
	}

	/* local output sets */
	const ArrayT<OutputSetT*>& element_sets = local_IO.ElementSets();

	/* load global geometry */
	if (fIO_map.HasValue(fComm.Rank())) ReadOutputGeometry(model_file, element_sets, format);

//cout << fComm.Rank() << ": IOManager_mpi::IOManager_mpi: constructing output sets" << endl;

	/* construct global output sets - all of them to preserve ID's */
	for (int i = 0; i < element_sets.Length(); i++)
	{
//cout << fComm.Rank() << ": IOManager_mpi::IOManager_mpi: set: " << i << endl;

		const OutputSetT& set = *(element_sets[i]);
		
//cout << "IOManager_mpi: set " << set.ID() << ": mode = " << set.Mode() << endl;
		
		/* output over element blocks */
		if (fIO_map[i] == fComm.Rank())
		{
//cout << fComm.Rank() << ": IOManager_mpi::IOManager_mpi: constructing set" << endl;

			/* output set ID */
			int IO_ID;
			
			/* regular output set */
			if (set.Mode() == OutputSetT::kElementBlock)
			{
				/* set block ID's */
				const ArrayT<StringT>& block_ID = set.BlockID();
			
				/* collect connectivities */
				ArrayT<const iArray2DT*> connect_list(block_ID.Length());
				for (int j = 0; j < block_ID.Length(); j++)
				{
					StringT block_name;
					block_name.Append(block_ID[j]);
			
					/* collect */
					connect_list[j] = fOutputGeometry->ElementGroupPointer(block_name);
				}

				/* construct output set */
				OutputSetT global_set(set.Geometry(), block_ID, connect_list,
					set.NodeOutputLabels(), set.ElementOutputLabels(), set.Changing());

//cout << fComm.Rank() << ": IOManager_mpi::IOManager_mpi: num nodes: " << global_set.NumNodes() << endl;
//cout << fComm.Rank() << ": IOManager_mpi::IOManager_mpi: num blocks: " << global_set.NumBlocks() << endl;
//cout << fComm.Rank() << ": IOManager_mpi::IOManager_mpi: num elements: " << global_set.NumElements() << endl;

				/* register */
				IO_ID = AddElementSet(global_set);
			}
			else /* construct free set */
			{
#ifndef __TAHOE_MPI__
cout << fComm.Rank() << ": skipping output set " << set.ID() << ": global free set requires MPI" << endl;
IO_ID = i;
#else /* __TAHOE_MPI__ */
			
//cout << fComm.Rank() << ": constructing free set here: " << i << endl;			
			
				/* collect number of elements from each processor */
				iArrayT elem_count(fComm.Size());
				elem_count[fComm.Rank()] = set.NumElements();
				if (MPI_Gather(elem_count.Pointer(fComm.Rank()), 1, MPI_INT, 
					elem_count.Pointer(), 1, MPI_INT, fComm.Rank(), fComm) != MPI_SUCCESS) throw ExceptionT::kMPIFail;

//cout << fComm.Rank() << ": counts:\n" << elem_count.wrap(5) << endl;			

				/* allocate space for incoming */
				const iArray2DT& my_connects = *(set.Connectivities(set.ID()));
				iArray2DT connects(elem_count.Sum(), my_connects.MinorDim());

				/* position in buffer */
				int offset = 0;
				for (int j = 0; j < fComm.Rank(); j++)
					offset += elem_count[j];
				iArrayT send(elem_count[fComm.Rank()], connects.Pointer(offset));	
				
				/* communicated size */
				elem_count *= connects.MinorDim();
				
				/* write my connects into send buffer with global node numbering */
				if (my_connects.Length() > 0)
				{				
					/* node map */
					const iArrayT& node_map = fPartition.NodeMap();
					for (int j = 0; j < my_connects.Length(); j++)
						send[j] = node_map[my_connects[j]];

//cout << fComm.Rank() << ": local:\n" << my_connects.wrap(5) << endl;
//cout << fComm.Rank() << ": global:\n" << send.wrap(5) << endl;
				}

				/* buffer shifts */
				iArrayT displ(fComm.Size());
				displ[0] = 0;
				for (int j = 1; j < displ.Length(); j++)
					displ[j] = displ[j-1] + elem_count[j-1];

				/* collect all */
				if (MPI_Gatherv(send.Pointer(), send.Length(), MPI_INT, 
					connects.Pointer(), elem_count.Pointer(), displ.Pointer(), MPI_INT,
					fComm.Rank(), fComm) != MPI_SUCCESS) throw ExceptionT::kMPIFail;

//cout << fComm.Rank() << ": incoming:\n" << connects.wrap(5) << endl;

				/* generate dummy block ID for "connectivities" of free set nodes */
				StringT dummy_ID = "900"; /* NOTE: same convention used in JoinOutputT::SetOutput */
				dummy_ID.Append(set.ID());
					
				/* add connectivities to the output model manager */
				fOutputGeometry->RegisterElementGroup(dummy_ID, connects, set.Geometry(), true);
									
				/* construct output set */
				OutputSetT global_set(set.Geometry(), fOutputGeometry->ElementGroup(dummy_ID), 
					set.NodeOutputLabels());

//cout << fComm.Rank() << ": IOManager_mpi::IOManager_mpi: num nodes: " << global_set.NumNodes() << endl;
//cout << fComm.Rank() << ": IOManager_mpi::IOManager_mpi: num blocks: " << global_set.NumBlocks() << endl;
//cout << fComm.Rank() << ": IOManager_mpi::IOManager_mpi: num elements: " << global_set.NumElements() << endl;

				/* register */
				IO_ID = AddElementSet(global_set);
#endif /* __TAHOE_MPI__ */
			}

			/* check */
			if (IO_ID != i)
			{
				cout << "\n IOManager_mpi::IOManager_mpi: expecting global I/O ID "
				     << IO_ID << " to be the\n" <<   "     same as the local I/O ID "
				     << i << endl;
				throw ExceptionT::kGeneralFail;
			}
		}
		else
		{
			/* construct a dummy set here */
			iArray2DT connects;
			ArrayT<const iArray2DT*> connects_list(set.BlockID().Length());
			connects_list = &connects;
			ArrayT<StringT> no_labels;
			OutputSetT dummy_set(set.Geometry(), set.BlockID(), connects_list, no_labels, no_labels, false);
					
			/* register */
			int IO_ID = AddElementSet(dummy_set);

			/* construct a free set */
			if (set.Mode() == OutputSetT::kFreeSet)
			{
#ifndef __TAHOE_MPI__
cout << fComm.Rank() << ": skipping output set " << set.ID() << ": global free set requires MPI" << endl;
#else /* __TAHOE_MPI__ */			
//cout << fComm.Rank() << ": sending free set" << endl;			
			
				/* collect number of elements from each processor */
				int* dummy;
				int count = set.NumElements();
				if (MPI_Gather(&count, 1, MPI_INT, dummy, 1, MPI_INT, fIO_map[i], fComm) 
					!= MPI_SUCCESS) throw ExceptionT::kMPIFail;

				/* local connects */
				const iArray2DT& my_connects = *(set.Connectivities(set.ID()));

				/* write my connects into send buffer with global node numbering */
				iArrayT send(my_connects.Length());
				if (my_connects.Length() > 0)
				{				
					/* node map */
					const iArrayT& node_map = fPartition.NodeMap();
					for (int j = 0; j < my_connects.Length(); j++)
						send[j] = node_map[my_connects[j]];

//cout << fComm.Rank() << ": local:\n" << my_connects.wrap(5) << endl;
//cout << fComm.Rank() << ": global:\n" << send.wrap(5) << endl;
				}
					
				/* gather to processor that will write */
				if (MPI_Gatherv(send.Pointer(), send.Length(), MPI_INT, 
					NULL, NULL, NULL, MPI_INT, fIO_map[i], fComm) != MPI_SUCCESS) throw ExceptionT::kMPIFail;
#endif
			}
		}
	}
	
//cout << fComm.Rank() << ": done constructing output sets" << endl;	

	/* distribute communication maps */
	SetCommunication(local_IO);
	
	/* debugging */
	//WriteMaps(cout);
}

/* destructor */
IOManager_mpi::~IOManager_mpi(void)
{ 
	delete fOutputGeometry;
	fOutputGeometry = NULL;
}

#ifdef __TAHOE_MPI__
/* distribute/assemble/write output */
void IOManager_mpi::WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values)
{
#pragma message("IOManager_mpi::WriteOutput: fix me")

//cout << fComm.Rank() << ": IOManager_mpi::WriteOutput" << endl;
	
	/* define message tag */
	int message_tag = 0;

	/* assembling here */
	if (fIO_map[ID] == fComm.Rank())
	{
//cout << fComm.Rank() << ": collecting ID " << ID << endl;

		/* global output set */
		const OutputSetT& set = *((fOutput->ElementSets())[ID]);

		/* assembly work space */
		const MapSetT& assembly_maps = fMapSets[ID];	
		dArray2DT all_n_values(set.NumNodes(), set.NumNodeValues());
		dArray2DT all_e_values(set.NumElements(), set.NumElementValues());

/*********************************
 ***** assemble nodal values *****
 *********************************/
//cout << fComm.Rank() << ": assembling nodal values" << endl;

		/* loop over source processors to assemble global nodal output */
		for (int i = 0; i < fComm.Size(); i++)
		{
//cout << fComm.Rank() << ": working processor " << i << endl;

			/* assembly map */
			const iArrayT& n_map = assembly_maps.NodeMap(i);
			if (i == fComm.Rank())
				all_n_values.Assemble(n_map, n_values);
			else if (n_map.Length() > 0)
			{
				/* incoming nodes buffer */
				dArray2DT n_values_in(n_map.Length(), all_n_values.MinorDim());		

//cout << fComm.Rank() << ": posting receive for " << n_values_in.MajorDim() << "x"
//     << n_values_in.MinorDim()  << " from " << i << endl;
		
				/* receive nodes */
				MPI_Status status;
				if (MPI_Recv(n_values_in.Pointer(), n_values_in.Length(), MPI_DOUBLE,
					i, message_tag, fComm, &status) != MPI_SUCCESS) throw ExceptionT::kMPIFail;

//cout << fComm.Rank() << ": posting received" << endl;
//cout << n_values_in << endl;

				/* assemble nodes */
				if (status.MPI_ERROR == MPI_SUCCESS)
					all_n_values.Assemble(n_map, n_values_in);
				else {
					WriteStatus(cout, "IOManager_mpi::WriteOutput", status);
					throw ExceptionT::kMPIFail;
				}
			}
		}
/*********************************
 ***** assemble nodal values *****
 *********************************/

		/* synchronize */
		if (MPI_Barrier(fComm) != MPI_SUCCESS) throw ExceptionT::kMPIFail;

/*********************************
 **** assemble element values ****
 *********************************/

//cout << fComm.Rank() << ": assembling element values" << endl;

		/* loop over source processors to assemble global nodal output */
		for (int i = 0; i < fComm.Size(); i++)
		{
			/* assembly map */
			const iArrayT& e_map = assembly_maps.ElementMap(i);
			if (i == fComm.Rank())
				all_e_values.Assemble(e_map, e_values);
			else if (e_map.Length() > 0)
			{
				/* incoming nodes buffer */
				dArray2DT e_values_in(e_map.Length(), all_e_values.MinorDim());		

//cout << fComm.Rank() << ": posting receive for " << e_values_in.MajorDim() << "x"
//     << e_values_in.MinorDim()  << " from " << i <<  endl;
		
				/* receive nodes */
				MPI_Status status;
				if (MPI_Recv(e_values_in.Pointer(), e_values_in.Length(), MPI_DOUBLE,
					i, message_tag, fComm, &status) !=
					MPI_SUCCESS) throw ExceptionT::kMPIFail;

//cout << fComm.Rank() << ": posting received:" << endl;
//cout << e_values_in << endl;

				/* assemble nodes */
				if (status.MPI_ERROR == MPI_SUCCESS)
					all_e_values.Assemble(e_map, e_values_in);
				else {
					WriteStatus(cout, "IOManager_mpi::WriteOutput", status);
					throw ExceptionT::kMPIFail;
				}
			}
/*********************************
 **** assemble element values ****
 *********************************/
		}

		/* inherited - do actual write through local IOManager */
		IOManager::WriteOutput(ID, all_n_values, all_e_values);
	}
	else /* assembling elsewhere */
	{
		/* send nodal values */
		if (fOutNodeCounts[ID] > 0) /* assumes ID is the same as the output set index */
		{
			/* check */
			if (fOutNodeCounts[ID] > n_values.MajorDim()) throw ExceptionT::kSizeMismatch;

			/* send only values for resident nodes (assume first in sequence) */
			dArray2DT out_n_values(fOutNodeCounts[ID], n_values.MinorDim(), n_values.Pointer());

//cout << fComm.Rank() << ": sending nodal data" << endl;
//cout << fComm.Rank() << ": sending ID " << ID << endl;
//cout << fComm.Rank() << ": posting send for " << out_n_values.MajorDim() << "x"
//     << out_n_values.MinorDim() << " to " << fIO_map[ID] << ':' << endl;
//cout << out_n_values << endl;

			/* send nodes */
			if (MPI_Send(out_n_values.Pointer(), out_n_values.Length(), MPI_DOUBLE,
				fIO_map[ID], message_tag, fComm) != MPI_SUCCESS) throw ExceptionT::kMPIFail;

//cout << fComm.Rank() << ": send received" << endl;
		}
		
		/* synchronize */
		if (MPI_Barrier(fComm) != MPI_SUCCESS) throw ExceptionT::kMPIFail;

		/* send element values */
		if (fElementCounts(fComm.Rank(), ID) > 0) /* assumes ID is the same as the output set index */
		{
			/* check */
			if (fElementCounts(fComm.Rank(), ID) > e_values.MajorDim()) throw ExceptionT::kSizeMismatch;

//cout << fComm.Rank() << ": sending element data" << endl;
//cout << fComm.Rank() << ": sending ID " << ID << endl;
//cout << fComm.Rank() << ": posting send for " << e_values.MajorDim() << "x"
//     << e_values.MinorDim() << " to " << fIO_map[ID] << ':' << endl;
//cout << e_values << endl;

			/* send nodes */
			if (MPI_Send(e_values.Pointer(), e_values.Length(), MPI_DOUBLE,
				fIO_map[ID], message_tag, fComm) != MPI_SUCCESS) throw ExceptionT::kMPIFail;

//cout << fComm.Rank() << ": send received" << endl;
		}
	}
	
	/* synchronize */
	if (MPI_Barrier(fComm) != MPI_SUCCESS) throw ExceptionT::kMPIFail;
}
#endif
	
/************************************************************************
* Private
************************************************************************/

/* write the assembly maps. Used for debugging */
void IOManager_mpi::WriteMaps(ostream& out) const
{
	out << "\n IOManager_mpi::WriteMaps: number of map sets: " << fMapSets.Length() << endl;
	for (int i = 0; i < fMapSets.Length(); i++)
	{
		const MapSetT& map_set = fMapSets[i];
		cout << "\n set: " << i << '\n';

		cout << " number of node maps: " << map_set.NumNodeMaps() << '\n';
		for (int j = 0; j < map_set.NumNodeMaps(); j++)
			cout << " map: " << j << '\n' << map_set.NodeMap(j).wrap(10) << '\n';

		cout << "\n number of element maps: " << map_set.NumElementMaps() << '\n';
		for (int j = 0; j < map_set.NumElementMaps(); j++)
			cout << " map: " << j << '\n' << map_set.ElementMap(j).wrap(10) << '\n';
			
		cout.flush();
	}
}

/* communicate output counts */
void IOManager_mpi::SetCommunication(const IOManager& local_IO)
#ifdef __TAHOE_MPI__
{
//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: start" << endl;

	/* local output sets */
	const ArrayT<OutputSetT*>& local_sets = local_IO.ElementSets();
	
	/* dimensions */
	int num_proc = fComm.Size();
	int num_sets = local_sets.Length();

	/* set local counts */
	fOutNodeCounts.Dimension(num_sets);

	/* resident nodes in each set */
	ArrayT<iArrayT> nodes(num_sets);

	/* count number of comminucated nodes/element per sent
	 * on this processor. for nodal data, only _i and _b nodes
	 * are communicated (_e data is calculated elsewhere).
	 * for element data, all elements are resident (_e elements
	 * are not stored). The elements are easy because ALL elements
	 * are resident and the element maps give the block global ID's */
	iArrayT elem_counts(num_sets);
	for (int i = 0; i < num_sets; i++)
	{
		OutputSetT& set = *(local_sets[i]);
	
		/* nodal output */
		if (set.NodeOutputLabels().Length() > 0)
		{
			/* send just resident nodes */
			GlobalSetNodes(set.NodesUsed(), nodes[i]);
			fOutNodeCounts[i] = nodes[i].Length();
		}
		else
			fOutNodeCounts[i] = 0;

		/* element output */
		if (set.ElementOutputLabels().Length() > 0)
			elem_counts[i] = set.NumElements(); // need to convert to global
		else
			elem_counts[i] = 0;
	}
		
	/* gather number of commumicated nodes in each set for each processor */
	fNodeCounts.Dimension(num_proc, num_sets);
	if (MPI_Allgather(fOutNodeCounts.Pointer(), num_sets, MPI_INT,
	                     fNodeCounts.Pointer(), num_sets, MPI_INT, fComm)
		!= MPI_SUCCESS) throw ExceptionT::kMPIFail;

//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: communicated nodes counts:\n" 
//     << fNodeCounts << endl;

	/* gather number of commumicated elements in each set for each processor */
	fElementCounts.Dimension(num_proc, num_sets);
	if (MPI_Allgather(elem_counts.Pointer(), num_sets, MPI_INT,
                   fElementCounts.Pointer(), num_sets, MPI_INT, fComm)
		!= MPI_SUCCESS) throw ExceptionT::kMPIFail;

//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: communicated element counts:\n" 
//     << fElementCounts << endl;

	/* allocate map sets */
	fMapSets.Dimension(num_sets);

	/* loop over sets */
	ArrayT<iArrayT> buffer(num_proc);
	for (int k = 0; k < num_sets; k++)
	{
//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: output set: " << k << endl;

		/* assembly maps for the output set */
		MapSetT& map_set = fMapSets[k];

		/* allocate assembly maps */	
		if (fIO_map[k] == fComm.Rank()) map_set.Dimension(num_proc, num_proc);

//###########################################################
//# set nodal communication maps ############################
//###########################################################
	
		/* data is incoming - post receives */
		if (fIO_map[k] == fComm.Rank())
		{
//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: writing set" << endl;

			/* output set spec from this processor */
			OutputSetT& global_set = *((fOutput->ElementSets())[k]);

//###########################################################
// using NON-blocking receives
#ifndef __SGI__
			/* allocate requests */
			int in_count = 0;
			for (int l = 0; l < num_proc; l++)
				if (fNodeCounts(l,k) > 0 && l != fComm.Rank()) in_count++;

			ArrayT<MPI_Request> r_requ(in_count);

			/* post receives */
			in_count = 0;
			for (int i = 0; i < num_proc; i++)
			{
				int num_incoming = fNodeCounts(i,k);
				if (num_incoming > 0 && i != fComm.Rank())
				{
					/* allocate receive buffer */
					buffer[i].Dimension(num_incoming);
				
					/* post non-blocking receives */
					if (MPI_Irecv(buffer[i].Pointer(), buffer[i].Length(),
						MPI_INT, i, k, fComm, r_requ.Pointer(in_count++)) !=
						MPI_SUCCESS) throw ExceptionT::kMPIFail;
				}
			}
#endif // __SGI__
// using NON-blocking receives
//###########################################################

			/* global nodes used by the set */
			const iArrayT& global_nodes_used = global_set.NodesUsed();

//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: setting inverse map" << endl;

			/* global to local map */
			int shift;
			iArrayT inv_global;
			SetInverseMap(global_nodes_used, inv_global, shift, -1);

//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: setting assembly map" << endl;

			/* process self */		
			SetAssemblyMap(inv_global, shift, nodes[k], map_set.NodeMap(fComm.Rank()));		

//###########################################################
// using NON-blocking receives
#ifndef __SGI__
			/* generate maps for incoming */
			for (int j = 0; j < in_count; j++)
			{
				/* grab completed receive */
				int index;
				MPI_Status status;
				if (MPI_Waitany(r_requ.Length(), r_requ.Pointer(), &index, &status) !=
					MPI_SUCCESS) throw ExceptionT::kMPIFail;

				/* process receive */
				if (status.MPI_ERROR == MPI_SUCCESS)
				{
//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: received from: " << status.MPI_SOURCE << endl;

					int source = status.MPI_SOURCE;
					SetAssemblyMap(inv_global, shift, buffer[source], map_set.NodeMap(source));

//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication:  set assembly map" << endl;
				}
				else
				{
					cout << "\n IOManager_mpi::SetCommunication: Waitany error setting node maps: "
					     << status.MPI_ERROR << " from " << status.MPI_SOURCE << endl;
					throw ExceptionT::kMPIFail;
				}				
			}
#endif // __SGI__
// using NON-blocking receives
//###########################################################

//###########################################################
// using blocking receives
#ifdef __SGI__
			/* generate maps for incoming */
			for (int j = 0; j < num_proc; j++)
			{
				int num_incoming = fNodeCounts(j,k);
				if (num_incoming > 0 && j != fComm.Rank())			
				{
					/* allocate receive buffer */
					buffer[j].Dimension(num_incoming);

					/* post non-blocking receives */
					MPI_Status status;
					if (MPI_Recv(buffer[j].Pointer(), buffer[j].Length(),
						MPI_INT, j, k, fComm, &status) !=
						MPI_SUCCESS) throw ExceptionT::kMPIFail;
				
					/* process receive */
					SetAssemblyMap(inv_global, shift, buffer[j], map_set.NodeMap(j));
				}
			}
#endif // __SGI__
// using blocking receives
//###########################################################
		}
		else /* set is outgoing - post sends */
		{
//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: sending set" << endl;

			int num_outgoing = nodes[k].Length();
			if (num_outgoing > 0)
			{
				/* post blocking send */
				if (MPI_Send(nodes[k].Pointer(), nodes[k].Length(),
					MPI_INT, fIO_map[k], k, fComm) != MPI_SUCCESS)
					throw ExceptionT::kMPIFail;			
			}
		}

//###########################################################
//# set nodal communication maps ############################
//###########################################################


//###########################################################
//# set element communication maps ##########################
//###########################################################

		/* data is incoming - post receives */
		if (fIO_map[k] == fComm.Rank())
		{
			/* set spec from this processor */
			OutputSetT& global_set = *((fOutput->ElementSets())[k]);

			/* block ID's in the set */
			const ArrayT<StringT>& block_ID = global_set.BlockID();

//###########################################################
// using NON-blocking receives
#ifndef __SGI__
			/* allocate requests */
			int in_count = 0;
			for (int l = 0; l < num_proc; l++)
				if (fElementCounts(l,k) > 0 && l != fComm.Rank()) in_count++;

			ArrayT<MPI_Request> r_requ(in_count);

			/* post receives */
			in_count = 0;
			for (int i = 0; i < num_proc; i++)
			{
				int num_incoming = fElementCounts(i,k);
				if (num_incoming > 0 && i != fComm.Rank())
				{
					/* allocate receive buffer: [block length list][concat-ed id list] */
					buffer[i].Dimension(block_ID.Length() + num_incoming);

					/* post non-blocking receives */
					if (MPI_Irecv(buffer[i].Pointer(), buffer[i].Length(),
						MPI_INT, i, k, fComm, r_requ.Pointer(in_count++)) !=
						MPI_SUCCESS) throw ExceptionT::kMPIFail;
				}
			}
#endif // __SGI__
// using NON-blocking receives
//###########################################################

			/* process my output set */
			if (elem_counts[k] > 0 ) {
				iArrayT& assem_map = map_set.ElementMap(fComm.Rank());
				assem_map.Dimension(elem_counts[k]);
				assem_map = -1;
				int offset = 0;
				iArrayT block_assem_map;
				for (int j = 0; j < block_ID.Length(); j++)
				{
					const iArrayT& elem_map = fPartition.ElementMap(block_ID[j]);
					block_assem_map.Set(elem_map.Length(), assem_map.Pointer(offset));
					BuildElementAssemblyMap(k, block_ID[j], elem_map, block_assem_map);
					offset += elem_map.Length();
				}
			}

//###########################################################
// using NON-blocking receives
#ifndef __SGI__
			/* generate maps for incoming */
			for (int j = 0; j < in_count; j++)
			{
				/* grab completed receive */
				int index;
				MPI_Status status;
				if (MPI_Waitany(r_requ.Length(), r_requ.Pointer(), &index, &status) !=
					MPI_SUCCESS) throw ExceptionT::kMPIFail;

				/* process receive */
				if (status.MPI_ERROR == MPI_SUCCESS)
				{
//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: received from: " << status.MPI_SOURCE 
//     << '\n' << buffer[status.MPI_SOURCE].wrap(5) << endl;

					int source = status.MPI_SOURCE;
					const iArrayT& rbuff = buffer[source];
					
					/* process a block at a time */
					iArrayT& assem_map = map_set.ElementMap(source);
					assem_map.Dimension(fElementCounts(source, k));
					assem_map = -1;					
					iArrayT elem_map, block_assem_map;
					int offset = 0; 
					for (int i = 0; i < block_ID.Length(); i++)
					{
						elem_map.Set(rbuff[i], rbuff.Pointer(offset + block_ID.Length())); /* skip over list of block dimensions */
						block_assem_map.Set(elem_map.Length(), assem_map.Pointer(offset));
						BuildElementAssemblyMap(k, block_ID[i], elem_map, block_assem_map);
						
						/* next */
						offset += elem_map.Length();
					}
				}
				else
				{
					cout << "\n IOManager_mpi::SetCommunication: Waitany error setting element maps: "
					     << status.MPI_ERROR << " from " << status.MPI_SOURCE << endl;
					throw ExceptionT::kMPIFail;
				}				
			}
#endif // __SGI__
// using NON-blocking receives
//###########################################################

//###########################################################
// using blocking receives
#ifdef __SGI__
			/* generate maps for incoming */
			for (int j = 0; j < num_proc; j++)
			{
				int num_incoming = fElementCounts(j,k);
				if (num_incoming > 0 && j != fComm.Rank())			
				{
					/* allocate receive buffer */
					buffer[j].Dimension(block_ID.Length() + num_incoming);

					/* post non-blocking receives */
					MPI_Status status;
					if (MPI_Recv(buffer[j].Pointer(), buffer[j].Length(),
						MPI_INT, j, k, fComm, &status) !=
						MPI_SUCCESS) throw ExceptionT::kMPIFail;

//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: receiving from " << buffer[j].Length() << " from " << j << '\n'
//     << buffer[j].wrap(10) << endl;

					/* process a block at a time */
					const iArrayT& rbuff = buffer[j];
					iArrayT& assem_map = map_set.ElementMap(j);
					assem_map.Dimension(fElementCounts(j,k));
					assem_map = -1;					
					iArrayT elem_map, block_assem_map;
					int offset = 0; 
					for (int i = 0; i < block_ID.Length(); i++)
					{
						elem_map.Set(rbuff[i], rbuff.Pointer(offset + block_ID.Length())); /* skip over list of block dimensions */
						block_assem_map.Set(elem_map.Length(), assem_map.Pointer(offset));
						BuildElementAssemblyMap(k, block_ID[i], elem_map, block_assem_map);
						
						/* next */
						offset += elem_map.Length();
					}
				}
			}
#endif // __SGI__
// using blocking receives
//###########################################################
		}
		else /* set is outgoing - post sends */
		{
			if (elem_counts[k] > 0)
			{
				/* output set spec from this processor */
				OutputSetT& local_set = *(local_sets[k]);

				/* block ID's in the set */
				const ArrayT<StringT>& block_ID = local_set.BlockID();
				
				/* set send buffer */
				iArrayT sbuff(block_ID.Length() + elem_counts[k]);
				int offset = block_ID.Length();
				for (int i = 0; i < block_ID.Length(); i++)
				{
					/* block local to block global element map */
					const iArrayT& map = fPartition.ElementMap(block_ID[i]); 
				
					/* check */
					if (map.Length() != local_set.NumBlockElements(block_ID[i])) {
						cout << "\n IOManager_mpi::SetCommunication: element block ID " << block_ID[i]
						     << " dimension in output set\n" <<   "     " << k << " (" 
						     << local_set.NumBlockElements(block_ID[i])
						     << ") does not match the element map in the partition data (" 
						     << map.Length() << ")"<< endl;
						throw ExceptionT::kSizeMismatch;
					}

					sbuff[i] = map.Length();
					sbuff.CopyPart(offset, map, 0, map.Length());
					
					/* next block */
					offset += map.Length();
				}

//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: sending to " << sbuff.Length() << " to " 
//     << fIO_map[k] << '\n' << sbuff.wrap(10) << endl;

				/* post blocking send */
				if (MPI_Send(sbuff.Pointer(), sbuff.Length(),
					MPI_INT, fIO_map[k], k, fComm) != MPI_SUCCESS)
					throw ExceptionT::kMPIFail;			
			}
		}
//###########################################################
//# set element communication maps ##########################
//###########################################################
		
		/* synchronize */
		if (MPI_Barrier(fComm) != MPI_SUCCESS) throw ExceptionT::kMPIFail;		
	}

//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: checking maps" << endl;
	
	/* check maps */
	CheckAssemblyMaps();

//cout << fComm.Rank() << ": IOManager_mpi::SetCommunication: done" << endl;
}
#else
{
#pragma unused(local_IO)
	cout << "\n IOManager_mpi::SetCommunication: invalid call" << endl;
	throw;
}
#endif

/* return the global node numbers of the set nodes residing
* on the current partition */
void IOManager_mpi::GlobalSetNodes(const iArrayT& local_set_nodes, iArrayT& nodes)
{
	/* assumes local nodes numbered sequentially through _i, _b, _e */
	const iArrayT& external_nodes = fPartition.Nodes_External();
	if (external_nodes.Length() > 0)
	{
		int cut_off = external_nodes[0];

		/* count non-external nodes */
		int count = 0;
		int   length = local_set_nodes.Length();
		int* p_local = local_set_nodes.Pointer();
		for (int i = 0; i < length; i++)
			if (*p_local++ < cut_off)
				count++;
				
		/* collect (global node numbers) */
		const iArrayT& to_global_nodes = fPartition.NodeMap();
		nodes.Dimension(count);
		int dex = 0;
		p_local = local_set_nodes.Pointer();
		for (int j = 0; j < length; j++)
		{
			if (*p_local < cut_off)
				nodes[dex++] = to_global_nodes[*p_local];
			p_local++;
		}
	}
	else /* no external nodes => all nodes resident */
	{
		int length = local_set_nodes.Length();
		nodes.Dimension(length);

		/* collect (global node numbers) */
		const iArrayT& to_global_nodes = fPartition.NodeMap();
		for (int j = 0; j < length; j++)
			nodes[j] = to_global_nodes[local_set_nodes[j]];
	}
}

/* determine map from local nodes into global array, such that:
*
*             global[lg_map[i]] = local[i]
*/
void IOManager_mpi::SetInverseMap(const iArrayT& global, iArrayT& inv_global,
	int& shift, int fill) const
{
	/* compressed number range */
	int max;
	global.MinMax(shift, max);
	int range = max - shift + 1;

	/* determine (all) used nodes */
	inv_global.Dimension(range);
	inv_global = fill;
	for (int i = 0; i < global.Length(); i++)
		inv_global[global[i] - shift] = i;
}

void IOManager_mpi::SetAssemblyMap(const iArrayT& inv_global, int shift, const iArrayT& local,
	iArrayT& lg_map) const
{
	/* set map */
	int n_map = local.Length();
	lg_map.Dimension(n_map);
	int dex = 0;
	int*  p = local.Pointer();
	for (int j = 0; j < n_map; j++)
	{
		int dex = inv_global[*p++ - shift];
		if (dex == -1) {
			cout << "\n IOManager_mpi::SetAssemblyMap: an error has been detected setting the assembly maps" << endl;
			throw ExceptionT::kGeneralFail;
		}
		lg_map[j] = dex;
	}	
}

/* determine the assembly map */
void IOManager_mpi::BuildElementAssemblyMap(int set, const StringT& block_ID, 
	const iArrayT& block_map, iArrayT& map) const
{
	/* must be dimensioned */
	if (block_map.Length() != map.Length()) {
		cout << "\n IOManager_mpi::BuildElementAssemblyMap: length of block map (" 
		     << block_map.Length() << ")\n" 
		     <<   "     does not match the length of the assembly map ("
		     << map.Length() << ")"<< endl;
		throw ExceptionT::kSizeMismatch;
	}

	/* global output set */
	const OutputSetT& output_set = *((fOutput->ElementSets())[set]);
	
	/* check */
	if (output_set.Mode() != OutputSetT::kElementBlock) {
		cout << "\n IOManager_mpi::BuildElementAssemblyMap: no element assembly map unless\n" 
		     <<   "     output set mode is " << OutputSetT::kElementBlock 
		     << ": " << output_set.Mode() << endl;
		throw ExceptionT::kGeneralFail;
	}	

	/* block ID's for the current set */
	const ArrayT<StringT>& ID_list = output_set.BlockID();
	
	/* compute global offset to this block */
	int offset = 0;
	for (int i = 0; i < ID_list.Length(); i++) {
		if (ID_list[i] == block_ID) break;
		offset += output_set.NumBlockElements(ID_list[i]);
	}

	/* set the map */
	for (int i = 0; i < block_map.Length(); i++)
		map[i] = block_map[i] + offset;
}

#ifdef __TAHOE_MPI__
/* clear all outstanding requests - returns 1 of all OK */
int IOManager_mpi::Clear(ArrayT<MPI_Request>& requests)
{
	int OK = 1;
	for (int i = 0; i < requests.Length(); i++)
		if (requests[i] != MPI_REQUEST_NULL)
		{
			cout << " IOManager_mpi::Clear: cancelling request" << '\n';
				
			/* cancel request */
			MPI_Cancel(&requests[i]);
			MPI_Status status;
			MPI_Wait(&requests[i], &status);
			int flag;
			MPI_Test_cancelled(&status, &flag);
			if (flag == true)
				cout << " IOManager_mpi::Clear: cancelling request: DONE" << endl;
			else	
			{
				cout << " IOManager_mpi::Clear: cancelling request: FAIL" << endl;
				OK = 0;
			}
		}
	return OK;
}

/* write status information */
void IOManager_mpi::WriteStatus(ostream& out, const char* caller, 
	const MPI_Status& status) const
{
	out << "\n " << caller << ": MPI status returned with error\n"
        <<   "     status.MPI_SOURCE: " << status.MPI_SOURCE << '\n'
        <<   "        status.MPI_TAG: " << status.MPI_TAG << '\n'
        <<   "      status.MPI_ERROR: " << status.MPI_ERROR << endl;
}
#endif

/* check that assembly maps are compact and complete */
void IOManager_mpi::CheckAssemblyMaps(void)
{
	/* global output sets */
	const ArrayT<OutputSetT*>& element_sets = fOutput->ElementSets();

	/* check */
	if (fIO_map.Length() != element_sets.Length()) {
		cout << "\n IOManager_mpi::CheckAssemblyMaps: length of the fIO_map (" 
		     << fIO_map.Length() << ") does not\n" 
		     <<   "     match the number of output sets (" << element_sets.Length() << ")" << endl;
		throw ExceptionT::kSizeMismatch;	
	}

	for (int i = 0; i < fIO_map.Length(); i++)
		if (fIO_map[i] == fComm.Rank())
		{
			/* output set data */
			OutputSetT& set = *(element_sets[i]);
			
			/* check node maps */
			if (set.NumNodeValues() > 0)
			{			
				const iArrayT& nodes_used = set.NodesUsed();
						
				/* assembly map */
				const MapSetT& map_set = fMapSets[i];
			
				/* check overall length */
				int node_count = 0;
				for (int j = 0; j < map_set.NumNodeMaps(); j++)
					node_count += map_set.NodeMap(j).Length();
				if (node_count != nodes_used.Length())
				{
					cout << "\n IOManager_mpi::CheckAssemblyMaps: node maps size error: " << node_count
					     << " should be " << nodes_used.Length() << " for set " << i << endl;
					throw ExceptionT::kGeneralFail;
				}

				/* check fill */
				iArrayT fill_check(nodes_used.Length());
				fill_check = 0;
			
				/* check for overlap */
				for (int k = 0; k < map_set.NumNodeMaps(); k++)
				{
					const iArrayT& node_assem_map = map_set.NodeMap(k);
					for (int j = 0; j < node_assem_map.Length(); j++)
					{
						int& check = fill_check[node_assem_map[j]];
						if (check != 0)
						{
							cout << "\n IOManager_mpi::CheckAssemblyMaps: duplicated fill for node "
							     << nodes_used[node_assem_map[j]] << "\n"
							     <<   "     in assembly map " << k << " for output set ID "
							     << set.ID() << endl;
							throw ExceptionT::kGeneralFail;
						}
						else
							check = 1;
					}
				}
			
				/* redundant check */
				if (fill_check.Count(0) != 0)
				{
					cout << "\n IOManager_mpi::CheckAssemblyMaps: node maps are incomplete" << endl;
					throw ExceptionT::kGeneralFail;
				}
			}
			
			/* check element maps */
			if (set.NumElementValues() > 0)
			{			
				/* assembly map */
				const MapSetT& map_set = fMapSets[i];
			
				/* check overall length */
				int element_count = 0;
				for (int j = 0; j < map_set.NumElementMaps(); j++)
					element_count += map_set.ElementMap(j).Length();

				/* sum elements may have redudant assembly, but there should be
				 * at least as many entries in the maps as there are elements */
				int num_elements = set.NumElements(); 
				if (element_count < num_elements)
				{
					cout << "\n IOManager_mpi::CheckAssemblyMaps: element maps size error: " << element_count
					     << " should be at least " << num_elements << " for set " << i << endl;
					throw ExceptionT::kGeneralFail;
				}

				/* check fill */
				iArrayT fill_check(num_elements);
				fill_check = 0;
			
				/* check for overlap */
				for (int k = 0; k < map_set.NumElementMaps(); k++)
				{
					const iArrayT& elem_assem_map = map_set.ElementMap(k);
					for (int j = 0; j < elem_assem_map.Length(); j++)
					{
						int& check = fill_check[elem_assem_map[j]];
						check = 1;
						
						/* NOTE: should not check element maps for duplicates because elements
						 *       are generally reproduced across processor boundaries. However,
						 *       the values for each of these elements should be identical. The
						 *       nArray2DT::Assemble used to collect element values in WriteOutput
						 *       overwrites, not accumulates, values in the global array */
					}
				}
			
				/* redundant check */
				if (fill_check.Count(0) != 0)
				{
					cout << "\n IOManager_mpi::CheckAssemblyMaps: element maps are incomplete" << endl;
					throw ExceptionT::kGeneralFail;
				}
			}
		}
}

/* load global geometry */
void IOManager_mpi::ReadOutputGeometry(const StringT& model_file,
	const ArrayT<OutputSetT*>& element_sets, IOBaseT::FileTypeT format)
{
	/* initialize model manager */
	fOutputGeometry = new ModelManagerT(cout);
	if (!fOutputGeometry) throw ExceptionT::kGeneralFail;
	if (!fOutputGeometry->Initialize(format, model_file, true)) {
		cout << "\n IOManager_mpi::ReadOutputGeometry: error initializing database: " 
		     << fOutputGeometry->DatabaseName() << endl;
		throw ExceptionT::kDatabaseFail;
	}
	
	/* set global coordinates */
	SetCoordinates(fOutputGeometry->Coordinates(), NULL);

	/* read connectivities needed for the local output sets */
	for (int i = 0; i < fIO_map.Length(); i++)
	{
		/* set info */
		const OutputSetT& output_set = *(element_sets[i]);

		/* free sets are constructed in place */
		if (output_set.Mode() == OutputSetT::kElementBlock)
		{
			if (fIO_map[i] == fComm.Rank())
			{
				/* element block ID's */
				const ArrayT<StringT>& block_ID = output_set.BlockID();
			
				/* read element block */
				for (int j = 0; j < block_ID.Length(); j++)
					fOutputGeometry->ReadConnectivity(block_ID[j]);
			}
		}
	}
}
