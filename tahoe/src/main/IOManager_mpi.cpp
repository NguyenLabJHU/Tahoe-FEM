/* $Id: IOManager_mpi.cpp,v 1.7 2002-01-07 00:56:22 paklein Exp $ */
/* created: paklein (03/14/2000) */

#include "IOManager_mpi.h"

#include "fstreamT.h"
#include "OutputBaseT.h"
#include "OutputSetT.h"
#include "PartitionT.h"
#include "ModelFileT.h"
#include "ExodusT.h"

/* constructor */
IOManager_mpi::IOManager_mpi(ifstreamT& in, const iArrayT& io_map,
	const IOManager& local_IO, const PartitionT& partition,
	const StringT& model_file, IOBaseT::FileTypeT format):
	IOManager(in, local_IO),
	fIO_map(io_map),
	fPartition(partition)
{
	if (io_map.Length() != local_IO.ElementSets().Length())
		throw eSizeMismatch;

	/* local output sets */
	const ArrayT<OutputSetT*>& element_sets = local_IO.ElementSets();

	/* count total number of element blocks 
	 * NOTE: repeated block ID's in separate sets are counted twice */
	int num_elem_blocks = 0;	
	bool do_warning_done = false; //TEMP - write warning about joining element data
	for (int i = 0; i < element_sets.Length(); i++)
	{
		/* count (possibly with redundance) the total number of blocks */
		num_elem_blocks += element_sets[i]->NumBlocks();
	
		//TEMP - write warning about joining element data
		if (!do_warning_done && element_sets[i]->ElementOutputLabels().Length() > 0) {
			cout << "\n IOManager_mpi::IOManager_mpi: runtime assembly of element output is\n" 
			     <<   "     not implemented. Use command line option \"-split_io\"."<< endl;
			do_warning_done = true; 
		}
	}
	fBlockID.Allocate(num_elem_blocks);
	fBlockID = -1;
	fConnectivities.Allocate(num_elem_blocks);

	/* load global geometry */
	if (fIO_map.HasValue(Rank())) ReadOutputGeometry(model_file, element_sets, format);

	/* construct global output sets - all of them to preserve ID's */
	for (int i = 0; i < element_sets.Length(); i++)
	{
		const OutputSetT& set = *(element_sets[i]);

		if (fIO_map[i] == Rank())
		{
			/* set block ID's */
			const iArrayT& block_ID = set.BlockID();
			
			/* collect connectivities */
			ArrayT<const iArray2DT*> connect_list(block_ID.Length());
			for (int j = 0; j < block_ID.Length(); j++)
			{
				int index;
				if (!fBlockID.HasValue(block_ID[j], index)) {
					cout << "\n IOManager_mpi::IOManager_mpi: block ID " << block_ID[j]
					     << " should be stored but was not found" << endl;
					throw eGeneralFail;
				}
			
				/* collect */
				connect_list[j] = fConnectivities.Pointer(index);
			}

			/* construct output set */
			OutputSetT global_set(set.ID(), set.Geometry(), block_ID, connect_list,
				set.NodeOutputLabels(), set.ElementOutputLabels(), set.Changing());

			/* register */
			int IO_ID = AddElementSet(global_set);

			/* check */
			if (IO_ID != i)
			{
				cout << "\n IOManager_mpi::IOManager_mpi: expecting global I/O ID "
				     << IO_ID << " to be the\n" <<   "     same as the local I/O ID "
				     << i << endl;
				throw eGeneralFail;
			}
		}
		else
		{
			/* dummy stuff */
			iArray2DT connects;
			ArrayT<const iArray2DT*> connects_list(1);
			connects_list[0] = &connects;
			ArrayT<StringT> no_labels;

			/* construct dummy set */
			OutputSetT dummy_set(set.ID(), set.Geometry(), set.BlockID(), connects_list,
				no_labels, no_labels, false);
				
			/* register */
			int IO_ID = AddElementSet(dummy_set);
		}
	}

	/* distribute communication maps */
	SetCommunication(local_IO);
}

/* destructor */
IOManager_mpi::~IOManager_mpi(void)
{
#ifdef __MPI__
	/* free any uncompleted requests */
	Clear(fRecvRequest);
	Clear(fSendRequest);
#endif
}

#ifdef __MPI__
/* distribute/assemble/write output */
void IOManager_mpi::WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values)
{
//cout << Rank() << ":IOManager_mpi::WriteOutput" << endl;

//NOTE - skipping distribution and assembly of element
//       output for now
#pragma unused(e_values)

	/* assembling here */
	if (fIO_map[ID] == Rank())
	{
//cout << Rank() << ":collecting ID " << ID << endl;

		/* global output set */
		const OutputSetT& set = *((fOutput->ElementSets())[ID]);
	
		dArray2DT all_n_values(set.NumNodes(), set.NumNodeValues());
		dArray2DT all_e_values;

		/* assemble global output */
		const MapSetT& assembly_maps = fMapSets[ID];
		for (int i = 0; i < Size(); i++)
		{
//cout << Rank() << ":working processor " << i << endl;

			/* assembly map */
			const iArrayT& n_map = assembly_maps.NodeMap(i);
		
			if (i == Rank())
				all_n_values.Assemble(n_map, n_values);
			else if (n_map.Length() > 0)
			{
				/* incoming nodes buffer */
				dArray2DT n_values_in(n_map.Length(), all_n_values.MinorDim());		

//cout << Rank() << ":posting receive for " << n_values_in.MajorDim() << "x"
//     << n_values_in.MinorDim()  << " from " << fIO_map[ID] << endl;
		
				/* receive nodes */
				MPI_Status status;
				if (MPI_Recv(n_values_in.Pointer(), n_values_in.Length(), MPI_DOUBLE,
					i, MPI_ANY_TAG, MPI_COMM_WORLD, &status) !=
					MPI_SUCCESS) throw eMPIFail;

//cout << Rank() << ":posting received" << endl;

				/* assemble nodes */
				if (status.MPI_ERROR == MPI_SUCCESS)
					all_n_values.Assemble(n_map, n_values_in);
				else
					throw eMPIFail;
			}
		}

		/* inherited - send for output */
		IOManager::WriteOutput(ID, all_n_values, all_e_values);
	}
	/* assembling elsewhere */
	else if (fOutNodeCounts[ID] > 0)
	{
		/* check */
		if (fOutNodeCounts[ID] > n_values.MajorDim()) throw eSizeMismatch;

		/* send only values for resident nodes (assume first in sequence) */
		dArray2DT out_n_values(fOutNodeCounts[ID], n_values.MinorDim(),
			n_values.Pointer());

//cout << Rank() << ":sending ID " << ID << endl;
//cout << Rank() << ":posting send for " << out_n_values.MajorDim() << "x"
//     << out_n_values.MinorDim() << " to " << fIO_map[ID] << endl;

		/* send nodes */
		int tag = 0;
		if (MPI_Send(out_n_values.Pointer(), out_n_values.Length(), MPI_DOUBLE,
			fIO_map[ID], tag, MPI_COMM_WORLD) != MPI_SUCCESS)
			throw eMPIFail;

//cout << Rank() << ":send received" << endl;			
	}
}
#endif
	
/************************************************************************
* Private
************************************************************************/

/* communicate output counts */
void IOManager_mpi::SetCommunication(const IOManager& local_IO)
#ifdef __MPI__
{
	/* local output sets */
	const ArrayT<OutputSetT*>& element_sets = local_IO.ElementSets();
	
	/* dimensions */
	int num_proc = Size();
	int num_sets = element_sets.Length();

	/* set local counts */
	fOutNodeCounts.Allocate(num_sets);
	ArrayT<iArrayT> nodes(num_sets);

	iArrayT elem_counts(num_sets);
	ArrayT<iArrayT> elements(num_sets);
	for (int i = 0; i < num_sets; i++)
	{
		OutputSetT& set = *(element_sets[i]);
	
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

//NOTE: since element groups are comprised of multiple element blocks,
//      it is not possible to map group-processor element numbers to group-global
//      numbers here. the element groups are the only ones who can do this.
//      DISABLE ELEMENT OUTPUT
	elem_counts = 0;
		
	/* gather */
	fNodeCounts.Allocate(num_proc, num_sets);
	fElementCounts.Allocate(num_proc, num_sets);
	if (MPI_Allgather(fOutNodeCounts.Pointer(), num_sets, MPI_INT,
	                     fNodeCounts.Pointer(), num_sets, MPI_INT, MPI_COMM_WORLD)
		!= MPI_SUCCESS) throw eMPIFail;

//	if (MPI_Allgather(elem_counts.Pointer(), num_sets, MPI_INT,
//	               fElementCounts.Pointer(), num_sets, MPI_INT, MPI_COMM_WORLD)
//		!= MPI_SUCCESS) throw eMPIFail;

	/* allocate map sets */
	fMapSets.Allocate(num_sets);

	/* loop over sets */
	ArrayT<iArrayT> set_nodes(num_proc);
	for (int k = 0; k < num_sets; k++)
	{
		/* data is incoming - post receives */
		if (fIO_map[k] == Rank())
		{
//###########################################################
// using NON-blocking receives
#ifndef __SGI__
			/* allocate requests */
			int in_count = 0;
			for (int l = 0; l < num_proc; l++)
				if (fNodeCounts(l,k) > 0 && l != Rank()) in_count++;

			ArrayT<MPI_Request> r_requ(in_count);

			/* post receives */
			in_count = 0;
			for (int i = 0; i < num_proc; i++)
			{
				int num_incoming = fNodeCounts(i,k);
				if (num_incoming > 0 && i != Rank())
				{
					/* allocate receive buffer */
					set_nodes[i].Allocate(num_incoming);
				
					/* post non-blocking receives */
					if (MPI_Irecv(set_nodes[i].Pointer(), set_nodes[i].Length(),
						MPI_INT, i, k, MPI_COMM_WORLD, r_requ.Pointer(in_count++)) !=
						MPI_SUCCESS) throw eMPIFail;
				}
			}
#endif // __SGI__
// using NON-blocking receives
//###########################################################

			/* global nodes used by the set */
			OutputSetT& set = *((fOutput->ElementSets())[k]);
			const iArrayT& global_nodes_used = set.NodesUsed();

			/* allocate assembly maps */
			MapSetT& map_set = fMapSets[k];
			map_set.Allocate(num_proc, 0); // no element maps

			/* global to local map */
			int shift;
			iArrayT inv_global;
			SetInverseMap(global_nodes_used, inv_global, shift, -1);

			/* process self */		
			SetAssemblyMap(inv_global, shift, nodes[k], map_set.NodeMap(Rank()));		

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
					MPI_SUCCESS) throw eMPIFail;

				/* process receive */
				if (status.MPI_ERROR == MPI_SUCCESS)
				{
					int source = status.MPI_SOURCE;
					SetAssemblyMap(inv_global, shift, set_nodes[source], map_set.NodeMap(source));
				}
				else
				{
					cout << "\n IOManager_mpi::SetCommunication: Waitany error: "
					     << status.MPI_ERROR << " from " << status.MPI_SOURCE << endl;
					throw eMPIFail;
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
				if (num_incoming > 0 && j != Rank())			
				{
					/* allocate receive buffer */
					set_nodes[j].Allocate(num_incoming);

					/* post non-blocking receives */
					MPI_Status status;
					if (MPI_Recv(set_nodes[j].Pointer(), set_nodes[j].Length(),
						MPI_INT, j, k, MPI_COMM_WORLD, &status) !=
						MPI_SUCCESS) throw eMPIFail;
				
					/* process receive */
					SetAssemblyMap(inv_global, shift, set_nodes[j], map_set.NodeMap(j));
				}
			}
#endif // __SGI__
// using blocking receives
//###########################################################
		}
		else
		{
			int num_outgoing = nodes[k].Length();
			if (num_outgoing > 0)
			{
				/* post blocking send */
				if (MPI_Send(nodes[k].Pointer(), nodes[k].Length(),
					MPI_INT, fIO_map[k], k, MPI_COMM_WORLD) != MPI_SUCCESS)
					throw eMPIFail;			
			}
		}
		
		/* synchronize */
		if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) throw eMPIFail;		
	}
	
	/* check maps */
	CheckAssemblyMaps();
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
	nodes.Allocate(count);
	int dex = 0;
	p_local = local_set_nodes.Pointer();
	for (int j = 0; j < length; j++)
	{
		if (*p_local < cut_off)
			nodes[dex++] = to_global_nodes[*p_local];
		p_local++;
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
	inv_global.Allocate(range);
	inv_global = fill;
	for (int i = 0; i < global.Length(); i++)
		inv_global[global[i] - shift] = i;
}

void IOManager_mpi::SetAssemblyMap(const iArrayT& inv_global, int shift, const iArrayT& local,
	iArrayT& lg_map) const
{
	/* set map */
	int n_map = local.Length();
	lg_map.Allocate(n_map);
	int dex = 0;
	int*  p = local.Pointer();
	for (int j = 0; j < n_map; j++)
	{
		int dex = inv_global[*p++ - shift];
		if (dex == -1) throw eGeneralFail;
		lg_map[j] = dex;
	}	
}

/* MPI information */
int IOManager_mpi::Rank(void) const
#ifdef __MPI__
{
	int rank;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) throw eMPIFail;
	return rank;
}
#else
{ return 0; }
#endif

int IOManager_mpi::Size(void) const
#ifdef __MPI__
{
	int size;
	if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) throw eMPIFail;
	return size;
}
#else
{ return 1; }
#endif

#ifdef __MPI__
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
#endif

/* check that assembly maps are compact and complete */
void IOManager_mpi::CheckAssemblyMaps(void)
{
	/* global output sets */
	const ArrayT<OutputSetT*>& element_sets = fOutput->ElementSets();

	/* check */
	if (fIO_map.Length() != element_sets.Length()) throw eSizeMismatch;

	for (int i = 0; i < fIO_map.Length(); i++)
		if (fIO_map[i] == Rank())
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
					throw eGeneralFail;
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
							throw eGeneralFail;
						}
						else
							check = 1;
					}
				}
			
				/* redundant check */
				if (fill_check.Count(0) != 0)
				{
					cout << "\n IOManager_mpi::CheckAssemblyMaps: node maps error" << endl;
					throw eGeneralFail;
				}
			}
			
			/* check element maps */
			//TEMP - not supported yet
		}
}

/* load global geometry */
void IOManager_mpi::ReadOutputGeometry(const StringT& model_file,
	const ArrayT<OutputSetT*>& element_sets, IOBaseT::FileTypeT format)
{
	switch (format)
	{
		case IOBaseT::kTahoeII:
		{
			/* database file */
			ModelFileT file;
			file.OpenRead(model_file);
			
			/* read coordinates */
			if (file.GetCoordinates(fCoordinates) != ModelFileT::kOK) throw eGeneralFail;
		
#if 0
			/* elements */
			iArrayT element_ID;
			if (file.GetElementSetID(element_ID) != ModelFileT::kOK) throw eGeneralFail;
#endif

			/* read connectivities needed for the local output sets */
			for (int i = 0; i < fIO_map.Length(); i++)
				if (fIO_map[i] == Rank())
				{
					/* set info */
					const OutputSetT& output_set = *(element_sets[i]);
					
					/* element block ID's */
					const iArrayT& block_ID = output_set.BlockID();
					for (int j = 0; j < block_ID.Length(); j++)
					{
						/* load if not already read */
						if (!fBlockID.HasValue(block_ID[j]))
						{
							/* find empty slot */
							int index = -1;
							if (!fBlockID.HasValue(-1, index)) {
								cout << "\n IOManager_mpi::ReadOutputGeometry: no more slots connectivities:\n" 
								     << fBlockID.wrap(5) << endl;
								throw eGeneralFail;
							}
							
							/* read and store */
							fBlockID[index] = block_ID[j];
							iArray2DT& connectivities = fConnectivities[index];
							if (file.GetElementSet(block_ID[j], connectivities) != ModelFileT::kOK)
								throw eGeneralFail;
				
							/* correct numbering offset */
							connectivities--;
						}
					}
				}
			break;
		}
		case IOBaseT::kExodusII:
		{
			/* open database */
			ExodusT exo(fLog);
			exo.OpenRead(model_file);

			/* read coordinates */
			int num_nodes = exo.NumNodes();
			int num_dim   = exo.NumDimensions();
			fCoordinates.Allocate(num_nodes, num_dim);
			exo.ReadCoordinates(fCoordinates);

#if 0
			/* element ID's */
			iArrayT element_ID(exo.NumElementBlocks());
			exo.ElementBlockID(element_ID);
#endif

			/* read connectivities needed for the local output sets */
			for (int i = 0; i < fIO_map.Length(); i++)
				if (fIO_map[i] == Rank())
				{
					/* set info */
					const OutputSetT& output_set = *(element_sets[i]);
					
					/* element block ID's */
					const iArrayT& block_ID = output_set.BlockID();
					for (int j = 0; j < block_ID.Length(); j++)
					{
						/* load if not already read */
						if (!fBlockID.HasValue(block_ID[j]))
						{
							/* find empty slot */
							int index = -1;
							if (!fBlockID.HasValue(-1, index)) {
								cout << "\n IOManager_mpi::ReadOutputGeometry: no more slots connectivities:\n" 
								     << fBlockID.wrap(5) << endl;
								throw eGeneralFail;
							}
							
							/* read and store */
							fBlockID[index] = block_ID[j];
							iArray2DT& connectivities = fConnectivities[index];
							int num_elems;
							int num_elem_nodes;
							exo.ReadElementBlockDims(block_ID[j], num_elems, num_elem_nodes);
							connectivities.Allocate(num_elems, num_elem_nodes);
							GeometryT::CodeT geometry_code;
							exo.ReadConnectivities(block_ID[j], geometry_code, connectivities);

							/* correct numbering offset */
							connectivities--;
						}
					}
				}
			break;
		}
		default:
			cout << "\n IOManager_mpi::ReadOutputGeometry: format must be "
			     << IOBaseT::kTahoeII << " or " << IOBaseT::kExodusII << endl;
			throw eGeneralFail;
	}
	
	/* set global coordinates */
	SetCoordinates(fCoordinates, NULL);
}
