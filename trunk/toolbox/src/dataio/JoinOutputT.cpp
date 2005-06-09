/* $Id: JoinOutputT.cpp,v 1.24 2005-06-09 16:18:38 paklein Exp $ */
/* created: paklein (03/24/2000) */
#include "JoinOutputT.h"

#include "ifstreamT.h"
#include "OutputBaseT.h"
#include "ModelManagerT.h"
#include "OutputSetT.h"
#include "StringT.h"
#include "iArray2DT.h"
#include "nVariArray2DT.h"
#include "VariArrayT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
JoinOutputT::JoinOutputT(const StringT& param_file, const StringT& model_file,
	IOBaseT::FileTypeT model_file_type, IOBaseT::FileTypeT results_file_type, 
	OutputBaseT* output, int size):
	fJobFile(param_file),
	fResultsFileType(results_file_type),
	fPartitions(size),
	fOutput(output)
{
	const char caller[] = "JoinOutputT::JoinOutputT";

	/* set model database manager */
	fModel = new ModelManagerT(cout);
	if (!fModel->Initialize(model_file_type, model_file, true))
		ExceptionT::DatabaseFail(caller, "error opening model file: %s",
			fModel->DatabaseName().Pointer());

	/* read partition data */
	for (int i = 0; i < fPartitions.Length(); i++)
	{
		/* file name */
		StringT file;
		file.Root(model_file);
		file.Append(".n", fPartitions.Length());
		file.Append(".part", i);

		/* open stream */
		ifstreamT part_in(file);
		if (!part_in.is_open())
			ExceptionT::GeneralFail(caller, "could not open decomposition file: %s",
				part_in.filename());

		/* read data */
		part_in >> fPartitions[i];
		
		/* set numbering scope */
		fPartitions[i].SetScope(PartitionT::kLocal);
	}
	
	/* set output */
	SetOutput();
	
	/* set assembly maps */
	SetMaps();
}
	
/* destructor */
JoinOutputT::~JoinOutputT(void)
{
	delete fModel;
	fModel = NULL;
}

/* do join */
void JoinOutputT::Join(void)
{
	/* assembly work space */
	dArray2DT part_n_values;
	dArray2DT part_e_values;
	nVariArray2DT<double> part_n_man(0, part_n_values, 0);
	nVariArray2DT<double> part_e_man(0, part_e_values, 0);	

	/* output sets data */
	const ArrayT<OutputSetT*>& element_sets = fOutput->ElementSets();
	for (int i = 0; i < element_sets.Length(); i++)
	{
		cout << "\n JoinOutputT: output set: " << i+1 << endl;
	
		/* set data */
		const OutputSetT& output_set = *(element_sets[i]);
		int io_ID = atoi(output_set.ID());
		const MapSetT& map_set = fMapSets[i];
	
		/* assembled values */
		dArray2DT all_n_values(output_set.NumNodes(), output_set.NumNodeValues());
		dArray2DT all_e_values;
		if (output_set.BlockID().Length() == 0 && output_set.NumElementValues() > 0) //outdated - this should not occur
			cout << "\n JoinOutputT::Join: skipping element output\n" << endl;
		else
			all_e_values.Dimension(output_set.NumElements(), output_set.NumElementValues());

		/* non-empty output */
		if (all_n_values.Length() > 0 || all_e_values.Length() > 0)
		{
			/* number of output steps */
			int num_steps = NumOutputSteps(i);
			cout << " JoinOutputT:      steps: " << num_steps;
			if (output_set.Changing()) cout << "+ (changing geometry)";
			cout << endl;
			int d_width = cout.precision() + kDoubleExtra;
			cout << setw(kIntWidth) << "step"
			     << setw(d_width)   << "time" << '\n';
			
			/* loop over steps */
			for (int j = 0; j < num_steps; j++)
			{
				/* output step index */
				int step_index = (output_set.Changing()) ? 0 : j;
			
				/* initialize */
				StringT filename;
				all_n_values = 0.0;
				all_e_values = 0.0;
			
				/* loop over partitions */
				double time = 0.0;
				bool found_time = false;
				for (int k = 0; k < fPartitions.Length(); k++)
				{
					/* file name */
					ResultFileName(k, io_ID, output_set.Changing(), j, filename);
					
					/* check if file is present */
					if (fstreamT::Exists(filename))
					{
						/* open the database file */
						ModelManagerT results(cout);
						if (!results.Initialize(fResultsFileType, filename, true)) {
							cout << "\n JoinOutputT::Join: error opening partial results file \""
							     << results.DatabaseName() << '\"' << endl;
							throw ExceptionT::kDatabaseFail;
						}

						/* get time */
						if (!found_time)
						{
							found_time = true;
							dArrayT steps;
							results.TimeSteps(steps);
							time = steps[step_index];
							cout << setw(kIntWidth) << j+1
							     << setw(d_width)   << time << endl;
						}
					
						/* assemble nodal values */
						if (all_n_values.Length() > 0)
						{
							/* dimensions */
							int num_nodes = results.NumNodes();
						
							/* assembly map: output_set_ID[partition_ID] */
							const iArrayT& node_map = map_set.NodeMap(k);

							/* very weak consistency check */
							if (node_map.Length() > num_nodes)
							{
								cout << "\n JoinOutputT::Join: assembly map of nodal values (" << node_map.Length() 
								     << ") is longer than\n" 
								     <<   "     is longer than the number of nodes (" << num_nodes
								     << ") in partial results file:\n"
								     <<   "     " << results.DatabaseName() << endl;
								throw ExceptionT::kSizeMismatch;
							}
							
							/* set work space */
							part_n_man.Dimension(num_nodes, all_n_values.MinorDim());

							/* read data */
							results.AllNodeVariables(step_index, part_n_values);

							/* assemble */
							all_n_values.Assemble(node_map, part_n_values);
						}					

						/* assemble element values */
						if (all_e_values.Length() > 0)
						{
							/* element map: global_output_block_ID[partition_output_block_ID] */
							const iArrayT& element_map = map_set.ElementMap(k);

							/* set work space */
							part_e_man.Dimension(element_map.Length(), all_e_values.MinorDim());

							/* block ID's in the set */
							const ArrayT<StringT>& block_ID = output_set.BlockID();
							
							/* read data by block - one block after the next */
							int row_offset = 0;
							for (int l = 0; l < block_ID.Length(); l++)
							{
								/* block name */
								StringT block_name;
								block_name.Append(block_ID[l]);
								
								/* block dimensions */
								int nel, nen;
								results.ElementGroupDimensions(block_name, nel, nen);

								/* weak check */
								if (nel > element_map.Length())
								{
									cout << "\n JoinOutputT::Join: number of elements (" << nel 
									     << ") in the partial results\n"
									     << "     file exceeds the number of elements (" << element_map.Length() 
									     << ") in the map for block ID " << block_ID[l] << " in output\n"
									     << "     set " << i << " in partition " << k << endl;
									throw ExceptionT::kSizeMismatch;
								}
								
								/* skip empty sets */
								if (nel > 0)
								{
									/* alias */
									dArray2DT block_values(nel, part_e_values.MinorDim(), part_e_values(row_offset));
									row_offset += nel;
								
									/* read variables */
									results.ElementVariables(step_index, block_name, block_values);
								}
							}
							
							/* assemble */
							all_e_values.Assemble(element_map, part_e_values);
						}
					}
				}

				/* write assembled data */
				fOutput->WriteOutput(time, i, all_n_values, all_e_values);
				
				/* look for next time step */
				ResultFileName(0, io_ID, output_set.Changing(), j+1, filename);
				if (fstreamT::Exists(filename)) num_steps++;
			}	
		}
		else cout << "\n JoinOutputT: output set: EMPTY" << endl;	
	}
}

/*************************************************************************
* Private
*************************************************************************/

/* set output */
void JoinOutputT::SetOutput(void)
{
	const char caller[] = "JoinOutputT::SetOutput";

	/* set coordinates */
	fOutput->SetCoordinates(fModel->Coordinates(), NULL); //what about the map?
	
	/* block ID's in io groups */
	StringT io_file;
	io_file.Root(fJobFile);
	io_file.Append(".io.ID");
	ifstreamT io('#', io_file);
	if (!io.is_open())
		ExceptionT::GeneralFail(caller, "error opening \"%s\"", io_file.Pointer());

/* see if file contains changing geometry information */
#pragma message("not needed after merging BRANCH_david_spatial_2")
io.clear_marker();
StringT tmp;
tmp.GetLineFromStream(io); // first comment line
tmp.GetLineFromStream(io); // second comment line
bool has_changing_info = (tmp.StringMatch("[changing]") != NULL);
io.set_marker('#');

	/* construct output sets for each ID */
	int ID;
	io >> ID;
	while (io.good())
	{
		int i_changing = -99;
		if (has_changing_info) //TEMP
			io >> i_changing;
		bool changing = (i_changing == 1);
		int num_ID = -99;
		io >> num_ID;
		ArrayT<StringT> block_ID(num_ID);
		for (int i = 0; i < block_ID.Length(); i++)
			io >> block_ID[i];
	
		/* get output labels */
		ArrayT<StringT> n_labels;
		ArrayT<StringT> e_labels;
		OutputLabels(ID, changing, n_labels, e_labels);

		/* no block ID's implies a "free set" otherwise the set
		 * is tied to element blocks in the geometry file */
		if (block_ID.Length() == 0)
		{
			cout << "\n JoinOutputT::SetOutput: configuring free set: " << ID << endl;

			/* geometry from each part */
			ArrayT<iArray2DT> part_elems(fPartitions.Length());
			
			/* read geometry */
			int num_elem_nodes = 0;
			int num_elem = 0;
			GeometryT::CodeT geometry_code;
			ArrayT<StringT> block_ID;
			for (int i = 0; i < part_elems.Length(); i++) 
			{
				/* look for part file */
				StringT data_file;
				ResultFileName(i, ID, changing, 0, data_file);
				if (fstreamT::Exists(data_file))
				{
					/* open database */
					ModelManagerT results(cout);
					if (!results.Initialize(fResultsFileType, data_file, true))
						ExceptionT::DatabaseFail(caller, "could not initialize file \"%s\"", data_file.Pointer());
					
					/* check */
					if (results.NumElementGroups() != 1)
						ExceptionT::DatabaseFail(caller, "expecting 1 not %d element groups in free output set", results.NumElementGroups());
				
					/* part geometry */
					if (block_ID.Length() == 0) {
						block_ID.Dimension(1);
						block_ID[0] = results.ElementGroupID(0);
						geometry_code = results.ElementGroupGeometry(block_ID[0]);
					}

					/* set connectivities */
					iArray2DT& connects = part_elems[i];
					connects = results.ElementGroup(block_ID[0]);
					
					/* file -> partition numbering */
					iArrayT nodes(results.NumNodes());
					results.AllNodeIDs(nodes);
					nodes--; /* id -> index */
					for (int k = 0; k < connects.Length(); k++) 
						connects[k] = nodes[connects[k]];

					/* partition -> global numbering - index decomp keeps global numbering 
					 * since the entire system is reproduced on all processors. */
					if (fPartitions[i].DecompType() != PartitionT::kIndex) {
						nodes.Alias(fPartitions[i].NodeMap());
						for (int k = 0; k < connects.Length(); k++) 
							connects[k] = nodes[connects[k]];
					}

					/* counts */					
					num_elem += connects.MajorDim();
					if (num_elem_nodes == 0) num_elem_nodes = connects.MinorDim();
				}
			}
			
			/* register free set */
			if (block_ID.Length() == 1) {
				/* concat connectivities (no check for redundancy for now) */
				iArray2DT all_connects(num_elem, num_elem_nodes);
				num_elem = 0;
				for (int i = 0; i < part_elems.Length(); i++) {
					all_connects.BlockRowCopyAt(part_elems[i], num_elem);
					num_elem += part_elems[i].MajorDim();
				}

				/* generate dummy block ID for "connectivities" of free set nodes */			
				StringT sID = "900"; /* NOTE: same convention used by IOManager_mpi::IOManager_mpi for run time assembly */
				sID.Append(block_ID[0]);

				/* store connectivities in the model manager */
				if (!fModel->RegisterElementGroup(sID, all_connects, geometry_code, true)) {
					cout << "\n JoinOutputT::SetOutput: error registering free set with the model manager: " << sID << endl;
					throw ExceptionT::kDatabaseFail;	
				}
	
				/* construct output set */
				const iArray2DT& connects = fModel->ElementGroup(sID);
				OutputSetT output_set(geometry_code, connects, n_labels, changing);

				/* register */
				fOutput->AddElementSet(output_set);
			}
			else
				cout << "\n JoinOutputT::SetOutput: ERROR configuring free set: " << ID 
				     << ": SKIPPING"<< endl;
		}
		else
		{
			try {
			
			/* collect element blocks */
			GeometryT::CodeT geometry_code;
			ArrayT<const iArray2DT*> connects_list(block_ID.Length());
			for (int i = 0; i < block_ID.Length(); i++)
			{
				/* block ID as string */
				StringT block_name;
				block_name.Append(block_ID[i]);

				/* geometry code */
				geometry_code = fModel->ElementGroupGeometry(block_name);
			
				/* load element group */
				const iArray2DT& connects = fModel->ElementGroup(block_name);
				connects_list[i] = &connects;
			}

			/* construct output set */
			OutputSetT output_set(geometry_code, block_ID, connects_list, n_labels, e_labels, changing);
	
			/* register */
			fOutput->AddElementSet(output_set);
			
			} /* end try */
			
			catch (ExceptionT::CodeT exc) {
				cout << "\n JoinOutputT::SetOutput: caught exception " << exc << " configuring\n"
				     <<   "     output set ID " << ID << " with block ID's:\n";
				for (int i = 0; i < block_ID.Length(); i++)
					cout << setw(kIntWidth) << i+1 << "    " << block_ID[i] << '\n';
				cout << endl;
				
				/* register dummy set keep IO ID's the same */
				ArrayT<const iArray2DT*> empty_connects_list;
				ArrayT<StringT> empty_list;
				OutputSetT output_set(GeometryT::kPoint, empty_list, empty_connects_list, empty_list, empty_list, false);
				fOutput->AddElementSet(output_set);
			}
		}
	
		/* next block */
		io >> ID;
	}
}

/* set assembly maps
 *    nodes: map partition local number -> output set local number 
 * elements: ??? */
void JoinOutputT::SetMaps(void)
{
	/* check */
	if (!fOutput) throw ExceptionT::kGeneralFail;

	/* output sets data */
	const ArrayT<OutputSetT*>& element_sets = fOutput->ElementSets();

	/* dimensions */
	int num_parts = fPartitions.Length();
	int num_sets  = element_sets.Length();

	/* global to set maps */
	fMapSets.Dimension(num_sets);
	iArrayT shift(num_sets);
	ArrayT<iArrayT> inv_global(num_sets);
	for (int i = 0; i < num_sets; i++)
	{
		/* set data */
		OutputSetT& output_set = *(element_sets[i]);

		/* global nodes used by the set */
		const iArrayT& global_nodes_used = output_set.NodesUsed();

		/* global to set map */
		SetInverseMap(global_nodes_used, inv_global[i], shift[i], -1);
		
		/* allocate maps set */
		MapSetT& map_set = fMapSets[i];
		int n_sets = (output_set.NumNodeValues()    > 0) ? num_parts : 0;
		int e_sets = (output_set.NumElementValues() > 0) ? num_parts : 0;
		map_set.Dimension(n_sets, e_sets);
	}

	/* resident partition for each node */
	iArrayT node_part_map;
	SetNodePartitionMap(node_part_map);

	/* construct nodal assembly maps (loops reversed) */
	for (int j = 0; j < num_parts; j++)
	{
		for (int i = 0; i < num_sets; i++)
		{
			/* set data */
			OutputSetT& output_set = *(element_sets[i]);
			MapSetT& map_set = fMapSets[i];
	
			/* non-empty nodal output */
			if (map_set.NumNodeMaps() > 0)
			{
				/* resident partition nodes in set */
				iArrayT nodes;
				PartitionSetNodes(j, node_part_map, output_set.NodesUsed(), nodes);

				/* set output assembly map */
				SetAssemblyMap(inv_global[i], shift[i], nodes, map_set.NodeMap(j));
			}			
		}
	}
	
	/* element assembly maps */
	for (int i = 0; i < num_sets; i++)
	{
		/* output set data */
		const OutputSetT& output_set = *(element_sets[i]);
		MapSetT& map_set = fMapSets[i];
		if (map_set.NumElementMaps() > 0)
		{
			const ArrayT<StringT>& block_ID = output_set.BlockID();
			if (block_ID.Length() > 0)
			{
				/* get block sizes from max element number 
				 * element may be duplicated across partitions */
				iArrayT block_size(block_ID.Length());
				block_size = 0;
				for (int j = 0; j < num_parts; j++)
					for (int k = 0; k < block_ID.Length(); k++)
					{
						const iArrayT& element_map = fPartitions[j].ElementMap(block_ID[k]);
						if (element_map.Length() > 0)
						{
							int max_element = element_map.Max();
						
							/* find max across parts */
							block_size[k] = (max_element > block_size[k]) ?  
								max_element : block_size[k];
						}
					}
				
				/* max number is one less than size */
				block_size += 1;

				/* check total against output size */
				if (block_size.Sum() != output_set.NumElements())
				{
					cout << "\n JoinOutputT::SetMaps: expecting " << output_set.NumElements() 
					     << " elements for output ID " << i << ",\n" 
					     << "     found " << block_size.Sum() 
					     << " counting by block in partitions" << endl;
					throw ExceptionT::kSizeMismatch;
				}
				
				/* loop over partitions */
				for (int n = 0; n < num_parts; n++)
				{
					/* element assembly map */
					iArrayT& element_map = map_set.ElementMap(n);
				
					/* map size - no duplicated elements within
					 * a partition */
					int num_elems = 0;
					for (int j = 0; j < block_ID.Length(); j++)
						num_elems += fPartitions[n].ElementMap(block_ID[j]).Length();
						
					/* allocate map */
					element_map.Dimension(num_elems);
					element_map = -1;
					
					/* fill map */
					int offset = 0;
					int dex = 0;
					for (int k = 0; k < block_ID.Length(); k++)
					{
						const iArrayT& part_element_map = fPartitions[n].ElementMap(block_ID[k]);
						for (int j = 0; j < part_element_map.Length(); j++)
							element_map[dex++] = part_element_map[j] + offset;
					
						/* numbering offset in next block */
						offset += block_size[k];
					}
				}
			}
		}
	}

	/* check maps (checks node maps only) */
	CheckAssemblyMaps();
}

/* resident partition for each node */
void JoinOutputT::SetNodePartitionMap(iArrayT& node_partition)
{
	/* initialize */
	node_partition.Dimension(fModel->NumNodes());
	node_partition = -1;
	
	for (int i = 0; i < fPartitions.Length(); i++)
	{
		/* nodes in the partition */
		iArrayT partition_nodes;
		fPartitions[i].PartitionNodes(partition_nodes, PartitionT::kGlobal);
		
		/* write to map */
		for (int j = 0; j < partition_nodes.Length(); j++)
		{
			int node = partition_nodes[j];
			if (node_partition[node] != -1)
			{
				cout << "\n JoinOutputT::SetNodePartitionMap: node already assigned "
				     << node << endl;
				throw ExceptionT::kGeneralFail;
			}
			else
				node_partition[node] = i;
		}
	}
	
	/* check map is complete */
	int count = node_partition.Count(-1);
	if (count != 0)
	{
		cout << "\n JoinOutputT::SetNodePartitionMap: " << count
		     << " nodes are unassigned" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/* determine map from local nodes into global array, such that:
*
*             global[lg_map[i]] = local[i]
*/
void JoinOutputT::SetInverseMap(const iArrayT& global, iArrayT& inv_global,
	int& shift, int fill) const
{
	if (global.Length() == 0) {
		shift = 0;
		inv_global.Dimension(0);
	}
	else {
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
}

/* return the global node numbers of the set nodes residing
* in the partition */
void JoinOutputT::PartitionSetNodes(int partition, const iArrayT& node_part_map,
	const iArrayT& set_nodes, iArrayT& nodes) const
{
	/* count */
	int count = 0;
	for (int i = 0; i < set_nodes.Length(); i++)
		if (node_part_map[set_nodes[i]] == partition) count++;

	/* allocate return space */
	nodes.Dimension(count);
	
	/* copy in */
	count = 0;
	for (int j = 0; j < set_nodes.Length(); j++)
	{
		int node = set_nodes[j];
		if (node_part_map[node] == partition)
			nodes[count++] = node;
	}
	
	/* need nodes in ascending order of local number */
	fPartitions[partition].SetNodeScope(PartitionT::kLocal, nodes);
	nodes.SortAscending();
	fPartitions[partition].SetNodeScope(PartitionT::kGlobal, nodes);
}

/* add partition contribution to set assembly map */
void JoinOutputT::SetAssemblyMap(const iArrayT& inv_global, int shift, const iArrayT& local,
	iArrayT& lg_map) const
{
	/* set map */
	int n_map = local.Length();
	lg_map.Dimension(n_map);
//	int dex = 0;
	const int*  p = local.Pointer();
	for (int j = 0; j < n_map; j++)
	{
		int dex = inv_global[*p++ - shift];
		if (dex == -1) throw ExceptionT::kGeneralFail;
		lg_map[j] = dex;
	}	
}

/* generate output file name */
void JoinOutputT::ResultFileName(int part, int group, bool changing, int print_step, StringT& name) const
{
	/* basic name */
	name.Root(fJobFile);
	name.Append(".p", part);
	name.Append(".io", group);
	
	/* changing geometry */
	if (changing) name.Append(".ps", print_step, 4);
	
	/* file format extension */
	if (fResultsFileType == IOBaseT::kTahoeResults)
		name.Append(".run");
	else if (fResultsFileType == IOBaseT::kExodusII)
		name.Append(".exo");
	else
		ExceptionT::GeneralFail("JoinOutputT::ResultFileName",
			"unsupported file format %d", fResultsFileType);

	//NOTE: not handling multiple time sequences ".sq" or changing
	//      geometry groups ".ps"
}

/* returns the number of output steps */
int JoinOutputT::NumOutputSteps(int group) const
{
	/* database  */
	int num_steps = -1;

	/* global output sets */
	const ArrayT<OutputSetT*>& element_sets = fOutput->ElementSets();

	/* some partitions could be missing */
	for (int i = 0; i < fPartitions.Length() && num_steps < 0; i++)
	{
		/* generate file name */
		StringT filename;
		ResultFileName(i, group, element_sets[0]->Changing(), 0, filename);
		
		/* file exists */
		if (fstreamT::Exists(filename)) {
			/* open the database file */
			ModelManagerT results(cout);
			if (results.Initialize(fResultsFileType, filename, true))
				num_steps = results.NumTimeSteps();
		}
	}
	return num_steps;
}

/* retrieve output labels */
void JoinOutputT::OutputLabels(int group, bool changing, ArrayT<StringT>& node_labels,
	ArrayT<StringT>& element_labels) const
{
	/* some partitions could be missing */
	bool found_file = false;
	for (int i = 0; i < fPartitions.Length() && !found_file; i++)
	{
		/* generate file name */
		StringT filename;
		ResultFileName(i, group, changing, 0, filename);
		
		/* check if found */
		if (fstreamT::Exists(filename))
		{
			/* database  */
			ModelManagerT model(cout);
			if (!model.Initialize(fResultsFileType, filename, true))
				ExceptionT::DatabaseFail("JoinOutputT::OutputLabels", "error opening database file \"%s\"",
					fModel->DatabaseName().Pointer());

			/* read labels */
			found_file = true;
			model.NodeLabels(node_labels);
			model.ElementLabels(element_labels);		
		}
	}
	
	/* no results files found */
	if (!found_file) {
		cout << "\n JoinOutputT::OutputLabels: could not find output labels for set " 
		     << group << endl;
	}
}

/* check that assembly maps are compact and complete */
void JoinOutputT::CheckAssemblyMaps(void)
{
	/* global output sets */
	const ArrayT<OutputSetT*>& element_sets = fOutput->ElementSets();

	for (int i = 0; i < element_sets.Length(); i++)
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
				cout << "\n JoinOutputT::CheckAssemblyMaps: node maps size error: " << node_count
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
						cout << "\n JoinOutputT::CheckAssemblyMaps: duplicated fill for node "
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
				cout << "\n JoinOutputT::CheckAssemblyMaps: node maps error" << endl;
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
					cout << "\n JoinOutputT::CheckAssemblyMaps: element maps size error: " << element_count
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
					cout << "\n JoinOutputT::CheckAssemblyMaps: element maps are incomplete" << endl;
					throw ExceptionT::kGeneralFail;
				}
			}
	}
}
