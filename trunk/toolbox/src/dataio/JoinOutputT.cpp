/* $Id: JoinOutputT.cpp,v 1.1 2002-01-09 12:36:51 paklein Exp $ */
/* created: paklein (03/24/2000)                                          */

#include "JoinOutputT.h"

#include "fstreamT.h"
#include "IOManager.h"
#include "OutputSetT.h"
#include "StringT.h"
#include "ExodusT.h"
#include "iArray2DT.h"
#include "nVariArray2DT.h"
#include "VariArrayT.h"
#include "dArrayT.h"

/* constructor */
JoinOutputT::JoinOutputT(ifstreamT& in, const StringT& model_file,
	const StringT& global_model_file, IOBaseT::FileTypeT file_type,
	int size):
	fJobFile(in.filename()),
	fModelFile(model_file),
	fGlobalModelFile(global_model_file),
	fFileType(file_type),
	fPartitions(size),
	fIO(NULL)
{
	/* must be ExodusII for now */
	if (fFileType != IOBaseT::kExodusII)
	{
		cout << "\n JoinOutputT::JoinOutputT: data format must be ExodusII, code "
		     << IOBaseT::kExodusII << endl;
		throw eGeneralFail;
	}

	/* read partition data */
	for (int i = 0; i < fPartitions.Length(); i++)
	{
		/* file name */
		StringT file;
		file.Root(model_file);
		file.Append(".n", fPartitions.Length());
		file.Append(".part", i);

		/* open stream */
		ifstreamT part_in(in.comment_marker(), file);
		if (!part_in.is_open())
		{
			cout << "\n JoinOutputT::JoinOutputT: could not open decomposition file: "
			     << part_in.filename() << endl;
			throw eGeneralFail;
		}
		
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
	delete fIO;
	fIO = NULL;
}

/* do join */
void JoinOutputT::Join(void)
{
	/* assembly work space */
	dArray2DT part_n_values;
	dArray2DT part_e_values;
	nVariArray2DT<double> part_n_man(0, part_n_values, 0);
	nVariArray2DT<double> part_e_man(0, part_e_values, 0);	
	dArrayT n_value;
	dArrayT e_value;
	VariArrayT<double> n_value_man(0, n_value);
	VariArrayT<double> e_value_man(0, e_value);

	/* output sets data */
	StringT exo_file;
	const ArrayT<OutputSetT*>& element_sets = fIO->ElementSets();
	for (int i = 0; i < element_sets.Length(); i++)
	{
		cout << "\n JoinOutputT: output set: " << i+1 << endl;
	
		/* set data */
		const OutputSetT& output_set = *(element_sets[i]);
		const MapSetT& map_set = fMapSets[i];
	
		/* assembled values */
		dArray2DT all_n_values(output_set.NumNodes(), output_set.NumNodeValues());
		dArray2DT all_e_values;
		if (output_set.BlockID().Length() == 0 && output_set.NumElementValues() > 0)
			cout << "\n JoinOutputT::Join: skipping element output\n" << endl;
		else
			all_e_values.Allocate(output_set.NumElements(), output_set.NumElementValues());

		/* non-empty output */
		if (all_n_values.Length() > 0 || all_e_values.Length() > 0)
		{
			/* number of output steps */
			int num_steps = NumOutputSteps(i);
			cout << " JoinOutputT:      steps: " << num_steps << endl;
			int d_width = cout.precision() + kDoubleExtra;
			cout << setw(kIntWidth) << "step"
			     << setw(d_width)   << "time" << '\n';
			
			/* loop over steps */
			for (int j = 0; j < num_steps; j++)
			{
				/* initialize */
				all_n_values = 0.0;
				all_e_values = 0.0;
			
				/* loop over partitions */
				double time;
				bool found_time = false;
				for (int k = 0; k < fPartitions.Length(); k++)
				{
					/* file name */
					ExodusT exo(cout);
					ResultFileName(k, i, exo_file);
					if (exo.OpenRead(exo_file))
					{
						/* get time */
						if (!found_time)
						{
							found_time = true;
							exo.ReadTime(j+1, time);
							cout << setw(kIntWidth) << j+1
							     << setw(d_width)   << time << endl;
						}
						
						/* assemble nodal values */
						if (all_n_values.Length() > 0)
						{
							/* assembly map: output_set_ID[partition_ID] */
							const iArrayT& node_map = map_set.NodeMap(k);

							/* very weak consistency check */
							if (node_map.Length() > exo.NumNodes())
							{
								cout << "\n JoinOutputT::Join: assembly map of nodal values (" << node_map.Length() 
								     << ") is longer than\n" 
								     <<   "     is longer than the number of nodes (" << exo.NumNodes()
								     << ") in partial results file:\n"
								     <<   "     " << exo_file << endl;
								throw eSizeMismatch;
							}

							/* set work space */
							part_n_man.Dimension(exo.NumNodes(), all_n_values.MinorDim());
							n_value_man.SetLength(exo.NumNodes(), false);
							
							/* loop over variables */
							for (int l = 0; l < all_n_values.MinorDim(); l++)
							{
								/* read value */
								exo.ReadNodalVariable(j+1, l+1, n_value);

								/* write into table */
								part_n_values.SetColumn(l, n_value);
							}

							/* assemble */
							all_n_values.Assemble(node_map, part_n_values);
						}

						/* assemble element values */
						if (all_e_values.Length() > 0)
						{
							/* element map: global_output_block_ID[partition_output_block_ID] */
							const iArrayT& element_map = map_set.ElementMap(k);
							part_e_man.Dimension(element_map.Length(), all_e_values.MinorDim());
							e_value_man.SetLength(element_map.Length(), false);

							/* dimension check */
							int nel, nen;
							exo.ReadElementBlockDims(i+1, nel, nen);
							if (nel != e_value.Length())
							{
								cout << "\n JoinOutputT::Join: number of elements (" << nel 
								     << ") in the partial results\n"
								     << "     file does not match the number of elements (" << e_value.Length() 
								     << ") expected for I/O\n"
								     << "     set " << i << " in partition " << k << endl;
								throw eSizeMismatch;
							}
						
							/* loop over variables */
							for (int ll = 0; ll < all_e_values.MinorDim(); ll++)
							{
								/* read value */
								exo.ReadElementVariable(j+1, i+1, ll+1, e_value);
						
								/* write into table */
								part_e_values.SetColumn(ll, e_value);
							}
						
							/* assemble */
							all_e_values.Assemble(element_map, part_e_values);
						}
					}
					else if (false) // skip
					{
						if (j == 0)
							cout << " JoinOutputT::Join: skipping: "
							     << exo_file << endl;
					}
				}

				/* write assembled data */
				fIO->SetOutputTime(time);
				fIO->WriteOutput(i, all_n_values, all_e_values);
			}	
		}
	}
}

/*************************************************************************
* Private
*************************************************************************/

/* set output */
void JoinOutputT::SetOutput(void)
{
	if (fFileType != IOBaseT::kExodusII) throw eGeneralFail;

	/* construct I/O */
	StringT program_name("tahoe");
	StringT nothing("none");
	fIO = new IOManager(cout, program_name, nothing, nothing, fJobFile, fFileType);
	if (!fIO) throw eOutOfMemory;

//TEMP - for now need a global io file, but in the future, should
//       assemble the global geometry from the partitioned output

	/* global output file */
	ExodusT global_model(cout);
	if (!global_model.OpenRead(fGlobalModelFile))
	{
		cout << "\n JoinOutputT::SetOutput: could not open file: " << fGlobalModelFile << endl;
		throw eGeneralFail;
	}
	
	/* read coordinates */
	fCoordinates.Allocate(global_model.NumNodes(), global_model.NumDimensions());
	global_model.ReadCoordinates(fCoordinates);
	fIO->SetCoordinates(fCoordinates, NULL);
	
	/* block ID's in io groups */
	StringT io_file;
	io_file.Root(fJobFile);
	io_file.Append(".io.ID");
	ifstreamT io('#', io_file);
	bool have_IO_ID = io.is_open();

	/* register output sets */
	iArrayT element_ID(global_model.NumElementBlocks());
	global_model.ElementBlockID(element_ID);
	fConnects.Allocate(element_ID.Length());
	for (int i = 0; i < element_ID.Length(); i++)
	{	
		/* read group data */
		int num_elems;
		int num_elem_nodes;
		global_model.ReadElementBlockDims(i+1, num_elems, num_elem_nodes);
		fConnects[i].Allocate(num_elems, num_elem_nodes);

		GeometryT::CodeT geometry_code;
		global_model.ReadConnectivities(i+1, geometry_code, fConnects[i]);
		fConnects[i]--;
	
		/* get output labels */
		ArrayT<StringT> n_labels;
		ArrayT<StringT> e_labels;
		OutputLabels(i, n_labels, e_labels);
		
		/* block ID's */
		iArrayT block_ID;
		if (have_IO_ID)
		{
			int ID; 
			io >> ID;
			
			/* unexpected output ID */
			if (ID != i)
			{
				if (e_labels.Length() > 0)
					cout << "\n JoinOutputT::SetOutput: output ID " << ID << " from \"" << io_file 
					     << "\" does not match " << i << '\n' << "     Cannot join element output" 
					     << endl;
				
				/* clear line */
				io.clear_line();
			}
			/* read block ID list */
			else
			{
				int num_ID = -99;
				io >> num_ID;
				block_ID.Allocate(num_ID);
				io >> block_ID;
			}
		}
		else if (e_labels.Length() > 0)
			cout << "\n JoinOutputT::SetOutput: no file of block ID per output ID: \"" << io_file 
			     << "\"\n" << "     Cannot join element output" << endl;

		/* construct output set */
		bool changing = false; // changing geometry not supported
		ArrayT<const iArray2DT*> connects_list(1);
		connects_list[0] = &fConnects[i];
		OutputSetT output_set(i+1, geometry_code, block_ID, connects_list, n_labels, e_labels, changing);
	
		/* register */
		fIO->AddElementSet(output_set);
	}
}

/* set assembly maps
 *    nodes: map partition local number -> output set local number 
 * elements: ??? */
void JoinOutputT::SetMaps(void)
{
	/* check */
	if (!fIO) throw eGeneralFail;

	/* output sets data */
	const ArrayT<OutputSetT*>& element_sets = fIO->ElementSets();

	/* dimensions */
	int num_parts = fPartitions.Length();
	int num_sets  = element_sets.Length();

	/* global to set maps */
	fMapSets.Allocate(num_sets);
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
		map_set.Allocate(n_sets, e_sets);
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
			const iArrayT& block_ID = output_set.BlockID();
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
					throw eSizeMismatch;
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
					element_map.Allocate(num_elems);
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
	node_partition.Allocate(fCoordinates.MajorDim());
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
				throw eGeneralFail;
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
		throw eGeneralFail;
	}
}

/* determine map from local nodes into global array, such that:
*
*             global[lg_map[i]] = local[i]
*/
void JoinOutputT::SetInverseMap(const iArrayT& global, iArrayT& inv_global,
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
	nodes.Allocate(count);
	
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

/* generate output file name */
void JoinOutputT::ResultFileName(int part, int group, StringT& name) const
{
	/* basic name */
	name.Root(fJobFile);
	name.Append(".p", part);
	name.Append(".io", group);
	
	/* file format extension */
	if (fFileType == IOBaseT::kExodusII)
		name.Append(".exo");
	else
		cout << "\n JoinOutputT::ResultFileName: no extension added for file type "
		     << fFileType << endl;
	
	//NOTE: not handling multiple time sequences ".sq" or changing
	//      geometry groups ".ps"
}

/* returns the number of output steps */
int JoinOutputT::NumOutputSteps(int group) const
{
	/* database  */
	int num_steps = -1;
	ExodusT exo(cout);
	StringT filename;
	
	/* some partitions could be missing */
	for (int i = 0; i < fPartitions.Length() && num_steps < 0; i++)
	{
		/* generate file name */
		ResultFileName(i, group, filename);
		if (exo.OpenRead(filename))
			num_steps = exo.NumTimeSteps();
	}
	return num_steps;
}

/* retrieve output labels */
void JoinOutputT::OutputLabels(int group, ArrayT<StringT>& node_labels,
	ArrayT<StringT>& element_labels) const
{
	/* database  */
	ExodusT exo(cout);
	StringT filename;
	
	/* some partitions could be missing */
	bool found_file = false;
	for (int i = 0; i < fPartitions.Length() && !found_file; i++)
	{
		/* generate file name */
		ResultFileName(i, group, filename);
		if (exo.OpenRead(filename))
		{
			found_file = true;
			exo.ReadNodeLabels(node_labels);
			exo.ReadElementLabels(element_labels);
		}
	}
}

/* check that assembly maps are compact and complete */
void JoinOutputT::CheckAssemblyMaps(void)
{
	/* global output sets */
	const ArrayT<OutputSetT*>& element_sets = fIO->ElementSets();

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
						cout << "\n JoinOutputT::CheckAssemblyMaps: duplicated fill for node "
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
				cout << "\n JoinOutputT::CheckAssemblyMaps: node maps error" << endl;
				throw eGeneralFail;
			}
		}
			
		/* check element maps */
		//TEMP - not supported yet
	}
}
