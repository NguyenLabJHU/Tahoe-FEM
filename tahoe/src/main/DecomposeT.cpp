/* $Id: DecomposeT.cpp,v 1.2 2004-09-28 15:35:37 paklein Exp $ */
#include "DecomposeT.h"

#include "ofstreamT.h"
#include "ifstreamT.h"

#include "ModelManagerT.h"
#include "ModelFileT.h"
#include "ExodusT.h"
#include "PartitionT.h"
#include "ParameterListT.h"
#include "FEManagerT.h"
#include "GraphT.h"

#include <time.h>

using namespace Tahoe;
#pragma message("clear out unneeded headers")

/* constructor */
DecomposeT::DecomposeT()
{
	// should take values of command line
}

/*****************************************
 * Public
 *****************************************/
 
/* returns true if the global output model file is not found */
bool DecomposeT::NeedModelFile(const StringT& model_file,
	IOBaseT::FileTypeT format) const
{
	switch (format)
	{
		case IOBaseT::kTahoeII:
		{
			ModelFileT file;
			if (file.OpenRead(model_file) == ModelFileT::kOK)
				return false;
			else
				return true;
		}	
		case IOBaseT::kExodusII:
		{
			ExodusT file(cout);
			if (file.OpenRead(model_file))
				return false;			
			else
				return true;
		}
		default:
		{
			ifstreamT test(model_file);
			if (test.is_open())
				return false;
			else
				return true;
		}
	}
}
 
/* returns 1 if a new decomposition is needed */
bool DecomposeT::NeedDecomposition(const StringT& model_file,
	int size) const
{
	/* model file root */
	StringT root;
	root.Root(model_file);

	for (int i = 0; i < size; i++)
	{
		/* partition data file name */
		StringT part_file = root;
		part_file.Append(".n", size);
		part_file.Append(".part", i);
	
		/* open partition file */
		ifstreamT part_in(PartitionT::CommentMarker(), part_file);
		if (part_in.is_open())
		{
			StringT version;
			part_in >> version;
		
			int num_parts, part_ID;
			part_in >> num_parts;
			part_in >> part_ID;

			if (!PartitionT::CheckVersion(version) ||
			     num_parts != size ||
			     part_ID != i) return true;
		}
		else
			return true;
	}
	return false;
}
 
 
void DecomposeT::CheckDecompose(const StringT& input_file, int size, int decomp_type, CommunicatorT& comm,
	const StringT& model_file, IOBaseT::FileTypeT format, const ArrayT<StringT>& commandlineoptions) const
{	
	/* dispatch */
	switch (decomp_type)
	{
		case PartitionT::kGraph:
			Decompose_graph(input_file, size, comm, model_file, format, commandlineoptions);
			break;

		case PartitionT::kAtom:
			Decompose_atom(input_file, size, model_file, format, commandlineoptions);
			break;
			
		case PartitionT::kSpatial:
			cout << "\n DecomposeT::Decompose: spatial decomposition not implemented yet" << endl;
			break;
						
		default:
			cout << "\n DecomposeT::Decompose: unrecognized method: " << decomp_type << endl;
	}
}

/**********************************************************************
 * Private
 **********************************************************************/
 void DecomposeT::Decompose_atom(const StringT& input_file, int size,
	const StringT& model_file, IOBaseT::FileTypeT format, const ArrayT<StringT>& commandlineoptions) const
{
#pragma unused(input_file)
	const char caller[] = "FEExecutionManagerT::Decompose_atom";

	/* files exist */
	bool need_decomp = NeedDecomposition(model_file, size);
	if (!need_decomp)
	{
		cout << "\n " << caller <<": decomposition files exist" << endl;
		return;
	}
	
	/* model manager */
	ModelManagerT model(cout);
	if (!model.Initialize(format, model_file, true))
		ExceptionT::BadInputValue(caller, "could not open model file: %s", (const char*) model_file);

	/* dimensions */
	int nnd = model.NumNodes();
	int nsd = model.NumDimensions();

	/* node-to-partition map */
	iArrayT part_map(nnd);
	
	/* number of nodes per processor - p_i = i/proc_size */
	int part_size = nnd/size;
	if (part_size < 1) ExceptionT::GeneralFail();
	
	/* labels nodes */
	int part = 0;
	int count = 0;
	for (int i = 0; i < nnd; i++)
	{
		part_map[i] = part;
		if (++count == part_size) {
			if (part < size - 1) {
				part++;
				count = 0;
			}
		}
	}

	/* get pointers to all blocks */
	const ArrayT<StringT>& IDs = model.ElementGroupIDs();
	ArrayT<const iArray2DT*> connects_1(IDs.Length());
	model.ElementGroupPointers(IDs, connects_1);
	ArrayT<const RaggedArray2DT<int>*> connects_2;

	/* set partition information and write partial geometry files*/
	for (int i = 0; i < size; i++)
	{
		/* partition data */
		PartitionT partition;

		/* mark nodes */
		partition.Set(size, i, part_map, connects_1, connects_2);
		
		/* set elements */
		const ArrayT<StringT>& elem_ID = model.ElementGroupIDs();
		partition.InitElementBlocks(elem_ID);
		for (int j = 0; j < elem_ID.Length(); j++)
		{
			const iArray2DT& elems = model.ElementGroup(elem_ID[j]);
			partition.SetElements(elem_ID[j], elems);
		}

		/* set to local scope */
		partition.SetScope(PartitionT::kLocal);
		
		/* set decomposition type */
		partition.SetDecompType(PartitionT::kAtom);
	
		/* output file name */
		StringT geom_file, suffix;
		suffix.Suffix(model_file);
		geom_file.Root(model_file);
		geom_file.Append(".n", size);
		geom_file.Append(".p", i);
		geom_file.Append(suffix);
				
		cout << "     Writing partial model file: " << geom_file << endl;
		try { EchoPartialGeometry(partition, model, geom_file, format); }
		catch (ExceptionT::CodeT error)
		{
			ExceptionT::Throw(error, caller, "exception writing file: %s", (const char*) geom_file);
		}
		
		/* partition information */
		StringT part_file;
		part_file.Root(model_file);
		part_file.Append(".n", size);
		part_file.Append(".part", i);

		ofstream part_out(part_file);
		part_out << "# data for partition: " << i << '\n';
		part_out << partition << '\n';
		part_out.close();
	}
}

void DecomposeT::Decompose_spatial(const StringT& input_file, int size,
	const StringT& model_file, IOBaseT::FileTypeT format) const
{
#pragma unused(input_file)
#pragma unused(size)
#pragma unused(model_file)
#pragma unused(format)
	cout << "\n DecomposeT::Decompose_spatial: not implemented" << endl;
}

/* graph-based decomposition */
void DecomposeT::Decompose_graph(const StringT& input_file, int size,
	CommunicatorT& comm, const StringT& model_file, IOBaseT::FileTypeT format, const ArrayT<StringT>& commandlineoptions) const
{
	const char caller[] = "DecomposeT::Decompose_graph";

	bool need_decomp = NeedDecomposition(model_file, size);
	if (need_decomp)
	{
		/* echo stream */
		StringT decomp_file;
		decomp_file.Root(input_file);
		decomp_file.Append(".out");
		ofstreamT decomp_out;
		decomp_out.open(decomp_file);
		
		/* generate validated parameter list */
		ParameterListT valid_list;
		FEManagerT::ParseInput(input_file, valid_list, true, false, false, commandlineoptions);
// DEBUG
//cout << "\n" << caller << " Prepare to construct global problem \n" << endl;
		/* construct global problem */
		FEManagerT global_FEman(input_file, decomp_out, comm, commandlineoptions, FEManagerT::kDecompose);
		try { 
			global_FEman.TakeParameterList(valid_list); 
		}
		catch (ExceptionT::CodeT code) {
			ExceptionT::Throw(code, caller, "exception \"%s\" constructing global problem",
				ExceptionT::ToString(code));
		}

		/* decompose */
		ArrayT<PartitionT> partition(size);
		if (need_decomp)
		{
			int method = 0;

/* use METIS by default */
#ifdef __METIS__
			method = 1;
			for (int i = 0; method = 0 && i < commandlineoptions.Length(); i++)
				if (commandlineoptions[i] == "-no_metis")
					method = 0;

			if (method == 1)
				cout << "\n using METIS. Disable with command-line option \"-no_metis\"" << endl;
#endif
			/* graph object */
			GraphT graph;	
			try {
				cout << "\n Decomposing: " << model_file << endl;
				DoDecompose_graph(partition, graph, true, method, global_FEman);
				cout << " Decomposing: " << model_file << ": DONE"<< endl;
			}
			catch (ExceptionT::CodeT code) {
				ExceptionT::Throw(code, caller, "exception \"%s\" during decomposition",
					ExceptionT::ToString(code));
			}
			
			/* write partition data out */
			for (int q = 0; q < partition.Length(); q++)
			{
				/* set to local scope */
				partition[q].SetScope(PartitionT::kLocal);

				/* set decomposition type */
				partition[q].SetDecompType(PartitionT::kGraph);

				StringT file_name;
				file_name.Root(model_file);
				file_name.Append(".n", partition.Length());
				file_name.Append(".part", q);
		
				ofstream out_q(file_name);
				out_q << "# data for partition: " << q << '\n';
				out_q << partition[q] << '\n';
				out_q.close();
			}
			
			/* write decomposition map */
			if (format == IOBaseT::kExodusII)
			{
				/* file name */
				StringT map_file;
				map_file.Root(model_file);
				map_file.Append(".n", size);
				map_file.Append(".decomp.exo");
				cout << " Node map file: " << map_file << endl;
			
				/* database */
				ExodusT exo(cout);
				StringT str;
				ArrayT<StringT> str_list;
				const dArray2DT& coords = global_FEman.ModelManager()->Coordinates();
				int nnd = coords.MajorDim();
				int nsd = coords.MinorDim();
				exo.Create(map_file, str, str_list, str_list, nsd, nnd,
					nnd, size, 0, 0);
			
				/* coordinates */
				exo.WriteCoordinates(coords);
				
				/* data in decomp file */
				int n_internal = 0;
				int n_border = 0;
				int n_external = 0;
				dArrayT part(nnd); part = -1;
				dArrayT inex(nnd); inex = 0.0;
				for (int i = 0; i < size; i++)
				{
					/* "owned" nodes */
					const iArrayT& nd_i = partition[i].Nodes_Internal();
					const iArrayT& nd_b = partition[i].Nodes_Border();
					iArray2DT connects(nd_i.Length() + nd_b.Length(), 1);
					connects.CopyPart(0, nd_i, 0, nd_i.Length());
					connects.CopyPart(nd_i.Length(), nd_b, 0, nd_b.Length());
					
					/* increment counts */
					n_internal += nd_i.Length();
					n_border += nd_b.Length();
					n_external += partition[i].Nodes_External().Length();
					
					/* convert to global numbering */
					if (partition[i].NumberScope() != PartitionT::kGlobal)
						partition[i].SetNodeScope(PartitionT::kGlobal, connects);
				
					/* label part */
					for (int j = 0; j < connects.Length(); j++)
						part[connects[j]] = i;
				
					/* label border nodes */
					int* pnd_b = connects.Pointer(nd_i.Length());
					for (int k = 0; k < nd_b.Length(); k++)
						inex[*pnd_b++] = 1;
				
					/* write to file */
					connects++;
					exo.WriteConnectivities(i+1, GeometryT::kPoint, connects);
				}
				
				/* degree of each node */
				iArrayT i_degree(nnd);
				int shift;
				graph.Degrees(i_degree, shift);
				if (shift != 0)				
					ExceptionT::GeneralFail(caller, "unexpected node number shift: %d", shift);

				dArrayT degree(nnd);
				for (int j = 0; j < nnd; j++)
					degree[j] = i_degree[j];

				/* part labels */
				ArrayT<StringT> labels(3);
				labels[0] = "part";
				labels[1] = "inex";
				labels[2] = "degree";
				exo.WriteLabels(labels, ExodusT::kNode);
				exo.WriteTime(1, 0.0);
				exo.WriteNodalVariable(1, 1, part);
				exo.WriteNodalVariable(1, 2, inex);
				exo.WriteNodalVariable(1, 3, degree);
				
				/* write statistics */
				cout << " Statistics:\n" 
				     << "     total number of nodes = " << coords.MajorDim() << '\n'
				     << "            internal nodes = " << n_internal << '\n'
				     << "              border nodes = " << n_border << '\n'
				     << "            external nodes = " << n_external << '\n';

				/* done */
				cout << " Node map file: " << map_file << ": DONE"<< endl;
			}
		}
		else
		{
			/* read partition information from stream */
			for (int i = 0; i < partition.Length(); i++)
			{
				StringT part_file;
				part_file.Root(model_file);
				part_file.Append(".n", size);
				part_file.Append(".part", i);
				ifstreamT part_in('#', part_file);
				
				part_in >> partition[i];
				partition[i].SetScope(PartitionT::kLocal);
			}
		}

/* when generating decomposition during parallel run, this was skipped and each
 * processor wrote it's own partial geometry file. When computing decomposition
 * explicitly ("-decomp"), all partial geometry files are written here. 
 */

		/* write partial geometry files */
		if (true)
		{
			/* model manager for the total geometry - can't use the one from global_FEman 
			 * because it may contain runtime-generated connectivities not in the original 
			 * model file */
			ModelManagerT model_ALL(cout);
			if (!model_ALL.Initialize(format, model_file, true))
				ExceptionT::GeneralFail(caller, "error opening file: %s", (const char*) model_file);
	
			/* write partial geometry files */
			for (int i = 0; i < partition.Length(); i++)
			{
				StringT partial_file, suffix;
				suffix.Suffix(model_file);
				partial_file.Root(model_file);
				partial_file.Append(".n", size);
				partial_file.Append(".p", i);
				partial_file.Append(suffix);
				
				if (NeedModelFile(partial_file, format))
				{			
					cout << "     Writing partial model file: " << partial_file << endl;
					try { 
						EchoPartialGeometry(partition[i], model_ALL, partial_file, format); 
					}
					catch (ExceptionT::CodeT error) {
						ExceptionT::Throw(error, caller, "exception writing file \"%s\"",
							partial_file.Pointer());
					}
				}
			}
		}
	}
	else
		cout << "\n " << caller << ": decomposition files exist" << endl;
}

/* domain decomposition */
void DecomposeT::DoDecompose_graph(ArrayT<PartitionT>& partition, GraphT& graphU,
	bool verbose, int method, const FEManagerT& feman) const
{
	const char caller[] = "DecomposeT::Decompose";

	//TEMP
	//TimeStamp("FEManagerT_mpi::Decompose");

	/* check */
	if (partition.Length() == 1)
		ExceptionT::GeneralFail(caller, "expecting more than 1 partition");
	
	IOBaseT::FileTypeT modelformat = feman.ModelFormat();
	/* geometry file must be ascii external */
	if (modelformat != IOBaseT::kTahoeII && modelformat != IOBaseT::kExodusII)
		ExceptionT::BadInputValue(caller, "expecting file format %d or %d, not %d",
			IOBaseT::kTahoeII, IOBaseT::kExodusII, modelformat);

	/* decomposition method */
	bool use_new_methods = false; //TEMP
	bool dual_graph = (feman.InterpolantDOFs() == 0);
	if (dual_graph && use_new_methods)
		DoDecompose_graph_2(partition, graphU, verbose, method, feman);
	else
		DoDecompose_graph_1(partition, graphU, verbose, method, feman);
	
	ModelManagerT* modelmanager = feman.ModelManager();
	const StringT& modelfile = modelmanager->DatabaseName();
	
	if (modelformat == IOBaseT::kTahoeII)
	{
		/* label element sets in partition data */
		for (int j = 0; j < partition.Length(); j++)
		{
			/* original model file */
			ModelFileT model_ALL;
			model_ALL.OpenRead(modelfile);
			
			/* set number of element sets */
			iArrayT elementID;
			if (model_ALL.GetElementSetID(elementID) != ModelFileT::kOK) ExceptionT::GeneralFail(caller);
			ArrayT<StringT> IDlist(elementID.Length());
			for (int i = 0; i < IDlist.Length(); i++)
				IDlist[i].Append(elementID[i]);
			partition[j].InitElementBlocks(IDlist);	
			for (int i = 0; i < elementID.Length(); i++)
			{
				/* get element set */
				iArray2DT set;
				if (model_ALL.GetElementSet(elementID[i], set) != ModelFileT::kOK)
					ExceptionT::GeneralFail(caller);
					
				/* correct node numbering offset */
				set--;	
				
				/* set partition */
				partition[j].SetElements(IDlist[i], set);
			}
		}
	}
	else if (modelformat == IOBaseT::kExodusII)
	{
		/* label element sets in partition data */
		for (int j = 0; j < partition.Length(); j++)
		{
			/* original model file */
			ExodusT model_ALL(cout);
			model_ALL.OpenRead(modelfile);
			
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
	else ExceptionT::GeneralFail(caller);
}

/* decomposition methods */
void DecomposeT::DoDecompose_graph_1(ArrayT<PartitionT>& partition, GraphT& graph, 
	bool verbose, int method, const FEManagerT& feman) const
{
	/* connectivities for partititioning */
	AutoArrayT<const iArray2DT*> connects_1;
	AutoArrayT<const RaggedArray2DT<int>*> connects_2;
	AutoArrayT<const iArray2DT*> equivalent_nodes;

	/* collect connectivies from all solver groups */
		feman.ConnectsU(connects_1,connects_2, equivalent_nodes);

	/* initialize graph */
	GraphT& graphU = graph;
	for (int r = 0; r < connects_1.Length(); r++)
		graphU.AddGroup(*(connects_1[r])); 

	for (int k = 0; k < connects_2.Length(); k++)
		graphU.AddGroup(*(connects_2[k]));

	for (int k = 0; k < equivalent_nodes.Length(); k++)
		graphU.AddEquivalentNodes(*(equivalent_nodes[k]));
		
	/* make graph */
	clock_t t0 = clock();
	if (verbose) cout << " DecomposeT::DoDecompose_graph_1: constructing graph" << endl;
	graphU.MakeGraph();
	clock_t t1 = clock();
	if (verbose)
		cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
		     << " sec: DecomposeT::DoDecompose_graph_1: construct graph" << endl;
	
	/* dual graph partitioning graph */
	int dual_graph = (feman.InterpolantDOFs() == 0) ? 1 : 0;
	AutoArrayT<const iArray2DT*> connectsX_1;
	GraphT graphX;
	if (dual_graph == 1)
	{
		if (verbose) cout << " DecomposeT::DoDecompose_graph_1: constructing dual graph" << endl;
		
		/* collect element groups */
			feman.ConnectsX(connectsX_1);

		/* initialize graph */
		for (int r = 0; r < connectsX_1.Length(); r++)
			graphX.AddGroup(*(connectsX_1[r]));
		
		/* make graph */
		clock_t t0 = clock();		
		graphX.MakeGraph();
		clock_t t1 = clock();
		if (verbose)
			cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
		     << " sec: DecomposeT::DoDecompose_graph_1: construct X graph" << endl;
	}
	
	/* generate partition */
	iArrayT config(1); //TEMP - will be scalar soon?
	config[0] = partition.Length();	
	iArrayT weight;
	feman.WeightNodalCost(weight);
	if (dual_graph == 1)
		graphX.Partition(config, weight, graphU, partition, true, method);
	else
		graphU.Partition(config, weight, partition, true, method);
}

void DecomposeT::DoDecompose_graph_2(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int
	method, const FEManagerT& feman) const
{
	/* connectivities for partititioning */
	AutoArrayT<const iArray2DT*> connects_1;
	AutoArrayT<const RaggedArray2DT<int>*> connects_2;
	AutoArrayT<const iArray2DT*> equivalent_nodes;

	/* collect element groups */
	feman.ConnectsU(connects_1, connects_2, equivalent_nodes);		
	
	/* dual graph partitioning graph */
	AutoArrayT<const iArray2DT*> connectsX_1;
	
	/* collect minimal connects */
	feman.ConnectsX(connectsX_1);

	/* initialize graph */
	GraphT& graphX = graph;
	for (int r = 0; r < connectsX_1.Length(); r++)
		graphX.AddGroup(*(connectsX_1[r]));
	for (int k = 0; k < equivalent_nodes.Length(); k++)
		graphX.AddEquivalentNodes(*(equivalent_nodes[k]));
		
	/* make graph */
	if (verbose) cout << " DecomposeT::DoDecompose_graph_2: constructing dual graph" << endl;
	clock_t t0 = clock();		
	graphX.MakeGraph();
	clock_t t1 = clock();
	if (verbose)
		cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
	     << " sec: DecomposeT::DoDecompose_graph_2: construct graph" << endl;
	
	/* generate partition */
	iArrayT config(1); //TEMP - will be scalar soon?
	config[0] = partition.Length();	
	iArrayT weight;
	feman.WeightNodalCost(weight);
	graphX.Partition(config, weight, connects_1, connects_2, partition, true, method);
}

/* write partial geometry files */
void DecomposeT::EchoPartialGeometry(const PartitionT& partition,
	ModelManagerT& model_ALL, const StringT& partial_file,
	IOBaseT::FileTypeT format) const
{
	switch (format)
	{
		case IOBaseT::kExodusII:
			EchoPartialGeometry_ExodusII(partition, model_ALL, partial_file);
			break;
	
		case IOBaseT::kTahoeII:
			EchoPartialGeometry_TahoeII(partition, model_ALL, partial_file);
			break;

		default:	
			ExceptionT::GeneralFail("FEExecutionManagerT::EchoPartialGeometry", 
				"unsupported file format: %d", format);
	}
}

void DecomposeT::EchoPartialGeometry_ExodusII(const PartitionT& partition,
	ModelManagerT& model_ALL, const StringT& partial_file) const
{
	/* partition */
	int part = partition.ID();

	/* collect file creation information */
	StringT title = "partition file: ";
	title.Append(part);
	ArrayT<StringT> nothing;
	const iArrayT& node_map = partition.NodeMap();
	int num_node = node_map.Length();
	int num_dim  = model_ALL.NumDimensions();
	int num_blks = model_ALL.NumElementGroups();
	int num_ns   = model_ALL.NumNodeSets();
	int num_ss   = model_ALL.NumSideSets();
	int num_elem = model_ALL.NumElements();
	
	/* partial model file */
	ExodusT model(cout);
	model.Create(partial_file, title, nothing, nothing, num_dim, num_node,
		num_elem, num_blks, num_ns, num_ss);
	
	/* coordinates */
	const dArray2DT& coords_ALL = model_ALL.Coordinates();
	dArray2DT coords(node_map.Length(), coords_ALL.MinorDim());		
	coords.RowCollect(node_map, coords_ALL);
	model.WriteCoordinates(coords);
	coords.Free();
		
	/* element sets */
	const ArrayT<StringT>& elem_ID = model_ALL.ElementGroupIDs();
	for (int j = 0; j < elem_ID.Length(); j++)
	{
		/* read global block */
		const iArray2DT& set_ALL = model_ALL.ElementGroup(elem_ID[j]);

		/* collect connectivities within the partition */
		const iArrayT& element_map = partition.ElementMap(elem_ID[j]);
		iArray2DT set(element_map.Length(), set_ALL.MinorDim());
		set.RowCollect(element_map, set_ALL);

		/* map to local scope */
		partition.SetNodeScope(PartitionT::kLocal, set);

		/* write to file */
		set++;
		GeometryT::CodeT geometry_code = model_ALL.ElementGroupGeometry(elem_ID[j]);
		int id = atoi(elem_ID[j]);
		model.WriteConnectivities(id, geometry_code, set);
	}
		
	/* node sets */
	const ArrayT<StringT>& node_ID = model_ALL.NodeSetIDs();
	for (int k = 0; k < node_ID.Length(); k++)
	{
		/* whole node set */
		const iArrayT& nodeset_ALL = model_ALL.NodeSet(node_ID[k]);
				
		/* partition node set */
		iArrayT local_indices;
		partition.ReturnPartitionNodes(nodeset_ALL, local_indices);
		
		/* non-empty set */
		if (1 || local_indices.Length() > 0)
		{
			iArrayT nodeset(local_indices.Length());
			nodeset.Collect(local_indices, nodeset_ALL);
			
			/* map to local numbering */
			partition.SetNodeScope(PartitionT::kLocal, nodeset);
				
			/* add */
			int id = atoi(node_ID[k]);
			nodeset++;
			model.WriteNodeSet(id, nodeset);
		}
	}

	/* side sets */		
	const ArrayT<StringT>& side_ID = model_ALL.SideSetIDs();
	for (int l = 0; l < side_ID.Length(); l++)
	{
		/* whole side set */
		iArray2DT sideset_ALL = model_ALL.SideSet(side_ID[l]);
		const StringT& element_set_ID = model_ALL.SideSetGroupID(side_ID[l]);

		/* partition side set */
		iArray2DT sideset;
		if (sideset_ALL.MajorDim() > 0)
		{
			iArrayT elements_ALL(sideset_ALL.MajorDim());
			sideset_ALL.ColumnCopy(0, elements_ALL);
			sideset_ALL.SetColumn(0, elements_ALL);
				
			iArrayT local_indices;
			partition.ReturnPartitionElements(element_set_ID, elements_ALL, local_indices);
						
			/* non-empty set */
			if (local_indices.Length() > 0)
			{
				sideset.Dimension(local_indices.Length(), sideset_ALL.MinorDim());
				sideset.RowCollect(local_indices, sideset_ALL);

				iArrayT elements(sideset.MajorDim());
				sideset.ColumnCopy(0, elements);
				partition.SetElementScope(PartitionT::kLocal, element_set_ID, elements);
				sideset.SetColumn(0, elements);
			}

			/* add */
			sideset++;
			int ss_id = atoi(side_ID[l]);
			int el_id = atoi(element_set_ID);
			model.WriteSideSet(ss_id, el_id, sideset);
		}			
	}
}

void DecomposeT::EchoPartialGeometry_TahoeII(const PartitionT& partition,
	ModelManagerT& model_ALL, const StringT& partial_file) const
{
	/* partition */
	int part = partition.ID();
	
	/* open model file */
	bool extern_file = true;
	ModelFileT model;
	model.OpenWrite(partial_file, extern_file);

	/* title */
	StringT title;
	title.Append("partition ", part);
	if (model.PutTitle(title) != ModelFileT::kOK) ExceptionT::GeneralFail();

	/* nodal coordinates */
	const iArrayT& node_map = partition.NodeMap();
	const dArray2DT& coords_ALL = model_ALL.Coordinates();
	dArray2DT coords(node_map.Length(), coords_ALL.MinorDim());		
	coords.RowCollect(node_map, coords_ALL);
	if (model.PutCoordinates(coords) != ModelFileT::kOK) ExceptionT::GeneralFail();
	coords.Free();	
		
	/* element sets */
	const ArrayT<StringT>& elem_ID = model_ALL.ElementGroupIDs();
	for (int j = 0; j < elem_ID.Length(); j++)
	{
		/* read global block */
		const iArray2DT& set_ALL = model_ALL.ElementGroup(elem_ID[j]);

		/* collect connectivities within the partition */
		const iArrayT& element_map = partition.ElementMap(elem_ID[j]);
		iArray2DT set(element_map.Length(), set_ALL.MinorDim());
		set.RowCollect(element_map, set_ALL);
			
		/* map to local scope */
		partition.SetNodeScope(PartitionT::kLocal, set);

		/* add */
		int id = atoi(elem_ID[j]);
		set++;
		if (model.PutElementSet(id, set) != ModelFileT::kOK) ExceptionT::GeneralFail();
	}
		
	/* node sets */
	const ArrayT<StringT>& node_ID = model_ALL.NodeSetIDs();
	for (int k = 0; k < node_ID.Length(); k++)
	{
		/* whole node set */
		const iArrayT& nodeset_ALL = model_ALL.NodeSet(node_ID[k]);
				
		/* partition node set */
		iArrayT local_indices;
		partition.ReturnPartitionNodes(nodeset_ALL, local_indices);
		iArrayT nodeset(local_indices.Length());
		nodeset.Collect(local_indices, nodeset_ALL);
			
		/* map to local numbering */
		partition.SetNodeScope(PartitionT::kLocal, nodeset);
			
		/* add */
		int id = atoi(node_ID[k]);
		nodeset++;
		if (model.PutNodeSet(id, nodeset) != ModelFileT::kOK) ExceptionT::GeneralFail();
	}

	/* side sets */		
	const ArrayT<StringT>& side_ID = model_ALL.SideSetIDs();
	for (int l = 0; l < side_ID.Length(); l++)
	{
		/* whole side set */
		iArray2DT sideset_ALL = model_ALL.SideSet(side_ID[l]);
		const StringT& element_set_ID = model_ALL.SideSetGroupID(side_ID[l]);

		/* partition side set */
		iArray2DT sideset;
		if (sideset_ALL.MajorDim() > 0)
		{
			iArrayT elements_ALL(sideset_ALL.MajorDim());
			sideset_ALL.ColumnCopy(0, elements_ALL);
			sideset_ALL.SetColumn(0, elements_ALL);
				
			iArrayT local_indices;
			partition.ReturnPartitionElements(element_set_ID, elements_ALL, local_indices);
			sideset.Dimension(local_indices.Length(), sideset_ALL.MinorDim());
			sideset.RowCollect(local_indices, sideset_ALL);
				
			/* map to (block) local numbering */
			if (sideset.MajorDim() > 0)
			{
				iArrayT elements(sideset.MajorDim());
				sideset.ColumnCopy(0, elements);
				partition.SetElementScope(PartitionT::kLocal, element_set_ID, elements);
				sideset.SetColumn(0, elements);
			}
		}
			
		/* add */
		sideset++;
		int ss_id = atoi(side_ID[l]);
		int el_id = atoi(element_set_ID);
		if (model.PutSideSet(ss_id, el_id, sideset) != ModelFileT::kOK) ExceptionT::GeneralFail();
	}

	/* close database */
	model.Close();
}