/* $Id: MergeResults.cpp,v 1.8 2003-02-25 15:05:37 sawimme Exp $ */
#include "MergeResults.h"
#include "ExceptionT.h"
#include "OutputSetT.h"

using namespace Tahoe;

MergeResults::MergeResults(ostream& out, istream& in, bool write):
	TranslateIOManager(out, in, write)
{

}

/* destructor */
MergeResults::~MergeResults(void)
{
	for (int i = 0; i < fInputs.Length(); i++)
		delete fInputs[i];
}

void MergeResults::Translate (const StringT& program, const StringT& version, const StringT& title)
{
	/* set sources */
	SetInput();
	SetOutput(program, version, title);

	/* collect comined coordinate list */
	dArray2DT coords;
	iArrayT nodes_used;
	CombinedCoordinates(coords, nodes_used);

	/* set coords - row is (global) id */
	iArrayT global_node_numbers(coords.MajorDim());
	global_node_numbers.SetValueToPosition();
	fOutput->SetCoordinates(coords, &global_node_numbers);
	
	/* collect list of connectivity ID's and node labels */
	AutoArrayT<StringT> union_elem_ids;
	AutoArrayT<StringT> union_node_labels;
	AutoArrayT<GeometryT::CodeT> union_elem_geom;
	AutoArrayT<const iArray2DT*> union_local_elem_connects;
	AutoArrayT<const iArrayT*> global_node_map;
	for (int i = 0; i < fInputs.Length(); i++)
	{
		/* element block ID's and geometries */
		const ArrayT<StringT>& ids = fInputs[i]->ElementGroupIDs();
		for (int j = 0; j < ids.Length(); j++)
			if (union_elem_ids.AppendUnique(ids[j]))
			{
				union_elem_geom.Append(fInputs[i]->ElementGroupGeometry(ids[j]));
				union_local_elem_connects.Append(&(fInputs[i]->ElementGroup(ids[j])));
				global_node_map.Append(fNodeMaps.Pointer(i));
			}
		
		/* node labels */
		ArrayT<StringT> labels;
		fInputs[i]->NodeLabels(labels);
		union_node_labels.AppendUnique(labels);
	}

	/* map connectivities to global numbering */
	ArrayT<iArray2DT> union_elem_connects(union_local_elem_connects.Length());
	for (int i = 0; i < union_elem_connects.Length(); i++)
	{
		const iArrayT& map = *(global_node_map[i]);
		const iArray2DT& local = *(union_local_elem_connects[i]);
		iArray2DT& global = union_elem_connects[i];
		global.Dimension(local);
		for (int j = 0; j < local.Length();  j++)
			global[j] = map[local[j]];
	}

	/* generate one output set per block ID */
	ArrayT<OutputSetT*> output_sets(union_elem_ids.Length());
	iArrayT output_ids(union_elem_ids.Length());
	for (int i = 0; i < output_sets.Length(); i++)
	{
		/* construct new output set */
		StringT ID;
		ID.Append(i);
		ArrayT<StringT> block_ID(1);
		block_ID[0] = union_elem_ids[i];
		ArrayT<const iArray2DT*> connectivities(1);
		connectivities[0] = union_elem_connects.Pointer(i);
		ArrayT<StringT> e_labels;
		OutputSetT* output_set = new OutputSetT(union_elem_geom[i], block_ID, connectivities, 
			union_node_labels, e_labels, false);
	
		/* register output */
		output_ids[i] = fOutput->AddElementSet(*output_set);
		output_sets[i] = output_set;
	}

	/* combine and write */
	dArrayT time_steps;
	fInputs[0]->TimeSteps(time_steps);
	for (int i = 0; i < time_steps.Length(); i++)
	{
		/* collect combined output data array */
		dArray2DT all_e_values;
		dArray2DT all_n_values(coords.MajorDim(), union_node_labels.Length());
		all_n_values = -1.0;
		for (int j = 0; j < fInputs.Length(); j++)
		{
			/* group to global node number map */
			const iArrayT& n_map = fNodeMaps[j];
		
			/* node values */
			ArrayT<StringT> labels;
			fInputs[j]->NodeLabels(labels);
			dArray2DT block_n_values(fInputs[j]->NumNodes() , fInputs[j]->NumNodeVariables());
			fInputs[j]->AllNodeVariables(i, block_n_values);
				
			/* assemble values into combined array */
			for (int s = 0; s < labels.Length(); s++)
			{
				int dex = -1;
				for (int r = 0; dex == -1 && r < union_node_labels.Length(); r++)
					if (labels[s] == union_node_labels[r])
						dex = r;
				if (dex == -1) throw ExceptionT::kOutOfRange;
					
				/* assemble */
				for (int r = 0; r < block_n_values.MajorDim(); r++)
					all_n_values(n_map[r], dex) = block_n_values(r, s);
			}
		}
		
		/* collect set values and write */	
		for (int j = 0; j < output_ids.Length(); j++)
		{
			/* nodes used by the set (local numbering) */
			const iArrayT& nodes_used = fOutput->NodesUsed(output_ids[j]);
		
			/* set values */
			dArray2DT e_values;
			dArray2DT n_values(nodes_used.Length(), all_n_values.MinorDim());
			n_values.RowCollect(nodes_used, all_n_values);
		
			/* write output */
			fOutput->WriteOutput(time_steps[i], output_ids[j], n_values, e_values);
		}
	}
	
	/* free memory */
	for (int i = 0; i < output_sets.Length(); i++)
		delete output_sets[i];
}

/**************** PROTECTED **********************/

void MergeResults::SetInput(void)
{
	int num_files = 0;
	cout << "\n Number of source files (> 1): ";
	fIn >> num_files;
	if (fEcho) fEchoOut << num_files << "\n";
	if (num_files < 2) ExceptionT::GeneralFail ("MergeResults::SetInput","Must have more than 1 file for merging.");
	
	fInputs.Dimension(num_files);
	fInputs = NULL;
	for (int i = 0; i < fInputs.Length(); i++)
	{
		StringT database;
		cout << "file " << i+1 << ": ";
		fIn >> database;
		if (fEcho) fEchoOut << database << "\n";
		database.ToNativePathName();
	
		/* try to guess file format */
		IOBaseT::FileTypeT file_type = IOBaseT::name_to_FileTypeT(database);
	
		/* new model manager */
		ModelManagerT* model = new ModelManagerT(cout);
		if (!model->Initialize(file_type, database, true))
		  ExceptionT::GeneralFail ("MergeResults::SetInput","Could not initialize results file %s", database.Pointer());

		fInputs[i] = model;
    }
    
    /* check that number of time steps is the same */
    int num_steps = -1;
	for (int i = 0; i < fInputs.Length(); i++)
	{
		int this_num_steps = fInputs[i]->NumTimeSteps();
		if (i > 0 && num_steps != this_num_steps)
 		        ExceptionT::SizeMismatch ("MergeResults::SetInput","num_steps=%i, this_num_steps=%i", num_steps, this_num_steps);
		else
			num_steps = this_num_steps;
	}
}

/**************** PRIVATE **********************/

/* generate combined coordinates list */
void MergeResults::CombinedCoordinates(dArray2DT& coords, iArrayT& nodes_used)
{
	/* find max id */	
	fNodeMaps.Dimension(fInputs.Length());
	int max_id = 0;
	for (int i = 0; i < fNodeMaps.Length(); i++)
	{
		/* read map */
		fInputs[i]->AllNodeIDs(fNodeMaps[i]);
	
		int this_max_id = fNodeMaps[i].Max();
		max_id = (this_max_id > max_id) ? this_max_id : max_id;
	}

	/* workspace */
	iArrayT all_ids(max_id+1);
	all_ids = -1;
	for (int i = 0; i < fNodeMaps.Length(); i++)
	{
		iArrayT& map = fNodeMaps[i];
		for (int j = 0; j < map.Length(); j++)
			all_ids[map[j]] = 1;
	}
	
	/* collect id list and construct assembly */
	int num_nodes = all_ids.Count(1);
	nodes_used.Dimension(num_nodes);
	int dex = 0;
	for (int i = 0; i < all_ids.Length(); i++)
	{
		if (all_ids[i] == 1)
			nodes_used[dex++] = i;
	}

	/* fill combined coordinate array */
	int nsd = fInputs[0]->NumDimensions();
	coords.Dimension(max_id+1, nsd);
	coords = 0.0;
	for (int i = 0; i < fInputs.Length(); i++)
	{
		const dArray2DT& coords_i = fInputs[i]->Coordinates();
		const iArrayT& map_i = fNodeMaps[i];
		for (int j = 0; j < coords_i.MajorDim(); j++)
			coords.SetRow(map_i[j], coords_i(j));
	}	
}
