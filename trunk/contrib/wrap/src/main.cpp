/* $Id: main.cpp,v 1.4 2004-12-27 06:05:40 paklein Exp $ */
#include "ModelManagerT.h"
#include "iArray2DT.h"
#include "OutputBaseT.h"
#include "OutputSetT.h"
#include "InverseMapT.h"

using namespace Tahoe;

int main (int argc, char* argv[])
{
#pragma unused(argc)
#pragma unused(argv)

	const char caller[] = "wrap";

	/* resolve input file name */
	StringT file;
	cout << " file path: ";
	cin >> file;
	file.ToNativePathName();
	IOBaseT::FileTypeT file_type = IOBaseT::name_to_FileTypeT(file);

	/* left-right node set pairs */
	int num_pairs = 0;
	cout << " number of left-right node set pairs: ";
	cin >> num_pairs;
	if (num_pairs < 1) return 0;
	ArrayT<StringT> left_ID(num_pairs), right_ID(num_pairs);
	for (int i = 0; i < num_pairs; i++)
	{
		cout << "  left(" << i+1 << "): ";
		cin >> left_ID[i];
		cout << " right(" << i+1 << "): ";
		cin >> right_ID[i];
		
		/* verify */
		cout << "(left, right) = (" << left_ID[i] << "," << right_ID[i] << ")\n";
		cout << " continue (y/n) ? ";
		char ans;
		cin >> ans;
		if (ans == 'n' || ans == 'N') i--;
	}

	/* open database */
	ModelManagerT model(cout);
	if (!model.Initialize(file_type, file, true))
		ExceptionT::GeneralFail(caller, "could not open file %s", file.Pointer());
	if (model.NumDimensions() != 2)
		ExceptionT::GeneralFail(caller, "dimension %dD must be 2D", model.NumDimensions());

	/* get coordinate bounds */
	const dArray2DT& coords = model.Coordinates();
	dArrayT x(coords.MajorDim());
	coords.ColumnCopy(0, x);
	double x_min, x_max;
	x.MinMax(x_min, x_max);
	cout << " coordinate bounds = {" << x_min << ", " << x_max << "}" << endl;
	
	/* find nodes on coordinate bounds */
	AutoArrayT<int> n_min, n_max;
	const double* px = coords.Pointer();
	int nsd = coords.MinorDim();
	for (int i = 0; i < coords.MajorDim(); i++)
	{
		double x = *px;
		if (fabs(x - x_min) < kSmall)
			n_min.Append(i);
		else if (fabs(x - x_max) < kSmall)
			n_max.Append(i);
		px += nsd;
	}
	if (n_min.Length() != n_max.Length())
		ExceptionT::GeneralFail(caller, "nodes at x_min %d must be nodes at x_max %d",
			n_min.Length(), n_max.Length());
	cout << " number of boundary nodes = " << 2*n_min.Length() << endl;
	

	/* mark pairs - loop over matching node set IDs */
	iArrayT map(coords.MajorDim());
	map = 0;
	ArrayT<iArrayT> n_min_match(left_ID.Length());
	for (int k = 0; k < left_ID.Length(); k++)
	{
		const iArrayT&  left_nodes = model.NodeSet(left_ID[k]);
		const iArrayT& right_nodes = model.NodeSet(right_ID[k]);
		iArrayT& match = n_min_match[k];
		if (left_nodes.Length() != right_nodes.Length()) 
			ExceptionT::GeneralFail(caller, "node set %s and %s have different length",
				left_ID[k].Pointer(), right_ID[k].Pointer());

		int n_edge = left_nodes.Length();
		match.Dimension(n_edge);
		for (int i = 0; i < n_edge; i++)
		{
			/* test to see if the node is on the boundaries */
			if (fabs(coords(left_nodes[i], 0) - x_min) > kSmall) ExceptionT::GeneralFail(caller, "node %d is not on the left edge", left_nodes[i]+1);
			if (fabs(coords(right_nodes[i], 0) - x_max) > kSmall) ExceptionT::GeneralFail(caller, "node %d is not on the right edge", right_nodes[i]+1);
		
			/* find match */
			double y_test = coords(left_nodes[i], 1);
			int j_match = -1;
			for (int j = 0; j_match == -1 && j < n_edge; j++)
				if (fabs(y_test - coords(right_nodes[j], 1)) < kSmall)
					j_match = j;
			
			/* match not found */
			if (j_match == -1) 
				ExceptionT::GeneralFail(caller, "no match found for node %d", left_nodes[i]+1);
			else /* update map */
			{
				match[i] = right_nodes[j_match];
				map[right_nodes[j_match]] = -1;
			}
		}
	}
	
	/* create renumbering map */
	dArray2DT coords_wrap(coords.MajorDim() - n_min.Length(), 2);
	int count = 0;
	for (int i = 0; i < map.Length(); i++)
		if (map[i] != -1) {
			coords_wrap.SetRow(count, coords(i));
			map[i] = count++;
		}
	iArrayT map_keepers(map);

	/* By here, map either contains -1 for all nodes on the right edge or node numbers for nodes not on the right edge */

	/* map node numbers of wrapper nodes */
	for (int k = 0; k < left_ID.Length(); k++) {
		const iArrayT& left_nodes = model.NodeSet(left_ID[k]);
		const iArrayT& match = n_min_match[k];
		for (int i = 0; i < left_nodes.Length(); i++)
			map[match[i]] = map[left_nodes[i]];
	}

	/*copied the leftnodes onto -1 map entries (i.e. right nodes)*/

	/* renumber connectivities */
	cout << " renumbering connectivities:\n";
	int num_groups = model.NumElementGroups();
	ArrayT<iArray2DT> connects(num_groups);

	iArrayT all_right_nodes;
	model.ManyNodeSets(right_ID,all_right_nodes);
	ArrayT<AutoArrayT<iArrayT> > split_connects(num_groups*2);
	ArrayT<AutoArrayT<int> > split_element_numbers(num_groups*2);
	const ArrayT<StringT>& elem_ID = model.ElementGroupIDs();
	for (int i = 0; i < num_groups; i++)
	{
		/* read element group */
		connects[i] = model.ElementGroup(elem_ID[i]);

		/* workspace */
		split_connects[i].Dimension(0);
		split_connects[i].SetHeadRoom(10);
		split_connects[i+num_groups].Dimension(0);
		split_connects[i+num_groups].SetHeadRoom(10);
		iArrayT tmp(connects[i].MinorDim());

		split_element_numbers[i].SetHeadRoom(25);
		split_element_numbers[i+num_groups].SetHeadRoom(25);

		/* map connectivites */
		for (int j = 0; j < connects[i].MajorDim(); j++) {
			int match_q = 0;
			for (int k = 0; match_q == 0 && k < connects[i].MinorDim(); k++) {
				for (int l = 0; match_q == 0 && l < all_right_nodes.Length(); l++) 
					match_q = connects[i](j,k) == all_right_nodes[l];
		  	}

			connects[i].RowCopy(j, tmp.Pointer());
			if (match_q == 1) {
				split_connects[i+num_groups].Append(tmp);
				split_element_numbers[i+num_groups].Append(j);
			}
			else {
				split_connects[i].Append(tmp);
				split_element_numbers[i].Append(j);				
			}
		}
	}

	/* update the number of groups */
	iArrayT update_group_index(num_groups);
	update_group_index = -1;
	int update_num_groups = num_groups;
	for (int i = num_groups; i < 2*num_groups; i++)
		if (split_connects[i].Length() > 0) {
			update_num_groups++;
			update_group_index[i - num_groups] = update_num_groups - 1; /* location of the split group */
		}

	/* copy into updated connectivities array */
	ArrayT<iArray2DT> update_connects(update_num_groups);
	int index =0;
	for (int i = 0; i < num_groups*2; i++) 
	{
		cout << '\t' << index << endl;
		if (split_connects[i].Length() > 0) {
			if (index > update_num_groups)
				ExceptionT::GeneralFail(caller, "index is greater than update_num_grops");
			
			/* map connectivities */
			update_connects[index].Dimension(split_connects[i].Length(), split_connects[i][0].Length());
			for (int j = 0; j < update_connects[index].MajorDim(); j++) {
				for (int k = 0; k < update_connects[index].MinorDim(); k++) {
	      			update_connects[index](j,k) = map[split_connects[i][j][k]];
				}
			}
			index++;
		}
	}
        
	/* assign elem_ID to new groups */
	ArrayT<StringT> update_elem_ID(update_num_groups);
	for (int i = 0; i < num_groups; i++)
		update_elem_ID[i] = elem_ID[i];
	for (int i = num_groups; i < update_num_groups; i++) {
		update_elem_ID[i] = "100";
		update_elem_ID[i].Append(elem_ID[i-num_groups]);
	}

	/* renumber/union node sets */
	const ArrayT<StringT>& node_ID = model.NodeSetIDs();
	ArrayT<iArrayT> node_sets(node_ID.Length());
	AutoArrayT<int> node_set_tmp;
	for (int i = 0; i < node_ID.Length(); i++)
	{
		/* collect nodes in wrapped mesh */
		node_set_tmp.Dimension(0);
		const iArrayT& nodes = model.NodeSet(node_ID[i]);
		for (int j = 0; j < nodes.Length(); j++)
		{
			int new_number = map_keepers[nodes[j]];
			if (new_number > -1)
				node_set_tmp.Append(new_number);
		}
	
		/* copy in */
		node_sets[i].Dimension(node_set_tmp.Length());
		node_set_tmp.CopyInto(node_sets[i]);
	}

	/* write wrapped geometry file */
	StringT program_name = "wrap";
	StringT version = "0.1";
	StringT title;
	StringT output_file = file;
	StringT ext;
	ext.Suffix(output_file);
	output_file.Root();
	output_file.Append(".wrap", ext);
	
	OutputBaseT* output = IOBaseT::NewOutput(program_name, version, title, output_file, file_type, cout);
	output->SetCoordinates(coords_wrap, NULL);
	ArrayT<StringT> labels;
	for (int i = 0; i < update_elem_ID.Length(); i++)
	{
		ArrayT<StringT> block_ID(1);
		block_ID[0] = update_elem_ID[i];

		ArrayT<const iArray2DT*> conn (1);
		conn[0] = update_connects.Pointer(i);

		/* block information */
		/*		OutputSetT set(model.ElementGroupGeometry(update_elem_ID[i]), block_ID, conn, labels, labels, false);
		 for now hard wire geometry code to quads*/
		OutputSetT set(GeometryT::kQuadrilateral, block_ID, conn, labels, labels, false);

		/* register */
		output->AddElementSet(set);
	}
	
	/* add node sets */
	for (int i = 0; i < node_ID.Length(); i++)
		output->AddNodeSet(node_sets[i], node_ID[i]);

	/* add side sets */
	const ArrayT<StringT>& side_ID = model.SideSetIDs();
	AutoArrayT<int> ss_1(25), ss_2(25);
	InverseMapT map_1, map_2;
	map_1.SetOutOfRange(InverseMapT::MinusOne);
	map_2.SetOutOfRange(InverseMapT::Throw);
	ArrayT<iArray2DT> ss_new(2*side_ID.Length());
	for (int i = 0; i < side_ID.Length(); i++) {
	
		/* read side set */
		const iArray2DT& ss = model.SideSet(side_ID[i]);
		
		/* associated element block */
		const StringT& elem_ID = model.SideSetGroupID(side_ID[i]);
		int elem_index = model.ElementGroupIndex(elem_ID);
		
		/* check if associated group is split */
		if (split_element_numbers[i+num_groups].Length() > 0) /* group is split */
		{
			/* set up */
			ss_1.Dimension(0);
			ss_2.Dimension(0);
			map_1.SetMap(split_element_numbers[i]);
			map_2.SetMap(split_element_numbers[i+num_groups]);
			
			/* split/renumber */
			for (int j = 0; j < ss.MajorDim(); j++) {
				int e_1 = map_1.Map(ss(j,0));
				if (e_1 > -1) {
					ss_1.Append(e_1);
					ss_1.Append(ss(j,1));
				} else {
					ss_2.Append(map_2.Map(ss(j,0)));
					ss_2.Append(ss(j,1));				
				}
			}
			
			/* add "old" set to output */
			iArray2DT ss_tmp;
			ss_tmp.Alias(ss_1.Length()/2, 2, ss_1.Pointer());
			ss_new[i] = ss_tmp;
			output->AddSideSet(ss_new[i], side_ID[i], elem_ID);

			/* add "new" set to output */
			ss_tmp.Alias(ss_2.Length()/2, 2, ss_2.Pointer());
			ss_new[i+side_ID.Length()] = ss_tmp;
			StringT ss_ID_new = "100";
			ss_ID_new.Append(side_ID[i]);
			int new_group_index = update_group_index[i];
			if (new_group_index == -1)
				ExceptionT::GeneralFail(caller, "could not find split group for side set \"%s\"",
					side_ID[i].Pointer());
			output->AddSideSet(ss_new[i+side_ID.Length()], ss_ID_new, update_elem_ID[new_group_index]);
		}
		else /* just add to output */
			output->AddSideSet(ss, side_ID[i], elem_ID);
	}

	output->WriteGeometry();
	cout << " wrote file: " << output_file << endl;

	/* generate 3D file? */
	double c = x_max - x_min;
	double r = c/acos(-1.0)/2.0;
	dArray2DT coords_wrap_3D(coords_wrap.MajorDim(), 3);	
	for (int i = 0; i < coords_wrap_3D.MajorDim(); i++)
	{
		/* copy y-coordinate */
		coords_wrap_3D(i,1) = coords_wrap(i,1);
	
		/* angle */
		double t = (coords_wrap(i,0) - x_min)/r;
		coords_wrap_3D(i,0) = r*cos(t);
		coords_wrap_3D(i,2) =-r*sin(t);
	}
	StringT output_file_3D = file;
	output_file_3D.Root();
	output_file_3D.Append(".3D", ext);
	OutputBaseT* output_3D = IOBaseT::NewOutput(program_name, version, title, output_file_3D, file_type, cout);
	output_3D->SetCoordinates(coords_wrap_3D, NULL);
	for (int i = 0; i < update_elem_ID.Length(); i++)
	{
		ArrayT<StringT> block_ID(1);
		block_ID[0] = update_elem_ID[i];

		ArrayT<const iArray2DT*> conn (1);
		conn[0] = update_connects.Pointer(i);

		/* block information */
		/*		OutputSetT set(model.ElementGroupGeometry(update_elem_ID[i]), block_ID, conn, labels, labels, false);
		 for now hard wire geometry code to quads*/
		OutputSetT set(GeometryT::kQuadrilateral, block_ID, conn, labels, labels, false);

		/* register */
		output_3D->AddElementSet(set);
	}

	/* add node sets */
	for (int i = 0; i < node_ID.Length(); i++)
		output_3D->AddNodeSet(node_sets[i], node_ID[i]);

	/* add side sets */
	ModelManagerT model_wrap(cout);
	if (!model_wrap.Initialize(file_type, output_file, true))
		ExceptionT::GeneralFail(caller, "could not open file %s", output_file.Pointer());
	const ArrayT<StringT>& side_ID_new = model_wrap.SideSetIDs();
	for (int i = 0; i < side_ID_new.Length(); i++)
	{
		/* read side set info */
		const iArray2DT& ss = model_wrap.SideSet(side_ID_new[i]);
		const StringT& elem_ID = model_wrap.SideSetGroupID(side_ID_new[i]);

		/* add to 3D output */
		output_3D->AddSideSet(ss, side_ID_new[i], elem_ID);	
	}

	output_3D->WriteGeometry();
	cout << " wrote file: " << output_file_3D << endl;

	return 0;
}
