/* $Id: main.cpp,v 1.2 2003-11-21 22:49:52 paklein Exp $ */
#include "ModelManagerT.h"
#include "iArray2DT.h"
#include "OutputBaseT.h"
#include "OutputSetT.h"

using namespace Tahoe;

int main (int argc, char* argv[])
{
#pragma unused(argc)
#pragma unused(argv)

	const char caller[] = "wrap";

#if 0
	//hard code the input
	StringT file = "/Volumes/Uster/USERS/tahoe/bin/conveyor/geometry/beam.CSE.0.geom";
	ArrayT<StringT> left_ID(2);
	left_ID[0] = "3";
	left_ID[1] = "5";
	ArrayT<StringT> right_ID(2);
	right_ID[0] = "4";
	right_ID[1] = "6";
#endif

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

	/* map node numbers of wrapper nodes */
	for (int k = 0; k < left_ID.Length(); k++)
	{
		const iArrayT& left_nodes = model.NodeSet(left_ID[k]);
		const iArrayT& match = n_min_match[k];
		for (int i = 0; i < left_nodes.Length(); i++)
			map[match[i]] = map[left_nodes[i]];
	}

	/* renumber connectivities */
	cout << " renumbering connectivities:\n";
	ArrayT<iArray2DT> connects(model.NumElementGroups());
	const ArrayT<StringT>& elem_ID = model.ElementGroupIDs();
	for (int i = 0; i < elem_ID.Length(); i++)
	{
		cout << '\t' << elem_ID[i] << endl;

		/* read element group */
		connects[i] = model.ElementGroup(elem_ID[i]);
		
		/* map */
		int* p = connects[i].Pointer();
		int len = connects[i].Length();
		for (int j = 0; j < len; j++)
		{
			*p = map[*p];
			p++;
		}
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
	
	OutputBaseT* output = IOBaseT::NewOutput(program_name, version, title, output_file, IOBaseT::kTahoeII, cout);
	output->SetCoordinates(coords_wrap, NULL);
	ArrayT<StringT> labels;
	for (int i = 0; i < elem_ID.Length(); i++)
	{
		ArrayT<StringT> block_ID(1);
		block_ID[0] = elem_ID[i];

		ArrayT<const iArray2DT*> conn (1);
		conn[0] = connects.Pointer(i);

		/* block information */
		OutputSetT set(model.ElementGroupGeometry(elem_ID[i]), block_ID, conn, labels, labels, false);

		/* register */
		output->AddElementSet(set);
	}
	
	for (int i = 0; i < node_ID.Length(); i++)
		output->AddNodeSet(node_sets[i], node_ID[i]);
	
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
	OutputBaseT* output_3D = IOBaseT::NewOutput(program_name, version, title, output_file_3D, IOBaseT::kTahoeII, cout);
	output_3D->SetCoordinates(coords_wrap_3D, NULL);
	for (int i = 0; i < elem_ID.Length(); i++)
	{
		ArrayT<StringT> block_ID(1);
		block_ID[0] = elem_ID[i];

		ArrayT<const iArray2DT*> conn (1);
		conn[0] = connects.Pointer(i);

		/* block information */
		OutputSetT set(model.ElementGroupGeometry(elem_ID[i]), block_ID, conn, labels, labels, false);

		/* register */
		output_3D->AddElementSet(set);
	}
	output_3D->WriteGeometry();
	cout << " wrote file: " << output_file_3D << endl;

	return 0;
}
