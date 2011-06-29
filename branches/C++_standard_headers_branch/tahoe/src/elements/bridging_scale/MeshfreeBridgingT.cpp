/* $Id: MeshfreeBridgingT.cpp,v 1.12 2005-04-28 23:54:50 paklein Exp $ */
#include "MeshfreeBridgingT.h"

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "PointInCellDataT.h"
#include "ShapeFunctionT.h"
#include "iGridManagerT.h"
#include "iNodeT.h"
#include "InverseMapT.h"
#include "MLSSolverT.h"
#include "VariLocalArrayT.h"
#include "VariArrayT.h"
#include "CommManagerT.h"
#include "OutputBaseT.h"
#include "OutputSetT.h"
#include "SolidElementT.h"

using namespace Tahoe;

/* constructor */
MeshfreeBridgingT::MeshfreeBridgingT(const ElementSupportT& support):
	BridgingScaleT(support),
	fMLS(NULL)
{
	SetName("meshfree_bridging");
}

/* destructor */
MeshfreeBridgingT::~MeshfreeBridgingT(void) { delete fMLS; }

/* initialize projection data */
void MeshfreeBridgingT::InitProjection(CommManagerT& comm, const iArrayT& points_used, const dArray2DT* init_coords, 
	const dArray2DT* curr_coords, PointInCellDataT& cell_data)
{
	const char caller[] = "MeshfreeBridgingT::InitProjection";

	/* collect point within each nodal neighborhood */
	BuildNodalNeighborhoods(comm, points_used, init_coords, curr_coords, cell_data);

	/* point coordinates */
	if (curr_coords && init_coords) ExceptionT::GeneralFail(caller, "cannot pass both init and curr coords");
	if (!curr_coords && !init_coords) ExceptionT::GeneralFail(caller, "must define init or curr coords");
	const dArray2DT& point_coordinates = (init_coords != NULL) ? *init_coords : *curr_coords;

	/* cell coordinates */
	const dArray2DT& cell_coordinates = (init_coords != NULL) ?
		ElementSupport().InitialCoordinates() : ElementSupport().CurrentCoordinates();
	LocalArrayT::TypeT coord_type = (init_coords != NULL) ? 
		LocalArrayT::kInitCoords : LocalArrayT::kCurrCoords;

	/* nodal neighbor data */
	InterpolationDataT& point_to_node = cell_data.PointToNode();
	const RaggedArray2DT<int>& nodal_neighbors = point_to_node.Neighbors();
	iArrayT neighbor_count(nodal_neighbors.MajorDim());
	nodal_neighbors.MinorDim(neighbor_count);
	RaggedArray2DT<double>& neighbor_weights = point_to_node.NeighborWeights();
	neighbor_weights.Configure(neighbor_count);
	neighbor_count.Free();
	
	/* point coordinates */
	dArray2DT neighbor_coords;
	nVariArray2DT<double> neighbor_coords_man(0, neighbor_coords, point_coordinates.MinorDim());

	/* support size parameters */
	dArray2DT neighbor_support;
	nVariArray2DT<double> neighbor_support_man(0, neighbor_support, 1);

	/* MLS fit weighting */
	dArrayT neighbor_volume;
	VariArrayT<double> neighbor_volume_man(0, neighbor_volume);

	/* compute weights for all cell nodes */
	dArrayT x_node;
	iArrayT neighbors;
	const iArrayT& cell_nodes = cell_data.CellNodes();
	iArrayT cell_nodes_OK(cell_nodes.Length());
	cell_nodes_OK = 0; /* not OK */
	for (int i = 0; i < cell_nodes.Length(); i++)
	{
		/* nodal neighbors */
		nodal_neighbors.RowAlias(i,neighbors);
		
		/* dimension work space */
		int nngh = neighbors.Length();
		neighbor_coords_man.SetMajorDimension(nngh, false);
		neighbor_support_man.SetMajorDimension(nngh, false);
		neighbor_volume_man.SetLength(nngh, false);
		neighbor_volume = 1.0;
		
		/* collect data */
		neighbor_coords.RowCollect(neighbors, point_coordinates);
		neighbor_support.RowCollect(neighbors, fSupportParams);
		cell_coordinates.RowAlias(cell_nodes[i], x_node);
		
		/* compute MLS fit */
		if (fMLS->SetField(neighbor_coords, neighbor_support, neighbor_volume, x_node, 0)) 
		{
			/* store MLS fit weights */
			neighbor_weights.SetRow(i, fMLS->phi());

			/* set flag */
			cell_nodes_OK[i] = 1;
		}
		else if (ElementSupport().Logging() != GlobalT::kSilent) /* report nodes that could not be fit */
		{
			/* write support size of the neighborhood nodes */
			bool write_support_size = true;
			if (write_support_size)
			{
				StringT file;
				file.Root(ElementSupport().InputFile());
				file.Append(".MLS", cell_nodes[i]+1);
				file.Append(".out");

				dArray2DT coords_used(neighbor_coords.MajorDim()+1, neighbor_coords.MinorDim());
				coords_used.SetRow(0, x_node);
				coords_used.BlockRowCopyAt(neighbor_coords, 1);

				iArrayT points_used(neighbor_coords.MajorDim()+1);
				points_used[0] = cell_nodes[i];
				points_used.CopyIn(1, neighbors);

				ArrayT<StringT> labels(1);
				labels[0] = "r";

				dArray2DT n_values(coords_used.MajorDim(), 1);
				n_values[0];
				n_values.BlockRowCopyAt(neighbor_support, 1);

				ElementSupport().WriteOutput(file, coords_used, points_used, n_values, labels);
			}
			
			bool write_node_connectivities = true;
			if (write_node_connectivities)
			{
				/* stream */
				ostream& out = ElementSupport().Output();
				int node = cell_nodes[i];
				out << "elements containing node: " << node+1 << '\n';
			
				/* the element group */
				const ContinuumElementT* element_group = cell_data.ContinuumElement();

				/* point in cell information */
				RaggedArray2DT<int>& point_in_cell = cell_data.PointInCell();
			
				/* elements containing the node */
				iArrayT connects_tmp, points_tmp;
				for (int i = 0; i < element_group->NumElements(); i++) {
					const iArrayT& connects = element_group->ElementCard(i).NodesX();
					if (connects.HasValue(node)) {
						point_in_cell.RowAlias(i, points_tmp);
						if (points_tmp.Length() > 0) {
							connects_tmp.Alias(connects);
							connects_tmp++;
							out << "element " << i+1 << ": " << connects_tmp.no_wrap() << '\n';
							connects_tmp--;
							
							points_tmp++;
							out << "points: " << points_tmp.no_wrap() << '\n';
							points_tmp--;
						}
					}
				}
			}
		}		
	}
	
	/* filter out nodes that could not be fit */
	int num_nodes_OK = cell_nodes_OK.Count(1);
	if (num_nodes_OK != cell_nodes.Length()) {

		/* collect nodes with good fits */
		iArrayT nodes_OK(num_nodes_OK);
		iArrayT nodes_OK_dim(num_nodes_OK);
		int count = 0;
		for (int i = 0; i < cell_nodes.Length(); i++)
			if (cell_nodes_OK[i] == 1) {
				nodes_OK[count] = cell_nodes[i];
				nodes_OK_dim[count] = nodal_neighbors.MinorDim(i);
				count++;	
			}

		/* filter down neighbor data */
		RaggedArray2DT<int>& nodal_neighbors = point_to_node.Neighbors();
		RaggedArray2DT<double>& neighbor_weights = point_to_node.NeighborWeights();
		RaggedArray2DT<int> nodal_neighbors_tmp = nodal_neighbors;
		RaggedArray2DT<double> neighbor_weights_tmp = neighbor_weights;
		nodal_neighbors.Configure(nodes_OK_dim);
		neighbor_weights.Configure(nodes_OK_dim);
		count = 0;
		for (int i = 0; i < cell_nodes.Length(); i++)
			if (cell_nodes_OK[i] == 1) {
				nodal_neighbors.SetRow(count, nodal_neighbors_tmp(i));
				neighbor_weights.SetRow(count, neighbor_weights_tmp(i));
				count++;
			}
		
		/* reset the map */
		point_to_node.Map().SetMap(nodes_OK);
	}
	
	/* verbose output */
	if (ElementSupport().Logging() == GlobalT::kVerbose)
	{
		/* stream */
		ostream& out = ElementSupport().Output();

		/* neighbor weights */
		out << "\n Nodal neighborhoods weights:\n";
		neighbor_weights.WriteNumbered(out);
	}	

	/* build neighborhoods for projecting points */
	BuildPointNeighborhoods(points_used, point_coordinates, cell_data);

	/* compute weights in point neighborhoods */
	InterpolationDataT& point_to_point = cell_data.PointToPoint();
	const RaggedArray2DT<int>& point_neighbors = point_to_point.Neighbors();
	neighbor_count.Dimension(point_neighbors.MajorDim());
	point_neighbors.MinorDim(neighbor_count);
	RaggedArray2DT<double>& point_neighbor_weights = point_to_point.NeighborWeights();
	point_neighbor_weights.Configure(neighbor_count);
	neighbor_count.Free();
	for (int i = 0; i < points_used.Length(); i++)
	{
		/* nodal neighbors */
		point_neighbors.RowAlias(i, neighbors);
		
		/* dimension work space */
		int nngh = neighbors.Length();
		neighbor_coords_man.SetMajorDimension(nngh, false);
		neighbor_support_man.SetMajorDimension(nngh, false);
		neighbor_volume_man.SetLength(nngh, false);
		neighbor_volume = 1.0;
		
		/* collect data */
		neighbor_coords.RowCollect(neighbors, point_coordinates);
		neighbor_support.RowCollect(neighbors, fSupportParams);
		point_coordinates.RowAlias(points_used[i], x_node);
		
		if (nngh < 2) /* some image points will be off the grid */
			point_neighbor_weights.SetRow(i, 0.0);
		else {
			/* compute MLS fit */
			if (!fMLS->SetField(neighbor_coords, neighbor_support, neighbor_volume, x_node, 0))
				ExceptionT::GeneralFail(caller, "could not compute MLS fit for point %d", points_used[i]+1);
		
			/* store weights */
			point_neighbor_weights.SetRow(i, fMLS->phi());
		}
	}

	/* compute node to node projection matrix */
	Compute_B_hatU_U(cell_data, cell_data.NodeToNode());
}

/* project the point values onto the mesh */
void MeshfreeBridgingT::ProjectField(const PointInCellDataT& cell_data,
	const dArray2DT& point_values, dArray2DT& projection)
{
	/* projected part of the mesh */
	const iArrayT& cell_nodes = cell_data.CellNodes();

	/* nodal neighbor data */
	const InterpolationDataT& point_to_node = cell_data.PointToNode();	
	const RaggedArray2DT<int>& nodal_neighbors = point_to_node.Neighbors();
	const RaggedArray2DT<double>& neighbor_weights = point_to_node.NeighborWeights();

	/* initialize return value */
	projection.Dimension(nodal_neighbors.MajorDim(), point_values.MinorDim());
	projection = 0.0;

	/* neighborhood values */
	LocalArrayT loc_values;
	VariLocalArrayT loc_values_man(0, loc_values, point_values.MinorDim());
	loc_values_man.SetNumberOfNodes(0); /* loc_values must be dimensioned before SetGlobal */
	loc_values.SetGlobal(point_values);

	/* calculate "projection" using weights */
	iArrayT neighbors;
	dArrayT weights;
	for (int i = 0; i < projection.MajorDim(); i++)
	{
		/* fetch neighbor data */
		nodal_neighbors.RowAlias(i, neighbors);
		neighbor_weights.RowAlias(i, weights);
		
		/* collect neighbor values */
		loc_values_man.SetNumberOfNodes(neighbors.Length());
		loc_values.SetLocal(neighbors);
		
		/* compute components of projection */
		for (int j = 0; j < projection.MinorDim(); j++)
			projection(i,j) = dArrayT::Dot(loc_values(j), weights);
	}
}

/* compute the coarse scale part of the source field */
void MeshfreeBridgingT::CoarseField(const PointInCellDataT& cell_data, const dArray2DT& field, 
	dArray2DT& coarse) const
{
	/* point neighbor data */
	const InterpolationDataT& point_to_point = cell_data.PointToPoint();	
	const RaggedArray2DT<int>& point_neighbors = point_to_point.Neighbors();
	const RaggedArray2DT<double>& point_neighbor_weights = point_to_point.NeighborWeights();

	/* initialize return value */
	coarse.Dimension(point_neighbors.MajorDim(), field.MinorDim());
	coarse = 0.0;

	/* neighborhood values */
	LocalArrayT loc_values;
	VariLocalArrayT loc_values_man(0, loc_values, field.MinorDim());
	loc_values_man.SetNumberOfNodes(0); /* loc_values must be dimensioned before SetGlobal */
	loc_values.SetGlobal(field);

	/* calculate "projection" using weights */
	iArrayT neighbors;
	dArrayT weights;
	for (int i = 0; i < coarse.MajorDim(); i++)
	{
		/* fetch neighbor data */
		point_neighbors.RowAlias(i, neighbors);
		point_neighbor_weights.RowAlias(i, weights);
		
		/* collect neighbor values */
		loc_values_man.SetNumberOfNodes(neighbors.Length());
		loc_values.SetLocal(neighbors);
		
		/* compute components of projection */
		for (int j = 0; j < coarse.MinorDim(); j++)
			coarse(i,j) = dArrayT::Dot(loc_values(j), weights);
	}
}

/* collect the cells without any free nodes */
void MeshfreeBridgingT::CollectProjectedCells(const PointInCellDataT& cell_data, iArrayT& cells) const
{
//TEMP - for now just assume a cell that contains any projecting points has no free nodes
	BridgingScaleT::CollectProjectedCells(cell_data, cells);
}

/* return list of projected nodes */
void MeshfreeBridgingT::CollectProjectedNodes(const PointInCellDataT& cell_data, iArrayT& nodes) const
{
	const InterpolationDataT& point_to_node = cell_data.PointToNode();
	const InverseMapT& driven_node_map = point_to_node.Map();
	driven_node_map.Forward(nodes);
}

/* compute \f$ B_{\hat{U}U} \f$ */
void MeshfreeBridgingT::Compute_B_hatU_U(const PointInCellDataT& projection, 
	InterpolationDataT& B_hatU_U) const
{
	const char caller[] = "MeshfreeBridgingT::Compute_B_hatU_U";

	/* N_Q_U */	
	const InverseMapT& Q_global_to_local = projection.GlobalToLocal();
	const iArray2DT& N_cell_connectivities = projection.CellConnectivities();
	const iArrayT& N_cell = projection.InterpolatingCell();
	const dArray2DT& N_cell_weights = projection.InterpolationWeights();

	/* B_hatU_Q */
	const InterpolationDataT& B_hatU_Q = projection.PointToNode();
	const RaggedArray2DT<int>& B_hatU_Q_neighbors = B_hatU_Q.Neighbors();
	const RaggedArray2DT<double>& B_hatU_Q_weights = B_hatU_Q.NeighborWeights();
	
	/* B_hatU_Q and B_hatU_U row maps */
	InverseMapT& B_hatU_U_row_map = B_hatU_U.Map();
	B_hatU_U_row_map = B_hatU_Q.Map(); /* copy */
	B_hatU_U_row_map.SetOutOfRange(InverseMapT::MinusOne);
	iArrayT nodes_Uhat;
	B_hatU_U_row_map.Forward(nodes_Uhat);
	
	/* coarse scale element group */
	const SolidElementT& solid = SolidElement();
	int nen = solid.NumElementNodes();

	/* compute B_hatU_U = B_hatU_Q x N_Q_U */
	AutoFill2DT<int> B_hatU_U_neighbors_tmp(nodes_Uhat.Length(), 1, 10, 10);
	AutoFill2DT<double> B_hatU_U_weights_tmp(nodes_Uhat.Length(), 1, 10, 10);
	iArrayT hatU_neighbors;
	for (int i = 0; i < nodes_Uhat.Length(); i++) /* loop over projected nodes */ 
	{		
		/* loop over points contributing to each projected node */
		B_hatU_Q_neighbors.RowAlias(i, hatU_neighbors);
		for (int j = 0; j < hatU_neighbors.Length(); j++) 
		{
			int j_loc = Q_global_to_local.Map(hatU_neighbors[j]);

			/* cell containing projection source point */
			int element = N_cell[j_loc];
			const iArrayT& nodes = SolidElement().ElementCard(element).NodesU();
			
			/* collect free cell nodes/constract over Q relating hatU and U */
			double B = B_hatU_Q_weights(i,j);
			for (int k = 0; k < nodes.Length(); k++)
			{
				int node = nodes[k];
				if (B_hatU_U_row_map.Map(node) == -1) /* not a hatU node */
				{
					int row_index = B_hatU_U_neighbors_tmp.PositionInRow(i, node);
					double BxN = B*N_cell_weights(j_loc,k);
					if (row_index == -1) /* new value */
					{
						B_hatU_U_neighbors_tmp.Append(i, node);
						B_hatU_U_weights_tmp.Append(i, BxN);
					}	
					else /* existing value - matrix multiply */
						B_hatU_U_weights_tmp(i, row_index) += BxN;
				}
			}
		}
	}

	/* B_hatU_U */
	B_hatU_U.Neighbors().Copy(B_hatU_U_neighbors_tmp);
	B_hatU_U.NeighborWeights().Copy(B_hatU_U_weights_tmp);
	B_hatU_U_row_map.SetOutOfRange(InverseMapT::Throw); /* allow no errors */
}

/* information about subordinate parameter lists */
void MeshfreeBridgingT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	BridgingScaleT::DefineSubs(sub_list);

	/* meshfree shape functions */
	sub_list.AddSub("RKPM");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MeshfreeBridgingT::NewSub(const StringT& name) const
{
	if (name == "RKPM")
		return fMeshFreeSupport.NewSub(name);
	else /* inherited */
		return BridgingScaleT::NewSub(name);
}
	
/* accept parameter list */
void MeshfreeBridgingT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	BridgingScaleT::TakeParameterList(list);

	/* construct the MLS solver */
	const ParameterListT& rkpm = list.GetList("RKPM");
	int completeness = rkpm.GetParameter("completeness");
	int cross_terms = rkpm.GetParameter("cross_terms");
	const ParameterListT& window = rkpm.GetListChoice(fMeshFreeSupport, "window_function_choice");		
	fMLS = MeshFreeSupportT::New_MLSSolverT(NumSD(), completeness, cross_terms, window);
	fMLS->Initialize();
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* determines points in the neighborhoods of nodes of each non-empty cell */
void MeshfreeBridgingT::BuildNodalNeighborhoods(CommManagerT& comm, const iArrayT& points_used, const dArray2DT* init_coords, 
	const dArray2DT* curr_coords, PointInCellDataT& cell_data)
{
	const char caller[] = "MeshfreeBridgingT::BuildNodalNeighborhoods";

	/* map all the points into cells */
//	MaptoCells(points_used, init_coords, curr_coords, cell_data);	

	/* init all data structures for interpolating values to the points */
	InitInterpolation(points_used, init_coords, curr_coords, cell_data);

	/* set map of node used to rows in neighbor data */
	cell_data.CollectCellNodes();
	const iArrayT& nodes_used = cell_data.CellNodes();
	InterpolationDataT& point_to_node = cell_data.PointToNode();
	InverseMapT& node_to_neighbor_data = point_to_node.Map();
	node_to_neighbor_data.SetMap(nodes_used);
	node_to_neighbor_data.SetOutOfRange(InverseMapT::MinusOne);

	/* point coordinates */
	if (curr_coords && init_coords) ExceptionT::GeneralFail(caller, "cannot pass both init and curr coords");
	if (!curr_coords && !init_coords) ExceptionT::GeneralFail(caller, "must define init or curr coords");
	const dArray2DT& point_coordinates = (init_coords != NULL) ? *init_coords : *curr_coords;

	/* cell coordinates */
	const dArray2DT& cell_coordinates = (init_coords != NULL) ?
		ElementSupport().InitialCoordinates() : ElementSupport().CurrentCoordinates();
	LocalArrayT::TypeT coord_type = (init_coords != NULL) ? 
		LocalArrayT::kInitCoords : LocalArrayT::kCurrCoords;

	/* compute (average) support size of each node */
	int nel = SolidElement().NumElements();	
	iArrayT counts(nodes_used.Length());
	counts = 0;
	dArray2DT support_param(nodes_used.Length(), 1);
	support_param = 0.0;
	const ParentDomainT& parent = ShapeFunction().ParentDomain();
	dArrayT centroid;
	LocalArrayT loc_cell_coords(coord_type, SolidElement().NumElementNodes(), NumSD());
	loc_cell_coords.SetGlobal(cell_coordinates);
	node_to_neighbor_data.SetOutOfRange(InverseMapT::MinusOne);
	for (int i = 0; i < nel; i++) {

		/* element nodes */
		const iArrayT& element_nodes = SolidElement().ElementCard(i).NodesX();
	
		/* gives domain (global) nodal coordinates */
		loc_cell_coords.SetLocal(element_nodes);

		/* centroid and radius */
		double radius = parent.AverageRadius(loc_cell_coords, centroid);

		/* contribution to support size of element nodes */
		for (int j = 0; j < element_nodes.Length(); j++)
		{
			int dex = node_to_neighbor_data.Map(element_nodes[j]);
			if (dex != -1) /* is in nodes_used */
			{
				counts[dex]++;
				support_param[dex] += radius;
			}
		}
	}
	node_to_neighbor_data.SetOutOfRange(InverseMapT::Throw);
	
	/* compute support size from average element size */
	for (int i = 0; i < counts.Length(); i++)
		if (counts[i] > 0)
			/* support parameters tend to be about twice the near-neighbor distance */
			support_param[i] /= counts[i]; 
		else
			ExceptionT::GeneralFail(caller, "could not compute suppose size for node %d",
				nodes_used[i]+1);

	/* configure search grid */
	iGridManagerT grid(10, 100, point_coordinates, &points_used);
	grid.Reset();

	/* support parameters */
	fSupportParams.Dimension(point_coordinates.MajorDim(), 1);
	fSupportParams = 0.0;
	dArrayT support_size(support_param.Length());
	fMLS->SphericalSupportSize(support_param, support_size);
	
	/* collect neighborhood nodes and set support size */
	InverseMapT& global_to_local = cell_data.GlobalToLocal();
	AutoFill2DT<int> auto_fill(nodes_used.Length(), 1, 10, 10);
	dArrayT x_node, x_point;
	dArrayT nodal_params;
	for (int i = 0; i < nodes_used.Length(); i++)
	{
		/* candidate points */
		cell_coordinates.RowAlias(nodes_used[i], x_node);
		const AutoArrayT<iNodeT>& hits = grid.HitsInRegion(x_node.Pointer(), support_size[i]);

		/* test all hits */
		for (int j = 0; j < hits.Length(); j++)
		{
			/* atom info */
			int point = hits[j].Tag();
			point_coordinates.RowAlias(point, x_point);

			/* add to neighbor list */
			nodal_params.Alias(1, support_param.Pointer(i));
			if (fMLS->Covers(x_node, x_point, nodal_params))
			{
				auto_fill.Append(i, point);
				
				/* take max support */
				fSupportParams[point] = (support_param[i] > fSupportParams[point]) ? 
					support_param[i] : fSupportParams[point];
			}
		}
	}
	
	/* distribute the support sizes */
	int id = comm.Init_AllGather(fSupportParams);
	comm.AllGather(id, fSupportParams);
	comm.Clear_AllGather(id);

	/* copy/compress contents */
	RaggedArray2DT<int>& nodal_neighbors = point_to_node.Neighbors();	
	nodal_neighbors.Copy(auto_fill);
	
	/* verbose output */
	if (ElementSupport().Logging() == GlobalT::kVerbose)
	{
		/* stream */
		ostream& out = ElementSupport().Output();

		/* neighbors */
		out << "\n Nodal neighborhoods:\n";
		iArrayT tmp(nodal_neighbors.Length(), nodal_neighbors.Pointer());
		tmp++;
		iArrayT neighbors;
		for (int i = 0; i < nodes_used.Length(); i++) {
			nodal_neighbors.RowAlias(i, neighbors);
			out << setw(kIntWidth) << nodes_used[i]+1 << ":" << neighbors.no_wrap() << '\n';
		}
		tmp--;
		
		/* support sizes for each source point */
		bool write_support_size = true;
		if (write_support_size) {
			StringT file;
			file.Root(ElementSupport().InputFile());
			file.Append(".support.out");

			dArray2DT coords_used(points_used.Length(), point_coordinates.MinorDim());
			coords_used.RowCollect(points_used, point_coordinates);

			ArrayT<StringT> labels(1);
			labels[0] = "r";

			dArray2DT n_values(points_used.Length(), 1);
			n_values.Collect(points_used, fSupportParams);
		
			ElementSupport().WriteOutput(file, coords_used, points_used, n_values, labels);			
		}

		/* number of neighbors for each projected node */
		bool write_num_neighbors = true;
		if (write_num_neighbors) {
			StringT file;
			file.Root(ElementSupport().InputFile());
			file.Append(".neighbors.out");

			/* nodal coordinates */
			const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
			dArray2DT coords_used(nodes_used.Length(), init_coords.MinorDim());
			coords_used.RowCollect(nodes_used, init_coords);

			ArrayT<StringT> labels(1);
			labels[0] = "n";

			dArray2DT n_values(nodes_used.Length(), 1);
			for (int i = 0; i < n_values.MajorDim(); i++)
				n_values[i] = double(nodal_neighbors.MinorDim(i));
		
			ElementSupport().WriteOutput(file, coords_used, nodes_used, n_values, labels);
		}
	}
}

/* determines points in the neighborhoods of nodes of each non-empty cell */
void MeshfreeBridgingT::BuildPointNeighborhoods(const iArrayT& points_used, const dArray2DT& point_coords,
	PointInCellDataT& cell_data)
{
	const char caller[] = "MeshfreeBridgingT::BuildPointNeighborhoods";

	/* set map of points used to rows in neighbor data */
	InterpolationDataT& point_to_point = cell_data.PointToPoint();
	InverseMapT& point_to_neighbor_data = point_to_point.Map();
	point_to_neighbor_data.SetMap(points_used);
	point_to_neighbor_data.SetOutOfRange(InverseMapT::MinusOne);

	/* set up search grid */
	iGridManagerT grid(10, 100, point_coords, &points_used);
	grid.Reset();

	/* collect neighborhood points */
	AutoFill2DT<int> auto_fill(points_used.Length(), 1, 10, 10);
	dArrayT nodal_params;
	dArrayT x_node, x_point;
	for (int i = 0; i < points_used.Length(); i++)
	{
	
		/* candidate points */
		int point = points_used[i];	
		fSupportParams.RowAlias(point, nodal_params);
		double support_size = fMLS->SphericalSupportSize(nodal_params);
		point_coords.RowAlias(point, x_node);
		const AutoArrayT<iNodeT>& hits = grid.HitsInRegion(x_node.Pointer(), support_size);

		/* test all hits */
		for (int j = 0; j < hits.Length(); j++)
		{
			/* point info */
			int hit = hits[j].Tag();
			point_coords.RowAlias(hit, x_point);

			/* add to neighbor list */
			if (fMLS->Covers(x_node, x_point, nodal_params))
				auto_fill.Append(i, hit);
		}
	}

	/* copy/compress contents */
	RaggedArray2DT<int>& point_neighbors = point_to_point.Neighbors();	
	point_neighbors.Copy(auto_fill);	
}
