/* $Id: MeshfreeBridgingT.cpp,v 1.2.12.1 2003-10-16 12:51:29 paklein Exp $ */
#include "MeshfreeBridgingT.h"

#include "PointInCellDataT.h"
#include "ShapeFunctionT.h"
#include "iGridManagerT.h"
#include "iNodeT.h"
#include "InverseMapT.h"
#include "MLSSolverT.h"
#include "VariLocalArrayT.h"
#include "VariArrayT.h"
#include "ofstreamT.h"

using namespace Tahoe;

/* constructor */
MeshfreeBridgingT::MeshfreeBridgingT(const ElementSupportT& support, const FieldT& field,
	const SolidElementT& solid):
	BridgingScaleT(support, field, solid),
	fMLS(NULL)
{
	/* MLS solver (using Gaussian window function) */
	dArrayT window_params(3);
	window_params[0] = 1.5; /* support size scaling    */
	window_params[1] = 0.4; /* sharpening factor       */
	window_params[2] = 2.0; /* neighbor cut-off factor */
	fMLS = new MLSSolverT(NumSD(), 1, MeshFreeT::kGaussian, window_params);
	fMLS->Initialize();
}

/* destructor */
MeshfreeBridgingT::~MeshfreeBridgingT(void) { delete fMLS; }

/* initialize projection data */
void MeshfreeBridgingT::InitProjection(const iArrayT& points_used, const dArray2DT* init_coords, 
	const dArray2DT* curr_coords, PointInCellDataT& cell_data)
{
	const char caller[] = "MeshfreeBridgingT::InitProjection";

	/* collect point within each nodal neighborhood */
	BuildNodalNeighborhoods(points_used, init_coords, curr_coords, cell_data);

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
	const RaggedArray2DT<int>& nodal_neighbors = cell_data.NodalNeighbors();
	iArrayT neighbor_count(nodal_neighbors.MajorDim());
	nodal_neighbors.MinorDim(neighbor_count);
	RaggedArray2DT<double>& neighbor_weights = cell_data.NodalNeighborWeights();
	neighbor_weights.Configure(neighbor_count);
	
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
		neighbor_support.RowCollect(neighbors, fSupport);
		cell_coordinates.RowAlias(cell_nodes[i], x_node);
		
		/* compute MLS fit */
		if (!fMLS->SetField(neighbor_coords, neighbor_support, neighbor_volume, x_node, 0))
			ExceptionT::GeneralFail(caller, "could not compute MLS fit for node %d", cell_nodes[i]+1);
	
		/* store weights */
		neighbor_weights.SetRow(i, fMLS->phi());
	}
	
	/* verbose output */
	if (ElementSupport().PrintInput())
	{
		/* stream */
		ostream& out = ElementSupport().Output();

		/* neighbor weights */
		out << "\n Nodal neighborhoods weights:\n";
		const RaggedArray2DT<double>& neighbors_weights = cell_data.NodalNeighborWeights();
		neighbors_weights.WriteNumbered(out);
	}	
}

/* project the point values onto the mesh */
void MeshfreeBridgingT::ProjectField(const PointInCellDataT& cell_data,
	const dArray2DT& point_values, dArray2DT& projection)
{
	/* projected part of the mesh */
	const iArrayT& cell_nodes = cell_data.CellNodes();

	/* nodal neighbor data */
	const RaggedArray2DT<int>& nodal_neighbors = cell_data.NodalNeighbors();
	const RaggedArray2DT<double>& neighbor_weights = cell_data.NodalNeighborWeights();

	/* initialize return value */
	projection.Dimension(cell_nodes.Length(), point_values.MinorDim());
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

/***********************************************************************
 * Protected
 ***********************************************************************/

/* determines points in the neighborhoods of nodes of each non-empty cell */
void MeshfreeBridgingT::BuildNodalNeighborhoods(const iArrayT& points_used, const dArray2DT* init_coords, 
	const dArray2DT* curr_coords, PointInCellDataT& cell_data)
{
	const char caller[] = "MeshfreeBridgingT::BuildNodalNeighborhoods";

	/* map all the points into cells */
	MaptoCells(points_used, init_coords, curr_coords, cell_data);	

	/* set map of node used to rows in neighbor data */
	cell_data.CollectCellNodes();
	const iArrayT& nodes_used = cell_data.CellNodes();
	InverseMapT node_to_neighbor_data = cell_data.NodeToNeighborData();
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
	iArrayT counts(nodes_used.Length());
	counts = 0;
	dArrayT support_size(nodes_used.Length());
	support_size = 0.0;
	const ParentDomainT& parent = ShapeFunction().ParentDomain();
	dArrayT centroid;
	LocalArrayT loc_cell_coords(coord_type, fSolid.NumElementNodes(), NumSD());
	loc_cell_coords.SetGlobal(cell_coordinates);
	node_to_neighbor_data.SetOutOfRange(InverseMapT::MinusOne);
	for (int i = 0; i < fSolid.NumElements(); i++) {

		/* element nodes */
		const iArrayT& element_nodes = fSolid.ElementCard(i).NodesX();
	
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
				support_size[dex] += radius;
			}
		}
	}
	node_to_neighbor_data.SetOutOfRange(InverseMapT::Throw);
	
	/* compute support size from average element size */
	for (int i = 0; i < counts.Length(); i++)
		if (counts[i] > 0)
			support_size[i] /= counts[i];
		else
			ExceptionT::GeneralFail(caller, "could not compute suppose size for node %d",
				nodes_used[i]+1);

	/* configure search grid */
	iGridManagerT grid(10, 100, point_coordinates, &points_used);
	grid.Reset();

	/* collect neighborhood nodes and set support size */
	fSupport.Dimension(point_coordinates.MajorDim(), 1);
	fSupport = 0.0;
	InverseMapT& global_to_local = cell_data.GlobalToLocal();
	AutoFill2DT<int> auto_fill(nodes_used.Length(), 1, 10, 10);
	dArrayT nodal_params(1);
	dArrayT x_node, x_point;
	for (int i = 0; i < nodes_used.Length(); i++)
	{
		/* candidate points */
		nodal_params[0] = support_size[i];
		cell_coordinates.RowAlias(nodes_used[i], x_node);
		const AutoArrayT<iNodeT>& hits = grid.HitsInRegion(x_node.Pointer(), nodal_params[0]);

		/* test all hits */
		for (int j = 0; j < hits.Length(); j++)
		{
			/* atom info */
			int point = hits[j].Tag();
			point_coordinates.RowAlias(point, x_point);

			/* add to neighbor list */
			if (fMLS->Covers(x_node, x_point, nodal_params))
			{
				auto_fill.Append(i, point);
				
				/* take max support */
				fSupport[point] = (support_size[i] > fSupport[point]) ? 
					support_size[i] : fSupport[point];
			}
		}
	}

	/* copy/compress contents */
	cell_data.NodalNeighbors().Copy(auto_fill);
	
	/* verbose output */
	if (ElementSupport().PrintInput())
	{
		/* stream */
		ostream& out = ElementSupport().Output();

		/* neighbors */
		out << "\n Nodal neighborhoods:\n";
		const RaggedArray2DT<int>& nodal_neighbors = cell_data.NodalNeighbors();
		iArrayT tmp(nodal_neighbors.Length(), nodal_neighbors.Pointer());
		tmp++;
		nodal_neighbors.WriteNumbered(out);
		tmp--;
	}
}
