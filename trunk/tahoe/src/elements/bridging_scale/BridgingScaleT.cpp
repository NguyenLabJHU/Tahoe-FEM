/* $Id: BridgingScaleT.cpp,v 1.52 2005-04-13 22:15:35 paklein Exp $ */
#include "BridgingScaleT.h"

#include <iostream.h>
#include <iomanip.h>

#include "SolidElementT.h"
#include "ShapeFunctionT.h"
#include "ofstreamT.h"
#include "iAutoArrayT.h"
#include "OutputSetT.h"
#include "AutoFill2DT.h"
#include "RaggedArray2DT.h"
#include "iGridManagerT.h"
#include "iNodeT.h"
#include "nArrayGroupT.h"
#include "PointInCellDataT.h"
#include "CCSMatrixT.h"

using namespace Tahoe;

/* constructor */
BridgingScaleT::BridgingScaleT(const ElementSupportT& support):
	ElementBaseT(support),
	fSolid(NULL),
	fElMatU(ElementMatrixT::kSymmetric),
	fLocInitCoords(LocalArrayT::kInitCoords),
	fLocDisp(LocalArrayT::kDisp),
	//fGlobalMass(support.Output(), 1, support.Communicator())
	fGlobalMass(NULL)
{
	SetName("bridging");
}

/* destructor */
BridgingScaleT::~BridgingScaleT(void) {
	delete fGlobalMass;
}

/* map coordinates into elements */
void BridgingScaleT::MaptoCells(const iArrayT& points_used, const dArray2DT* init_coords, 
	const dArray2DT* curr_coords, PointInCellDataT& cell_data)
{
	const char caller[] = "BridgingScaleT::MaptoCells";

	/* global to local map */
	InverseMapT& global_to_local = cell_data.GlobalToLocal();
	global_to_local.SetMap(points_used);

	/* map data */
	cell_data.SetContinuumElement(SolidElement());
	RaggedArray2DT<int>& point_in_cell = cell_data.PointInCell();
	RaggedArray2DT<double>& point_in_cell_coords = cell_data.PointInCellCoords();

	/* point coordinates */
	if (curr_coords && init_coords) ExceptionT::GeneralFail(caller, "cannot pass both init and curr coords");
	if (!curr_coords && !init_coords) ExceptionT::GeneralFail(caller, "must define init or curr coords");
	const dArray2DT& point_coordinates = (init_coords != NULL) ? *init_coords : *curr_coords;

	/* cell coordinates */
	const dArray2DT& cell_coordinates = (init_coords != NULL) ?
		ElementSupport().InitialCoordinates() :
		ElementSupport().CurrentCoordinates();
	LocalArrayT::TypeT coord_type = (init_coords != NULL) ? 
		LocalArrayT::kInitCoords : 
		LocalArrayT::kCurrCoords;

	/* stream */
	ostream& out = ElementSupport().Output();

	/* configure search grid */
	iGridManagerT grid(10, 100, point_coordinates, &points_used);
	grid.Reset();

	/* verbose output */
	if (ElementSupport().Logging() == GlobalT::kVerbose) {
		grid.WriteStatistics(out);
		grid.DumpGrid(out);
	}

	/* track cell containing each point, so only one cell is associated with each point */
	iArrayT found_in_cell(points_used.Length());
	found_in_cell = -1;

	/* check all cells for points */
	int nel = SolidElement().NumElements();
	const ParentDomainT& parent = ShapeFunction().ParentDomain();
	AutoFill2DT<int> auto_fill(nel, 1, 10, 10);
	dArrayT x_atom, centroid;
	LocalArrayT loc_cell_coords(coord_type, SolidElement().NumElementNodes(), NumSD());
	loc_cell_coords.SetGlobal(cell_coordinates);
	for (int i = 0; i < nel; i++) {
	
		/* gives domain (global) nodal coordinates */
		loc_cell_coords.SetLocal(SolidElement().ElementCard(i).NodesX());

		/* centroid and radius */
		double radius = parent.AverageRadius(loc_cell_coords, centroid);

		/* candidate points */
		const AutoArrayT<iNodeT>& hits = grid.HitsInRegion(centroid.Pointer(), 1.01*radius);

		/* check if points are within the element domain */
		for (int j = 0; j < hits.Length(); j++)
		{
			int global = hits[j].Tag();
			int local = global_to_local.Map(global);
			
			/* not mapped yet */
			if (found_in_cell[local] == -1)
			{
				x_atom.Alias(NumSD(), hits[j].Coords());
				if (parent.PointInDomain(loc_cell_coords, x_atom)) 
				{
					found_in_cell[local] = i;
					auto_fill.Append(i, global);
				}
			}
		}
	}

	/* map should return -1 if point not in map */
	global_to_local.SetOutOfRange(InverseMapT::MinusOne);

	/* copy/compress contents */
	point_in_cell.Copy(auto_fill);
	auto_fill.Free();
	found_in_cell.Free();
	
	/* verbose output */
	if (ElementSupport().Logging() == GlobalT::kVerbose) {
		iArrayT tmp(point_in_cell.Length(), point_in_cell.Pointer());
		out << "\n Particles in cells:\n"
		    << setw(kIntWidth) << "no." << '\n';
		tmp++;
		point_in_cell.WriteNumbered(out);
		tmp--;
		out.flush();
	}

	if (point_in_cell.Length() > points_used.Length()) {
		cout << '\n' << caller << ": WARNING: number of particles in cells " 
		     << point_in_cell.Length() << " exceeds\n" 
		     <<   "     the total number of particles " << points_used.Length() << endl;
	}

	/* map points in every cell to parent domain coordinates */
	dArrayT mapped(NumSD()), point;
	AutoFill2DT<double> inverse(point_in_cell.MajorDim(), 1, 10, 10);
	for (int i = 0; i < point_in_cell.MajorDim(); i++) 
	{
		/* cell coordinates */
		loc_cell_coords.SetLocal(SolidElement().ElementCard(i).NodesX()); 

		/* run through list and map to parent domain */
		int* particles = point_in_cell(i);
		for (int j = 0; j < point_in_cell.MinorDim(i); j++) 
		{
			point_coordinates.RowAlias(particles[j], point);
			if (parent.MapToParentDomain(loc_cell_coords, point, mapped))
				inverse.Append(i, mapped);
			else
				ExceptionT::GeneralFail(caller, "mapping to parent domain failed");
		}
	}

	/* copy compress coordinate list */
	point_in_cell_coords.Copy(inverse);
	inverse.Free();
	
	/* verbose output */
	if (ElementSupport().Logging() == GlobalT::kVerbose) {

		int nsd = NumSD();
		out << "\n Mapped coordinates of particles in elements:\n";
		
		/* run though element and dump coordinates of particles in parent domain */
		dArray2DT inv_coords;
		for (int i = 0; i < point_in_cell.MajorDim(); i++) {
		
			int np = point_in_cell.MinorDim(i);
			out << "element = " << i+1 << '\n'
			    << "  count	= " << np << '\n';

			/* shallow copy of inverse coords */
			inv_coords.Set(np, nsd, point_in_cell_coords(i));

			/* write coordinates */
			int* p_list = point_in_cell(i);
			for (int j = 0; j < np; j++)
			{
				out << setw(kIntWidth) << p_list[j]+1 << ": ";
				inv_coords.PrintRow(j, out);
			}
		}
	}
}

const ShapeFunctionT& BridgingScaleT::ShapeFunction(void) const {
	return SolidElement().ShapeFunction();
}

/* initialize interpolation data */
void BridgingScaleT::InitInterpolation(const iArrayT& points_used, const dArray2DT* init_coords, 
	const dArray2DT* curr_coords, PointInCellDataT& cell_data)
{
	/* initialize the point-in-cell data */
	MaptoCells(points_used, init_coords, curr_coords, cell_data);
	
	/* dimension return value */
	dArray2DT& weights = cell_data.InterpolationWeights();
	weights.Dimension(points_used.Length(), SolidElement().NumElementNodes());
	iArrayT& cell = cell_data.InterpolatingCell();
	cell.Dimension(points_used.Length());
	cell = -1;
	
	/* global to local map */
	const InverseMapT& global_to_local = cell_data.GlobalToLocal();

	/* cell shape functions */
	int nsd = NumSD();
	const ParentDomainT& parent = ShapeFunction().ParentDomain();

	/* point in cells data */
	const RaggedArray2DT<int>& point_in_cell = cell_data.PointInCell();
	const RaggedArray2DT<double>& point_in_cell_coords = cell_data.PointInCellCoords();

	/* run through cells */
	dArrayT Na;
	dArrayT point_coords;
	dArray2DT mapped_coords;

	for (int i = 0; i < point_in_cell.MajorDim(); i++)
	{
		int np = point_in_cell.MinorDim(i);
		if (np > 0)
		{
			const int* points = point_in_cell(i);
			
			/* mapped coordinates of points in this cell */
			mapped_coords.Alias(np, nsd, point_in_cell_coords(i));
			
			/* run through points */
			for (int j = 0; j < np; j++)
			{
				int dex = global_to_local.Map(points[j]);
				if (cell[dex] == -1) /* not set yet */
				{
					cell[dex] = i;
					
					/* fetch point data */
					mapped_coords.RowAlias(j, point_coords);
					weights.RowAlias(dex, Na);
		
					/* evaluate shape functions */
					parent.EvaluateShapeFunctions(point_coords, Na);
				}
			}
		}
	}
}

void BridgingScaleT::InterpolateField(const StringT& field, int order, const PointInCellDataT& cell_data,
	dArray2DT& point_values) const
{
	int nen = SolidElement().NumElementNodes();

	/* get the field */
	const FieldT* the_field = ElementSupport().Field(field);
	LocalArrayT loc_field(LocalArrayT::kDisp, nen, the_field->NumDOF());
	loc_field.SetGlobal((*the_field)[order]); /* change order so can accomodate any field */
	
	/* interpolation data */
	const dArray2DT& weights = cell_data.InterpolationWeights();
	const iArrayT& cell = cell_data.InterpolatingCell();
	
	/* dimension return value */
	point_values.Dimension(weights.MajorDim(), the_field->NumDOF());
	
	/* loop over points */
	for (int i = 0; i < point_values.MajorDim(); i++)
	{
		/* element nodes */
		int element = cell[i];
		const iArrayT& nodes = SolidElement().ElementCard(element).NodesU();

		/* collect local values */
		loc_field.SetLocal(nodes);
		
		/* interpolate */
		for (int j = 0; j < point_values.MinorDim(); j++)
			point_values(i,j) = weights.DotRow(i, loc_field(j));
	}
}

/* compute the projection matrix */
void BridgingScaleT::InitProjection(CommManagerT& comm, const iArrayT& points_used, const dArray2DT* init_coords, 
	const dArray2DT* curr_coords, PointInCellDataT& cell_data)
{
#pragma unused(comm)

	/* compute interpolation data */
	InitInterpolation(points_used, init_coords, curr_coords, cell_data);

	/* collect nodes in non-empty cells and generate cell connectivities 
	 * in local numbering*/
	cell_data.GenerateCellConnectivities();

	/* projected part of the mesh */
	const iArrayT& cell_nodes = cell_data.CellNodes();
	const iArray2DT& cell_connects = cell_data.CellConnectivities();

	/* cell connectivities are (matrix equations) - 1 */
	iArrayT tmp_shift;
	tmp_shift.Alias(cell_connects);
	tmp_shift++;

	/* configure the matrix */
	if (!fGlobalMass) fGlobalMass = new CCSMatrixT(ElementSupport().Output(), 1, ElementSupport().Communicator());
	int num_projected_nodes = cell_nodes.Length();
	fGlobalMass->AddEquationSet(cell_connects);
	fGlobalMass->Initialize(num_projected_nodes, num_projected_nodes, 1);
	fGlobalMass->Clear();

	/* points in cell data */
	const RaggedArray2DT<int>& point_in_cell = cell_data.PointInCell();
	const dArray2DT& weights = cell_data.InterpolationWeights();
	const InverseMapT& global_to_local = cell_data.GlobalToLocal();

	/* loop over mesh */
	int cell_dex = 0;
	iArrayT cell_eq;
	dArrayT Na;
	for (int i = 0; i < point_in_cell.MajorDim(); i++)
	{
		int np = point_in_cell.MinorDim(i);
		if (np > 0)
		{
			fElMatU = 0.0;
			const int* points = point_in_cell(i);
			for (int j = 0; j < np; j++)
			{
				/* fetch interpolation weights */
				int point_dex = global_to_local.Map(points[j]);
				weights.RowAlias(point_dex, Na);
				//atommass = mdmass[point_dex];
				
				/* element mass */
				fElMatU.Outer(Na, Na, 1.0, dMatrixT::kAccumulate);
				//fElMatU *= atommass;	// need to multiply by atomic mass, i.e. M = N^{T}M_{A}N
			}

			/* equations of cell in projector */
			cell_connects.RowAlias(cell_dex++, cell_eq);

			/* assemble into matrix */
			fGlobalMass->Assemble(fElMatU, cell_eq);
		}
	}

	/* shift back */
	tmp_shift--;
	
//TEMP - write mass matrix to file
#if 0
ostream& out = ElementSupport().Output();
out << "\n weights =\n" << weights << endl;
out << "\n mass matrix =\n" << fGlobalMass << endl;
#endif
//TEMP
}

/* project the point values onto the mesh */
void BridgingScaleT::ProjectField(const PointInCellDataT& cell_data,
	const dArray2DT& point_values, dArray2DT& projection)
{
	/* projected part of the mesh */
	const iArrayT& cell_nodes = cell_data.CellNodes();
	const iArray2DT& cell_connects = cell_data.CellConnectivities();

	/* points in cell data */
	const RaggedArray2DT<int>& point_in_cell = cell_data.PointInCell();
	const dArray2DT& weights = cell_data.InterpolationWeights();
	const InverseMapT& global_to_local = cell_data.GlobalToLocal();
	
	/* initialize return value */
	projection.Dimension(cell_nodes.Length(), point_values.MinorDim());
	projection = 0.0;
	
	/* loop over mesh */
	int cell_dex = 0;
	iArrayT cell_eq;
	dArrayT Na, point_value;
	dMatrixT Nd(cell_connects.MinorDim(), point_values.MinorDim());
	for (int i = 0; i < point_in_cell.MajorDim(); i++)
	{
		int np = point_in_cell.MinorDim(i);
		if (np > 0)
		{
			const int* points = point_in_cell(i);
			Nd = 0.0;
			for (int j = 0; j < np; j++)
			{
				int point = points[j];
			
				/* fetch interpolation weights */
				int point_dex = global_to_local.Map(point);
				weights.RowAlias(point_dex, Na);

				/* source values of the point */
				point_values.RowAlias(point, point_value);
			
				/* rhs during projection - calculating part of w */
				Nd.Outer(Na, point_value, 1.0, dMatrixT::kAccumulate);
			}

			/* equations of cell in projector */
			cell_connects.RowAlias(cell_dex++, cell_eq);

			/* assemble */
			for (int j = 0; j < Nd.Cols(); j++)
				projection.Accumulate(j, cell_eq, Nd(j));
		}
	}
	
//TEMP - write mass matrix to file
#if 0
ostream& out = ElementSupport().Output();
out << "\n values =\n" << values << endl;
out << "\n residual =\n" << projection << endl;
#endif
//TEMP

	/* calculate projection - requires global matrix that supports 
	 * multiple solves - projection = w after operations within this loop */
	dArrayT u_tmp(projection.MajorDim());
	for (int i = 0; i < projection.MinorDim(); i++)
	{
		projection.ColumnCopy(i, u_tmp);
		fGlobalMass->Solve(u_tmp);
		projection.SetColumn(i, u_tmp);
	}
	u_tmp.Free();

	/* initialize return values */
	fFineScale.Dimension(point_values);
	fFineScale = 0.0;
	fCoarseScale.Dimension(point_values);
	fCoarseScale = 0.0;

	cell_dex = 0;
	iArrayT cell_connect;
	dArray2DT cell_projection(cell_connects.MinorDim(), projection.MinorDim());
	for (int i = 0; i < point_in_cell.MajorDim(); i++)
	{
		int np = point_in_cell.MinorDim(i);
		if (np > 0)
		{
			/* gather cell information */
			cell_connects.RowAlias(cell_dex++, cell_connect);
			cell_projection.RowCollect(cell_connect, projection);
		
			const int* points = point_in_cell(i);
			for (int j = 0; j < np; j++)
			{
				int point = points[j];
			
				/* fetch interpolation weights */
				int point_dex = global_to_local.Map(point);
				weights.RowAlias(point_dex, Na);

				/* interpolate to point */
				for (int k = 0; k < projection.MinorDim(); k++)
				{
					/* interpolate coarse scale to point (gives Nw) */
					fCoarseScale(point, k) = cell_projection.DotColumn(k, Na);

					/* error = source - projection = q-Nw*/
					fFineScale(point, k) = point_values(point, k) - fCoarseScale(point, k);
				}
			}
		}
	}
}

/* compute the coarse scale part of the source field */
void BridgingScaleT::CoarseField(const PointInCellDataT& cell_data, const dArray2DT& field, dArray2DT& coarse) const
{
	const char caller[] = "BridgingScaleT::CoarseField";

	/* projected part of the mesh */
	const iArrayT& cell_nodes = cell_data.CellNodes();
	const iArray2DT& cell_connects = cell_data.CellConnectivities();
	
	/* points in cell data */
	const RaggedArray2DT<int>& point_in_cell = cell_data.PointInCell();
	const dArray2DT& weights = cell_data.InterpolationWeights();
	const InverseMapT& global_to_local = cell_data.GlobalToLocal();
	
	/* initialize return value */
	coarse.Dimension(cell_nodes.Length(), field.MinorDim());
	coarse = 0.0;
	
	/* loop over mesh */
	int cell_dex = 0;
	iArrayT cell_eq;
	dArrayT Na, point_value;
	dMatrixT Nd(cell_connects.MinorDim(), field.MinorDim());
	//double atommass;
	for (int i = 0; i < point_in_cell.MajorDim(); i++)
	{
		int np = point_in_cell.MinorDim(i);
		if (np > 0)
		{
			const int* points = point_in_cell(i);
			Nd = 0.0;
			for (int j = 0; j < np; j++)
			{
				int point = points[j];
			
				/* fetch interpolation weights */
				int point_dex = global_to_local.Map(point);
				weights.RowAlias(point_dex, Na);
				//atommass = mdmass[point_dex];

				/* source values of the point */
				field.RowAlias(point, point_value);
			
				/* rhs during projection - calculating part of w */
				Nd.Outer(Na, point_value, 1.0, dMatrixT::kAccumulate);
				//Nd *= atommass;	// need to multiply by atomic mass, i.e. w = M^{-1}N^{T}M_{A}q
			}

			/* equations of cell in projector */
			cell_connects.RowAlias(cell_dex++, cell_eq);

			/* assemble */
			for (int j = 0; j < Nd.Cols(); j++)
				coarse.Accumulate(j, cell_eq, Nd(j));
		}
	}

	/* need non-const global matrix to solve */
	BridgingScaleT* non_const_this = const_cast<BridgingScaleT*>(this);

	/* calculate projection - requires global matrix that supports 
	 * multiple solves - projection = w after operations within this loop */
	dArrayT u_tmp(coarse.MajorDim());
	for (int i = 0; i < coarse.MinorDim(); i++)
	{
		coarse.ColumnCopy(i, u_tmp);
		non_const_this->fGlobalMass->Solve(u_tmp);
		coarse.SetColumn(i, u_tmp);
	}
}

/* collect the cells without any free nodes */
void BridgingScaleT::CollectProjectedCells(const PointInCellDataT& cell_data, iArrayT& cells) const
{
	/* mark cells */
	const RaggedArray2DT<int>& point_in_cell = cell_data.PointInCell();
	iArrayT projected_cell(point_in_cell.MajorDim());
	for (int i = 0; i < projected_cell.Length(); i++)
		if (point_in_cell.MinorDim(i) > 0)
			projected_cell[i] = 1;
		else
			projected_cell[i] = 0;

	/* collect filled cells */
	cells.Dimension(projected_cell.Count(1));
	int index = 0;
	for (int i = 0; i < projected_cell.Length(); i++)
		if (projected_cell[i])
			cells[index++] = i;
}

/* return list of projected nodes */
void BridgingScaleT::CollectProjectedNodes(const PointInCellDataT& cell_data, iArrayT& nodes) const {
	nodes = cell_data.CellNodes();
}

/* write projection-interpolation matrix from projection_data into cell_data */
void BridgingScaleT::ComputeProjectionInterpolation(
	const PointInCellDataT& cell_data,
	const PointInCellDataT& projection_data,
	const iArrayT& projection_source,
	const iArrayT& projection_dest) const
{
#pragma unused(cell_data)
#pragma unused(projection_data)
#pragma unused(projection_source)
#pragma unused(projection_dest)

	ExceptionT::GeneralFail("BridgingScaleT::ComputeProjectionInterpolation",
		"not implemented. Requires explicit matrix inverse");
}

/* Project point values onto mesh, write into displacement field.  Used to compute initial
   displacements from point values to mesh. */
void BridgingScaleT::InitialProject(const StringT& field, const PointInCellDataT& cell_data,
	const dArray2DT& point_values, dArray2DT& projection, dArray2DT& projectedu)
{
#pragma unused(field)
	/* projected part of the mesh */
	const iArrayT& cell_nodes = cell_data.CellNodes();
	const iArray2DT& cell_connects = cell_data.CellConnectivities();

	/* points in cell data */
	const RaggedArray2DT<int>& point_in_cell = cell_data.PointInCell();
	const dArray2DT& weights = cell_data.InterpolationWeights();
	const InverseMapT& global_to_local = cell_data.GlobalToLocal();
	
	/* initialize return value */
	projection.Dimension(cell_nodes.Length(), point_values.MinorDim());
	projection = 0.0;
	
	/* loop over mesh */
	int cell_dex = 0;
	iArrayT cell_eq;
	dArrayT Na, point_value;
	dMatrixT Nd(cell_connects.MinorDim(), point_values.MinorDim());
	//double atommass;
	for (int i = 0; i < point_in_cell.MajorDim(); i++)
	{
		int np = point_in_cell.MinorDim(i);
		if (np > 0)
		{
			const int* points = point_in_cell(i);
			Nd = 0.0;
			for (int j = 0; j < np; j++)
			{
				int point = points[j];
		
				/* fetch interpolation weights */
				int point_dex = global_to_local.Map(point);
				weights.RowAlias(point_dex, Na);
				//atommass = mdmass[point_dex];
				
				/* source values of the point */
				point_values.RowAlias(point, point_value);

				/* rhs during projection - calculating part of w */
				Nd.Outer(Na, point_value, 1.0, dMatrixT::kAccumulate);
				//Nd *= atommass;	// need to multiply by atomic mass, i.e. w = M^{-1}N^{T}M_{A}q
			}

			/* equations of cell in projector */
			cell_connects.RowAlias(cell_dex++, cell_eq);

			/* assemble */
			for (int j = 0; j < Nd.Cols(); j++)
				projection.Accumulate(j, cell_eq, Nd(j));
		}
	}

//TEMP - write mass matrix to file
#if 0
ostream& out = ElementSupport().Output();
out << "\n values =\n" << values << endl;
out << "\n residual =\n" << projection << endl;
#endif
//TEMP

	/* calculate projection - requires global matrix that supports 
	 * multiple solves - projection = w after operations within this loop */
	dArrayT u_tmp(projection.MajorDim());
	for (int i = 0; i < projection.MinorDim(); i++)
	{
		projection.ColumnCopy(i, u_tmp);
		fGlobalMass->Solve(u_tmp);
		projection.SetColumn(i, u_tmp);
	}	

#if 0
	ofstream project, fenodes1, fenodes2;
	project.open("project.dat");
	fenodes1.open("fenodes1.dat");
	fenodes2.open("fenodes2.dat");
	project.precision(13);
	for (int i = 0; i < projection.MajorDim(); i++)
	{
		project << i+1 << " " << 1 << " " << 1 << " " << projection(i,0) << endl;
		project << i+1 << " " << 2 << " " << 1 << " " << projection(i,1) << endl;
		fenodes1 << "*set" << endl;
		fenodes1 << 1 << endl;
		fenodes1 << cell_nodes[i]+1 << endl;
		fenodes2 << i+1 << " " << 1 << endl;
	}
	
	project.close();
	fenodes1.close();
	fenodes2.close();
#endif
	
	u_tmp.Free();

	fFineScale.Dimension(point_values);
	fFineScale = 0.0;
	fCoarseScale.Dimension(point_values);
	fCoarseScale = 0.0;
	projectedu.Dimension(point_values);
	projectedu = 0.0;

	cell_dex = 0;
	iArrayT cell_connect;
	dArray2DT cell_projection(cell_connects.MinorDim(), projection.MinorDim());
	
	for (int i = 0; i < point_in_cell.MajorDim(); i++)
	{
		int np = point_in_cell.MinorDim(i);
		if (np > 0)
		{
			/* gather cell information */
			cell_connects.RowAlias(cell_dex++, cell_connect);
			cell_projection.RowCollect(cell_connect, projection);
		
			const int* points = point_in_cell(i);
			for (int j = 0; j < np; j++)
			{
				int point = points[j];
		
				/* fetch interpolation weights */
				int point_dex = global_to_local.Map(point);
				weights.RowAlias(point_dex, Na);

				/* interpolate to point */
				for (int k = 0; k < projection.MinorDim(); k++)
				{
					/* interpolate coarse scale to point (gives Nw) */
					fCoarseScale(point, k) = cell_projection.DotColumn(k, Na);

					/* error = source - projection = q-Nw*/
					fFineScale(point, k) = point_values(point, k) - fCoarseScale(point, k);
					
					/* calculate totalu as function of fine scale and projected u (Nw) */
					projectedu(point,k) = fCoarseScale(point,k) + fFineScale(point,k);
				}
			}
		}
	}
}

/* calculate the fine scale part of MD solution as well as total solution u - same as project Field
 * except for those changes */
void BridgingScaleT::BridgingFields(const StringT& field, const PointInCellDataT& cell_data,
	const dArray2DT& mddisp, const dArray2DT& fedisp, dArray2DT& projection, dArray2DT& totalu)
{
#pragma unused(field)

	/* projected part of the mesh */
	const iArrayT& cell_nodes = cell_data.CellNodes();
	const iArray2DT& cell_connects = cell_data.CellConnectivities();
	
	/* points in cell data */
	const RaggedArray2DT<int>& point_in_cell = cell_data.PointInCell();
	const dArray2DT& weights = cell_data.InterpolationWeights();
	const InverseMapT& global_to_local = cell_data.GlobalToLocal();
	
	/* initialize return value */
	projection.Dimension(cell_nodes.Length(), mddisp.MinorDim());
	projection = 0.0;
	
	/* loop over mesh */
	int cell_dex = 0;
	iArrayT cell_eq;
	dArrayT Na, point_value;
	dMatrixT Nd(cell_connects.MinorDim(), mddisp.MinorDim());
	//double atommass;
	for (int i = 0; i < point_in_cell.MajorDim(); i++)
	{
		int np = point_in_cell.MinorDim(i);
		if (np > 0)
		{
			const int* points = point_in_cell(i);
			Nd = 0.0;
			for (int j = 0; j < np; j++)
			{
				int point = points[j];
			
				/* fetch interpolation weights */
				int point_dex = global_to_local.Map(point);
				weights.RowAlias(point_dex, Na);
				//atommass = mdmass[point_dex];

				/* source values of the point */
				mddisp.RowAlias(point, point_value);
			
				/* rhs during projection - calculating part of w */
				Nd.Outer(Na, point_value, 1.0, dMatrixT::kAccumulate);
				//Nd *= atommass;	// need to multiply by atomic mass, i.e. w = M^{-1}N^{T}M_{A}q
			}

			/* equations of cell in projector */
			cell_connects.RowAlias(cell_dex++, cell_eq);

			/* assemble */
			for (int j = 0; j < Nd.Cols(); j++)
				projection.Accumulate(j, cell_eq, Nd(j));
		}
	}
	
//TEMP - write mass matrix to file
#if 0
ostream& out = ElementSupport().Output();
out << "\n values =\n" << values << endl;
out << "\n residual =\n" << projection << endl;
#endif
//TEMP

	/* calculate projection - requires global matrix that supports 
	 * multiple solves - projection = w after operations within this loop */
	dArrayT u_tmp(projection.MajorDim());
	for (int i = 0; i < projection.MinorDim(); i++)
	{
		projection.ColumnCopy(i, u_tmp);
		fGlobalMass->Solve(u_tmp);
		projection.SetColumn(i, u_tmp);
	}
	u_tmp.Free();

	/* initialize return values */
	fFineScale.Dimension(mddisp);
	fFineScale = 0.0;
	totalu.Dimension(mddisp);
	totalu = 0.0;
	fCoarseScale.Dimension(mddisp);
	fCoarseScale = 0.0;
	dArray2DT coarse_scale;
	dArray2DT coarse(cell_connects.MinorDim(), projection.MinorDim());
	coarse_scale.Dimension(mddisp);
	coarse_scale = 0.0;

	cell_dex = 0;
	iArrayT cell_connect;
	dArray2DT cell_projection(cell_connects.MinorDim(), projection.MinorDim());

	/* element group information */
	const ContinuumElementT* continuum = cell_data.ContinuumElement();
	const iArrayT& cell = cell_data.InterpolatingCell();

	for (int i = 0; i < point_in_cell.MajorDim(); i++)
	{
		int np = point_in_cell.MinorDim(i);
		if (np > 0)
		{
			/* gather cell information */
			cell_connects.RowAlias(cell_dex++, cell_connect);
			cell_projection.RowCollect(cell_connect, projection);
				
			/* element info */
			const int* points = point_in_cell(i);
			int off = global_to_local.Map(points[0]);
			const ElementCardT& element_card = continuum->ElementCard(cell[off]);
			const iArrayT& fenodes = element_card.NodesU();
			coarse.RowCollect(fenodes, fedisp);  
		
			for (int j = 0; j < np; j++)
			{
				int point = points[j];
				
				/* fetch interpolation weights */
				int point_dex = global_to_local.Map(point);
				weights.RowAlias(point_dex, Na);

				/* interpolate to point */
				for (int k = 0; k < projection.MinorDim(); k++)
				{
					/* interpolate projected coarse scale to point (gives Nw) */
					fCoarseScale(point, k) = cell_projection.DotColumn(k, Na);
					
					/* interpolate actual FEM solution to atoms (point) */
					coarse_scale(point,k) = coarse.DotColumn(k,Na);
					
					/* error = source - projection = q-Nw*/
					fFineScale(point, k) = mddisp(point, k) - fCoarseScale(point, k);
					
					/* compute total displacement u = FEM + fine scale */
					totalu(point,k) = coarse_scale(point,k) + fFineScale(point,k);
				}
			}
		}
	}
}

/* writing output */
void BridgingScaleT::RegisterOutput(void)
{
	/* collect variable labels */
	int ndof = NumDOF();
	if (ndof > 3) throw;
	ArrayT<StringT> n_labels(2*ndof);
	const char* coarse_labels[] = {"FE_X", "FE_Y", "FE_Z"};
	const char* fine_labels[] = {"fine_X", "fine_Y", "fine_Z"};
	int dex = 0;
	for (int i = 0; i < ndof; i++) n_labels[dex++] = coarse_labels[i];
	for (int i = 0; i < ndof; i++) n_labels[dex++] = fine_labels[i];

	/* register output at solid nodes */
#if 0
	OutputSetT output_set_solid(GeometryT::kPoint, fSolidNodesUsed, n_labels);
	fSolidOutputID = ElementSupport().RegisterOutput(output_set_solid);
#endif

	/* register output at particles */
//	OutputSetT output_set_particle(GeometryT::kPoint, fParticlesUsed, n_labels);
//	fParticleOutputID = ElementSupport().RegisterOutput(output_set_particle);
}

//NOTE - this function is/was identical to CSEBaseT::WriteOutput
void BridgingScaleT::WriteOutput(void)
{
	/* calculate output values */
	dArray2DT n_values; // dimension: [number of output nodes] x [number of values]
//	ComputeOutput(n_counts, n_values);
//  rows in n_values correspond to the nodes listed in fSolidNodesUsed

	/* send to output */
//	ElementSupport().WriteOutput(fOutputID, n_values);
}

/* describe the parameters needed by the interface */
void BridgingScaleT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

	/* coarse scale element group */
	list.AddParameter(ParameterT::Integer, "solid_element_group");
}
	
/* accept parameter list */
void BridgingScaleT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "BridgingScaleT::TakeParameterList";

	/* resolve the solid element group */
	int group = list.GetParameter("solid_element_group");
	group--;
	ElementBaseT& element = ElementSupport().ElementGroup(group);
#ifdef __NO_RTTI__
	fSolid = element.dynamic_cast_SolidElementT();
#else
	fSolid = dynamic_cast<const SolidElementT*>(&element);
#endif
	if (!fSolid) ExceptionT::GeneralFail(caller, "could not resolve element group %d", group+1);

	/* solid element class needs to store the internal force vector */
	SolidElementT* solid = const_cast<SolidElementT*>(fSolid);
	solid->SetStoreInternalForce(true);

	/* inherited */
	ElementBaseT::TakeParameterList(list);

	/* dimension workspace */
	fElMatU.Dimension(ShapeFunction().ParentDomain().NumNodes());
	fDOFvec.Dimension(NumDOF());
	fConnect.Dimension(fSolid->NumElements(), fSolid->NumElementNodes());
	fWtempU.Dimension(ShapeFunction().ParentDomain().NumNodes(), NumDOF());
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* initialize local arrays */
void BridgingScaleT::SetLocalArrays(void)
{
	/* dimension */
	fLocInitCoords.Dimension(NumElementNodes(), NumSD());
	fLocDisp.Dimension(NumElementNodes(), NumDOF());

	/* set source */
	ElementSupport().RegisterCoordinates(fLocInitCoords);
	Field().RegisterLocal(fLocDisp);	
}

/* write all current element information to the stream */
void BridgingScaleT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	ElementBaseT::CurrElementInfo(out);
	dArray2DT temp;
	temp.Dimension(fLocInitCoords.NumberOfNodes(), fLocInitCoords.MinorDim());
	
	out <<   " initial coords:\n";
	temp.Dimension(fLocInitCoords.NumberOfNodes(), fLocInitCoords.MinorDim());
	fLocInitCoords.ReturnTranspose(temp);
	temp.WriteNumbered(out);

	out <<   " displacements:\n";
	temp.Dimension(fLocDisp.NumberOfNodes(), fLocDisp.MinorDim());
	fLocDisp.ReturnTranspose(temp);
	temp.WriteNumbered(out);
}
