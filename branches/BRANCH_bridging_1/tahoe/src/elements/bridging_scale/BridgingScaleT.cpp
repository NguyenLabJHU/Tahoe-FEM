/* $Id: BridgingScaleT.cpp,v 1.30.2.2 2003-02-10 09:25:37 paklein Exp $ */
#include "BridgingScaleT.h"

#include <iostream.h>
#include <iomanip.h>

#include "ShapeFunctionT.h"
//#include "RodT.h"
#include "fstreamT.h"
#include "iAutoArrayT.h"
#include "OutputSetT.h"
#include "AutoFill2DT.h"
#include "RaggedArray2DT.h"
#include "iGridManagerT.h"
#include "iNodeT.h"
#include "nArrayGroupT.h"
#include "PointInCellDataT.h"

using namespace Tahoe;

/* constructor */
BridgingScaleT::BridgingScaleT(const ElementSupportT& support, 
	const FieldT& field,
//	const RodT& particle,
	const SolidElementT& solid):
	ElementBaseT(support, field),
//	fParticle(particle),
	fSolid(solid),
	fElMatU(ShapeFunction().ParentDomain().NumNodes(), ElementMatrixT::kSymmetric),
	fLocInitCoords(LocalArrayT::kInitCoords),
	fLocDisp(LocalArrayT::kDisp),
	fDOFvec(NumDOF()),
	fGlobalMass(support.Output(),1),
	fConnect(solid.NumElements(), solid.NumElementNodes()),
	fWtempU(ShapeFunction().ParentDomain().NumNodes(), NumDOF())
{

}

/* map coordinates into elements */
void BridgingScaleT::MaptoCells(const iArrayT& points_used, const dArray2DT* init_coords, 
	const dArray2DT* curr_coords, PointInCellDataT& cell_data)
{
	const char caller[] = "BridgingScaleT::MaptoCells";

	/* map data */
	cell_data.SetContinuumElement(fSolid);
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

	/* check all cells for points */
	const ParentDomainT& parent = ShapeFunction().ParentDomain();
	AutoFill2DT<int> auto_fill(fSolid.NumElements(), 1, 10, 10);
	dArrayT x_atom, centroid;
	LocalArrayT loc_cell_coords(coord_type, fSolid.NumElementNodes(), NumSD());
	loc_cell_coords.SetGlobal(cell_coordinates);
	for (int i = 0; i < fSolid.NumElements(); i++) {
	
		/* gives domain (global) nodal coordinates */
		loc_cell_coords.SetLocal(fSolid.ElementCard(i).NodesX());

		/* centroid and radius */
		double radius = parent.AverageRadius(loc_cell_coords, centroid);

		/* candidate points */
		const AutoArrayT<iNodeT>& hits = grid.HitsInRegion(centroid.Pointer(), 1.01*radius);

		/* check if points are within the element domain */
		for (int j = 0; j < hits.Length(); j++)
		{
			x_atom.Set(NumSD(), hits[j].Coords());
			if (parent.PointInDomain(loc_cell_coords, x_atom)) 
				auto_fill.Append(i, hits[j].Tag());
		}
	}

	/* copy/compress contents */
	point_in_cell.Copy(auto_fill);
	auto_fill.Free();
	
	/* verbose output */
	if (ElementSupport().PrintInput()) {
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
		loc_cell_coords.SetLocal(fSolid.ElementCard(i).NodesX()); 

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
	if (ElementSupport().PrintInput()) {

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

void BridgingScaleT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* inherited */
	ElementBaseT::Equations(eq_1, eq_2);
}

/* initialize/finalize step */
void BridgingScaleT::InitStep(void)
{
	/* inherited */
	ElementBaseT::InitStep();
}

/* initialize/finalize step */
void BridgingScaleT::CloseStep(void)
{
	/* inherited */
	ElementBaseT::CloseStep();
	
	/* Do Bridging scale calculations after MD displacements
	   computed using RodT */
	CoarseFineFields();
}

/* resets to the last converged solution */
void BridgingScaleT::ResetStep(void)
{
	/* inherited */
	ElementBaseT::ResetStep();
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

#if 0
	/* register output at particles */
	OutputSetT output_set_particle(GeometryT::kPoint, fParticlesUsed, n_labels);
	fParticleOutputID = ElementSupport().RegisterOutput(output_set_particle);
#endif
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

/* print element group data */
void BridgingScaleT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

//DEV
//	out << " Particle group number . . . . . . . . . . . . . = " << ElementSupport().ElementGroupNumber(&fParticle) + 1 << '\n';
	out << " Continuum group number. . . . . . . . . . . . . = " << ElementSupport().ElementGroupNumber(&fSolid) + 1 << '\n';
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

/***********************************************************************
* Private
***********************************************************************/

void BridgingScaleT::CoarseFineFields(void)
{
  /* computes coarse and fine scale displacement fields */
  const ParentDomainT& parent = ShapeFunction().ParentDomain();
  const FieldT& field = Field();
  iArrayT atoms, elemconnect(parent.NumNodes());
  const dArray2DT& displacements = field[0];

//  fFineScaleU.Dimension(fParticlesUsed.MajorDim(), displacements.MinorDim());

  dArrayT map, shape(parent.NumNodes()), temp(NumSD()), disp;
  dMatrixT tempmass(parent.NumNodes()), Nd(parent.NumNodes(), NumSD()), wglobalu(fTotalNodes, NumSD());
  /* set global matrices for bridging scale calculations */
  fGlobalMass.AddEquationSet(fConnect);
  fGlobalMass.Initialize(fTotalNodes,fTotalNodes,1); // hard wired for serial calculations
  fGlobalMass.Clear();
  wglobalu = 0.0;

  for (int i = 0; i < fParticlesInCell.MajorDim(); i++)
  {
      fElMatU = 0.0, fWtempU = 0.0;
      fParticlesInCell.RowAlias(i,atoms);
      fInverseMapInCell.RowAlias(i,map);
      fConnect.RowAlias(i,elemconnect);
      for (int j = 0; j < fParticlesInCell.MinorDim(i); j++)
      {
	  // still need to access individual atomic masses!!!
	  displacements.RowAlias(atoms[j], disp);
	  temp.CopyPart(0, map, NumSD()*j, NumSD());
	  parent.EvaluateShapeFunctions(temp,shape);
	  tempmass.Outer(shape,shape);
	  fElMatU += tempmass;
	  Nd.Outer(shape, disp);
	  fWtempU += Nd;
      }
      
      /* Assemble local equations to global equations for LHS */
      fGlobalMass.Assemble(fElMatU, elemconnect);

      /* Assemble local equations to global equations for RHS */
      for (int i = 0; i < parent.NumNodes(); i++)
      {
	for (int j = 0; j < NumSD(); j++)
	    wglobalu(elemconnect[i] - 1,j) += fWtempU(i,j);
      }
  }
  
  /* Solve for global coarse scale (FEM) displacements */
  if (NumSD() == 1)
    {
      wglobalu.ColumnAlias(0, fUx);
      fGlobalMass.Solve(fUx);
      cout << "Ux = \n" << fUx << endl;
    }
  else if (NumSD() == 2)
    {
      wglobalu.ColumnAlias(0, fUx);
      wglobalu.ColumnAlias(1, fUy);
      fGlobalMass.Solve(fUx);
      fGlobalMass.Solve(fUy);
      cout << "Ux = \n" << fUx << endl;
      cout << "Uy = \n" << fUy << endl;
    }

  /* Compute resulting fine scale fields = MD - FE */
  dArrayT ux(parent.NumNodes()), uy(parent.NumNodes());
  iArrayT atomconn, conn;
  double nux, nuy;

  for (int i = 0; i < fParticlesInCell.MajorDim(); i++)
  {
    fParticlesInCell.RowAlias(i,atoms);
    fInverseMapInCell.RowAlias(i,map);
    fConnect.RowCopy(i,conn);
    fAtomConnect.RowCopy(i,atomconn);
    atomconn -= 1;
    conn -= 1; // so that array indices start from 0
    ux.Collect(conn, fUx);
    for (int j = 0; j < fParticlesInCell.MinorDim(i); j++)
      {
	displacements.RowAlias(atoms[j],disp);
	temp.CopyPart(0, map, NumSD()*j, NumSD());
	parent.EvaluateShapeFunctions(temp,shape);
	nux = dArrayT::Dot(shape,ux);
	int which = atomconn[j];
	double fux = disp[0] - nux;
	fFineScaleU(which,0) = fux;
	if (NumSD() == 2)
	  {
	    uy.Collect(conn, fUy);
	    nuy = dArrayT::Dot(shape,uy);
	    double fuy = disp[1] - nuy;
	    fFineScaleU(which,1) = fuy; 
	  }
      }
  }
  cout << "MD U = \n" << displacements << endl;
  cout << "fine scale U = \n" << fFineScaleU << endl; 
}
