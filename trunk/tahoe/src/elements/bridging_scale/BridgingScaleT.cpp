/* $Id: BridgingScaleT.cpp,v 1.19 2002-08-15 23:41:26 paklein Exp $ */
#include "BridgingScaleT.h"

#include <iostream.h>
#include <iomanip.h>

#include "ShapeFunctionT.h"
#include "RodT.h"
#include "fstreamT.h"
#include "iAutoArrayT.h"
#include "OutputSetT.h"
#include "AutoFill2DT.h"
#include "RaggedArray2DT.h"
#include "iGridManagerT.h"
#include "iNodeT.h"
#include "nArrayGroupT.h"

using namespace Tahoe;

/* constructor */
BridgingScaleT::BridgingScaleT(const ElementSupportT& support, 
	const FieldT& field,
	const RodT& particle,
	const ElasticT& solid):
	ElementBaseT(support, field),
	fParticle(particle),
	fSolid(solid),
	fElMat(ShapeFunction().ParentDomain().NumNodes(), ElementMatrixT::kSymmetric),
	fLocInitCoords(LocalArrayT::kInitCoords),
	fLocDisp(LocalArrayT::kDisp),
	fDOFvec(NumDOF()),
	fGlobal(support.Output(),1),
	fMass(ShapeFunction().ParentDomain().NumNodes()),
	fConnect(solid.NumElements(), solid.NumElementNodes()),
	//fWtemp(support.NumNodes(),ShapeFunction().ParentDomain().NumNodes()),
	//fW(support.NumNodes(),ShapeFunction().ParentDomain().NumNodes())
	fWtempU(ShapeFunction().ParentDomain().NumNodes(), NumDOF()),
	fWtempV(ShapeFunction().ParentDomain().NumNodes(), NumDOF()),
	fWtempA(ShapeFunction().ParentDomain().NumNodes(), NumDOF()),
	//fWU(ShapeFunction().ParentDomain().NumNodes()),
	fWV(ShapeFunction().ParentDomain().NumNodes()),
	fWA(ShapeFunction().ParentDomain().NumNodes()),
	//fError(ShapeFunction().ParentDomain().NumNodes()),
  	//fFineScaleU(ShapeFunction().ParentDomain().NumNodes()),
	//fCoarseScaleU(ShapeFunction().ParentDomain().NumNodes()),
	//fTotalU(ShapeFunction().ParentDomain().NumNodes()),
	fMassInv(ShapeFunction().ParentDomain().NumNodes())
	//	fProjection(support.NumNodes())
{

}

/* destructor */
BridgingScaleT::~BridgingScaleT(void)
{	

}

/* allocates space and reads connectivity data */
void BridgingScaleT::Initialize(void)
{
	/* inherited */
	ElementBaseT::Initialize();

	/* stream */
	ostream& out = ElementSupport().Output();
	
	/* current coordinates of all nodes/points in the calculation */
	const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();

	/* distinguish FEM nodes vs. atoms to separate their respective coordinates */
	iArrayT atoms_used, nodes_used;
	fParticle.NodesUsed(atoms_used);
	fSolid.NodesUsed(nodes_used);
	iArrayT node1(nodes_used.Max() + 1); // +1 because starts from 0
	node1 = -1;
	for (int i = 0; i < nodes_used.Length(); i++)
	    node1[nodes_used[i]] = i + 1;

#if 0
	dArray2DT atom_coords, node_coords;
	atom_coords.RowCollect(atoms_used, curr_coords);
	node_coords.RowCollect(nodes_used, curr_coords);
#endif
	atoms_used++;
	out << " Particles used:\n" << atoms_used.wrap(5) << '\n';
	atoms_used--;
	nodes_used++;
	out << " Nodes used:\n" << nodes_used.wrap(5) << '\n';
	nodes_used--;

	/* now take an atom, check against every element */
	const ParentDomainT& parent = ShapeFunction().ParentDomain();
	int totalatoms = atoms_used.Length();
	fTotalNodes = nodes_used.Length();
	AutoFill2DT<int> auto_fill(fSolid.NumElements(), 10, 10);
	dArrayT x_atom, centroid;
	LocalArrayT cell_coords(LocalArrayT::kCurrCoords, fSolid.NumElementNodes(), NumSD());
	cell_coords.SetGlobal(curr_coords); // Sets address of cell_coords

	iGridManagerT grid(10, 100, curr_coords, &atoms_used);
	grid.Reset();
	for (int i = 0; i < fSolid.NumElements(); i++) {
	
			/* gives domain (global) nodal coordinates */
	                cell_coords.SetLocal(fSolid.ElementCard(i).NodesX());

			/* get global connectivity matrix */
		        fConnect.SetRow(i,fSolid.ElementCard(i).NodesX());

			/* inverse map from global node numbers to local node numbers */
			for (int j = 0; j < fSolid.NumElementNodes(); j++)
			{
			    int x = fConnect(i,j);
			    fConnect(i,j) = node1[x];
			}

			/* centroid and radius */
			double radius = parent.AverageRadius(cell_coords, centroid);

			/* candidate particles */
			const AutoArrayT<iNodeT>& hits = grid.HitsInRegion(centroid.Pointer(), 1.01*radius);

			for (int j = 0; j < hits.Length(); j++)
			{
				x_atom.Set(NumSD(), hits[j].Coords());
				if (parent.PointInDomain(cell_coords, x_atom)) 
				        auto_fill.Append(i, hits[j].Tag());
			}
	}

	/* copy/compress contents */
	fParticlesInCell.Copy(auto_fill);
	iArrayT tmp(fParticlesInCell.Length(), fParticlesInCell.Pointer());
	out << " Particles in cells:\n"
	    << setw(kIntWidth) << "no." << '\n';
	tmp++;
	fParticlesInCell.WriteNumbered(out);
	tmp--;
	out.flush();

	if (fParticlesInCell.Length() > atoms_used.Length()) {
		cout << "\n BridgingScaleT::Initialize: WARNING: number of particles in cells " 
		     << fParticlesInCell.Length() << " exceeds\n" 
		     <<   "     the total number of particles " << atoms_used.Length() << endl;
	}

	// (2) compute the inverse map using list fParticlesInCell
	LocalArrayT cell_coord(LocalArrayT::kCurrCoords, fSolid.NumElementNodes(), NumSD());
	cell_coord.SetGlobal(curr_coords); // Sets address of cell_coords
	iArrayT atom_nums;
	dArrayT mapped(NumSD()), point(NumSD());
	dArray2DT atom_coords(fParticlesInCell.MaxMinorDim(), NumSD());
	AutoFill2DT<double> inverse(fParticlesInCell.MajorDim(),10,10);
	for (int i = 0; i < fParticlesInCell.MajorDim(); i++) {

	                fParticlesInCell.RowAlias(i,atom_nums);
	                atom_coords.RowCollect(atom_nums,curr_coords);
			/* gives domain (global) nodal coordinates */
	                cell_coord.SetLocal(fSolid.ElementCard(i).NodesX()); 
			
	                for (int j = 0; j < fParticlesInCell.MinorDim(i); j++) 
			{
			    //need a ColumnAlias function in nArray2DT similar to RowAlias!
			    atom_coords.RowAlias(j, point);
			    if (parent.MapToParentDomain(cell_coord,point,mapped))
				inverse.Append(i,mapped);   
			}
	}
        fInverseMapInCell.Copy(inverse);
	dArrayT temp(fInverseMapInCell.Length(), fInverseMapInCell.Pointer());
	out << " Inverse maps in cell:\n"
	    << setw(kDoubleWidth) << "no." << '\n';
	fInverseMapInCell.WriteNumbered(out);
	out.flush();	

	/* store solid nodes used */
	fSolidNodesUsed.Dimension(fTotalNodes, 1);
	fSolidNodesUsed.SetColumn(0, nodes_used);

	/* store particles used */
	fParticlesUsed.Dimension(atoms_used.Length(), 1);
	fParticlesUsed.SetColumn(0, atoms_used);
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
	ComputeCoarseScaleFields();
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

	/* group ID */	
	StringT set_ID;
	set_ID.Append(ElementSupport().ElementGroupNumber(this) + 1);

	/* register output at solid nodes */
	OutputSetT output_set_solid(set_ID, GeometryT::kPoint, fSolidNodesUsed, n_labels);
	fSolidOutputID = ElementSupport().RegisterOutput(output_set_solid);

	/* register output at particles */
	OutputSetT output_set_particle(set_ID, GeometryT::kPoint, fParticlesUsed, n_labels);
	fParticleOutputID = ElementSupport().RegisterOutput(output_set_particle);
}

//NOTE - this function is/was identical to CSEBaseT::WriteOutput
void BridgingScaleT::WriteOutput(IOBaseT::OutputModeT mode)
{
//TEMP - not handling general output modes yet
	if (mode != IOBaseT::kAtInc)
	{
		cout << "\n BridgingScaleT::WriteOutput: only handling \"at increment\"\n"
		     <<   "     print mode. SKIPPING." << endl;
		return;
	}

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
	fLocInitCoords.Allocate(NumElementNodes(), NumSD());
	fLocDisp.Allocate(NumElementNodes(), NumDOF());

	/* set source */
	ElementSupport().RegisterCoordinates(fLocInitCoords);
	Field().RegisterLocal(fLocDisp);	
}

/* print element group data */
void BridgingScaleT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

	out << " Particle group number . . . . . . . . . . . . . = " << ElementSupport().ElementGroupNumber(&fParticle) + 1 << '\n';
	out << " Continuum group number. . . . . . . . . . . . . = " << ElementSupport().ElementGroupNumber(&fSolid) + 1 << '\n';
}

/* write all current element information to the stream */
void BridgingScaleT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	ElementBaseT::CurrElementInfo(out);
	dArray2DT temp;
	temp.Allocate(fLocInitCoords.NumberOfNodes(), fLocInitCoords.MinorDim());
	
	out <<   " initial coords:\n";
	temp.Allocate(fLocInitCoords.NumberOfNodes(), fLocInitCoords.MinorDim());
	fLocInitCoords.ReturnTranspose(temp);
	temp.WriteNumbered(out);

	out <<   " displacements:\n";
	temp.Allocate(fLocDisp.NumberOfNodes(), fLocDisp.MinorDim());
	fLocDisp.ReturnTranspose(temp);
	temp.WriteNumbered(out);
}

/***********************************************************************
* Private
***********************************************************************/

void BridgingScaleT::ComputeCoarseScaleFields(void)
{
  /* computes coarse scale fields */
  const ParentDomainT& parent = ShapeFunction().ParentDomain();
  const FieldT& field = Field();
  iArrayT atoms, elemconnect(parent.NumNodes());
  const dArray2DT& displacements = field[0];
  dArrayT map, shape(parent.NumNodes()), temp(NumSD()), disp, vel, acc;
  dMatrixT tempmass(parent.NumNodes()); 
  dMatrixT Nd(parent.NumNodes(), NumSD()), Nv(parent.NumNodes(), NumSD()), Na(parent.NumNodes(), NumSD());
  /* set global matrices for bridging scale calculations */
  fGlobal.AddEquationSet(fConnect);
  fGlobal.Initialize(fTotalNodes,fTotalNodes,1); // hard wired for serial calculations
  fGlobal.Clear();

  for (int i = 0; i < fParticlesInCell.MajorDim(); i++)
  {
      fElMat = 0.0, fWtempU = 0.0, fWtempV = 0.0, fWtempA = 0.0;
      fParticlesInCell.RowAlias(i,atoms);
      fInverseMapInCell.RowAlias(i,map);
      fConnect.RowAlias(i,elemconnect);
      for (int j = 0; j < fParticlesInCell.MinorDim(i); j++)
      {
	  // still need to access individual atomic masses!!!
	  displacements.RowAlias(atoms[j],disp);
//	  velocities.RowAlias(atoms[j],vel);
//	  accelerations.RowAlias(atoms[j],acc);
	  temp.CopyPart(0, map, NumSD()*j, NumSD());
	  parent.EvaluateShapeFunctions(temp,shape);
	  tempmass.Outer(shape,shape);
	  fElMat += tempmass;
	  Nd.Outer(shape, disp);
//	  Nv.Outer(shape, vel);
//	  Na.Outer(shape, acc);
	  fWtempU += Nd;
//	  fWtempV += Nv;
//	  fWtempA += Na;
      }
      
      /* Assemble local equations to global equations for LHS */
      fGlobal.Assemble(fElMat, elemconnect);

      //ComputeU(displacements, velocities, accelerations);
  }
  cout << "Global = \n" << fGlobal << endl;
}

void BridgingScaleT::ComputeFineScaleFields(void)
{
  /* computes fine scale field solutions after coarse scale fields are known */
  const ParentDomainT& parent = ShapeFunction().ParentDomain();
  const FieldT& field = Field();
  const dArray2DT& displacements = field[0];
  const dArray2DT& velocities = field[1];
  const dArray2DT& accelerations = field[2];
  dArrayT map, shape1(parent.NumNodes()), shape2(parent.NumNodes());
  dArrayT map1(NumSD()), map2(NumSD());
  iArrayT disp1;
  dArray2DT disp(fParticlesInCell.MaxMinorDim(),NumDOF()), vel(fParticlesInCell.MaxMinorDim(),NumDOF());
  dArray2DT acc(fParticlesInCell.MaxMinorDim(),NumDOF());
  double mult;
  
	/* arrays with changing dimensions */
  	dArrayT ErrorU, FineScaleU, CoarseScaleU, TotalU, ErrorV, FineScaleV;
  	dArrayT CoarseScaleV, TotalV, ErrorA, FineScaleA, CoarseScaleA, TotalA;
	nArrayGroupT<double> dArrayT_group(25);
  	dArrayT_group.Register(ErrorU);
  	dArrayT_group.Register(FineScaleU);
  	dArrayT_group.Register(CoarseScaleU);
  	dArrayT_group.Register(TotalU);
  	dArrayT_group.Register(ErrorV);
  	dArrayT_group.Register(FineScaleV);
  	dArrayT_group.Register(CoarseScaleV);
  	dArrayT_group.Register(TotalV);
  	dArrayT_group.Register(ErrorA);
  	dArrayT_group.Register(FineScaleA);
  	dArrayT_group.Register(CoarseScaleA);
  	dArrayT_group.Register(TotalA);
  
  for (int i = 0; i < fParticlesInCell.MajorDim(); i++)
  {
      /* resize array group */
      dArrayT_group.Dimension(fParticlesInCell.MinorDim(i), false);
  
      fInverseMapInCell.RowAlias(i,map);
      fParticlesInCell.RowAlias(i,disp1);
      disp.RowCollect(disp1,displacements);
      vel.RowCollect(disp1,velocities);
      acc.RowCollect(disp1,accelerations);
      dMatrixT Projection(fParticlesInCell.MinorDim(i));
      ErrorU = 0.0, FineScaleU = 0.0, CoarseScaleU = 0.0, TotalU = 0.0, Projection = 0.0;
      ErrorA = 0.0, FineScaleA = 0.0, CoarseScaleA = 0.0, TotalA = 0.0;
      ErrorV = 0.0, FineScaleV = 0.0, CoarseScaleV = 0.0, TotalV = 0.0;

      for (int j = 0; j < fParticlesInCell.MinorDim(i); j++)
      {
	map1 = map[j];
	parent.EvaluateShapeFunctions(map1,shape1);
	ErrorU[j] = dArrayT::Dot(shape1,fWU);
	ErrorV[j] = dArrayT::Dot(shape1,fWV);
	ErrorA[j] = dArrayT::Dot(shape1,fWA);
	FineScaleU[j] = disp(j,0) - ErrorU[j];
	FineScaleV[j] = vel(j,0) - ErrorV[j];
	FineScaleA[j] = acc(j,0) - ErrorA[j];
	for (int k = 0; k < fParticlesInCell.MinorDim(i); k++)
	{
	  // still need to access individual atomic masses
	  map2 = map[k];
	  parent.EvaluateShapeFunctions(map2,shape2);
	  mult = fMassInv.MultmBn(shape1,shape2);
	  Projection(j,k) += mult;
	}
      }
      Projection.Multx(disp,CoarseScaleU);
      Projection.Multx(vel,CoarseScaleV);
      Projection.Multx(acc,CoarseScaleA);
      TotalU += CoarseScaleU;
      TotalU += FineScaleU;
      TotalV += CoarseScaleV;
      TotalV += FineScaleV;
      TotalA += CoarseScaleA;
      TotalA += FineScaleA;
  }
  cout << "FineScaleU = \n" << disp << endl;
  cout << "FineScaleV = \n" << TotalU << endl;
  cout << "FineScaleA = \n" << CoarseScaleU << endl;
}
