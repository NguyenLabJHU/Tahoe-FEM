/* $Id: BridgingScaleT.cpp,v 1.12 2002-08-08 23:24:11 hspark Exp $ */
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

using namespace Tahoe;

/* constructor */
BridgingScaleT::BridgingScaleT(const ElementSupportT& support, 
	const FieldT& field,
	const RodT& particle,
	const ElasticT& solid):
	ElementBaseT(support, field),
	fParticle(particle),
	fSolid(solid),
	fLocInitCoords(LocalArrayT::kInitCoords),
	fLocDisp(LocalArrayT::kDisp),
	fDOFvec(NumDOF()),
	fMass(ShapeFunction().ParentDomain().NumNodes()),
	//fWtemp(support.NumNodes(),ShapeFunction().ParentDomain().NumNodes()),
	//fW(support.NumNodes(),ShapeFunction().ParentDomain().NumNodes())
	fWtempU(ShapeFunction().ParentDomain().NumNodes()),
	fWtempV(ShapeFunction().ParentDomain().NumNodes()),
	fWtempA(ShapeFunction().ParentDomain().NumNodes()),
	fWU(ShapeFunction().ParentDomain().NumNodes()),
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
	int totalnodes = nodes_used.Length();
	AutoFill2DT<int> auto_fill(fSolid.NumElements(), 10, 10);
	dArrayT x_atom, centroid;
	LocalArrayT cell_coords(LocalArrayT::kCurrCoords, fSolid.NumElementNodes(), NumSD());
	cell_coords.SetGlobal(curr_coords); // Sets address of cell_coords
	// SetGlobal sets source for SetLocal which copies the relevant data 
	iGridManagerT grid(10, 100, curr_coords, &atoms_used);
	grid.Reset();
	for (int i = 0; i < fSolid.NumElements(); i++) {
	
			/* gives domain (global) nodal coordinates */
	                cell_coords.SetLocal(fSolid.ElementCard(i).NodesX()); 
			/* centroid and radius */
			double radius = parent.AverageRadius(cell_coords, centroid);
			/* candidate particles */
			const AutoArrayT<iNodeT>& hits = grid.HitsInRegion(centroid.Pointer(), 1.01*radius);
			//hits contains atom #'s of atoms within search region.  Atom #'s converted
			//to coordinates below by x_atom.Set()
			//cout << hits.Length() << endl; // returns all atoms currently...Exhaustive search?
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
	dArray2DT atom_coords(NumSD(),fParticlesInCell.MaxMinorDim());
	AutoFill2DT<double> inverse(fParticlesInCell.MajorDim(),10,10);
	for (int i = 0; i < fParticlesInCell.MajorDim(); i++) {

	                fParticlesInCell.RowAlias(i,atom_nums);
	                atom_coords.RowCollect(atom_nums,curr_coords);
			/* gives domain (global) nodal coordinates */
	                cell_coord.SetLocal(fSolid.ElementCard(i).NodesX()); 

	                for (int j = 0; j < fParticlesInCell.MinorDim(i); j++) 
			{
			    //need a ColumnAlias function in nArray2DT similar to RowAlias!
			    point[0] = atom_coords(0,j);
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
	ComputeMass();
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
#if 0
//NOTE: could loop over each output mode and register
//      it with the output separately. for now just register
//      "kAtInc"
	
	/* nodal output */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, fNodalOutputCodes, n_counts);

	/* element output */
	iArrayT e_counts;
	SetElementOutputCodes(IOBaseT::kAtInc, fElementOutputCodes, e_counts);
	ArrayT<StringT> block_ID(fBlockData.Length());
	for (int i = 0; i < block_ID.Length(); i++)
		block_ID[i] = fBlockData[i].ID();

	/* collect variable labels */
	ArrayT<StringT> n_labels(n_counts.Sum());
	ArrayT<StringT> e_labels(e_counts.Sum());
	GenerateOutputLabels(n_counts, n_labels, e_counts, e_labels);

	/* set output specifier */
	StringT set_ID;
	set_ID.Append(ElementSupport().ElementGroupNumber(this) + 1);
	OutputSetT output_set(set_ID, fGeometryCode, block_ID, fConnectivities,
		n_labels, e_labels, false);
		
	/* register and get output ID */
	fOutputID = ElementSupport().RegisterOutput(output_set);
#endif
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

	/* map output flags to count of values */
	iArrayT n_counts;
//	SetNodalOutputCodes(mode, fNodalOutputCodes, n_counts);
	iArrayT e_counts;
//	SetElementOutputCodes(mode, fElementOutputCodes, e_counts);

	/* calculate output values */
	dArray2DT n_values;
	dArray2DT e_values;
//	ComputeOutput(n_counts, n_values, e_counts, e_values);

	/* send to output */
//	ElementSupport().WriteOutput(fOutputID, n_values, e_values);
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

void BridgingScaleT::ComputeMass(void)
{
  /* computes the coarse scale mass matrix once inverse map has
   * been performed */
  const ParentDomainT& parent = ShapeFunction().ParentDomain();
  const FieldT& field = Field();
  iArrayT atoms;
  const dArray2DT& displacements = field[0];
  const dArray2DT& velocities = field[1];
  const dArray2DT& accelerations = field[2];
  dArrayT map, shape(parent.NumNodes()), shape2(parent.NumNodes()), shape3(parent.NumNodes());
  dArrayT temp(NumSD()), disp, vel, acc;
  double dp, ve, ac;

  dMatrixT tempmass(parent.NumNodes());
  for (int i = 0; i < fParticlesInCell.MajorDim(); i++)
  {
      fMass = 0.0, fWtempU = 0.0, fWtempV = 0.0, fWtempA = 0.0;
      fParticlesInCell.RowAlias(i,atoms);
      fInverseMapInCell.RowAlias(i,map);
      for (int j = 0; j < fParticlesInCell.MinorDim(i); j++)
      {
	  // still need to access individual atomic masses
	  displacements.RowAlias(atoms[j],disp);
	  velocities.RowAlias(atoms[j],vel);
	  accelerations.RowAlias(atoms[j],acc);
	  dp = disp[0];
	  ve = vel[0];
	  ac = acc[0];
	  temp = map[j];
	  parent.EvaluateShapeFunctions(temp,shape);
	  parent.EvaluateShapeFunctions(temp,shape2);
	  parent.EvaluateShapeFunctions(temp,shape3);
	  tempmass.Outer(shape,shape);
	  fMass += tempmass;
	  shape *= dp;	  
	  shape2 *= ve;
	  shape3 *= ac;
	  fWtempU += shape;
	  fWtempV += shape2;
	  fWtempA += shape3;
	  //fW.SetRow(atoms[j],shape); -> multiD implementation
      }
      fMassInv.Inverse(fMass);
      fMassInv.Multx(fWtempU,fWU);
      fMassInv.Multx(fWtempV,fWV);
      fMassInv.Multx(fWtempA,fWA);
      ComputeU(displacements, velocities, accelerations);
  }
}

void BridgingScaleT::ComputeU(const dArray2DT& field1, const dArray2DT& field2, const dArray2DT& field3)
{
  /* compute the coarse scale (FEM) solution by projecting the fine scale
     (MD) solution onto a finite dimensional basis space */
  const ParentDomainT& parent = ShapeFunction().ParentDomain();
  dArrayT map, shape1(parent.NumNodes()), shape2(parent.NumNodes());
  dArrayT map1(NumSD()), map2(NumSD());
  iArrayT disp1;
  dArray2DT disp(fParticlesInCell.MaxMinorDim(),NumDOF()), vel(fParticlesInCell.MaxMinorDim(),NumDOF());
  dArray2DT acc(fParticlesInCell.MaxMinorDim(),NumDOF());
  double mult;

  for (int i = 0; i < fParticlesInCell.MajorDim(); i++)
  {
      fInverseMapInCell.RowAlias(i,map);
      fParticlesInCell.RowAlias(i,disp1);
      disp.RowCollect(disp1,field1);
      vel.RowCollect(disp1,field2);
      acc.RowCollect(disp1,field3);
      dMatrixT Projection(fParticlesInCell.MinorDim(i));
      dArrayT ErrorU(fParticlesInCell.MinorDim(i)), FineScaleU(fParticlesInCell.MinorDim(i));
      dArrayT CoarseScaleU(fParticlesInCell.MinorDim(i)), TotalU(fParticlesInCell.MinorDim(i));
      dArrayT ErrorV(fParticlesInCell.MinorDim(i)), FineScaleV(fParticlesInCell.MinorDim(i));
      dArrayT CoarseScaleV(fParticlesInCell.MinorDim(i)), TotalV(fParticlesInCell.MinorDim(i));
      dArrayT ErrorA(fParticlesInCell.MinorDim(i)), FineScaleA(fParticlesInCell.MinorDim(i));
      dArrayT CoarseScaleA(fParticlesInCell.MinorDim(i)), TotalA(fParticlesInCell.MinorDim(i));
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
      cout << "MD = \n" << disp << endl;
      cout << "TotalU = \n" << TotalU << endl;
      cout << "FEM = \n" << CoarseScaleU << endl;
      cout << "Fine Scale = \n" << FineScaleU << endl;
      cout << "MD = \n" << vel << endl;
      cout << "TotalV = \n" << TotalV << endl;
      cout << "FEM = \n" << CoarseScaleV << endl;
      cout << "Fine Scale = \n" << FineScaleV << endl;
      cout << "MD = \n" << acc << endl;
      cout << "TotalA = \n" << TotalA << endl;
      cout << "FEM = \n" << CoarseScaleA << endl;
      cout << "Fine Scale = \n" << FineScaleA << endl;
  }
}
