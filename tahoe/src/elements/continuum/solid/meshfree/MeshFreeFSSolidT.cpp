/* $Id: MeshFreeFSSolidT.cpp,v 1.6 2002-02-20 09:45:48 paklein Exp $ */
/* created: paklein (09/16/1998)                                          */
/* large deformation elasticity with MLS shapefunctions for the           */
/* field (displacement) representation                                    */

#include "MeshFreeFSSolidT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "fstreamT.h"
#include "Constants.h"
#include "ExceptionCodes.h"
#include "MeshFreeShapeFunctionT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"

//TEMP
#include "MaterialListT.h"
#include "StructuralMaterialT.h"

//TEMP - parallel non-repeatability testing
#ifdef __MPI__
#include "mpi.h"
#endif

/* parameters */
const double Pi = acos(-1.0);

/* constructor */
MeshFreeFSSolidT::MeshFreeFSSolidT(FEManagerT& fe_manager):
	TotalLagrangianT(fe_manager),
	MeshFreeFractureSupportT(fFEManager.Input()),
	fB_wrap(10, fB),
	fGradNa_wrap(10, fGradNa),
	fStressStiff_wrap(10, fStressStiff)
{

}

/* data initialization */
void MeshFreeFSSolidT::Initialize(void)
{
	/* inherited */
	TotalLagrangianT::Initialize();

	/* free memory associated with "other" eqnos */
	fEqnos.Free(); // is this OK ? can't be freed earlier b/c of
	               // base class initializations

	/* register dynamic local arrays */
	fLocGroup.Register(fLocDisp);     // ContinuumElementT
	fLocGroup.Register(fLocVel);      // ContinuumElementT
	fLocGroup.Register(fLocAcc);      // ContinuumElementT
	fLocGroup.Register(fLocLastDisp); // TotalLagrangianT

	/* register other variable length workspace */
	fNEEArray.Register(fRHS);         // ElementBaseT
	fNEEArray.Register(fNEEvec);      // ContinuumElementT
	fNEEArray.Register(fTemp2);       // TotalLagrangianT
	fNEEMatrix.Register(fLHS);        // ElementBaseT

	/* dimension */
	fDNa_x_wrap.SetWard(10, fDNa_x, fNumElemNodes);

	/* set MLS data base (connectivities must be set 1st) */
	fMFShapes->SetSupportSize();

	/* exchange nodal parameters (only Dmax for now) */
	iArrayT nodes_in;
	fFEManager.IncomingNodes(nodes_in);
	if (nodes_in.Length() > 0)
	{
		/* send all */
		const dArray2DT& nodal_params = fMFShapes->NodalParameters();
		fFEManager.SendExternalData(nodal_params);

		/* receive */
		dArray2DT all_params_in(nodes_in.Length(), nodal_params.MinorDim());
		fFEManager.RecvExternalData(all_params_in);
	
		/* set values */
		fMFShapes->SetNodalParameters(nodes_in, all_params_in);

		/* skip MLS fit at external nodes */
		fMFShapes->SetSkipNodes(nodes_in);
	}

	/* set nodal neighborhoods */
	fMFShapes->SetNeighborData();

	/* initialize support data */
	iArrayT surface_nodes;
	if (fAutoBorder) SurfaceNodes(surface_nodes);
	MeshFreeFractureSupportT::InitSupport(fFEManager.Input(), fFEManager.Output(),
		fElementCards, surface_nodes, NumDOF(), fNodes->NumNodes(),
		fFEManager.ModelManager());

	/* final MLS initializations */
	fMFShapes->SetExactNodes(fAllFENodes);
	fMFShapes->WriteStatistics(fFEManager.Output());
	
	//TEMP - only works for one material right now, else would have to check
	//       for the material active within the integration cell (element)
	if (HasActiveCracks() && fMaterialList->Length() != 1)
	{	
		cout << "\n MeshFreeFSSolidT::Initialize: can only have 1 material in the group\n";
		cout <<   "     with active cracks" << endl;
		throw eBadInputValue;
	}
	
	/* check for localizing materials */
	if (FractureCriterion() == MeshFreeFractureSupportT::kAcoustic &&
	   !fMaterialList->HasLocalizingMaterials())
	{
		cout << "\n MeshFreeFSSolidT::Initialize: failure criterion requires\n"
		     <<   "     localizing materials: " << MeshFreeFractureSupportT::kAcoustic
		     << endl;
		throw eBadInputValue;
	}

//TEMP - write nodal parameters
#if 0
ostream& out = FEManager().Output();
cout << "\n MeshFreeFSSolidT::Initialize: writing d_max to file :"
<< FEManager().Input().filename() << endl;
out << "\n MeshFreeFSSolidT::Initialize: d_max:\n";
const dArrayT& Dmax = fMFShapes->Dmax();
const iArrayT* node_map = FEManager().NodeMap();
int d_width = out.precision() + kDoubleExtra;
out << setw(kIntWidth) << "node"
<< setw(  d_width) << "d_max" << '\n';
if (node_map)
	for (int i = 0; i < Dmax.Length(); i++)
		out << setw(kIntWidth) << (*node_map)[i] + 1
		    << setw(  d_width) << Dmax[i] << '\n';
else
	for (int j = 0; j < Dmax.Length(); j++)
		out << setw(kIntWidth) << j + 1
		    << setw(  d_width) << Dmax[j] << '\n';
#endif

//TEMP - U connectivities
#if 0
out << "\n MeshFreeFSSolidT::Initialize: connects U: 0, 1,...\n";
fElemNodesEX->WriteNumbered(out);
out.flush();
#endif

//TEMP - trace node
#if 0
int rank, size;
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
int hit_node = -1;
if (size == 1)
	hit_node = 13607 - 1;
else if (size == 2 && rank == 0)
	hit_node = 5420 - 1;
if (hit_node > 0) TraceNode(fFEManager.Output(), hit_node, *this);
#endif
}

/* append element equations numbers to the list */
void MeshFreeFSSolidT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)

	/* get local equations numbers */
	fNodes->SetLocalEqnos(*fElemNodesEX, fElemEqnosEX);

	/* add to list */
	eq_2.Append(&fElemEqnosEX);
	
	/* update active cells */
	int num_active = MarkActiveCells(fElementCards);
	if (num_active != fNumElements)
	{
		/* collect inactive */
		int num_inactive = fNumElements - num_active;
		iArrayT skip_elements(num_inactive);
		int count = 0;
		for (int i = 0; i < fNumElements; i++)
			if (fElementCards[i].Flag() != 1)
				skip_elements[count++] = i;
	
		/* send to MLS */
		fMFShapes->SetSkipElements(skip_elements);
	}
}

/* appends group connectivities to the array */
void MeshFreeFSSolidT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_1)

	/* integration cell field connects */
	connects_2.Append(fElemNodesEX);

	/* nodal field connects */
	connects_2.Append(&(fMFShapes->NodeNeighbors()));
}

/* write output */
void MeshFreeFSSolidT::WriteOutput(IOBaseT::OutputModeT mode)
{
	/* inherited */
	TotalLagrangianT::WriteOutput(mode);

//TEMP - crack path
	if (mode == IOBaseT::kAtInc)
	{
		ostream& out = fFEManager.Output();

		/* time stamp */
		out << "\n time = " << fFEManager.Time() << '\n';
		
		/* inherited */
		MeshFreeFractureSupportT::WriteOutput(out);
	}
}

/* returns true if the internal force has been changed since
* the last time step */
GlobalT::RelaxCodeT MeshFreeFSSolidT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = TotalLagrangianT::RelaxSystem();
	if (HasActiveCracks())
	{
		//TEMP - need to replace material/element interface for evaluation
		//       of stresses/material properties at the sampling points. This
		//       includes the current element pointer, gradient operators,
		//       state variables, etc. Not implemented
		cout << "\n MeshFreeFSSolidT::RelaxSystem: crack growth not available" << endl;
		throw eGeneralFail;
		fElementCards.Current(0);
		
		/* check for crack growth */
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[0];
		StructuralMaterialT* pmat = (StructuralMaterialT*) pcont_mat;
		bool verbose = false;
	 	if (CheckGrowth(pmat, &fLocDisp, verbose))
	 	{
	 		relax = GlobalT::MaxPrecedence(relax, GlobalT::kRelax);
			//TEMP - currently, neighborlists are not reset when cutting
			//    	 facets are inserted because it would require reconfiguring
			//       all of the storage. Therefore, just call for relaxation
			//       since this is needed, though not optimal
			
			/* write new facets to output stream */
			ostream& out = fFEManager.Output();
			const dArray2DT& facets = Facets();
			const ArrayT<int>& reset_facets = ResetFacets();
			out << "\n MeshFreeFSSolidT::RelaxSystem:\n";
			out << "               time: " << fFEManager.Time() << '\n';
			out << " new cutting facets: " << reset_facets.Length() << '\n';
			for (int i = 0; i < reset_facets.Length(); i++)
				facets.PrintRow(reset_facets[i], out);
			out.flush();	
		}
	}
	else if (CheckGrowth(NULL, &fLocDisp, false)) {
		cout << "\n MeshFreeFSSolidT::RelaxSystem: unexpected crack growth" << endl;
		throw eGeneralFail;
	}
		
	return relax;
}

/* returns 1 if DOF's are interpolants of the nodal values */
int MeshFreeFSSolidT::InterpolantDOFs(void) const
{
	/* unknowns are not nodal displacements */
	return 0;
}

/* retrieve nodal DOF's */
void MeshFreeFSSolidT::NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const
{
	/* inherited */
	GetNodalField(fNodes->Displacements(), nodes, DOFs);
}

/* weight the computational effort of every node */
void MeshFreeFSSolidT::WeightNodalCost(iArrayT& weight) const
{
	/* inherited */
	WeightNodes(weight);
}

/* initialize/finalize time increment */
void MeshFreeFSSolidT::InitStep(void)
{
	/* inherited */
	MeshFreeFractureSupportT::InitStep();
}

void MeshFreeFSSolidT::CloseStep(void)
{
	/* inherited */
	MeshFreeFractureSupportT::CloseStep();
}

void MeshFreeFSSolidT::ResetStep(void)
{
	/* inherited */
	MeshFreeFractureSupportT::ResetStep();
}

/***********************************************************************
* Protected
***********************************************************************/

/* print element group data */
void MeshFreeFSSolidT::PrintControlData(ostream& out) const
{
	/* inherited */
	TotalLagrangianT::PrintControlData(out);
	MeshFreeFractureSupportT::PrintControlData(out);
}

/* initialization functions */
void MeshFreeFSSolidT::SetShape(void)
{
	/* only support single list of integration cells for now */
	if (fConnectivities.Length() > 1) {
		cout << "\n MeshFreeFSSolidT::SetShape: multiple element blocks within an"
		     <<   "     element group not supported. Number of blocks: " 
		     << fConnectivities.Length() << endl;
		throw eGeneralFail;
	}

	/* constructors */
	fMFShapes = new MeshFreeShapeFunctionT(fGeometryCode, fNumIP,
		fLocInitCoords, fNodes->InitialCoordinates(), *fConnectivities[0], fOffGridNodes,
		fElementCards.Position(), fFEManager.Input());
	if (!fMFShapes) throw eOutOfMemory;

	/* echo parameters */
	fMFShapes->WriteParameters(fFEManager.Output());
	
	/* initialize (set internal database) */
	fMFShapes->Initialize();
	
	/* set base class pointer */
	fShapes = fMFShapes;
}

/* current element operations */
bool MeshFreeFSSolidT::NextElement(void)
{
	/* inherited (skip inactive cells) */
	bool OK = ElasticT::NextElement();
	while (OK && CurrentElement().Flag() != 1)
		OK = ElasticT::NextElement();

	/* configure for current element */
	if (OK)
	{
		/* current number of element neighbors */
		int nen = SetElementNodes(fElementCards.Position());

		/* resize */
		fStressStiff_wrap.SetDimensions(nen, nen);
		fB_wrap.SetDimensions(fB.Rows(), fNumSD*nen);
		fGradNa_wrap.SetDimensions(fGradNa.Rows(), nen);
		fDNa_x_wrap.Dimension(fNumSD, nen);

		/* reset base class dimensions */
		fNumElemEqnos = nen*fNumDOF;
	}

	return OK;
}

/* driver for nodal value calculations */
void MeshFreeFSSolidT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
	/* set nodal displacements data */
	if (n_codes[iNodalDisp] == fNumDOF) SetNodalField(fNodes->Displacements());

	/* inherited */
	TotalLagrangianT::ComputeOutput(n_codes, n_values, e_codes, e_values);

	/* free work space memory */
	if (n_codes[iNodalDisp] == fNumDOF) FreeNodalField();
}

/***********************************************************************
* Private
***********************************************************************/

/* write displacement field and gradients */
void MeshFreeFSSolidT::WriteField(void)
{
	cout << "\n MeshFreeFSSolidT::WriteField: writing full field" << endl;
	
	const dArray2DT& DOFs = fNodes->Displacements();
	
	/* reconstruct displacement field and all derivatives */
	dArray2DT u;
	dArray2DT Du;
	iArrayT nodes;
	fMFShapes->NodalField(DOFs, u, Du, nodes);

	/* write data */
	ifstreamT& in = fFEManager.Input();
	
	/* output filenames */
	StringT s_u, s_Du;
	s_u.Root(in.filename());
	s_Du.Root(in.filename());
	
	s_u.Append(".u.", fFEManager.StepNumber());
	s_Du.Append(".Du.", fFEManager.StepNumber());
	
	/* open output streams */
	ofstreamT out_u(s_u), out_Du(s_Du);

	/* write */
	for (int i = 0; i < nodes.Length(); i++)
	{
		out_u << setw(kIntWidth) << nodes[i] + 1;
		out_Du << setw(kIntWidth) << nodes[i] + 1;

		u.PrintRow(i, out_u);		
		Du.PrintRow(i, out_Du);		
	}	
	
	/* close */
	out_u.close();
	out_Du.close();
}
