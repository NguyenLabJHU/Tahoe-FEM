/* $Id: MeshFreeFSSolidAxiT.cpp,v 1.2 2004-06-17 07:41:26 paklein Exp $ */
/* created: paklein (09/16/1998) */
#include "MeshFreeFSSolidAxiT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "MeshFreeShapeFunctionT.h"
#include "ModelManagerT.h"
#include "CommManagerT.h"

//TEMP
#include "MaterialListT.h"
#include "SolidMaterialT.h"

using namespace Tahoe;

/* parameters */
const double Pi = acos(-1.0);

/* constructor */
MeshFreeFSSolidAxiT::MeshFreeFSSolidAxiT(const ElementSupportT& support, const FieldT& field):
	TotalLagrangianAxiT(support, field),
	MeshFreeFractureSupportT(ElementSupport().Input()),
	fB_wrap(10, fB),
	fGradNa_wrap(10, fGradNa),
	fStressStiff_wrap(10, fStressStiff)
{
	/* check */
	if (AutoBorder() && ElementSupport().Size() > 1)
		ExceptionT::BadInputValue("MeshFreeFSSolidAxiT::MeshFreeFSSolidAxiT", "auto-border not support in parallel");
}

/* data initialization */
void MeshFreeFSSolidAxiT::Initialize(void)
{
	/* inherited */
	TotalLagrangianAxiT::Initialize();

	/* free memory associated with "other" eqnos */
	fEqnos.Free(); // is this OK ? can't be freed earlier b/c of
	               // base class initializations

	/* register dynamic local arrays */
	fLocGroup.Register(fLocDisp);     // ContinuumElementT
	fLocGroup.Register(fLocVel);      // ContinuumElementT
	fLocGroup.Register(fLocAcc);      // ContinuumElementT
	fLocGroup.Register(fLocLastDisp); // TotalLagrangianAxiT

	/* register other variable length workspace */
	fNEEArray.Register(fRHS);         // ElementBaseT
	fNEEArray.Register(fNEEvec);      // ContinuumElementT
	fNEEArray.Register(fTemp2);       // TotalLagrangianAxiT
	fNEEMatrix.Register(fLHS);        // ElementBaseT

	/* dimension */
	fDNa_x_wrap.SetWard(10, fDNa_x, MeshFreeElementSupportT::NumElementNodes());

	/* set MLS data base (connectivities must be set 1st) */
	fMFShapes->SetSupportSize();

	/* exchange nodal parameters (only Dmax for now) */
	const ArrayT<int>* p_nodes_in = ElementSupport().ExternalNodes();
	if (p_nodes_in)
	{
		/* skip MLS fit at external nodes */
		iArrayT nodes_in;
		nodes_in.Alias(*p_nodes_in);
		fMFShapes->SetSkipNodes(nodes_in);
		
		/* exchange */
		CommManagerT& comm = ElementSupport().CommManager();

		/* send all */
		dArray2DT& nodal_params = fMFShapes->NodalParameters();

		/* initialize the exchange */
		int id = comm.Init_AllGather(nodal_params);
		
		/* do the exchange */
		comm.AllGather(id, nodal_params);
		
		/* clear the communication */
		comm.Clear_AllGather(id);
	}

	/* set nodal neighborhoods */
	fMFShapes->SetNeighborData();

	/* initialize support data */
	iArrayT surface_nodes;
	if (fAutoBorder) {
		ArrayT<StringT> IDs;
		ElementBlockIDs(IDs);
		ElementSupport().Model().SurfaceNodes(IDs, surface_nodes,
			&(ShapeFunction().ParentDomain().Geometry()));
	}
	MeshFreeFractureSupportT::InitSupport(
		ElementSupport().Input(), 
		ElementSupport().Output(),
		fElementCards, 
		surface_nodes, 
		NumDOF(), 
		ElementSupport().NumNodes(),
		&ElementSupport().Model());

	/* final MLS initializations */
	fMFShapes->SetExactNodes(fAllFENodes);
	fMFShapes->WriteStatistics(ElementSupport().Output());
	
	//TEMP - only works for one material right now, else would have to check
	//       for the material active within the integration cell (element)
	if (HasActiveCracks() && fMaterialList->Length() != 1)
	{	
		cout << "\n MeshFreeFSSolidAxiT::Initialize: can only have 1 material in the group\n";
		cout <<   "     with active cracks" << endl;
		throw ExceptionT::kBadInputValue;
	}	
}

/* append element equations numbers to the list */
void MeshFreeFSSolidAxiT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)

	/* get local equations numbers */
	Field().SetLocalEqnos(*fElemNodesEX, fElemEqnosEX);

	/* add to list */
	eq_2.Append(&fElemEqnosEX);
	
	/* update active cells */
	int num_active = MarkActiveCells(fElementCards);
	if (num_active != NumElements())
	{
		/* collect inactive */
		int num_inactive = NumElements() - num_active;
		iArrayT skip_elements(num_inactive);
		int count = 0;
		for (int i = 0; i < NumElements(); i++)
			if (fElementCards[i].Flag() != 1)
				skip_elements[count++] = i;
	
		/* send to MLS */
		fMFShapes->SetSkipElements(skip_elements);
	}
}

/* appends group connectivities to the array */
void MeshFreeFSSolidAxiT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_1)

	/* integration cell field connects */
	connects_2.Append(fElemNodesEX);

	/* nodal field connects */
	connects_2.Append(&(fMFShapes->NodeNeighbors()));
}

/* write output */
void MeshFreeFSSolidAxiT::WriteOutput(void)
{
	/* inherited */
	TotalLagrangianAxiT::WriteOutput();

//TEMP - crack path
	ostream& out = ElementSupport().Output();
	out << "\n time = " << ElementSupport().Time() << '\n';
	MeshFreeFractureSupportT::WriteOutput(out);
}

/* returns true if the internal force has been changed since
* the last time step */
GlobalT::RelaxCodeT MeshFreeFSSolidAxiT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = TotalLagrangianAxiT::RelaxSystem();
	if (HasActiveCracks())
	{
		//TEMP - need to replace material/element interface for evaluation
		//       of stresses/material properties at the sampling points. This
		//       includes the current element pointer, gradient operators,
		//       state variables, etc. Not implemented
		cout << "\n MeshFreeFSSolidAxiT::RelaxSystem: crack growth not available" << endl;
		throw ExceptionT::kGeneralFail;
		fElementCards.Current(0);
		
		/* check for crack growth */
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[0];
		SolidMaterialT* pmat = (SolidMaterialT*) pcont_mat;
		bool verbose = false;
	 	if (CheckGrowth(pmat, &fLocDisp, verbose))
	 	{
	 		relax = GlobalT::MaxPrecedence(relax, GlobalT::kRelax);
			//TEMP - currently, neighborlists are not reset when cutting
			//    	 facets are inserted because it would require reconfiguring
			//       all of the storage. Therefore, just call for relaxation
			//       since this is needed, though not optimal
			
			/* write new facets to output stream */
			ostream& out = ElementSupport().Output();
			const dArray2DT& facets = Facets();
			const ArrayT<int>& reset_facets = ResetFacets();
			out << "\n MeshFreeFSSolidAxiT::RelaxSystem:\n";
			out << "               time: " << ElementSupport().Time() << '\n';
			out << " new cutting facets: " << reset_facets.Length() << '\n';
			for (int i = 0; i < reset_facets.Length(); i++)
				facets.PrintRow(reset_facets[i], out);
			out.flush();	
		}
	}
	else if (CheckGrowth(NULL, &fLocDisp, false)) {
		cout << "\n MeshFreeFSSolidAxiT::RelaxSystem: unexpected crack growth" << endl;
		throw ExceptionT::kGeneralFail;
	}
		
	return relax;
}

/* returns 1 if DOF's are interpolants of the nodal values */
int MeshFreeFSSolidAxiT::InterpolantDOFs(void) const
{
	/* unknowns are not nodal displacements */
	return 0;
}

/* retrieve nodal DOF's */
void MeshFreeFSSolidAxiT::NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const
{
	/* inherited */
	GetNodalField(Field()[0], nodes, DOFs); /* displacements */
}

/* weight the computational effort of every node */
void MeshFreeFSSolidAxiT::WeightNodalCost(iArrayT& weight) const
{
	/* inherited */
	WeightNodes(weight);
}

/* initialize/finalize time increment */
void MeshFreeFSSolidAxiT::InitStep(void)
{
	/* inherited */
	TotalLagrangianAxiT::InitStep();
	MeshFreeFractureSupportT::InitStep();
}

void MeshFreeFSSolidAxiT::CloseStep(void)
{
	/* inherited */
	TotalLagrangianAxiT::CloseStep();
	MeshFreeFractureSupportT::CloseStep();
}

GlobalT::RelaxCodeT MeshFreeFSSolidAxiT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = TotalLagrangianAxiT::ResetStep();
	MeshFreeFractureSupportT::ResetStep();
	return relax;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* print element group data */
void MeshFreeFSSolidAxiT::PrintControlData(ostream& out) const
{
	/* inherited */
	TotalLagrangianAxiT::PrintControlData(out);
	MeshFreeFractureSupportT::PrintControlData(out);
}

/* initialization functions */
void MeshFreeFSSolidAxiT::SetShape(void)
{
	//TEMP - quick and dirty attempt to run with multiple element blocks
	const iArray2DT* mf_connect = fConnectivities[0];
	if (fConnectivities.Length() > 1) {

		/* dimension */
		fConnectsAll.Dimension(NumElements(), TotalLagrangianAxiT::NumElementNodes());

		/* copy connectivities from blocks in */
		int count = 0;
		for (int i = 0; i < fConnectivities.Length(); i++)
		{
			fConnectsAll.BlockRowCopyAt(*fConnectivities[i], count);
			count += fConnectivities[i]->MajorDim();
		}
		mf_connect = &fConnectsAll;

		//TEMP write warning
		cout << "\n MeshFreeFSSolidAxiT::SetShape: WARNING multiple element blocks within the element group" << endl;
	}

	/* constructors */
	fMFShapes = new MeshFreeShapeFunctionT(GeometryCode(), NumIP(),
		fLocInitCoords, ElementSupport().InitialCoordinates(), *mf_connect, fOffGridNodes,
		fElementCards.Position(), ElementSupport().Input());
	if (!fMFShapes) throw ExceptionT::kOutOfMemory;

	/* echo parameters */
	fMFShapes->WriteParameters(ElementSupport().Output());
	
	/* initialize (set internal database) */
	fMFShapes->Initialize();
	
	/* set base class pointer */
	fShapes = fMFShapes;
}

/* current element operations */
bool MeshFreeFSSolidAxiT::NextElement(void)
{
	/* inherited (skip inactive cells) */
	bool OK = SolidElementT::NextElement();
	while (OK && CurrentElement().Flag() != 1)
		OK = SolidElementT::NextElement();

	/* configure for current element */
	if (OK)
	{
		/* current number of element neighbors */
		int nen = SetElementNodes(fElementCards.Position());

		/* resize */
		fStressStiff_wrap.SetDimensions(nen, nen);
		fB_wrap.SetDimensions(fB.Rows(), NumSD()*nen);
		fGradNa_wrap.SetDimensions(fGradNa.Rows(), nen);
		fDNa_x_wrap.Dimension(NumSD(), nen);
	}

	return OK;
}

/* driver for nodal value calculations */
void MeshFreeFSSolidAxiT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
	/* set nodal displacements data */
	if (n_codes[iNodalDisp] == NumDOF()) SetNodalField(Field()[0]);

	/* inherited */
	TotalLagrangianAxiT::ComputeOutput(n_codes, n_values, e_codes, e_values);

	/* free work space memory */
	if (n_codes[iNodalDisp] == NumDOF()) FreeNodalField();
}

/***********************************************************************
 * Private 
 ***********************************************************************/

/* write displacement field and gradients */
void MeshFreeFSSolidAxiT::WriteField(void)
{
	cout << "\n MeshFreeFSSolidAxiT::WriteField: writing full field" << endl;
	
	const dArray2DT& DOFs = Field()[0]; /* displacements */
	
	/* reconstruct displacement field and all derivatives */
	dArray2DT u;
	dArray2DT Du;
	iArrayT nodes;
	fMFShapes->NodalField(DOFs, u, Du, nodes);

	/* write data */
	ifstreamT& in = ElementSupport().Input();
	
	/* output filenames */
	StringT s_u, s_Du;
	s_u.Root(in.filename());
	s_Du.Root(in.filename());
	
	s_u.Append(".u.", ElementSupport().StepNumber());
	s_Du.Append(".Du.", ElementSupport().StepNumber());
	
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
