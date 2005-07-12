/* $Id: MeshFreeFSSolidT.cpp,v 1.21.8.2 2005-07-08 06:19:15 paklein Exp $ */
/* created: paklein (09/16/1998) */
#include "MeshFreeFSSolidT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "MeshFreeShapeFunctionT.h"
#include "ModelManagerT.h"
#include "CommManagerT.h"
#include "ParameterContainerT.h"

#include "MeshFreeFractureSupportT.h"
#include "MeshFreeSupport2DT.h"
#include "MeshFreeSupport3DT.h"

//TEMP
#include "MaterialListT.h"
#include "SolidMaterialT.h"

using namespace Tahoe;

/* parameters */
const double Pi = acos(-1.0);

/* constructor */
MeshFreeFSSolidT::MeshFreeFSSolidT(const ElementSupportT& support):
	TotalLagrangianT(support),
	fAutoBorder(false),
	fB_wrap(10, fB),
	fGradNa_wrap(10, fGradNa),
	fStressStiff_wrap(10, fStressStiff),
	fMFShapes(NULL),
	fMFFractureSupport(NULL),
	fMeshfreeParameters(NULL)		
{
	SetName("large_strain_meshfree");
}

MeshFreeFSSolidT::~MeshFreeFSSolidT(void)
{
	delete fMFFractureSupport;
}

/* append element equations numbers to the list */
void MeshFreeFSSolidT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)

	/* element arrays */
	const RaggedArray2DT<int>& element_nodes = fMFFractureSupport->ElementNodes();
	RaggedArray2DT<int>& element_equations = fMFFractureSupport->ElementEquations();

	/* get local equations numbers */
	Field().SetLocalEqnos(element_nodes, element_equations);

	/* add to list */
	eq_2.Append(&element_equations);
	
	/* update active cells */
	int num_active = fMFFractureSupport->MarkActiveCells(fElementCards);
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
void MeshFreeFSSolidT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_1)

	/* element arrays */
	const RaggedArray2DT<int>& element_nodes = fMFFractureSupport->ElementNodes();

	/* integration cell field connects */
	connects_2.Append(&element_nodes);

	/* nodal field connects */
	connects_2.Append(&(fMFShapes->NodeNeighbors()));
}

/* write output */
void MeshFreeFSSolidT::WriteOutput(void)
{
	/* inherited */
	TotalLagrangianT::WriteOutput();

//TEMP - crack path
	ostream& out = ElementSupport().Output();
	out << "\n time = " << ElementSupport().Time() << '\n';
	fMFFractureSupport->WriteOutput(out);
}

/* returns true if the internal force has been changed since
* the last time step */
GlobalT::RelaxCodeT MeshFreeFSSolidT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = TotalLagrangianT::RelaxSystem();
	if (fMFFractureSupport->HasActiveCracks())
	{
		//TEMP - need to replace material/element interface for evaluation
		//       of stresses/material properties at the sampling points. This
		//       includes the current element pointer, gradient operators,
		//       state variables, etc. Not implemented
		cout << "\n MeshFreeFSSolidT::RelaxSystem: crack growth not available" << endl;
		throw ExceptionT::kGeneralFail;
		fElementCards.Current(0);
		
		/* check for crack growth */
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[0];
		SolidMaterialT* pmat = (SolidMaterialT*) pcont_mat;
		bool verbose = false;
	 	if (fMFFractureSupport->CheckGrowth(pmat, &fLocDisp, verbose))
	 	{
	 		relax = GlobalT::MaxPrecedence(relax, GlobalT::kRelax);
			//TEMP - currently, neighborlists are not reset when cutting
			//    	 facets are inserted because it would require reconfiguring
			//       all of the storage. Therefore, just call for relaxation
			//       since this is needed, though not optimal
			
			/* write new facets to output stream */
			ostream& out = ElementSupport().Output();
			const dArray2DT& facets = fMFFractureSupport->Facets();
			const ArrayT<int>& reset_facets = fMFFractureSupport->ResetFacets();
			out << "\n MeshFreeFSSolidT::RelaxSystem:\n";
			out << "               time: " << ElementSupport().Time() << '\n';
			out << " new cutting facets: " << reset_facets.Length() << '\n';
			for (int i = 0; i < reset_facets.Length(); i++)
				facets.PrintRow(reset_facets[i], out);
			out.flush();	
		}
	}
	else if (fMFFractureSupport->CheckGrowth(NULL, &fLocDisp, false)) {
		cout << "\n MeshFreeFSSolidT::RelaxSystem: unexpected crack growth" << endl;
		throw ExceptionT::kGeneralFail;
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
void MeshFreeFSSolidT::NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const {
	fMFFractureSupport->GetNodalField(Field()[0], nodes, DOFs); /* displacements */
}

/* weight the computational effort of every node */
void MeshFreeFSSolidT::WeightNodalCost(iArrayT& weight) const {
	fMFFractureSupport->WeightNodes(weight);
}

/* initialize/finalize time increment */
void MeshFreeFSSolidT::InitStep(void)
{
	/* inherited */
	TotalLagrangianT::InitStep();
	fMFFractureSupport->InitStep();
}

void MeshFreeFSSolidT::CloseStep(void)
{
	/* inherited */
	TotalLagrangianT::CloseStep();
	fMFFractureSupport->CloseStep();
}

GlobalT::RelaxCodeT MeshFreeFSSolidT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = TotalLagrangianT::ResetStep();
	fMFFractureSupport->ResetStep();

	return relax;
}

MeshFreeSupportT& MeshFreeFSSolidT::MeshFreeSupport(void) const
{
	if (!fMFShapes) ExceptionT::GeneralFail("MeshFreeFSSolidT::MeshFreeSupport", "shape functions not set");
	return fMFShapes->MeshFreeSupport();
}

/* describe the parameters needed by the interface */
void MeshFreeFSSolidT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	TotalLagrangianT::DefineParameters(list);

	/* shape function storage handled by meshless classes */
	list.RemoveParameter("store_shapefunctions");
	
	ParameterT auto_border(fAutoBorder, "auto_border");
	auto_border.SetDefault(false);
	list.AddParameter(auto_border);
}

/* information about subordinate parameter lists */
void MeshFreeFSSolidT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	TotalLagrangianT::DefineSubs(sub_list);
	
	/* parameters for the meshfree support */
	sub_list.AddSub("meshfree_support_choice", ParameterListT::Once, true);

	/* element support */
	sub_list.AddSub("meshfree_fracture_support");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MeshFreeFSSolidT::NewSub(const StringT& name) const
{
	if (name == "meshfree_support_choice")
	{
		ParameterContainerT* mf_choice = new ParameterContainerT(name);
		mf_choice->SetSubSource(this);
		mf_choice->SetListOrder(ParameterListT::Choice);
		
		mf_choice->AddSub("meshfree_support_2D");
		mf_choice->AddSub("meshfree_support_3D");
		
		return mf_choice;
	}
	else if (name == "meshfree_support_2D")
		return new MeshFreeSupport2DT;
	else if (name == "meshfree_support_3D")
		return new MeshFreeSupport3DT;
	else if (name == "meshfree_fracture_support")
		return new MeshFreeFractureSupportT;
	else /* inherited */
		return TotalLagrangianT::NewSub(name);
}

/* accept parameter list */
void MeshFreeFSSolidT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MeshFreeFSSolidT::TakeParameterList";

	/* construct meshfree support before calling inherited method because
	 * support class needed to construct shape functions */
	fMFFractureSupport = new MeshFreeFractureSupportT;
	fMFFractureSupport->TakeParameterList(list.GetList("meshfree_fracture_support"));

	/* get parameters needed to construct shape functions */
	fMeshfreeParameters = list.ListChoice(*this, "meshfree_support_choice");

	/* inherited */
	TotalLagrangianT::TakeParameterList(list);

	/* make field at border nodes nodally exact */
	fAutoBorder = list.GetParameter("auto_border");
	if (fAutoBorder && ElementSupport().Size() > 1)
		ExceptionT::BadInputValue(caller, "auto-border not support in parallel");

	/* free memory associated with "other" eqnos */
	fEqnos.Free(); // is this OK ? can't be freed earlier b/c of
	               // base class initializations

	/* register dynamic local arrays */
	fMFFractureSupport->Register(fLocDisp);     // ContinuumElementT
	fMFFractureSupport->Register(fLocVel);      // ContinuumElementT
	fMFFractureSupport->Register(fLocAcc);      // ContinuumElementT
	fMFFractureSupport->Register(fLocLastDisp); // TotalLagrangianT

	/* register other variable length workspace */
	fMFFractureSupport->Register(fRHS);         // ElementBaseT
	fMFFractureSupport->Register(fNEEvec);      // ContinuumElementT
	fMFFractureSupport->Register(fTemp2);       // TotalLagrangianT
	fMFFractureSupport->Register(fLHS);         // ElementBaseT

	/* dimension */
	fDNa_x_wrap.SetWard(10, fDNa_x, fMFFractureSupport->NumElementNodes());

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
		ElementSupport().ModelManager().SurfaceNodes(IDs, surface_nodes,
			&(ShapeFunction().ParentDomain().Geometry()));
	}

	/* initialize meshfree support class */
	fMFFractureSupport->InitSupport(ElementSupport().Output(),
		fElementCards, 
		surface_nodes,
		NumDOF(), 
		ElementSupport().NumNodes(),
		&ElementSupport().ModelManager());

	/* final MLS initializations */
	fMFShapes->SetExactNodes(fMFFractureSupport->InterpolantNodes());
	fMFShapes->WriteStatistics(ElementSupport().Output());
	
	//TEMP - only works for one material right now, else would have to check
	//       for the material active within the integration cell (element)
	if (fMFFractureSupport->HasActiveCracks() && fMaterialList->Length() != 1)
		ExceptionT::BadInputValue(caller, "can only have 1 material in the group with active cracks");

//TEMP - needs rethinking
#if 0
	/* check for localizing materials */
	if (FractureCriterion() == MeshFreeFractureSupportT::kAcoustic &&
	   !fMaterialList->HasLocalizingMaterials())
	   ExceptionT::BadInputValue(caller, "failure criterion requires localizing materials: %d",
	   	MeshFreeFractureSupportT::kAcoustic);
#endif
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* initialization functions */
void MeshFreeFSSolidT::SetShape(void)
{
	const char caller[] = "MeshFreeFSSolidT::SetShape";

	/* only support single list of integration cells for now */
	if (fConnectivities.Length() > 1)
		ExceptionT::GeneralFail(caller, "multiple (%d) element blocks not supported",
			fConnectivities.Length());

	//TEMP - quick and dirty attempt to run with multiple element blocks
	const iArray2DT* mf_connect = fConnectivities[0];
	if (fConnectivities.Length() > 1) {

		/* dimension */
		fConnectsAll.Dimension(NumElements(), TotalLagrangianT::NumElementNodes());

		/* copy connectivities from blocks in */
		int count = 0;
		for (int i = 0; i < fConnectivities.Length(); i++)
		{
			fConnectsAll.BlockRowCopyAt(*fConnectivities[i], count);
			count += fConnectivities[i]->MajorDim();
		}
		mf_connect = &fConnectsAll;

		//TEMP write warning
		cout << "\n MeshFreeFSSolidT::SetShape: WARNING multiple element blocks within the element group" << endl;
	}

	/* construct */
	if (!fMeshfreeParameters) ExceptionT::GeneralFail(caller, "shape function parameters not set");
	fMFShapes = new MeshFreeShapeFunctionT(GeometryCode(), NumIP(),
		fLocInitCoords, ElementSupport().InitialCoordinates(), *mf_connect, fMFFractureSupport->OffGridNodes(),
		fElementCards.Position(), *fMeshfreeParameters);

	/* echo parameters */
	fMFShapes->WriteParameters(ElementSupport().Output());
	
	/* initialize (set internal database) */
	fMFShapes->Initialize();
	
	/* set base class pointer */
	fShapes = fMFShapes;	

	/* set support class */
	fMFFractureSupport->SetShape(fMFShapes);
}

/* current element operations */
bool MeshFreeFSSolidT::NextElement(void)
{
	/* inherited (skip inactive cells) */
	bool OK = SolidElementT::NextElement();
	while (OK && CurrentElement().Flag() != 1)
		OK = SolidElementT::NextElement();

	/* configure for current element */
	if (OK)
	{
		/* current number of element neighbors */
		int nen = fMFFractureSupport->SetElementNodes(fElementCards.Position());

		/* resize */
		fStressStiff_wrap.SetDimensions(nen, nen);
		fB_wrap.SetDimensions(fB.Rows(), NumSD()*nen);
		fGradNa_wrap.SetDimensions(fGradNa.Rows(), nen);
		fDNa_x_wrap.Dimension(NumSD(), nen);
	}

	return OK;
}

/* driver for nodal value calculations */
void MeshFreeFSSolidT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
	/* set nodal displacements data */
	if (n_codes[iNodalDisp] == NumDOF()) 
		fMFFractureSupport->SetNodalField(Field()[0]);

	/* inherited */
	TotalLagrangianT::ComputeOutput(n_codes, n_values, e_codes, e_values);

	/* free work space memory */
	if (n_codes[iNodalDisp] == NumDOF()) 
		fMFFractureSupport->FreeNodalField();
}

/***********************************************************************
* Private
***********************************************************************/

/* write displacement field and gradients */
void MeshFreeFSSolidT::WriteField(void)
{
	cout << "\n MeshFreeFSSolidT::WriteField: writing full field" << endl;
		
	const dArray2DT& DOFs = Field()[0]; /* displacements */
	
	/* reconstruct displacement field and all derivatives */
	dArray2DT u;
	dArray2DT Du;
	iArrayT nodes;
	fMFShapes->NodalField(DOFs, u, Du, nodes);

	/* write data */
	const StringT& input_file = ElementSupport().InputFile();
	
	/* output filenames */
	StringT s_u, s_Du;
	s_u.Root(input_file);
	s_Du.Root(input_file);
	
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
