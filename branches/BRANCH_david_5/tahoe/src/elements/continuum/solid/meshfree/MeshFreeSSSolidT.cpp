/* $Id: MeshFreeSSSolidT.cpp,v 1.21 2005-07-08 23:38:50 paklein Exp $ */
/* created: paklein (09/11/1998) */
#include "MeshFreeSSSolidT.h"

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
MeshFreeSSSolidT::MeshFreeSSSolidT(const ElementSupportT& support):
	SmallStrainT(support),
	fAutoBorder(false),
	fB_wrap(10, fB),
	fMFShapes(NULL),
	fMFFractureSupport(NULL),
	fMeshfreeParameters(NULL)	
{
	SetName("small_strain_meshfree");
}

MeshFreeSSSolidT::~MeshFreeSSSolidT(void)
{
	delete fMFFractureSupport;
}

/* append element equations numbers to the list */
void MeshFreeSSSolidT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
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
void MeshFreeSSSolidT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
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
void MeshFreeSSSolidT::WriteOutput(void)
{
	/* inherited */
	SmallStrainT::WriteOutput();

//TEMP - crack path
	ostream& out = ElementSupport().Output();
	out << "\n time = " << ElementSupport().Time() << '\n';
	fMFFractureSupport->WriteOutput(out);
}

/* returns true if the internal force has been changed since
* the last time step */
GlobalT::RelaxCodeT MeshFreeSSSolidT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = SmallStrainT::RelaxSystem();
	if (fMFFractureSupport->HasActiveCracks())
	{
		//TEMP - need to replace material/element interface for evaluation
		//       of stresses/material properties at the sampling points. This
		//       includes the current element pointer, gradient operators,
		//       state variables, etc. Not implemented
		cout << "\n MeshFreeFDSolidT::RelaxSystem: crack growth not available" << endl;
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
			out << "\n MeshFreeSSSolidT::RelaxSystem:\n";
			out << "               time: " << ElementSupport().Time() << '\n';
			out << " new cutting facets: " << reset_facets.Length() << '\n';
			for (int i = 0; i < reset_facets.Length(); i++)
				facets.PrintRow(reset_facets[i], out);
			out.flush();	
		}
	}
	else if (fMFFractureSupport->CheckGrowth(NULL, &fLocDisp, false)) {
		cout << "\n MeshFreeSSSolidT::RelaxSystem: unexpected crack growth" << endl;
		throw ExceptionT::kGeneralFail;
	}

	return relax;
}

/* returns 1 if DOF's are interpolants of the nodal values */
int MeshFreeSSSolidT::InterpolantDOFs(void) const
{
	/* unknowns are not nodal displacements */
	return 0;
}

/* weight the computational effort of every node */
void MeshFreeSSSolidT::WeightNodalCost(iArrayT& weight) const
{
	fMFFractureSupport->WeightNodes(weight);
}

/* retrieve nodal DOF's */
void MeshFreeSSSolidT::NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const {
	const FieldT& field = Field();
	fMFFractureSupport->GetNodalField(field[0], nodes, DOFs);
}

/* initialize/finalize time increment */
void MeshFreeSSSolidT::InitStep(void)
{
	/* inherited */
	SmallStrainT::InitStep();
	fMFFractureSupport->InitStep();
}

void MeshFreeSSSolidT::CloseStep(void)
{
	/* inherited */
	SmallStrainT::CloseStep();
	fMFFractureSupport->CloseStep();
}

GlobalT::RelaxCodeT MeshFreeSSSolidT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = SmallStrainT::ResetStep();
	fMFFractureSupport->ResetStep();

	return relax;
}

MeshFreeSupportT& MeshFreeSSSolidT::MeshFreeSupport(void) const
{
	if (!fMFShapes) ExceptionT::GeneralFail("MeshFreeSSSolidT::MeshFreeSupport", "shape functions not set");
	return fMFShapes->MeshFreeSupport();
}

/* describe the parameters needed by the interface */
void MeshFreeSSSolidT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SmallStrainT::DefineParameters(list);

	/* shape function storage handled by meshless classes */
	list.RemoveParameter("store_shapefunctions");
	
	ParameterT auto_border(fAutoBorder, "auto_border");
	auto_border.SetDefault(false);
	list.AddParameter(auto_border);
}

/* information about subordinate parameter lists */
void MeshFreeSSSolidT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SmallStrainT::DefineSubs(sub_list);
	
	/* parameters for the meshfree support */
	sub_list.AddSub("meshfree_support_choice", ParameterListT::Once, true);

	/* element support */
	sub_list.AddSub("meshfree_fracture_support");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MeshFreeSSSolidT::NewSub(const StringT& name) const
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
		return SmallStrainT::NewSub(name);
}

/* accept parameter list */
void MeshFreeSSSolidT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MeshFreeSSSolidT::TakeParameterList";

	/* construct meshfree support before calling inherited method because
	 * support class needed to construct shape functions */
	fMFFractureSupport = new MeshFreeFractureSupportT;
	fMFFractureSupport->TakeParameterList(list.GetList("meshfree_fracture_support"));

	/* get parameters needed to construct shape functions */
	fMeshfreeParameters = list.ListChoice(*this, "meshfree_support_choice");

	/* inherited */
	SmallStrainT::TakeParameterList(list);

	/* override any strain-displacement options */
	fStrainDispOpt = kStandardB;

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

	/* register other variable length workspace */
	fMFFractureSupport->Register(fRHS);         // ElementBaseT
	fMFFractureSupport->Register(fNEEvec);      // ContinuumElementT
	fMFFractureSupport->Register(fLHS);         // ElementBaseT

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

	/* output nodal shape function information */
	if (ElementSupport().Logging() == GlobalT::kVerbose)
	{
		/* output file root */
		StringT root;
		root.Root(ElementSupport().InputFile());
		ofstreamT out;

		/* nodal neighbors */
		StringT neighbor_file = root;
		neighbor_file.Append(".", Name(), ".nodal_neighbors");
		out.open(neighbor_file);
		fMFShapes->MeshFreeSupport().WriteNodalNeighbors(out);
		out.close();

		/* nodal shape functions */
		StringT shape_file = root;
		shape_file.Append(".", Name(), ".nodal_phi");
		out.open(shape_file);
		fMFShapes->MeshFreeSupport().WriteNodalShapes(out);
		out.close();
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* initialization functions */
void MeshFreeSSSolidT::SetShape(void)
{
	const char caller[] = "MeshFreeSSSolidT::SetShape";

	/* only support single list of integration cells for now */
	if (fConnectivities.Length() > 1)
		ExceptionT::GeneralFail(caller, "multiple (%d) element blocks not supported",
			fConnectivities.Length());

	/* constructors */
	if (!fMeshfreeParameters) ExceptionT::GeneralFail(caller, "shape function parameters not set");	
	fMFShapes = new MeshFreeShapeFunctionT(GeometryCode(), NumIP(),
		fLocInitCoords, ElementSupport().InitialCoordinates(), *fConnectivities[0], 
		fMFFractureSupport->OffGridNodes(),
		fElementCards.Position(), *fMeshfreeParameters);
	
	/* echo parameters */
	fMFShapes->WriteParameters(ElementSupport().Output());
	
	/* initialize */
	fMFShapes->Initialize();

	/* set base class pointer */
	fShapes = fMFShapes;

	/* set support class */
	fMFFractureSupport->SetShape(fMFShapes);
}

/* current element operations */
bool MeshFreeSSSolidT::NextElement(void)
{
	/* inherited (skip inactive cells) */
	bool OK = SmallStrainT::NextElement();
	while (OK && CurrentElement().Flag() != 1)
		OK = SmallStrainT::NextElement();
	
	/* configure for current element */
	if (OK)
	{
		/* current number of element neighbors */
		int nen = fMFFractureSupport->SetElementNodes(fElementCards.Position());
		
		/* resize */
		fB_wrap.SetDimensions(fB.Rows(), NumSD()*nen);
	}

	return OK;
}

/* driver for nodal value calculations */
void MeshFreeSSSolidT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
	/* set nodal displacements data */
	if (n_codes[iNodalDisp] == NumDOF()) 
		fMFFractureSupport->SetNodalField(Field()[0]);

	/* inherited */
	SmallStrainT::ComputeOutput(n_codes, n_values, e_codes, e_values);

	/* free work space memory */
	if (n_codes[iNodalDisp] == NumDOF()) 
		fMFFractureSupport->FreeNodalField();
}

/***********************************************************************
* Private
***********************************************************************/

/* write displacement field and gradients */
void MeshFreeSSSolidT::WriteField(void)
{
	cout << "\n MeshFreeFDSolidT::WriteField: writing full field" << endl;
	
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
