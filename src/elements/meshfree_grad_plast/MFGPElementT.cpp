/* $Id: MFGPElementT.cpp,v 1.8 2005-08-24 22:17:15 kyonten Exp $ */
#include "MFGPElementT.h"

/* materials lists */
#include "MFGPMatSupportT.h"
#include "MFGPSSSolidMatT.h"
#include "MFGPMaterialT.h"
#include "MFGPSSSolidMatList2DT.h"
#include "MFGPSSSolidMatList3DT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "CommManagerT.h"
#include "ParameterContainerT.h"
#include "ModelManagerT.h"

#include "MeshFreeFractureSupportT.h"
#include "D3MeshFreeSupport2DT.h"
#include "D3MeshFreeShapeFunctionT.h"

using namespace Tahoe;

/* initialize static data */
const int MFGPElementT::NumNodalOutputCodes = 8;
static const char* NodalOutputNames[] = {
	"coordinates",
	"displacements",
	"lambdas",
	"stress",
	"strain",
	"lap_strain",
	"lap_lambda",
	"material_output"};

const int MFGPElementT::NumElementOutputCodes = 6;
static const char* ElementOutputNames[] = {
	"lambdas",
	"stress",
	"strain",
	"lap_strain",
	"lap_lambda",
	"material_output"};

/* constructor */
MFGPElementT::MFGPElementT(const ElementSupportT& support):
	MFGPAssemblyT(support), //pass the displacement field to the base class
	fAutoBorder(false),
	fB1_wrap(10, fB1),
	fB3_wrap(10, fB3),
	fB4_wrap(10, fB4),
	fPsiLam_wrap(10, fPsiLam),
	fKulambda_wrap(10, fKulambda),
	fKulambda_temp_wrap(10, fKulambda_temp),
	fKlambdau_wrap(10, fKlambdau),
	fKlambdau_temp_wrap(10, fKlambdau_temp),
	fMFFractureSupport_displ(NULL),
	fMFFractureSupport_plast(NULL),
	fMeshfreeParameters(NULL),
	fMassType(kConsistentMass),
	fMFGPMatSupport(NULL)
{
	SetName("mfgp_element");
}

/* destructor */
MFGPElementT::~MFGPElementT(void) 
{  
	delete fMFFractureSupport_displ;
	delete fMFFractureSupport_plast;
	delete fMFGPMatSupport;
}

/* accumulate the residual force on the specified node */
void MFGPElementT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	const char caller[] = "MFGPElementT::AddNodalForce";
	// not implemented
}

/* dummy function */
double MFGPElementT::InternalEnergy (void)
{
	//not implemented
	return 0.0;
}

void MFGPElementT::SendOutput(int kincode)
{
//not implemented
}

/* contribution to the nodal residual forces */
const dArray2DT& MFGPElementT::InternalForce(int group)
{
	const char caller[] = "MFGPElementT::InternalForce";

	/* check */
	if (group != Group())
		ExceptionT::GeneralFail(caller, "expecting solver group %d not %d", 
			Group(), group);
			
	/* must be storing force */
	if (!fStoreInternalForce)
		ExceptionT::GeneralFail(caller, "internal force not being stored");
			
	return fForce;
}

/* relate local and global equation numbers */
void MFGPElementT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)

	/* doing monolithic solution */
	if (fDispl->Group() == fPlast->Group())
	{
		/* element arrays */
		const RaggedArray2DT<int>& connects_displ = fMFFractureSupport_displ->ElementNodes();
		RaggedArray2DT<int>& displ_equ = fMFFractureSupport_displ->ElementEquations();
		
		const RaggedArray2DT<int>& connects_plast = fMFFractureSupport_plast->ElementNodes();
		RaggedArray2DT<int>& plast_equ = fMFFractureSupport_plast->ElementEquations();
		
		/* get local equations numbers */
		fDispl->SetLocalEqnos(connects_displ, displ_equ);
		fPlast->SetLocalEqnos(connects_plast, plast_equ);
		
		/* combine the two equation sets via row-by-row method */
		fMonolithicEquations.Combine(displ_equ, plast_equ);
		
		/* add to list */
		eq_2.Append(&fMonolithicEquations);
	}
	else
	{
#pragma message("correct initialization for staggered solution")
	
		/* ElementBaseT handles equation array for displacements */
		if (ElementSupport().CurrentGroup() == fDispl->Group())
			ElementBaseT::Equations(eq_1, eq_2);

		/* plastic multiplier equation */
		if (ElementSupport().CurrentGroup() == fPlast->Group())
		{
			/* collect local equation numbers */
			//fPlast.SetLocalEqnos(fConnectivities_plast, fEqnos_plast);
		
			//eq_1.Append(&fEqnos_plast);
		}
	}
	
	/* update active cells: displacement */
	int num_active_u = fMFFractureSupport_displ->MarkActiveCells(fElementCards_displ);
	if (num_active_u != NumElements())
	{
		/* collect inactive */
		int num_inactive_u = NumElements() - num_active_u;
		iArrayT skip_elements(num_inactive_u);
		int count = 0;
		for (int i = 0; i < NumElements(); i++)
			if (fElementCards_displ[i].Flag() != 1)
				skip_elements[count++] = i;
	
		/* send to MLS */
		fShapes_displ->SetSkipElements(skip_elements);
	}
	
	/* update active cells: plastic multiplier */
	int num_active_l = fMFFractureSupport_plast->MarkActiveCells(fElementCards_plast);
	if (num_active_l != NumElements())
	{
		/* collect inactive */
		int num_inactive_l = NumElements() - num_active_l;
		iArrayT skip_elements(num_inactive_l);
		int count = 0;
		for (int i = 0; i < NumElements(); i++)
			if (fElementCards_plast[i].Flag() != 1)
				skip_elements[count++] = i;
	
		/* send to MLS */
		fShapes_plast->SetSkipElements(skip_elements);
	}
}

/* appends group connectivities to the array */
void MFGPElementT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_1)

	/* element arrays */
	const RaggedArray2DT<int>& elem_nodes_displ = fMFFractureSupport_displ->ElementNodes();
	//const RaggedArray2DT<int>& elem_nodes_plast = fMFFractureSupport_plast->ElementNodes();
	
	/* mesh-free neighbor sets */
	connects_2.Append(&elem_nodes_displ); 

	/* nodal field connects */
	connects_2.Append(&(fShapes_displ->NodeNeighbors())); 
	//connects_2.Append(&(fShapes_plast->NodeNeighbors())); 	
}

/* weight the computational effort of every node */
void MFGPElementT::WeightNodalCost(iArrayT& weight) const
{
	fMFFractureSupport_displ->WeightNodes(weight);
	fMFFractureSupport_plast->WeightNodes(weight);
}

/* retrieve nodal DOF's */
void MFGPElementT::NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const 
{
	if (DOFs.MinorDim() == 1) { // doesn't support 1D
		const FieldT& field_plast = *fPlast;
		fMFFractureSupport_plast->GetNodalField(field_plast[0], nodes, DOFs);
	}
	else {
		const FieldT& field_displ = *fDispl;
		fMFFractureSupport_displ->GetNodalField(field_displ[0], nodes, DOFs);
	}
}

D3MeshFreeSupportT& MFGPElementT::D3MeshFreeSupport(void) const
{
	if (!fShapes_displ) ExceptionT::GeneralFail("MFGPElementT::D3MeshFreeSupport", "shape functions not set");
	return fShapes_displ->D3MeshFreeSupport();
	
	if (!fShapes_plast) ExceptionT::GeneralFail("MFGPElementT::D3MeshFreeSupport", "shape functions not set");
	return fShapes_plast->D3MeshFreeSupport();
}

/* initialize/finalize step */
void MFGPElementT::InitStep(void)
{
	/* inherited */
	MFGPAssemblyT::InitStep();
	
	fMFFractureSupport_displ->InitStep();
	fMFFractureSupport_plast->InitStep();
}

/* initialize/finalize step */
void MFGPElementT::CloseStep(void)
{
	/* inherited */
	MFGPAssemblyT::CloseStep();
	
	fMFFractureSupport_displ->CloseStep();
	fMFFractureSupport_plast->CloseStep();  
}

/* resets to the last converged solution */
GlobalT::RelaxCodeT MFGPElementT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = MFGPAssemblyT::ResetStep();

	fMFFractureSupport_displ->ResetStep();
	fMFFractureSupport_plast->ResetStep();
	return relax;
}

void MFGPElementT::WriteOutput(void)
{
	/* inherited */
	MFGPAssemblyT::WriteOutput();
	
	//TEMP - crack path
	ostream& out = ElementSupport().Output();
	out << "\n time = " << ElementSupport().Time() << '\n';
	fMFFractureSupport_displ->WriteOutput(out);
}	

/* element level reconfiguration for the current time increment */
GlobalT::RelaxCodeT MFGPElementT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = MFGPAssemblyT::RelaxSystem();
	return relax;
}

/* extract the list of material parameters */
void MFGPElementT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
	const char caller[] = "MFGPElementT::CollectMaterialInfo";
	
	/* initialize */
	mat_params.Clear();
	
	/* collected material parameters */
	int num_blocks = all_params.NumLists("mfgp_element_element_block");
	for (int i = 0; i < num_blocks; i++) {

		/* block information */	
		const ParameterListT& block = all_params.GetList("mfgp_element_element_block", i);
		
		/* resolve material list name */
		if (i == 0) {
			const ParameterListT& mat_list_params = block.GetListChoice(*this, "mfgp_element_material_choice");
			mat_params.SetName(mat_list_params.Name());
		}
		
		/* collect material parameters */
		const ParameterListT& mat_list = block.GetList(mat_params.Name());
		const ArrayT<ParameterListT>& mat = mat_list.Lists();
		mat_params.AddList(mat[0]);
	}
}

/* extrapolate the integration point stresses and internal variables, 
   check the yield condition on the nodes and pass the flag whether 
   the nodes are elastically or plastically loaded
*/ 
void MFGPElementT::CheckNodalYield()
{
	/* dimensions */
	int nen = NumElementNodes();
	int nnd = ElementSupport().NumNodes();
	int nstrs = dSymMatrixT::NumValues(NumSD());
	
	/* number of output values */
	int n_out, n_stressdata, n_matdata, nsd;
	if (nstrs != 4) {
	    nsd = NumSD();
		n_stressdata = dSymMatrixT::NumValues(nsd);
	}
	else {
	    nsd = dSymMatrixT::k3D_plane;
		n_stressdata = dSymMatrixT::NumValues(nsd);
	}
		
	n_matdata = (*fMFGPMatList)[0]->NumOutputVariables();  
	n_out = n_stressdata + n_matdata;

	/* reset averaging workspace */
	ElementSupport().ResetAverage(n_out);

	/* nodal work arrays */
	dArray2DT nodal_space(nen, n_out);
	dArray2DT nodal_all(nen, n_out);
	dArray2DT nodalstress, matdat;

	/* ip values */
	dSymMatrixT cauchy((nstrs != 4) ? NumSD() : dSymMatrixT::k3D_plane), nstr_tmp;
	dArrayT ipmat((*fMFGPMatList)[0]->NumOutputVariables());

	/* set shallow copies */
	double* pall = nodal_space.Pointer();
	nodalstress.Alias(nen, nstrs, pall); pall += nodalstress.Length();
	matdat.Alias(nen, (*fMFGPMatList)[0]->NumOutputVariables(), pall); pall += matdat.Length();
        
	Top();
	while (NextElement())
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			/* initialize */
			nodal_space = 0.0;

			/* global shape function values */
			SetGlobalShape();

			/* integrate */
			fShapes_displ->TopIP();
			while (fShapes_displ->NextIP())
			{
				/* get Cauchy stress */
				const dSymMatrixT& stress = fCurrMaterial->s_ij();
				cauchy.Translate(stress);

				/* stresses */       
				fShapes_displ->Extrapolate(cauchy, nodalstress);

				/* compute material output */
				fCurrMaterial->ComputeOutput(ipmat);
				
				//cout << endl << "ip_int_var = " << endl << ipmat << endl;

				/* store nodal data */
				fShapes_displ->Extrapolate(ipmat, matdat);
			} // while (fShapes_displ->NextIP())

			/* copy in the cols */
			int colcount = 0;
			nodal_all.BlockColumnCopyAt(nodalstress, colcount); colcount += nodalstress.MinorDim();
			nodal_all.BlockColumnCopyAt(matdat     , colcount); colcount += matdat.MinorDim();

			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_all);
		}

	/* get nodally averaged values of stresses and internal variables */
	dArray2DT n_values;
	ElementSupport().OutputUsedAverage(n_values); // nodal values arranged in increasing node order
	int n_nodes = n_values.MajorDim();
	ArrayT<dArrayT> stress_data(n_nodes);
	ArrayT<dArrayT> mat_data(n_nodes);
	for (int i = 0; i < stress_data.Length(); i++) {
		stress_data[i].Dimension(n_stressdata);
		mat_data[i].Dimension(n_matdata);
	}
	
	for (int i = 0; i < n_nodes; i++) {
		for (int j = 0; j < n_out; j++) { 
			if (j < n_stressdata)
				stress_data[i] = n_values(i, j); 
			else
				mat_data[i] = n_values(i,j);
		}
	}
	
	/* calculate nodal yield function and set up yield flags */
	double yield;
	const double tol = 1.0e-10;
	iArrayT yield_flags(n_nodes);
	yield_flags = 0;
	for (int i = 0; i < n_nodes; i++) {
		ComputeNodalYield(stress_data[i], mat_data[i], yield);
		//cout << endl << "stress = " << endl << stress_data[i] << endl;
		//cout << endl << "int variable =" << endl << mat_data[i] << endl;
		//cout << "ff = " << endl << yield << endl;
		if (yield > tol) {
			yield_flags[i] = 1;
			cout << "nodal yield condition satisfied! " << endl;
		}
	}
	
	#if __option(extended_errorcheck)
		if (fNodalYieldFlags.Length() != yield_flags.Length())
		ExceptionT::SizeMismatch("MFGPElementT::CheckNodalYield");
	#endif
	
	fNodalYieldFlags = yield_flags; 
}

/* calculate yield condition at the nodes */ 
const double& MFGPElementT::ComputeNodalYield(const dArrayT& Sig, const dArrayT& qn, 
                double& ff)
{
  double kTemp1, kTemp2, kTemp3, kTemp4;
  double fc, fchi, ffriction, fpress;
  dMatrixT devstress(3,3);
  devstress = 0.0;
  if (Sig.Length() == 4) 
  {
    fpress  = Sig[0]+Sig[1]+Sig[3];
    fpress /=3.;
  	devstress(0,0) = Sig[0] - fpress;
  	devstress(1,1) = Sig[1] - fpress;
  	devstress(2,2) = Sig[3] - fpress;
  	devstress(0,1) = devstress(1,0) = Sig[2];
  }
  if (Sig.Length() == 6)
  {
  	fpress  = Sig[0]+Sig[1]+Sig[2];
  	fpress /=3.;
  	devstress(0,0) = Sig[0] - fpress;
  	devstress(1,1) = Sig[1] - fpress;
  	devstress(2,2) = Sig[2] - fpress;
  	devstress(1,2) = devstress(2,1) = Sig[3];
  	devstress(0,2) = devstress(2,0) = Sig[4];
  	devstress(0,1) = devstress(1,0) = Sig[5];
  }
  
  fc = qn[1];
  ffriction = qn[2];
  fchi = qn[0];
  ff = dMatrixT::Dot(devstress,devstress);
  ff /= 2.;
  kTemp2  = (fc - ffriction*fpress);
  kTemp1  = kTemp2;
  kTemp1 *= kTemp2;
  ff  -= kTemp1;
  kTemp3  = (fc - ffriction*fchi);
  kTemp4  = kTemp3;
  kTemp4 *= kTemp3;
  ff  += kTemp4;
  return ff;
}

/* describe the parameters needed by the interface */
void MFGPElementT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	MFGPAssemblyT::DefineParameters(list);

	/* mass type */
	ParameterT mass_type(ParameterT::Enumeration, "mass_type");
	mass_type.AddEnumeration("automatic", kAutomaticMass);
	mass_type.AddEnumeration("no_mass", kNoMass);
    mass_type.AddEnumeration("consistent_mass", kConsistentMass);
    mass_type.AddEnumeration("lumped_mass", kLumpedMass);
    mass_type.SetDefault(fMassType);
	list.AddParameter(mass_type);
	
	ParameterT auto_border(fAutoBorder, "auto_border");
	auto_border.SetDefault(false);
	list.AddParameter(auto_border);
}

/* information about subordinate parameter lists */
void MFGPElementT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	MFGPAssemblyT::DefineSubs(sub_list);
	
	/* nodal output codes (optional) */
	sub_list.AddSub("mfgp_element_nodal_output", ParameterListT::ZeroOrOnce);
	
	sub_list.AddSub("mfgp_element_element_output", ParameterListT::ZeroOrOnce);
	
	/* element block/material specification */
	sub_list.AddSub("mfgp_element_element_block", ParameterListT::OnePlus);
	
	/* parameters for the meshfree support */
	sub_list.AddSub("mfgp_support_choice", ParameterListT::Once, true);
	
	/* element support */
	sub_list.AddSub("meshfree_fracture_support");
}

/* return the description of the given inline subordinate parameter list */
void MFGPElementT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "mfgp_element_material_choice")
	{
		order = ParameterListT::Choice;
		
		/* list of choices */
		sub_lists.AddSub("mfgp_material_2D");
		sub_lists.AddSub("mfgp_material_3D");
	}
	else /* inherited */
		MFGPAssemblyT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MFGPElementT::NewSub(const StringT& name) const
{
	if (name == "mfgp_support_choice")
	{
		ParameterContainerT* mf_choice = new ParameterContainerT(name);
		mf_choice->SetSubSource(this);
		mf_choice->SetListOrder(ParameterListT::Choice);
		
		/* mfgp_support_choice available only for 2D */
		mf_choice->AddSub("D3_meshfree_support_2D");
		
		return mf_choice;
	}
		
	else if (name == "D3_meshfree_support_2D")
		return new D3MeshFreeSupport2DT;
		
	else if (name == "meshfree_fracture_support")
		return new MeshFreeFractureSupportT;
		
	else if (name == "mfgp_element_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("mfgp_element_material_choice", ParameterListT::Once, true);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	
	else if (name == "mfgp_element_nodal_output")
	{
		ParameterContainerT* node_output = new ParameterContainerT(name);
		
		/* wave speed sampling direction */
		ParameterContainerT wave_direction("wave_direction");
		wave_direction.SetListOrder(ParameterListT::Choice);
		wave_direction.AddSub("Vector_2");
		wave_direction.AddSub("Vector_3");
		node_output->AddSub(wave_direction, ParameterListT::ZeroOrOnce);
		
		/* all false by default */
		for (int i = 0; i < NumNodalOutputCodes; i++) 
		{
			ParameterT output(ParameterT::Integer, NodalOutputNames[i]);
			output.SetDefault(1);
			node_output->AddParameter(output, ParameterListT::ZeroOrOnce);
		}

		return node_output;
	}
	else if (name == "mfgp_element_element_output")
	{
		ParameterContainerT* element_output = new ParameterContainerT(name);
		
		/* all false by default */
		for (int i = 0; i < NumElementOutputCodes; i++) 
		{
			ParameterT output(ParameterT::Integer, ElementOutputNames[i]);
			output.SetDefault(1);
			element_output->AddParameter(output, ParameterListT::ZeroOrOnce);
		}

		return element_output;	
	}
	
	else /* inherited */
		return MFGPAssemblyT::NewSub(name);	
}

/* accept parameter list */
void MFGPElementT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MFGPElementT::TakeParameterList";
	
	/* construct meshfree support before calling inherited method because
	 * support class needed to construct shape functions */
	 /* displacement */
	fMFFractureSupport_displ = new MeshFreeFractureSupportT;
	fMFFractureSupport_displ->TakeParameterList(list.GetList("meshfree_fracture_support"));
	
	 /* plastic multiplier */
	fMFFractureSupport_plast = new MeshFreeFractureSupportT;
	fMFFractureSupport_plast->TakeParameterList(list.GetList("meshfree_fracture_support"));
	
	/* get parameters needed to construct mfgp shape functions */
	fMeshfreeParameters = list.ListChoice(*this, "mfgp_support_choice");
	
	/* set mass type before calling MFGPElementT::TangentType */
	int i_mass = list.GetParameter("mass_type");
	fMassType = int2MassTypeT(i_mass);
	
	/* inherited */
	MFGPAssemblyT::TakeParameterList(list);
	
	/* allocate work space */
	fB1.Dimension(dSymMatrixT::NumValues(NumSD()), NumSD()*NumElementNodes());
	fB3.Dimension(dSymMatrixT::NumValues(NumSD()), NumSD()*NumElementNodes());
	fB4.Dimension(1, NumElementNodes());
	fPsiLam.Dimension(1, NumElementNodes());

	fCUU1.Dimension(dSymMatrixT::NumValues(NumSD()));
	fCUU2.Dimension(dSymMatrixT::NumValues(NumSD()));
	fCULam1.Dimension(dSymMatrixT::NumValues(NumSD()),1);
	fCULam2.Dimension(dSymMatrixT::NumValues(NumSD()),1);
	fCLamU1.Dimension(1, dSymMatrixT::NumValues(NumSD()));
	fCLamU2.Dimension(1, dSymMatrixT::NumValues(NumSD()));
	fCLamLam1.Dimension(1,1);
	fCLamLam2.Dimension(1,1);
	fGradU.Dimension(NumSD());
	
	/* nodal output codes */
	fNodalOutputCodes.Dimension(NumNodalOutputCodes);
	fNodalOutputCodes = IOBaseT::kAtNever;
	qNoExtrap = false;	
	const ParameterListT* node_output = list.List("mfgp_element_nodal_output");
	if (node_output) 
	{
		/* set flags */
		for (int i = 0; i < NumNodalOutputCodes; i++)
		{
			/* look for entry */
			const ParameterT* nodal_value = node_output->Parameter(NodalOutputNames[i]);
			if (nodal_value) {
				int do_write = *nodal_value;

				/* Additional smoothing flags */
	    		if (!qNoExtrap && do_write == 2) {
	    			qNoExtrap = true;
	    			fNodalOutputCodes[i] = IOBaseT::kAtInc;
	    		}
	    		else if (do_write == 1)
	    			fNodalOutputCodes[i] = IOBaseT::kAtInc;
			}
		}
	}

	/* element output codes */
	fElementOutputCodes.Dimension(NumElementOutputCodes);
	fElementOutputCodes = IOBaseT::kAtNever;
	const ParameterListT* element_output = list.List("mfgp_element_element_output");
	if (element_output)
		for (int i = 0; i < NumElementOutputCodes; i++)
		{
			/* look for entry */
			const ParameterT* element_value = element_output->Parameter(ElementOutputNames[i]);
			if (element_value) {
				int do_write = *element_value;
				if (do_write == 1)
					fElementOutputCodes[i] = IOBaseT::kAtInc;
			}
		}
		
	/* allocate strain list */
	fStrain_List.Dimension(fNumIP_plast);
	for (int i = 0; i < fNumIP_plast; i++)
		fStrain_List[i].Dimension(NumSD());
	
	/* allocate "last" strain list */
	fStrain_last_List.Dimension(fNumIP_plast);
	for (int i = 0; i < fNumIP_plast; i++)
		fStrain_last_List[i].Dimension(NumSD());
	
	/* allocate laplacian of strain list */
	fLapStrain_List.Dimension(fNumIP_plast);
	for (int i = 0; i < fNumIP_plast; i++)
		fLapStrain_List[i].Dimension(NumSD());
	
	/* allocate laplacian of "last" strain list */
	fLapStrain_last_List.Dimension(fNumIP_plast);
	for (int i = 0; i < fNumIP_plast; i++)
		fLapStrain_last_List[i].Dimension(NumSD());
	
	/* allocate lambda list */
	int dum = 1;
	fLambda_List.Dimension(fNumIP_plast);
	for (int i = 0; i < fNumIP_plast; i++)
		fLambda_List[i].Dimension(dum);
	
	/* allocate "last" lambda list */
	fLambda_last_List.Dimension(fNumIP_plast);
	for (int i = 0; i < fNumIP_plast; i++)
		fLambda_last_List[i].Dimension(dum);
	
	/* allocate laplacian of lambda list */
	fLapLambda_List.Dimension(fNumIP_plast);
	for (int i = 0; i < fNumIP_plast; i++)
		fLapLambda_List[i].Dimension(dum);
	
	/* allocate laplacian of "last" lambda list */
	fLapLambda_last_List.Dimension(fNumIP_plast);
	for (int i = 0; i < fNumIP_plast; i++)
		fLapLambda_last_List[i].Dimension(dum);
	
	/* make field at border nodes nodally exact */
	fAutoBorder = list.GetParameter("auto_border");
	if (fAutoBorder && ElementSupport().Size() > 1)
		ExceptionT::BadInputValue(caller, "auto-border not support in parallel");
		
	/* free memory associated with "other" eqnos */
	fEqnos.Free(); // is this OK ? can't be freed earlier b/c of
                   // base class initializations
                   
	/* register dynamic local arrays: displacement */
	fMFFractureSupport_displ->Register(u); 
	fMFFractureSupport_displ->Register(u_n);
	//fMFFractureSupport_displ->Register(Du);
	fMFFractureSupport_displ->Register(DDu);
	
	/* register dynamic local arrays: plastic multiplier */
	fMFFractureSupport_plast->Register(lambda); 
	fMFFractureSupport_plast->Register(lambda_n); 
	
	/* register other variable length workspace */
	/* dynamic square matrices */ 
	fMFFractureSupport_displ->Register(fLHS);         // ElementBaseT
	fMFFractureSupport_displ->Register(fKuu);
	fMFFractureSupport_displ->Register(fKuu_temp);  
	fMFFractureSupport_plast->Register(fKlambdalambda);
	fMFFractureSupport_plast->Register(fKlambdalambda_temp);
	
	/* dynamic vectors */
	fMFFractureSupport_displ->Register(fRHS);         // ElementBaseT
	fMFFractureSupport_displ->Register(fNEEvec);     
	fMFFractureSupport_displ->Register(del_u);        
	fMFFractureSupport_plast->Register(del_lambda);   
	fMFFractureSupport_displ->Register(fFu_int);
	fMFFractureSupport_plast->Register(fFlambda);
	fMFFractureSupport_displ->Register(fFu_int_temp);
	fMFFractureSupport_plast->Register(fFlambda_temp);

	/* set MLS data base (connectivities must be set 1st) */
	fShapes_displ->SetSupportSize();
	fShapes_plast->SetSupportSize();
	
	/* exchange nodal parameters (only Dmax for now) */
	/* displacement */
	const ArrayT<int>* p_nodes_in = ElementSupport().ExternalNodes();
	if (p_nodes_in)
	{
		/* skip MLS fit at external nodes */
		iArrayT nodes_in;
		nodes_in.Alias(*p_nodes_in);
		fShapes_displ->SetSkipNodes(nodes_in);
		
		/* exchange */
		CommManagerT& comm = ElementSupport().CommManager();
		
		/* send all */
		dArray2DT& nodal_params = fShapes_displ->NodalParameters();

		/* initialize the exchange */
		int id = comm.Init_AllGather(nodal_params);
		
		/* do the exchange */
		comm.AllGather(id, nodal_params);
		
		/* clear the communication */
		comm.Clear_AllGather(id);
	}
	
	/* plastic multiplier */
	const ArrayT<int>* p_nodes_in_l = ElementSupport().ExternalNodes();
	if (p_nodes_in_l)
	{
		/* skip MLS fit at external nodes */
		iArrayT nodes_in;
		nodes_in.Alias(*p_nodes_in_l);
		fShapes_plast->SetSkipNodes(nodes_in);
		
		/* exchange */
		CommManagerT& comm = ElementSupport().CommManager();
		
		/* send all */
		dArray2DT& nodal_params = fShapes_plast->NodalParameters();

		/* initialize the exchange */
		int id = comm.Init_AllGather(nodal_params);
		
		/* do the exchange */
		comm.AllGather(id, nodal_params);
		
		/* clear the communication */
		comm.Clear_AllGather(id);
	}

	/* set nodal neighborhoods */
	fShapes_displ->SetNeighborData();
	fShapes_plast->SetNeighborData(); //

	/* initialize support data */
	iArrayT surface_nodes;
	if (fAutoBorder) {
		ArrayT<StringT> IDs;
		ElementBlockIDs(IDs);
		ElementSupport().ModelManager().SurfaceNodes(IDs, surface_nodes,
			&(ShapeFunction().ParentDomain().Geometry()));
	}
	
	/* initialize meshfree support class: displacement */
	fMFFractureSupport_displ->InitSupport(ElementSupport().Output(),
		fElementCards_displ, 
		surface_nodes,
		fDispl->NumDOF(), 
		ElementSupport().NumNodes(),
		&ElementSupport().ModelManager());
		
	/* initialize meshfree support class: plastic multiplier */
	fMFFractureSupport_plast->InitSupport(ElementSupport().Output(),
		fElementCards_plast, 
		surface_nodes, 
		fPlast->NumDOF(), 
		ElementSupport().NumNodes(),
		&ElementSupport().ModelManager());
     
	/* final MLS initializations: displacement */
	fShapes_displ->SetExactNodes(fMFFractureSupport_displ->InterpolantNodes());
	fShapes_displ->WriteStatistics(ElementSupport().Output());
	
	/* final MLS initializations: plastic multiplier */
	fShapes_plast->SetExactNodes(fMFFractureSupport_plast->InterpolantNodes());
	fShapes_plast->WriteStatistics(ElementSupport().Output());

	//TEMP - only works for one material right now, else would have to check
	//       for the material active within the integration cell (element)
	if (fMFFractureSupport_displ->HasActiveCracks())
		ExceptionT::BadInputValue(caller, "can only have 1 material in the group with active cracks");
		
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
		fShapes_displ->D3MeshFreeSupport().WriteNodalNeighbors(out);
		out.close();

		/* nodal shape functions */
		StringT shape_file = root;
		shape_file.Append(".", Name(), ".nodal_phi");
		out.open(shape_file);
		fShapes_displ->D3MeshFreeSupport().WriteNodalShapes(out);
		out.close();
	}
} //MFGPElementT::TakeParameterList(const ParameterListT& list)

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct output labels array */
void MFGPElementT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* initialize */
	counts.Dimension(flags.Length());
	counts = 0;

	/* set output flags */
	if (flags[iNodalCoord] == mode)
		counts[iNodalCoord] = NumSD();
	if (flags[iNodalDisp] == mode)
		counts[iNodalDisp] = fDispl->NumDOF();
	if (flags[iNodalLambda] == mode)
		counts[iNodalLambda] = fPlast->NumDOF();
	if (flags[iNodalLapLambda] == mode)
		counts[iNodalLapLambda] = fPlast->NumDOF();
	if (flags[iNodalStress] == mode)
		counts[iNodalStress] = fB1.Rows(); //
	if (flags[iNodalStrain] == mode)
		counts[iNodalStrain] = fB1.Rows(); //
	if (flags[iNodalLapStrain] == mode)
		counts[iNodalLapStrain] = fB1.Rows(); //	
	if (flags[iMaterialData] == mode)
		counts[iMaterialData] = (*fMFGPMatList)[0]->NumOutputVariables();
}

void MFGPElementT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* initialize */
	counts.Dimension(flags.Length());
	counts = 0;
	int ndof_plast = fPlast->NumDOF();

	/* set output flags */
	if (fElementOutputCodes[iIPLambda] == mode) counts[iIPLambda] = ndof_plast*NumIP();
	if (fElementOutputCodes[iIPStress] == mode) counts[iIPStress] = 2*fB1.Rows()*NumIP();
	if (fElementOutputCodes[iIPStrain] == mode) counts[iIPStrain] = 2*fB1.Rows()*NumIP();
	if (fElementOutputCodes[iIPLapStrain] == mode) counts[iIPLapStrain] = 2*fB1.Rows()*NumIP();
	if (fElementOutputCodes[iIPLapLambda] == mode) counts[iIPLapLambda] = ndof_plast*NumIP();
	if (fElementOutputCodes[iIPMaterialData] == mode) 
		counts[iIPMaterialData] = (*fMFGPMatList)[0]->NumOutputVariables()*NumIP();
}

/* construct a new material support and return a pointer */
MFGPMatSupportT* MFGPElementT::NewMFGPMatSupport(MFGPMatSupportT* p) const
{
	if (!p) p = new MFGPMatSupportT(fDispl->NumDOF(), fPlast->NumDOF(), NumIP(), NumIP());

	/* inherited initializations */
	MFGPAssemblyT::NewMFGPMatSupport(p);
	
	/* MFGPElementT sources */
	p->SetLinearStrain(&fStrain_List);
	p->SetLinearStrain_last(&fStrain_last_List);
	p->SetLapLinearStrain(&fLapStrain_List);
	p->SetLapLinearStrain_last(&fLapStrain_last_List);
	
	p->SetLambdaPM(&fLambda_List);
	p->SetLambdaPM_last(&fLambda_last_List);
	p->SetLapLambdaPM(&fLapLambda_List);
	p->SetLapLambdaPM_last(&fLapLambda_last_List);

	return p;
}

/* return a pointer to a new material list */
MFGPMatListT* MFGPElementT::NewMFGPMatList(const StringT& name, int size)
{
	/* resolve dimension */
	int nsd = -1;
	if (name == "mfgp_material_2D") nsd = 2;
	else if (name == "mfgp_material_3D") nsd = 3;
	
	/* no match */
	if (nsd == -1) return NULL;

	/* full list */
	if (size > 0)
	{
		/* material support */
		if (!fMFGPMatSupport) {
			fMFGPMatSupport = TB_DYNAMIC_CAST(MFGPMatSupportT*, NewMFGPMatSupport());
			if (!fMFGPMatSupport)
				ExceptionT::GeneralFail("MFGPElementT::NewMFGPMatList");
		}

		if (nsd == 2)
			return new MFGPSSSolidMatList2DT(size, *fMFGPMatSupport);
		else if (nsd == 3)
			return new MFGPSSSolidMatList3DT(size, *fMFGPMatSupport);
	}
	else
	{
		if (nsd == 2)
			return new MFGPSSSolidMatList2DT;
		else if (nsd == 3)
			return new MFGPSSSolidMatList3DT;
	}
	
	/* no match */
	return NULL;
}

/* form of tangent matrix */
GlobalT::SystemTypeT MFGPElementT::TangentType(void) const
{
	return  MFGPAssemblyT::TangentType();; 
}

/* first derivative of the displacement shape function: [nsd] x [nnd] */  
void MFGPElementT::Set_B1(const dArray2DT& DNa, dMatrixT& B1 )
{
#if __option(extended_errorcheck)
	if (B1.Rows() != dSymMatrixT::NumValues(DNa.MajorDim()) ||
	    B1.Cols() != DNa.Length())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = DNa.MinorDim();
	double* pB1 = B1.Pointer();

	/* 1D */
	if (DNa.MajorDim() == 1)
	{
		const double* pNax = DNa(0);
		for (int i = 0; i < nnd; i++)
			*pB1++ = *pNax++;
	}
	/* 2D */
	else if (DNa.MajorDim() == 2)
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		for (int i = 0; i < nnd; i++)
		{
			
			*pB1++ = *pNax;
			*pB1++ = 0.0;
			*pB1++ = *pNay;

			*pB1++ = 0.0;
			*pB1++ = *pNay++;
			*pB1++ = *pNax++;
		}
	}
	/* 3D */
	else		
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		const double* pNaz = DNa(2);
		for (int i = 0; i < nnd; i++)
		{
			
			*pB1++ = *pNax;
			*pB1++ = 0.0;
			*pB1++ = 0.0;
			*pB1++ = 0.0;
			*pB1++ = *pNaz;
			*pB1++ = *pNay;

			*pB1++ = 0.0;
			*pB1++ = *pNay;
			*pB1++ = 0.0;
			*pB1++ = *pNaz;
			*pB1++ = 0.0;
			*pB1++ = *pNax;

			*pB1++ = 0.0;
			*pB1++ = 0.0;
			*pB1++ = *pNaz++;
			*pB1++ = *pNay++;
			*pB1++ = *pNax++;
			*pB1++ = 0.0;
		}
	}
}

/* laplacian of the displacement shape function: [nsd*nsd] x [nnd] */  
void MFGPElementT::Set_B3(const dArray2DT& DDDNa, dMatrixT& B3)
{
#if __option(extended_errorcheck)
//	if (B3.Rows() != dSymMatrixT::NumValues(sqrt(DDDNa.MajorDim())) ||
//	    B3.Cols() != sqrt(DDDNa.MajorDim())*DDDNa.MinorDim())
//	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = DDDNa.MinorDim();
	int nsd = 3;
	if (DDDNa.MajorDim() == 4)
		nsd = 2;
	else if (DDDNa.MajorDim() == 1)
		nsd = 1;
	double* pB3 = B3.Pointer();

	/* 1D */
	if (nsd == 1)
	{
		const double* pNaxxx = DDDNa(0);
		for (int i = 0; i < nnd; i++)
			*pB3++ = *pNaxxx++;
	}
	/* 2D */
	else if (nsd == 2)
	{
		const double* pNaxxx = DDDNa(0);
		const double* pNayyx = DDDNa(1);
		const double* pNaxxy = DDDNa(2);
		const double* pNayyy = DDDNa(3);
		for (int i = 0; i < nnd; i++)
		{
			
			*pB3++ = *pNaxxx + (*pNayyx);
			*pB3++ = 0.0;
			*pB3++ = *pNayyy + (*pNaxxy);

			*pB3++ = 0.0;
			*pB3++ = *pNayyy + (*pNaxxy);
			*pB3++ = *pNaxxx + (*pNayyx);
		}
	}
	/* 3D */
	else		
	{
		const double* pNaxxx = DDDNa(0);
		const double* pNayyx = DDDNa(1);
		const double* pNazzx = DDDNa(2);
		const double* pNaxxy = DDDNa(3); 
		const double* pNayyy = DDDNa(4); 
		const double* pNazzy = DDDNa(5);
		const double* pNaxxz = DDDNa(6);
		const double* pNayyz = DDDNa(7);
		const double* pNazzz = DDDNa(8);
		
		for (int i = 0; i < nnd; i++)
		{
			
			*pB3++ = *pNaxxx + (*pNayyx) + (*pNazzx);
			*pB3++ = 0.0;
			*pB3++ = 0.0;
			*pB3++ = 0.0;
			*pB3++ = *pNazzz + (*pNaxxz) + (*pNayyz);
			*pB3++ = *pNayyy + (*pNaxxy) + (*pNazzy);

			*pB3++ = 0.0;
			*pB3++ = *pNayyy + (*pNaxxy) + (*pNazzy);
			*pB3++ = 0.0;
			*pB3++ = *pNazzz + (*pNaxxz) + (*pNayyz);
			*pB3++ = 0.0;
			*pB3++ = *pNaxxx + (*pNayyx) + (*pNazzx);

			*pB3++ = 0.0;
			*pB3++ = 0.0;
			*pB3++ = *pNazzz++ + (*pNaxxz++) + (*pNayyz++);
			*pB3++ = *pNayyy++ + (*pNaxxy++) + (*pNazzy++);
			*pB3++ = *pNaxxx++ + (*pNayyx++) + (*pNazzx++);
			*pB3++ = 0.0;
		}
	}
}

/* shape function of plastic multiplier: [1]x[nnd] */ 
void MFGPElementT::Set_PsiLam(const double* Na, dMatrixT& PsiLam) 
{
	int nnd = PsiLam.Cols();
	double* pphi = PsiLam.Pointer();
	for (int i = 0; i < nnd; i++)
	  	*pphi++ = *Na++; //Na[i]		
}

/* laplacian of the shape function of plastic multiplier: [1]x[nnd] */
void MFGPElementT::Set_B4(const dArray2DT& DDNa, dMatrixT& B4)  
{
#if __option(extended_errorcheck)
	if (B4.Cols() != DDNa.MinorDim())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = DDNa.MinorDim();
	double* pB4 = B4.Pointer();

	/* 1D */
	if (DDNa.MajorDim() == 1)
	{
		const double* pNaxx = DDNa(0);
		for (int i = 0; i < nnd; i++)
			*pB4++ = *pNaxx++;
	}
	/* 2D */
	else if (DDNa.MajorDim() == 2)
	{
		const double* pNaxx = DDNa(0);
		const double* pNayy = DDNa(1);
		for (int i = 0; i < nnd; i++)
			*pB4++ = *pNaxx + (*pNayy);
	}
	/* 3D */
	else		
	{
		const double* pNaxx = DDNa(0);
		const double* pNayy = DDNa(1);
		const double* pNazz = DDNa(2);
		
		for (int i = 0; i < nnd; i++)
			*pB4++ = *pNaxx + (*pNayy) + (*pNazz);
	}
}

/* calculate the internal force contribution ("-k*d") */
void MFGPElementT::FormKd(double constK)
{
	/* residual for displacement */
	const double* Det_d    = fShapes_displ->IPDets();
	const double* Weight_d = fShapes_displ->IPWeights();	

	fShapes_displ->TopIP();
	while (fShapes_displ->NextIP())
	{
		/* strain displacement matrix */
		Set_B1(fShapes_displ->Derivatives_U(), fB1);

		/* B1^T * Cauchy stress */
		fB1.MultTx(fCurrMaterial->s_ij(), /* fNEEvec */ fFu_int_temp); // Fu_int: [nsd*nnd]
		
		/* accumulate */
		fFu_int.AddScaled(constK*(*Weight_d++)*(*Det_d++), fFu_int_temp);
	}
	
	/* residual for plastic multiplier */
	const double* Det_p    = fShapes_plast->IPDets();
	const double* Weight_p = fShapes_plast->IPWeights();	

	fShapes_plast->TopIP();
	while (fShapes_plast->NextIP())
	{
		/* Flambda_int */
		int ip = fShapes_plast->CurrIP();
		fShapes_displ->SetIP(ip);
		double scale = constK*Det_p[ip]*Weight_p[ip];
		
		/* plastic multiplier shape functions */
		Set_PsiLam(fShapes_plast->IPShapeU(ip), fPsiLam);
		
		fFlambda_temp = fPsiLam[0];
		fFlambda_temp *= fCurrMaterial->YieldF(); // Flambda_int: [nnd]
		
		/* accumulate */
		fFlambda.AddScaled(scale, fFlambda_temp);
	}
}

/* form the element stiffness matrices */
void MFGPElementT::FormStiffness(double constK)
{
	/* matrix format */
    dMatrixT::SymmetryFlagT format =
        (fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
        dMatrixT::kWhole :
        dMatrixT::kUpperOnly;

	/********DEBUG*******/
	bool print = false; 
	int pos = fElementCards.Position(); 
	if (pos == 1&&0)  
	  print = true; 
	/*******************/
	
	/* integrate element stiffness */
	/* tangent for displacement */
	const double* Det_d    = fShapes_displ->IPDets();
	const double* Weight_d = fShapes_displ->IPWeights();	
	fShapes_displ->TopIP();
	while ( fShapes_displ->NextIP() )
	{
		int ip = fShapes_displ->CurrIP();
		fShapes_plast->SetIP(ip);
		double scale = constK*Det_d[ip]*Weight_d[ip];
	
		/* strain displacement matrix */
		Set_B1(fShapes_displ->Derivatives_U(ip), fB1);
		Set_B3(fShapes_displ->DDDerivatives_U(ip), fB3);
		Set_B4(fShapes_plast->DDerivatives_U(ip), fB4);
		Set_PsiLam(fShapes_plast->IPShapeU(ip), fPsiLam);

		/* get C matrices */
		fCUU1.SetToScaled(scale, fCurrMaterial->c_UU1_ijkl());
		fCUU2.SetToScaled(scale, fCurrMaterial->c_UU2_ijkl());
		fCULam1.SetToScaled(scale, fCurrMaterial->c_ULam1_ij());
		fCULam2.SetToScaled(scale, fCurrMaterial->c_ULam2_ij());
		
		if (print) {
			cout << "\nCUU1: "<<fCurrMaterial->c_UU1_ijkl() << endl;
			cout << "\nCUU2: "<<fCurrMaterial->c_UU2_ijkl() << endl;
			cout << "\nCULam1: "<<fCurrMaterial->c_ULam1_ij() << endl;
			cout << "\nCULam2: "<<fCurrMaterial->c_ULam2_ij() << endl;
		} 
							
		/* form Kuu */
        fKuu.MultQTBQ(fB1, fCUU1, format, dMatrixT::kAccumulate);
        fKuu_temp.MultATBC(fB1, fCUU2, fB3, format, dMatrixT::kAccumulate);
        fKuu += fKuu_temp; // Kuu: [nsd*nnd]x[nsd*nnd]
	
		/* form Kulambda */
		fKulambda.MultATBC(fB1, fCULam1, fPsiLam, format, dMatrixT::kAccumulate);
        fKulambda_temp.MultATBC(fB1, fCULam2, fB4, format, dMatrixT::kAccumulate);
        fKulambda += fKulambda_temp;        // Kulambda: [nsd*nnd]x[nnd] 	
	}
	
	/* tangent for plastic multiplier */
	const double* Det_p    = fShapes_plast->IPDets();
	const double* Weight_p = fShapes_plast->IPWeights();	
	fShapes_plast->TopIP();
	while ( fShapes_plast->NextIP() )
	{
		int ip = fShapes_plast->CurrIP();
		fShapes_displ->SetIP(ip);
		double scale = constK*Det_p[ip]*Weight_p[ip];
	
		/* strain displacement matrix */
		Set_B1(fShapes_displ->Derivatives_U(ip), fB1);
		Set_B3(fShapes_displ->DDDerivatives_U(ip), fB3);
		Set_B4(fShapes_plast->DDerivatives_U(ip), fB4);
		Set_PsiLam(fShapes_plast->IPShapeU(ip), fPsiLam);

		/* get C matrices */
		fCLamU1.SetToScaled(scale, fCurrMaterial->c_LamU1_ij());
		fCLamU2.SetToScaled(scale, fCurrMaterial->c_LamU2_ij());
		fCLamLam1.SetToScaled(scale, fCurrMaterial->c_LamLam1());
		fCLamLam2.SetToScaled(scale, fCurrMaterial->c_LamLam2());
		if (print) {
			cout << "\nCLamU1: "<<fCurrMaterial->c_LamU1_ij() << endl;
			cout << "\nCLamU2: "<<fCurrMaterial->c_LamU2_ij() << endl;
			cout << "\nCLamLam1: "<<fCurrMaterial->c_LamLam1() << endl;
			cout << "\nCLamLam2: "<<fCurrMaterial->c_LamLam2() << endl;
		} 
							
		/* form Klambdau */
		fKlambdau.MultATBC(fPsiLam, fCLamU1, fB1, format, dMatrixT::kAccumulate);
        fKlambdau_temp.MultATBC(fPsiLam, fCLamU2, fB3, format, dMatrixT::kAccumulate);
        fKlambdau += fKlambdau_temp;        // Klambdau :[nnd]x[nsd*nnd]
		
		/* form Klambdalambda */
		fKlambdalambda.MultQTBQ(fPsiLam, fCLamLam1, format, dMatrixT::kAccumulate);
        fKlambdalambda_temp.MultATBC(fPsiLam, fCLamLam2, fB4, format, dMatrixT::kAccumulate);
        fKlambdalambda += fKlambdalambda_temp;        //Klambdalambda: [nnd]x[nnd]   
	}
}

/* set the correct shape functions */
void MFGPElementT::SetShape(void)
{
	const char caller[] = "MFGPElementT::SetShape";
	
	/* displacement field */
	/* constructors */
	if (!fMeshfreeParameters) ExceptionT::GeneralFail(caller, "shape function parameters not set");		
	fShapes_displ = new D3MeshFreeShapeFunctionT(fGeometryCode_displ, fNumIP_displ,
				fInitCoords_displ, ElementSupport().InitialCoordinates(), *(fConnectivities_displ[0]), 
				fMFFractureSupport_displ->OffGridNodes(),
				fElementCards_displ.Position(), *fMeshfreeParameters);
	
	/* echo parameters */
	fShapes_displ->WriteParameters(ElementSupport().Output());
	
	/* initialize */
	fShapes_displ->Initialize();
	
	/* set support class */
	fMFFractureSupport_displ->SetShape(fShapes_displ);
	
	/* plastic multiplier field */
	/* constructors */
	fShapes_plast = new D3MeshFreeShapeFunctionT(fGeometryCode_plast, fNumIP_plast,
				fInitCoords_plast, ElementSupport().InitialCoordinates(), *(fConnectivities_plast[0]), 
				fMFFractureSupport_plast->OffGridNodes(),
				fElementCards_plast.Position(), *fMeshfreeParameters);
	
	/* echo parameters */
	fShapes_plast->WriteParameters(ElementSupport().Output());
	
	/* initialize */
	fShapes_plast->Initialize();
	
	/* set support class */
	fMFFractureSupport_plast->SetShape(fShapes_plast);
}

/* form global shape function derivatives */
void MFGPElementT::SetGlobalShape(void)
{
	/* inherited */
	MFGPAssemblyT::SetGlobalShape();
	
	/* compute the measures of strain/deformation over the element */
	/* loop over integration points */
	for (int i = 0; i < fNumIP_displ; i++)
	{
		/* deformation gradient */
		ShapeFunction().GradU(u, fGradU, i);
	
		/* symmetric part */
		fStrain_List[i].Symmetrize(fGradU);
		
		/* "last" deformation gradient */
		ShapeFunction().GradU(u_n, fGradU, i);
		/* symmetric part */
		fStrain_last_List[i].Symmetrize(fGradU);
		
		/* laplacian of strains */
    	dMatrixT gradgradgrad_u;
    	ShapeFunction().GradGradGradU(u, gradgradgrad_u, i);
    	/** Note: other option is to use B3 which will make the code more compact **/
    	if (NumSD() == 2) {
    		fLapStrain_List[i](0,0)=gradgradgrad_u(0,0)+gradgradgrad_u(0,1);
    		fLapStrain_List[i](1,1)=gradgradgrad_u(1,2)+gradgradgrad_u(1,3);
    		fLapStrain_List[i](0,1)=gradgradgrad_u(0,2)+gradgradgrad_u(0,3);
    		fLapStrain_List[i](0,1)+=gradgradgrad_u(1,0)+gradgradgrad_u(1,1);
    		fLapStrain_List[i](0,1)*=0.5;
    		fLapStrain_List[i](1,0)=fLapStrain_List[i](0,1);
    	}
    	else if (NumSD() == 3) {
    		fLapStrain_List[i](0,0)=gradgradgrad_u(0,0)+gradgradgrad_u(0,1)+gradgradgrad_u(0,2);
    		fLapStrain_List[i](1,1)=gradgradgrad_u(1,3)+gradgradgrad_u(1,4)+gradgradgrad_u(1,5);
    		fLapStrain_List[i](2,2)=gradgradgrad_u(2,6)+gradgradgrad_u(2,7)+gradgradgrad_u(2,8);
    		fLapStrain_List[i](1,2)=gradgradgrad_u(1,6)+gradgradgrad_u(1,7)+gradgradgrad_u(1,8);
    		fLapStrain_List[i](1,2)+=gradgradgrad_u(2,3)+gradgradgrad_u(2,4)+gradgradgrad_u(2,5);
    		fLapStrain_List[i](1,2)*=0.5;
    		fLapStrain_List[i](0,2)=gradgradgrad_u(0,6)+gradgradgrad_u(0,7)+gradgradgrad_u(0,8);
    		fLapStrain_List[i](0,2)+=gradgradgrad_u(2,0)+gradgradgrad_u(2,1)+gradgradgrad_u(2,2);
    		fLapStrain_List[i](0,2)*=0.5;
    		fLapStrain_List[i](0,1)=gradgradgrad_u(0,3)+gradgradgrad_u(0,4)+gradgradgrad_u(0,5);
    		fLapStrain_List[i](0,1)+=gradgradgrad_u(1,0)+gradgradgrad_u(1,1)+gradgradgrad_u(1,2);
    		fLapStrain_List[i](0,1)*=0.5;
    		fLapStrain_List[i](2,1)=fLapStrain_List[i](1,2);
    		fLapStrain_List[i](2,0)=fLapStrain_List[i](0,2);
    		fLapStrain_List[i](1,0)=fLapStrain_List[i](0,1);
    	}
    	else {
    		cout << "\n MFGPElementT::SetGlobalShape: invalid nsd " << endl;
			throw ExceptionT::kBadInputValue;
    	}
    	
    	/* "last" laplacian of strains */
    	dMatrixT gradgradgrad_u_n;
    	ShapeFunction().GradGradGradU(u_n, gradgradgrad_u_n, i);
    	/** Note: other option is to use B3 which will make the code more compact **/
    	if (NumSD() == 2) {
    		fLapStrain_last_List[i](0,0)=gradgradgrad_u_n(0,0)+gradgradgrad_u_n(0,1);
    		fLapStrain_last_List[i](1,1)=gradgradgrad_u_n(1,2)+gradgradgrad_u_n(1,3);
    		fLapStrain_last_List[i](0,1)=gradgradgrad_u_n(0,2)+gradgradgrad_u_n(0,3);
    		fLapStrain_last_List[i](0,1)+=gradgradgrad_u_n(1,0)+gradgradgrad_u_n(1,1);
    		fLapStrain_last_List[i](0,1)*=0.5;
    		fLapStrain_last_List[i](1,0)=fLapStrain_last_List[i](0,1);
    	}
    	else if (NumSD() == 3) {
    		fLapStrain_last_List[i](0,0)=gradgradgrad_u_n(0,0)+gradgradgrad_u_n(0,1)+gradgradgrad_u_n(0,2);
    		fLapStrain_last_List[i](1,1)=gradgradgrad_u_n(1,3)+gradgradgrad_u_n(1,4)+gradgradgrad_u_n(1,5);
    		fLapStrain_last_List[i](2,2)=gradgradgrad_u_n(2,6)+gradgradgrad_u_n(2,7)+gradgradgrad_u_n(2,8);
    		fLapStrain_last_List[i](1,2)=gradgradgrad_u_n(1,6)+gradgradgrad_u_n(1,7)+gradgradgrad_u_n(1,8);
    		fLapStrain_last_List[i](1,2)+=gradgradgrad_u_n(2,3)+gradgradgrad_u_n(2,4)+gradgradgrad_u_n(2,5);
    		fLapStrain_last_List[i](1,2)*=0.5;
    		fLapStrain_last_List[i](0,2)=gradgradgrad_u_n(0,6)+gradgradgrad_u_n(0,7)+gradgradgrad_u_n(0,8);
    		fLapStrain_last_List[i](0,2)+=gradgradgrad_u_n(2,0)+gradgradgrad_u_n(2,1)+gradgradgrad_u_n(2,2);
    		fLapStrain_last_List[i](0,2)*=0.5;
    		fLapStrain_last_List[i](0,1)=gradgradgrad_u_n(0,3)+gradgradgrad_u_n(0,4)+gradgradgrad_u_n(0,5);
    		fLapStrain_last_List[i](0,1)+=gradgradgrad_u_n(1,0)+gradgradgrad_u_n(1,1)+gradgradgrad_u_n(1,2);
    		fLapStrain_last_List[i](0,1)*=0.5;
    		fLapStrain_last_List[i](2,1)=fLapStrain_last_List[i](1,2);
    		fLapStrain_last_List[i](2,0)=fLapStrain_last_List[i](0,2);
    		fLapStrain_last_List[i](1,0)=fLapStrain_last_List[i](0,1);
    	}
    	else {
    		cout << "\n MFGPElementT::SetGlobalShape: invalid nsd " << endl;
			throw ExceptionT::kBadInputValue;
    	}
	} //for (int i = 0; i < NumIP(); i++)
	
	/* compute the plastic multiplier and its laplacian over the element */
	/* loop over integration points */
	for (int i = 0; i < fNumIP_plast; i++)
	{
		/* plastic multiplier */
		IP_Interpolate(lambda, fLambda_List[i], i);
		
		/* "last" plastic multiplier */
		IP_Interpolate(lambda_n, fLambda_last_List[i], i);
		
		/* laplacian of plastic multiplier */
		IP_ComputeLaplacian(lambda, fLapLambda_List[i], i);
		
		/* "last" laplacian of plastic multiplier */
		IP_ComputeLaplacian(lambda_n, fLapLambda_last_List[i], i);
	}
}

/* current element operations */
bool MFGPElementT::NextElement(void)
{
	/* inherited (skip inactive cells) */
	bool OK = MFGPAssemblyT::NextElement();
	while (OK && CurrentElement().Flag() != 1)
		OK = MFGPAssemblyT::NextElement();
	
	/* configure for current element */
	if (OK)
	{
		/* get material pointer */
		MFGPMaterialT* pcont_mat = (*fMFGPMatList)[CurrentElement().MaterialNumber()];
		
		/* cast is safe since class contructs materials list */
		fCurrMaterial = (MFGPMaterialT*) pcont_mat;
		
		/* current number of element neighbors */
		int nen_displ = fMFFractureSupport_displ->SetElementNodes(fElementCards_displ.Position());
		int nen_plast = fMFFractureSupport_plast->SetElementNodes(fElementCards_plast.Position());
		
		/* resize off-diagonal, non-square matrices */
		int dof_displ = fDispl->NumDOF();
		int dof_plast = fPlast->NumDOF();
		fKulambda_wrap.SetDimensions(nen_displ*dof_displ, nen_plast*dof_plast);
		fKulambda_temp_wrap.SetDimensions(nen_displ*dof_displ, nen_plast*dof_plast);
		fKlambdau_wrap.SetDimensions(nen_plast*dof_plast, nen_displ*dof_displ);
		fKlambdau_temp_wrap.SetDimensions(nen_plast*dof_plast, nen_displ*dof_displ);
		
		/* resize B matrices */
		fB1_wrap.SetDimensions(fB1.Rows(), NumSD()*nen_displ);
		fB3_wrap.SetDimensions(fB3.Rows(), NumSD()*nen_displ);
		fB4_wrap.SetDimensions(fB4.Rows(), nen_plast);
		fPsiLam_wrap.SetDimensions(fPsiLam.Rows(), nen_plast);
	}

	return OK;
}

/* driver for output value calculations
 * extrapolate the integration point stresses and strains and extrapolate */
void MFGPElementT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
	/* set nodal displacements data */
	if (n_codes[iNodalDisp] == fDispl->NumDOF()) {
		const FieldT& field_displ = *fDispl; 
		fMFFractureSupport_displ->SetNodalField(field_displ[0]);
	}
		
	/* set nodal plastic multipliers data */
	if (n_codes[iNodalLambda] == fPlast->NumDOF()) {
		const FieldT& field_plast = *fPlast; 
		fMFFractureSupport_plast->SetNodalField(field_plast[0]);
	}
	
	/* number of output values */
	int n_out = n_codes.Sum();
	int e_out = e_codes.Sum();
        
	int n_extrap;
	n_extrap = n_out;
	       
	/* nothing to output */
	if (n_out == 0 && e_out == 0) return;

	/* dimensions */
	int nsd = NumSD();
	int nen = NumElementNodes();
	int nnd = ElementSupport().NumNodes();
	int nstrs = fB1.Rows();

	/* reset averaging workspace */
	ElementSupport().ResetAverage(n_extrap);

	/* allocate element results space */
	e_values.Dimension(NumElements(), e_out);

	/* nodal work arrays */
	dArray2DT nodal_space(nen, n_extrap);
	dArray2DT nodal_all(nen, n_extrap);
	dArray2DT coords, disp, nodal_lambda, nodal_lap_lambda;
	dArray2DT nodalstress, nodalstrain, nodallapstrain, matdat;

	/* ip values */
	dArrayT iplambda(n_codes[iNodalLambda]);
	dArrayT iplaplambda(n_codes[iNodalLapLambda]);
	dSymMatrixT cauchy((nstrs != 4) ? nsd : dSymMatrixT::k3D_plane), nstr_tmp;
	dSymMatrixT strain(cauchy.Rows()), nstr_tmp2;
	dSymMatrixT strain_tmp(strain.Rows());
	dSymMatrixT lapstrain(strain.Rows()), nstr_tmp3;
	dSymMatrixT lapstrain_tmp(strain.Rows());
	dArrayT ipmat(n_codes[iMaterialData]);

	/* set shallow copies */
	double* pall = nodal_space.Pointer();
	coords.Alias(nen, n_codes[iNodalCoord], pall); pall += coords.Length();
	disp.Alias(nen, n_codes[iNodalDisp], pall)   ; pall += disp.Length();
	nodal_lambda.Alias(nen, n_codes[iNodalLambda], pall); pall += nodal_lambda.Length();
	nodal_lap_lambda.Alias(nen, n_codes[iNodalLapLambda], pall); pall += nodal_lap_lambda.Length();
	nodalstress.Alias(nen, n_codes[iNodalStress], pall); pall += nodalstress.Length();
	nodalstrain.Alias(nen, n_codes[iNodalStrain], pall); pall += nodalstrain.Length();
	nodallapstrain.Alias(nen, n_codes[iNodalLapStrain], pall); pall += nodallapstrain.Length();
	matdat.Alias(nen, n_codes[iMaterialData], pall)    ; pall += matdat.Length();
        
	/* element work arrays */
	dArrayT element_values(e_values.MinorDim());
	pall = element_values.Pointer();
	dArrayT ip_coords(nsd);

	dArray2DT ip_lambda;
	if (e_codes[iIPLambda]) {
		ip_lambda.Alias(NumIP(), e_codes[iIPLambda]/NumIP(), pall);
		pall += ip_lambda.Length();
		iplambda.Dimension(ip_lambda.MinorDim());
	}
	
	dArray2DT ip_lap_lambda;
	if (e_codes[iIPLapLambda]) {
		ip_lap_lambda.Alias(NumIP(), e_codes[iIPLapLambda]/NumIP(), pall);
		pall += ip_lap_lambda.Length();
		iplaplambda.Dimension(ip_lap_lambda.MinorDim());
	}
	
	dArray2DT ip_stress;
	if (e_codes[iIPStress]) {
		ip_stress.Alias(NumIP(), e_codes[iIPStress]/NumIP(), pall);
		pall += ip_stress.Length();
	}
	
	dArray2DT ip_strain;
	if (e_codes[iIPStrain]) {
		ip_strain.Alias(NumIP(), e_codes[iIPStrain]/NumIP(), pall);
		pall += ip_strain.Length();
	}
	
	dArray2DT ip_lap_strain;
	if (e_codes[iIPLapStrain]) {
		ip_lap_strain.Alias(NumIP(), e_codes[iIPLapStrain]/NumIP(), pall);
		pall += ip_lap_strain.Length();
	}
	
	dArray2DT ip_material_data;
	if (e_codes[iIPMaterialData]) {
		ip_material_data.Alias(NumIP(), e_codes[iIPMaterialData]/NumIP(), pall);
		pall += ip_material_data.Length();
		ipmat.Dimension(ip_material_data.MinorDim());
	}

	bool is_axi = Axisymmetric();
	double Pi2 = 2.0*acos(-1.0);

	Top();
	while (NextElement())
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			/* initialize */
			nodal_space = 0.0;

			/* global shape function values */
			SetGlobalShape();

			/* coordinates, displacements, and plastic multipliers all at once */
			if (n_codes[iNodalCoord]) fInitCoords_displ.ReturnTranspose(coords);
			if (n_codes[iNodalDisp]) NodalDOFs(CurrentElement().NodesX(), disp);
			if (n_codes[iNodalLambda]) NodalDOFs(CurrentElement().NodesX(), nodal_lambda);
 			
			/* initialize element values */
			const double* j = fShapes_displ->IPDets();
			const double* w = fShapes_displ->IPWeights();

			/* integrate */
			dArray2DT Na_X_ip_w;
			fShapes_displ->TopIP();
			while (fShapes_displ->NextIP())
			{
				/* density may change with integration point */
				double density = fCurrMaterial->Density();

				/* element integration weight */
				double ip_w = (*j++)*(*w++);
				if (is_axi) {
					fShapes_displ->IPCoords(ip_coords);
					ip_w *= Pi2*ip_coords[0]; /* radius is the x1 coordinate */
				}

				if (qNoExtrap) {
					Na_X_ip_w.Dimension(nen,1);
					for (int k = 0; k < nen; k++)
						Na_X_ip_w(k,0) = 1.;
				}
				
				/* get Cauchy stress */
				const dSymMatrixT& stress = fCurrMaterial->s_ij();
				cauchy.Translate(stress);

				/* stresses */
				if (n_codes[iNodalStress]) {        
					if (qNoExtrap)
						for (int k = 0; k < nen; k++)
							nodalstress.AddToRowScaled(k,Na_X_ip_w(k,0),cauchy);
					else
						fShapes_displ->Extrapolate(cauchy, nodalstress);
				}

				if (e_codes[iIPStress]) {
					double* row = ip_stress(fShapes_displ->CurrIP());
					nstr_tmp.Set(nsd, row);
					nstr_tmp = cauchy;
				}
				
				/* get strain and its laplacian */
				fCurrMaterial->Strain(strain_tmp); fCurrMaterial->LapStrain(lapstrain_tmp);
				strain.Translate(strain_tmp); lapstrain.Translate(lapstrain_tmp);

				/* strains */
				if (n_codes[iNodalStrain]) {        
					if (qNoExtrap)
						for (int k = 0; k < nen; k++)
							nodalstrain.AddToRowScaled(k,Na_X_ip_w(k,0),strain);
					else
						fShapes_displ->Extrapolate(strain, nodalstrain);
				}
				
				if (n_codes[iNodalLapStrain]) {        
					if (qNoExtrap)
						for (int k = 0; k < nen; k++)
							nodallapstrain.AddToRowScaled(k,Na_X_ip_w(k,0),lapstrain);
					else
						fShapes_displ->Extrapolate(lapstrain, nodallapstrain);
				}

				if (e_codes[iIPStrain]) {
					double* row = ip_strain(fShapes_displ->CurrIP());
					nstr_tmp2.Set(nsd, row);
					nstr_tmp2 = strain;
				}
				
				if (e_codes[iIPLapStrain]) {
					double* row = ip_lap_strain(fShapes_displ->CurrIP());
					nstr_tmp3.Set(nsd, row);
					nstr_tmp3 = lapstrain;
				}

				/* get laplacian of plastic multiplier */
				if (n_codes[iNodalLapLambda] || e_codes[iIPLapLambda] || e_codes[iIPLambda]) 
				{
					/* store nodal values */
					if (n_codes[iNodalLapLambda]) {
						fCurrMaterial->LapPlasticMultiplier(iplaplambda);
						if (qNoExtrap)
							for (int k = 0; k < nen; k++)
								nodal_lap_lambda.AddToRowScaled(k,Na_X_ip_w(k,0), iplaplambda);
						else 
							fShapes_displ->Extrapolate(iplaplambda, nodal_lap_lambda);
					}
					
					/* store element values */
					if (e_codes[iIPLapLambda]) { 
						fCurrMaterial->LapPlasticMultiplier(iplaplambda);
						ip_lap_lambda.SetRow(fShapes_displ->CurrIP(), iplaplambda);
					}
					
					if (e_codes[iIPLambda]) fCurrMaterial->PlasticMultiplier(iplaplambda); 
					IP_Interpolate(lambda, iplambda, fShapes_displ->CurrIP());
				} 
				
				/* material stuff */
				if (n_codes[iMaterialData] || e_codes[iIPMaterialData])
				{
					/* compute material output */
					fCurrMaterial->ComputeOutput(ipmat);

					/* store nodal data */
					if (n_codes[iMaterialData]) {
						if (qNoExtrap)
							for (int k = 0; k < nen; k++)
								matdat.AddToRowScaled(k,Na_X_ip_w(k,0),ipmat);
						else 
							fShapes_displ->Extrapolate(ipmat, matdat);
					}

					/* store element data */
					if (e_codes[iIPMaterialData]) ip_material_data.SetRow(fShapes_displ->CurrIP(), ipmat);
					
				}
			} // while (fShapes->NextIP())

			/* copy in the cols */
			int colcount = 0;
			nodal_all.BlockColumnCopyAt(disp       , colcount); colcount += disp.MinorDim();
			nodal_all.BlockColumnCopyAt(nodal_lambda       , colcount); colcount += nodal_lambda.MinorDim();
			nodal_all.BlockColumnCopyAt(coords     , colcount); colcount += coords.MinorDim();
			if (qNoExtrap) {
				double nip(fShapes_displ->NumIP());
				nodal_lap_lambda /= nip;
				nodalstress /= nip;
				matdat /= nip;
			}
			nodal_all.BlockColumnCopyAt(nodal_lap_lambda     , colcount); colcount += nodal_lap_lambda.MinorDim();
			nodal_all.BlockColumnCopyAt(nodalstress, colcount); colcount += nodalstress.MinorDim();
			nodal_all.BlockColumnCopyAt(nodalstrain, colcount); colcount += nodalstrain.MinorDim();
			nodal_all.BlockColumnCopyAt(nodallapstrain, colcount); colcount += nodallapstrain.MinorDim();
			nodal_all.BlockColumnCopyAt(matdat     , colcount); colcount += matdat.MinorDim();

			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_all);

			/* store results */
			e_values.SetRow(CurrElementNumber(), element_values);
		} // while (NextElement())

	/* get nodally averaged values */
	const OutputSetT& output_set = ElementSupport().OutputSet(fOutputID);
	const iArrayT& nodes_used = output_set.NodesUsed();
	dArray2DT extrap_values(nodes_used.Length(), n_extrap);
	extrap_values.RowCollect(nodes_used, ElementSupport().OutputAverage());

	int tmpDim = extrap_values.MajorDim();
	n_values.Dimension(tmpDim,n_out);
	n_values.BlockColumnCopyAt(extrap_values,0);
	
	/* free work space memory */
	if (n_codes[iNodalDisp] == fDispl->NumDOF()) 
		fMFFractureSupport_displ->FreeNodalField();
	if (n_codes[iNodalLambda] == fPlast->NumDOF()) 
		fMFFractureSupport_plast->FreeNodalField();
}

/* construct output labels array */
void MFGPElementT::GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels, 
	const iArrayT& e_codes, ArrayT<StringT>& e_labels) const
{
	const char caller[] = "MFGPElementT::GenerateOutputLabels";

	/* allocate */
	n_labels.Dimension(n_codes.Sum());

	int count = 0;
	if (n_codes[iNodalDisp]) {
		/* labels from the displacement field */
		const ArrayT<StringT>& labels = fDispl->Labels();
		for (int i = 0; i < labels.Length(); i++)
			n_labels[count++] = labels[i];
	}
	
	if (n_codes[iNodalLambda]) {
		/* labels from the plastic multiplier field */
		const ArrayT<StringT>& labels = fPlast->Labels();
		for (int i = 0; i < labels.Length(); i++)
			n_labels[count++] = labels[i];
	}
	
	if (n_codes[iNodalCoord]) {
		const char* xlabels[] = {"x1", "x2", "x3"};
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = xlabels[i];
	}
	
	if (n_codes[iNodalLapLambda]) {
		const char* lap_lambda_labels[] = {"lap_lambda"};
		const char** lap_lambdalabels = NULL;
		lap_lambdalabels = lap_lambda_labels;

		n_labels[count++] = lap_lambdalabels[0];
	}

	if (n_codes[iNodalStress]) {
		const char* slabels1D[] = {"s11"};
		const char* slabels2D[] = {"s11", "s22", "s12"};
		const char* slabels2D_axi[] = {"srr", "szz", "srz", "stt"};
		const char* slabels3D[] = {"s11", "s22", "s33", "s23", "s13", "s12"};
		int nstrs = fB1.Rows();
		const char** slabels = NULL;
		if (nstrs == 1)
			slabels = slabels1D;
		else if (nstrs == 3)
			slabels = slabels2D;
		else if (nstrs == 4)
			slabels = slabels2D_axi;
		else if (nstrs == 6)
			slabels = slabels3D;
		else
			ExceptionT::GeneralFail(caller);

		for (int i = 0; i < nstrs; i++)
			n_labels[count++] = slabels[i];
	}
	
	if (n_codes[iNodalStrain]) {
		const char* epslabels1D[] = {"eps11"};
		const char* epslabels2D[] = {"eps11", "eps22", "eps12"};
		const char* epslabels2D_axi[] = {"epsrr", "epszz", "epsrz", "epstt"};
		const char* epslabels3D[] = {"eps11", "eps22", "eps33", "eps23", "eps13", "eps12"};
		int nstrs = fB1.Rows();
		const char** epslabels = NULL;
		if (nstrs == 1)
			epslabels = epslabels1D;
		else if (nstrs == 3)
			epslabels = epslabels2D;
		else if (nstrs == 4)
			epslabels = epslabels2D_axi;
		else if (nstrs == 6)
			epslabels = epslabels3D;
		else
			ExceptionT::GeneralFail(caller);

		for (int i = 0; i < nstrs; i++)
			n_labels[count++] = epslabels[i];
	}
	
	if (n_codes[iNodalLapStrain]) {
		const char* lap_epslabels1D[] = {"lap_eps11"};
		const char* lap_epslabels2D[] = {"lap_eps11", "lap_eps22", "lap_eps12"};
		const char* lap_epslabels2D_axi[] = {"lap_epsrr", "lap_epszz", "lap_epsrz", "lap_epstt"};
		const char* lap_epslabels3D[] = {"lap_eps11", "lap_eps22", "lap_eps33", "lap_eps23", "lap_eps13", "lap_eps12"};
		int nstrs = fB1.Rows();
		const char** lap_epslabels = NULL;
		if (nstrs == 1)
			lap_epslabels = lap_epslabels1D;
		else if (nstrs == 3)
			lap_epslabels = lap_epslabels2D;
		else if (nstrs == 4)
			lap_epslabels = lap_epslabels2D_axi;
		else if (nstrs == 6)
			lap_epslabels = lap_epslabels3D;
		else
			ExceptionT::GeneralFail(caller);

		for (int i = 0; i < nstrs; i++)
			n_labels[count++] = lap_epslabels[i];
	}

	/* material output labels */
	if (n_codes[iMaterialData]) {
		ArrayT<StringT> matlabels;
		(*fMFGPMatList)[0]->OutputLabels(matlabels);	
		
		for (int i = 0; i < matlabels.Length(); i++)
			n_labels[count++] = matlabels[i];
	}

	/* allocate */
	e_labels.Dimension(e_codes.Sum());
	count = 0;
	if (n_codes[iIPLambda]) {
		const char* lambda_labels[] = {"lambda"};
		const char** lambdalabels = NULL;
		lambdalabels = lambda_labels;
		
		/* over integration points */
		for (int j = 0; j < NumIP(); j++)
		{
			StringT ip_label;
			ip_label.Append("ip", j+1);
			
			/* over lap lambda components */
			e_labels[count++].Clear();
			e_labels[count++].Append(ip_label, ".", lambdalabels[0]);
		}		
	}
	
	if (n_codes[iIPLapLambda]) {
		const char* lap_lambda_labels[] = {"lap_lambda"};
		const char** lap_lambdalabels = NULL;
		lap_lambdalabels = lap_lambda_labels;
		
		/* over integration points */
		for (int j = 0; j < NumIP(); j++)
		{
			StringT ip_label;
			ip_label.Append("ip", j+1);
			
			/* over lap lambda components */
			e_labels[count++].Clear();
			e_labels[count++].Append(ip_label, ".", lap_lambdalabels[0]);
		}		
	}
	
	if (e_codes[iIPStress]) {
		const char* slabels1D[] = {"s11", "e11"};
		const char* slabels2D[] = {"s11", "s22", "s12"};
		const char* slabels2D_axi[] = {"srr", "szz", "srz", "stt"};
		const char* slabels3D[] = {"s11", "s22", "s33", "s23", "s13", "s12"};
		int nstrs = fB1.Rows();
		const char** slabels = NULL;
		if (nstrs == 1)
			slabels = slabels1D;
		else if (nstrs == 3)
			slabels = slabels2D;
		else if (nstrs == 4)
			slabels = slabels2D_axi;
		else if (nstrs == 6)
			slabels = slabels3D;
		else
			ExceptionT::GeneralFail(caller);

		/* over integration points */
		for (int j = 0; j < NumIP(); j++)
		{
			StringT ip_label;
			ip_label.Append("ip", j+1);
			
			/* over stress/strain components */
			for (int i = 0; i < 2*nstrs; i++) {
				e_labels[count].Clear();
				e_labels[count].Append(ip_label, ".", slabels[i]);
				count++;
			}
		}		
	}
	
	if (e_codes[iIPStrain]) {
		const char* epslabels1D[] = {"eps11", "eps11"};
		const char* epslabels2D[] = {"eps11", "eps22", "eps12"};
		const char* epslabels2D_axi[] = {"epsrr", "epszz", "epsrz", "epstt"};
		const char* epslabels3D[] = {"eps11", "eps22", "eps33", "eps23", "eps13", "eps12"};
		int nstrs = fB1.Rows();
		const char** epslabels = NULL;
		if (nstrs == 1)
			epslabels = epslabels1D;
		else if (nstrs == 3)
			epslabels = epslabels2D;
		else if (nstrs == 4)
			epslabels = epslabels2D_axi;
		else if (nstrs == 6)
			epslabels = epslabels3D;
		else
			ExceptionT::GeneralFail(caller);

		/* over integration points */
		for (int j = 0; j < NumIP(); j++)
		{
			StringT ip_label;
			ip_label.Append("ip", j+1);
			
			/* over stress/strain components */
			for (int i = 0; i < 2*nstrs; i++) {
				e_labels[count].Clear();
				e_labels[count].Append(ip_label, ".", epslabels[i]);
				count++;
			}
		}		
	}
	
	if (e_codes[iIPLapStrain]) {
		const char* lap_epslabels1D[] = {"lap_eps11", "lap_eps11"};
		const char* lap_epslabels2D[] = {"lap_eps11", "lap_eps22", "lap_eps12"};
		const char* lap_epslabels2D_axi[] = {"lap_epsrr", "lap_epszz", "lap_epsrz", "lap_epstt"};
		const char* lap_epslabels3D[] = {"lap_eps11", "lap_eps22", "lap_eps33", "lap_eps23", "lap_eps13", "lap_eps12"};
		int nstrs = fB1.Rows();
		const char** lap_epslabels = NULL;
		if (nstrs == 1)
			lap_epslabels = lap_epslabels1D;
		else if (nstrs == 3)
			lap_epslabels = lap_epslabels2D;
		else if (nstrs == 4)
			lap_epslabels = lap_epslabels2D_axi;
		else if (nstrs == 6)
			lap_epslabels = lap_epslabels3D;
		else
			ExceptionT::GeneralFail(caller);

		/* over integration points */
		for (int j = 0; j < NumIP(); j++)
		{
			StringT ip_label;
			ip_label.Append("ip", j+1);
			
			/* over stress/strain components */
			for (int i = 0; i < 2*nstrs; i++) {
				e_labels[count].Clear();
				e_labels[count].Append(ip_label, ".", lap_epslabels[i]);
				count++;
			}
		}		
	}
	
	/* material output labels */
	if (e_codes[iIPMaterialData])
	{
		ArrayT<StringT> matlabels;
		(*fMFGPMatList)[0]->OutputLabels(matlabels);	

		/* over integration points */
		for (int j = 0; j < NumIP(); j++)
		{
			StringT ip_label;
			ip_label.Append("ip", j+1);
			
			/* over stress components */
			for (int i = 0; i < matlabels.Length(); i++)
			{
				e_labels[count].Clear();
				e_labels[count].Append(ip_label, ".", matlabels[i]);
				count++;
			}
		}		
	}
}

