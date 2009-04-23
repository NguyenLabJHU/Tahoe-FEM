/* $Id: FSFiber_OptNS.cpp,v 1.1 2009-04-23 03:03:42 thao Exp $ */

#include "FSFiber_OptNS.h"

#include "ShapeFunctionT.h"
#include "SolidMaterialT.h"

/* materials lists */
#include "FSFiberOptimize_MatListT.h"
#include "FSFiberOptimize_MatT.h"

/*material support*/
#include "FSFiberMatSupportT.h"

#include "ModelManagerT.h"
#include "LocalArrayT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "DomainIntegrationT.h"

using namespace Tahoe;

static const int flag = 1e+9;
static const double tol = 1.0e-9;

/* constructor */
FSFiber_OptNS::FSFiber_OptNS(const ElementSupportT& support):
	FSFiber_Optimize_Dual(support)
{
	SetName("uplag_fiber_opt_ns");
	
}

void FSFiber_OptNS::InitialCondition(void)
{
	const char caller[] = "FSFiber_OptNS::InitialCondition";

	/*inherited*/
	UpLagFiberCompT::InitialCondition();
	
	ModelManagerT& model = fFiberSupport->ModelManager();
	const int nnd = model.NumNodes();

	/*assumes all elements in group have the same parameters*/
	ContinuumMaterialT* pcont_mat = (*fMaterialList)[0];
	/* cast is safe since class contructs materials list */
	if(!fOptimize_Mat)
		ExceptionT::GeneralFail(caller, "Unable to cast material");

	const ArrayT<StringT>& labels = fOptimize_Mat->ParamLabels();	
	const int num_params = labels.Length();
	
	fParamLabels.Dimension(num_params);
	fParamLabels = labels;
	
	cout <<"\nSolves the dual problem, and computes the cost functional and gradients"
		 <<"\nfor optimization using DAKOTA's optpp_q_newton algorithm.  The order "
		 <<"\nthat the parameters are written to the .out file are: \n"<<endl;
	
	for (int i = 0; i < num_params; i++)
		cout << fParamLabels[i]<<endl;
	

	/*count the number of parameters and deal with block params*/	
	int nfparams = fField_Params.Length();

	/*get all of the element block IDs*/
	int nblocks = fMaterialList->Length();
	fblock_index.Dimension(nblocks);
	for(int i = 0; i < nblocks; i++)
	{
		fblock_index[i] = i;
	}
	
	fGradients.Dimension(num_params+ nfparams*(nblocks-1));
	fConstraint_Grad.Dimension(num_params);
	fParamVals.Dimension(num_params);

	
	fparam_map.Dimension(num_params);
	fparam_map = 0;
	
	for (int i = 0; i < fField_Params.Length(); i++)
	{
		for (int j = 0; j < num_params; j++)
		{
			bool same = (fField_Params[i] == fParamLabels[j]);
			if (same)
				fparam_map[j] = i+1;
		}
	}	
	
//	cout << "\nfparam_map: "<<fparam_map;
			
	/*initialize the cost function and gradients*/
	fCostFunction = 0.0;
	fconstraint = 0.0;
	fGradients = 0.0;
		
	foutput_step = 0;
	
	/*set up surface information*/
	fData.Dimension(nnd, NumDOF());
	fData_Coords.Dimension(nnd,NumDOF());
	fData = flag;	
	fData_Coords = 0.0;
}

void FSFiber_OptNS::InitStep(void)
{
	const char caller[] = "FSFiber_OptNS::InitStep";

	/*inherited*/
	UpLagFiberCompT::InitStep();

	fDataFile = fDataFileRoot;
//	cout << "\nfDataFileRoot: "<<fDataFileRoot;
	fDataFile.Append(".");
	fDataFile.Append(foutput_step,4);	
//	cout << "\nstep: "<<foutput_step;
	fDataInput.open(fDataFile);
	
	int numnodes, ndof;
	double time;
	fDataInput >> numnodes;
	fDataInput >> ndof;
	fDataInput >> time;
	foutput = (fabs(ElementSupport().Time()-time) < 1.0e-4*time);
	fNumNodes = numnodes;

	if (foutput)
	{
		for (int ind = 0; ind < ndof; ind++)
			fDataInput >> fweight_cost[ind];
//		cout << "\nfweight_cost: "<< fweight_cost;
		
		const ModelManagerT& model = fFiberSupport->ModelManager();
		if (ndof > NumDOF())
			ExceptionT::SizeMismatch(caller, "number of dof %d of datafile lager than expected by model.", ndof);

//		cout << "\nnumnodes: "<<numnodes;
		
		for (int n = 0 ; n < numnodes; n++)
		{
			int node;
			fDataInput >> node;
			node--;
//			cout << "\nnode: "<<node;
			for (int m = 0; m < ndof; m++)
			{
				double temp;
				fDataInput >> fData_Coords(node,m);
			}
			for (int m = 0; m < ndof; m++)
				fDataInput >> fData(node,m);
//			cout << "\ndatacoords: "<<fData_Coords(node, 0)<<"\t"<<fData_Coords(node,1)<<"\t"<<fData_Coords(node,2);
		}
	}
}

/* implementation of the ParameterInterfaceT interface */
void FSFiber_OptNS::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	UpLagFiberCompT::DefineParameters(list);

	list.SetDescription("Solves adjoint problem.  Must define FSFiber_Optimize_Primal first.");
	ParameterT primal_element_group(ParameterT::Integer, "primal_element_group");
	list.AddParameter(primal_element_group);

	ParameterT datafile(ParameterT::String, "data_input_file_root");
	datafile.SetDescription("Experimental data. Suffix .[XXX] will be added to filename, where XXX is the time step number.");
	list.AddParameter(datafile);
	
	ParameterT inputfile(ParameterT::String, "dakota_output_file");
	inputfile.SetDescription("Reports cost functional and gradients. Cost_Functional: sum_time(int_body( sum_nsd(U_i-U^m_i)^2 ) ) )");
	list.AddParameter(inputfile);

	
}

/* information about subordinate parameter lists */
void FSFiber_OptNS::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidElementT::DefineSubs(sub_list);	
	
	/* element block/material specification */
	sub_list.AddSub("fiber_opto_element_block", ParameterListT::OnePlus);
	sub_list.AddSub("fiber_orientations", ParameterListT::OnePlus);
	sub_list.AddSub("data_nodeset");
	sub_list.AddSub("field_params_ID_list", ParameterListT::Any);
}


/* return the description of the given inline subordinate parameter list */
ParameterInterfaceT* FSFiber_OptNS::NewSub(const StringT& name) const
{
	if (name == "fiber_opto_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("fiber_opto_material", ParameterListT::Once);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	else if(name == "data_nodeset")
	{
		ParameterContainerT* ns = new ParameterContainerT(name);
		ns->AddSub("nodeset_ID_list");
		return(ns);
	}
	else /* inherited */
		return UpLagFiberCompT::NewSub(name);
}

void FSFiber_OptNS::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FSFiber_OptNS::TakeParameterList";

	/* resolve primal problem element group before calling inherited */
	int primal_element_group = list.GetParameter("primal_element_group");
	primal_element_group--;
	ElementBaseT& element = ElementSupport().ElementGroup(primal_element_group);	
	fPrimal_Element = TB_DYNAMIC_CAST(UpLagFiberCompT*, &element);		

	if(!fPrimal_Element)
		ExceptionT::GeneralFail(caller, "Element group %d for primal problem not found", primal_element_group+1);


	/* inherited */
	UpLagFiberCompT::TakeParameterList(list);
		
	/* what's needed */
	bool need_F = false;
	bool need_F_last = false;
	for (int i = 0; i < fMaterialList->Length(); i++)
	{
		need_F = need_F || Needs_F(i);		
		need_F_last = need_F_last || Needs_F_last(i);
	}	

	/* allocate deformation gradient list */
	if (need_F)
	{
		fDualGrad_List.Dimension(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fDualGrad_List[i].Dimension(NumSD());
	}
	

	fip_residual.Dimension(NumSD());
	fip_data.Dimension(NumSD());

	fmat.Dimension(NumSD());
	fmat2.Dimension(NumSD());
	fGradU.Dimension(NumSD());
	fweight_cost.Dimension(NumSD());

	/*construct output file*/
	fOutFile = list.GetParameter("dakota_output_file");

	/*read data file*/
	fDataFileRoot = list.GetParameter ("data_input_file_root");
	foutput = false;
	
	/*read in for surface information and dimension work spaces*/
	const ParameterListT& surface_params = list.GetList("data_nodeset");
	StringListT::Extract(surface_params.GetList("nodeset_ID_list"),  fnodeset_IDs);
	int numlist = fnodeset_IDs.Length();
	fdata_nodeset.Dimension(numlist); 
	feqnos.Dimension(NumSD());
	fforce.Dimension(NumSD());
	
	int num_field_params = list.NumLists("field_params_ID_list");
	if (num_field_params)
	{
		const ParameterListT& field_params = list.GetList("field_params_ID_list");
		StringListT::Extract(field_params, fField_Params);
	}
	/*count the number of  facets and facet nodes*/
	int nsd = NumSD();
	double numnodes=0;
	
	/* nodes on*/
	ModelManagerT& model = ElementSupport().ModelManager();
	for (int i = 0; i< numlist; i++)
	{
		const StringT& ID = fnodeset_IDs[i];
		const iArrayT& set = model.NodeSet(ID);
		numnodes += set.Length();
	}
	fdata_nodeset.Dimension(numnodes);

	int dex = 0;
	for (int i=0; i<numlist; i++)
	{
		const StringT& ID = fnodeset_IDs[i];
		const iArrayT& set = model.NodeSet(ID);
		int nn = set.Length();
		for (int j = 0; j<nn; j++)
			fdata_nodeset[dex++] = set[j];
	}
	
}
	

/***********************************************************************
 * Protected
 ***********************************************************************/
/* calculate the cost function */
void FSFiber_OptNS::Compute_Cost(void)
{
	const char caller[] = "FSFiber_OptNS::Compute_Cost";

	int nsd = NumSD();
	int ndof = NumDOF();

	double tot = 0.0;
	/*get field data*/
	const FieldT& primal_field = fPrimal_Element->Field();
	const dArray2DT& displacement = primal_field[0];
	const dArray2DT& coordinates = ElementSupport().InitialCoordinates();

	for (int i = 0; i< fdata_nodeset.Length(); i++)
	{
		int node = fdata_nodeset[i];
//		cout << "\nnode: "<<node;
		for (int j = 0; j<NumSD() && ElementSupport().StepNumber() > 0; j++)
		{
			if(fabs(((fData_Coords(node,j) - coordinates(node,j))/fData_Coords(node,j)) > tol && fabs(fData_Coords(node,j))>tol) )
			{
//				cout <<"\nMismatch coords for node "<<node;
					ExceptionT::GeneralFail(caller);	
			}
			double diff = 0.0;
			if (fData(node,j) <flag)
			{
//				cout <<"\nNo data for node "<<node;
//					ExceptionT::GeneralFail(caller);	
				diff = displacement(node,j) - fData(node, j);
			}				
			tot += 0.5*diff*diff*fweight_cost[j];
//			cout << "\ndisp: "<<displacement(node,j);
//			cout << "\ndiff: "<<diff;
//			cout << "\ndata: "<<fData(node, j);
//			cout << "\ntot: "<<tot;
		}
	}
	fCostFunction += tot/fNumNodes;
}


/* calculate the cost function */
void FSFiber_OptNS::Compute_Gradients(void)
{
	const char caller[] = "FSFiber_Optimize_Dual::Compute_Gradients";
	int nsd = NumSD();
	
	Top();
	while (NextElement())
	{
		/* set shape function derivatives */
		SetGlobalShape();
	
		const double* det    = fCurrShapes->IPDets();
		const double* weight = fCurrShapes->IPWeights();

		fCurrShapes->TopIP();
		while ( fCurrShapes->NextIP() )
		{
			/*interpolate displacement and target_data to calculate residual at integration point*/		
			const dSymMatrixT& dual_grad = DualGradient();
//			cout << "\ndualgrad: "<<dual_grad;
			const dArray2DT& stress_grad = fOptimize_Mat->ds_ij_dlambda_q();
			double temp = (*det++)*(*weight++);
			double product=0;
			
			int nparams = fParamLabels.Length();
			int nblocks = fblock_index.Length();
			int nfp = fField_Params.Length();
			
			int blocknum = CurrentElement().MaterialNumber();
			
//			cout << "\n nparams: "<< nparams;
			for (int i = 0; i < nparams; i++)
			{				
				if (NumSD() == 2){
					product = dual_grad[0]*stress_grad(0,i) 
							+ dual_grad[1]*stress_grad(1,i)
							+ 2.0*dual_grad[2]*stress_grad(2,i);
				}
				else if (NumSD() == 3) {
					product = dual_grad[0]*stress_grad(0,i) 
							+ dual_grad[1]*stress_grad(1,i) 
							+ dual_grad[2]*stress_grad(2,i) 
							+ 2.0*dual_grad[3]*stress_grad(3,i) 
							+ 2.0*dual_grad[4]*stress_grad(4,i) 
							+ 2.0*dual_grad[5]*stress_grad(5,i);
				}
				else 
					ExceptionT::GeneralFail(caller, "nsd = %d not supported", NumSD());
				if (fparam_map[i] && fblock_index[blocknum] > 0)
				{
//					cout << "\ni: "<<i;
//					cout << "\nblocknum: "<<blocknum;
//					cout << "\nindex: "<<fblock_index[blocknum];
//					cout << "\nparam_map: "<<fparam_map[i];

					int dex =  nparams +  fblock_index[blocknum]*(fparam_map[i]-1);
					fGradients[dex] += (product)*temp;
				}
				else
					fGradients[i] += (product)*temp;

			} /*for nparams*/
		}/*while ip*/
	}/*while elem*/
}


void FSFiber_OptNS::RHSDriver(void)
{
	/* inherited */
	UpLagFiberCompT::RHSDriver();

	/* element contribution */
	if (foutput)
		ApplyNodalForce();
}

/* calculate the internal force contribution ("-k*d") */
void FSFiber_OptNS::FormKd(double constK)
{
	const char caller[] = "FSFiber_SurfOpt::FormKd";
	
	/* current element info */
	int nen = NumElementNodes();
	int nsd = NumSD();

	/* integration */
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* initialize */
	fRHS = 0.0;
	fCurrShapes->TopIP();
	while ( fCurrShapes->NextIP() && foutput)
	{
		double scale = constK*(*Det++)*(*Weight++);

		/* strain displacement matrix */
		Set_B(fCurrShapes->Derivatives_U(), fB);

	/* MATERIAL STIFFNESS */
		/* get D matrix */
		const dMatrixT& modulus = fCurrMaterial->c_ijkl();

		const dSymMatrixT& dual_strain = fDualGrad_List[CurrIP()];
		fmat.A_ijkl_B_kl(modulus, dual_strain);
		fB.MultTx(fmat, fNEEvec);
		fRHS.AddScaled(scale, fNEEvec);		


	/* GEOMETRIC STIFFNESS */

		/* B^T * Cauchy stress */
		const dSymMatrixT& cauchy = fCurrMaterial->s_ij();

		cauchy.ToMatrix(fCauchyStress);
		
		/* integration constants */		
		fCauchyStress *= scale;

		/*gradient of W*/
		fCurrShapes->GradU(fLocDisp, fGradU);
		
		fmat2.MultAB(fGradU, fCauchyStress);
		
		/* get shape function gradients matrix */
		fCurrShapes->GradNa(fGradNa);
	
		for (int a = 0; a < nen; a++)
		{
			for (int i = 0; i < nsd; i++)			
			{
				double p = a*nsd+i;
				fRHS[p] += fGradNa(0,a)*fmat2(i,0)+fGradNa(1,a)*fmat2(i,1)+fGradNa(2,a)*fmat2(i,2);						 
			}
		}
	}/*while*/
//	cout << "\nfRHS: "<<fRHS;	
}

/* compute contribution to RHS from traction BC's */
void FSFiber_OptNS::ApplyNodalForce(void)
{
	const char caller[] = "FSFiber_OptNS::ApplyNodalForce";

	int nsd = NumSD();
	int ndof = NumDOF();

	/*collect equation nos*/
	const iArray2DT& eqnos = Field().Equations();

	/*get field data*/
	const FieldT& primal_field = fPrimal_Element->Field();
	const dArray2DT& displacement = primal_field[0];
	const dArray2DT& coordinates = ElementSupport().InitialCoordinates();
	fforce = 0.0;
	for (int i = 0; i< fdata_nodeset.Length(); i++)
	{
		int node = fdata_nodeset[i];
//			cout << "\nData_Coords: "<<fData_Coords(node,0)<<"\t"<<fData_Coords(node,1)<<"\t"<<fData_Coords(node,2);
//			cout << "\ncoords: "<<coordinates(node,0)<<"\t"<<coordinates(node,1)<<"\t"<<coordinates(node,2);
		for (int j = 0; j<NumSD(); j++)
		{
			if(fabs((fData_Coords(node,j) - coordinates(node,j))/fData_Coords(node,j)) > tol && fabs(fData_Coords(node,j))>tol)
			{
				cout <<"\nMismatch coords for node "<<node+1;
					ExceptionT::GeneralFail(caller);	
			}
			double diff = 0.0;
			if (fData(node,j) < flag)
			{
//				cout <<"\nNo data for node "<<node;
//					ExceptionT::GeneralFail(caller);	
				diff = displacement(node,j) - fData(node, j);				
			}
			
//			cout << "\nnode: "<<node;
//			cout << "\ndiff: "<<diff;
			fforce[j] = -diff/fNumNodes*fweight_cost[j];
		}
		eqnos.RowAlias(node, feqnos);
		ElementSupport().AssembleRHS(Group(), fforce, feqnos);
	}
}


/* compute the measures of strain/deformation over the element */
/* called before internal force vector and stiffness matrix eval*/
void FSFiber_OptNS::SetGlobalShape(void)
{
	const char caller[] = "FSFiber_Optimize_Dual::SetGlobalShape";

	/* inherited */
	SolidElementT::SetGlobalShape();

	/* shape function wrt current config */
	SetLocalX(fLocCurrCoords);
	fCurrShapes->SetDerivatives();

	/* what needs to get computed */
	int material_number = CurrentElement().MaterialNumber();
	bool needs_F = Needs_F(material_number);
	bool needs_F_last = Needs_F_last(material_number);
	
	/* material dependent local arrays.  Obtain element values (in local node ordering) from global fields*/
	if (needs_F) {
		SetLocalU(fLocPrimalDisp);	
		SetLocalU(fLocDisp);
	}
	if (needs_F_last) SetLocalU(fLocPrimalDisp_last);	

	/* loop over integration points */
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* deformation gradient */
		if (needs_F)
		{
			dMatrixT& mat = fF_List[ip];
			
			/* displacement gradient */
			fShapes->GradU(fLocPrimalDisp, mat, ip);
			mat.PlusIdentity();

			/* calculate dual strains */
			fCurrShapes->GradU(fLocDisp, fGradU, ip);
			fDualGrad_List[ip].Symmetrize(fGradU);
		}

		/* "last" deformation gradient */
		if (needs_F_last)
		{
			dMatrixT& mat = fF_last_List[ip];

			/* displacement gradient */
			fShapes->GradU(fLocPrimalDisp_last, mat, ip);

			/* add identity */
			mat.PlusIdentity();

		}
	}
	
}
