/* $Id: FSFiber_OptSurf.cpp,v 1.1 2009-04-23 03:03:43 thao Exp $ */

#include "FSFiber_OptSurf.h"

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
FSFiber_OptSurf::FSFiber_OptSurf(const ElementSupportT& support):
	FSFiber_Optimize_Dual(support),
	fLocSurfDisp(LocalArrayT::kDisp),
	fLocSurfCoords(LocalArrayT::kInitCoords)		
{
	SetName("uplag_fiber_optsurf");
	
}

void FSFiber_OptSurf::InitialCondition(void)
{
	const char caller[] = "FSFiber_OptSurf::InitialCondition";

	/*inherited*/
	UpLagFiberCompT::InitialCondition();
	
	ModelManagerT& model = fFiberSupport->ModelManager();
	const int nnd = model.NumNodes();

	/*assumes all elements in group have the same parameters*/
	ContinuumMaterialT* pcont_mat = (*fMaterialList)[0];
	/* cast is safe since class contructs materials list */
//	fOptimize_Mat =  TB_DYNAMIC_CAST(SSOptimize_MatT*, pcont_mat);
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
	ElementBlockIDs(fblock_IDs);
	int nblocks = fblock_IDs.Length();
	int maxID = -1;
	for (int i = 0; i< nblocks; i++)
	{
		int temp = atoi(fblock_IDs[i]);
		if (temp > maxID)
			maxID = temp;
	}
	fblock_index.Dimension(maxID);
	fblock_index = -1;
	
	for(int i = 0; i < nblocks; i++)
	{
		int temp = atoi(fblock_IDs[i])-1;
		fblock_index[temp] = i;
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
	
	int numlist = fside_set_IDs.Length();
	int nsd = NumSD();
	
	iArray2DT nd_tmp, eq_tmp;

	for (int i = 0; i< numlist; i++)
	{
		const StringT& side_ID = fside_set_IDs[i];
		iArray2DT& side_set = fdata_sides[i];
		iArray2DT& facets_global = fdata_global_nodes[i];    /* number facets x num_facet nodes (global numbering) */
		iArray2DT& facets_local = fdata_local_nodes[i];     /* number facets x num_facet nodes (global numbering) */
		iArray2DT& facets_eqnos = fdata_eqnos[i];			/* number of facets x num facet nodes*nsd */
		
		int numfacets = side_set.MajorDim();
		/*get element geometry*/
		StringT elemID;   /*element group ID*/
		elemID = model.SideSetGroupID(side_ID);

		/*get the geometry of the element */
		GeometryT::CodeT geometry_code = model.ElementGroupGeometry(elemID);
		GeometryBaseT* geometry = GeometryT::New(geometry_code, NumElementNodes());
		
		/*get connectivities*/
		const iArray2DT& connectivities = model.ElementGroup(elemID);

		iArrayT local_nodes(fNumFacetNodes);  /*for each facet*/
		iArrayT global_nodes(fNumFacetNodes); /*for each facet*/
		iArrayT equation(fNumFacetNodes); /*for each facet*/
		nd_tmp.Set(1,fNumFacetNodes, global_nodes.Pointer());
		for (int j = 0; j < numfacets; j++)
		{
			int elem = side_set(j,0);
			int side = side_set(j,1);
			/* gets nodes on faces in local numbering*/
			geometry->NodesOnFacet(side, local_nodes);
			facets_local.SetRow(j,local_nodes);

			/*gets nodes on faces in global numbering*/
			global_nodes.Collect(local_nodes, connectivities(elem));
			facets_global.SetRow(j,global_nodes);

			eq_tmp.Set(1,nsd*fNumFacetNodes, facets_eqnos(j));
			Field().SetLocalEqnos(nd_tmp, eq_tmp);
		}
/*
		cout << "\nsideset: "<<side_set;
		cout << "\nlocal nodes: "<<facets_local;
		cout << "\nglobal nodes: "<<facets_global;
*/
	}
	
}

void FSFiber_OptSurf::InitStep(void)
{
	const char caller[] = "FSFiber_OptSurf::InitStep";

	/*inherited*/
//	cout << caller << endl;
	
	UpLagFiberCompT::InitStep();
	
	foutput = ElementSupport().WriteOutput();

	/*only read in data for time >0*/
	fDataFile = fDataFileRoot;
	if (foutput)
	{
//		cout << "\nfoutput_step: "<<foutput_step;
		fDataFile.Append(".");
		fDataFile.Append(foutput_step,4);
	
		fDataInput.open(fDataFile);
	
		int numnodes, ndof;
		double time;
	
		fDataInput >> numnodes;
		fDataInput >> ndof;
		fDataInput >> time;
		for (int ind = 0; ind < ndof; ind++)
			fDataInput >> fweight_cost[ind];

/*		cout << "\nstep number: "<<ElementSupport().StepNumber();
		cout << "\nRoot: "<<fDataFileRoot;
		cout << "\nfDataInput: "<<fDataFile;
	   cout << "\ntime: "<<time;
	   cout << "\nanalysis time: "<< ElementSupport().Time();
*/
		const ModelManagerT& model = fFiberSupport->ModelManager();
		/*check time stamp*/
		if (fabs(ElementSupport().Time()-time) > 1.0e-4*time)
			ExceptionT::GeneralFail(caller, "time stamp of data %f and simulation %f do not match", time, ElementSupport().Time());

		if (ndof > NumDOF())
			ExceptionT::SizeMismatch(caller, "number of dof %d of datafile lager than expected by model.", ndof);

		for (int n = 0 ; n < numnodes; n++)
		{
			int node;
			fDataInput >> node;
			node--;
			for (int m = 0; m < ndof; m++)
				fDataInput >> fData_Coords(node,m);
			for (int m = 0; m < ndof; m++)
				fDataInput >> fData(node,m);
		}
	}
}

/* implementation of the ParameterInterfaceT interface */
void FSFiber_OptSurf::DefineParameters(ParameterListT& list) const
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
void FSFiber_OptSurf::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidElementT::DefineSubs(sub_list);	
	
	/* element block/material specification */
	sub_list.AddSub("fiber_opto_element_block", ParameterListT::OnePlus);
	sub_list.AddSub("fiber_orientations", ParameterListT::OnePlus);
	sub_list.AddSub("data_surface");
	sub_list.AddSub("field_params_ID_list", ParameterListT::Any);
}


/* return the description of the given inline subordinate parameter list */
ParameterInterfaceT* FSFiber_OptSurf::NewSub(const StringT& name) const
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
	else if(name == "data_surface")
	{
		ParameterContainerT* surf = new ParameterContainerT(name);
		surf->AddSub("side_set_ID_list");
		return(surf);
	}
	else /* inherited */
		return UpLagFiberCompT::NewSub(name);
}

void FSFiber_OptSurf::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FSFiber_OptSurf::TakeParameterList";

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

	/*construct output file*/
	fOutFile = list.GetParameter("dakota_output_file");

	/*read data file*/
	fDataFileRoot = list.GetParameter ("data_input_file_root");
	
	
	/*read in for surface information and dimension work spaces*/
	const ParameterListT& surface_params = list.GetList("data_surface");
	StringListT::Extract(surface_params.GetList("side_set_ID_list"),  fside_set_IDs);
	int numlist = fside_set_IDs.Length();
	fdata_sides.Dimension(numlist); 
	fdata_global_nodes.Dimension(numlist);
	fdata_local_nodes.Dimension(numlist);
	fdata_eqnos.Dimension(numlist);
	
	int num_field_params = list.NumLists("field_params_ID_list");
	if (num_field_params)
	{
		const ParameterListT& field_params = list.GetList("field_params_ID_list");
		StringListT::Extract(field_params, fField_Params);
	}
	/*count the number of  facets and facet nodes*/
	int nsd = NumSD();
	int countnodes=0;
	int countfacets=0;
	
	/* nodes on element facets */
	iArrayT num_facet_nodes;
	fShapes->NumNodesOnFacets(num_facet_nodes);
	
	ModelManagerT& model = ElementSupport().ModelManager();
	for (int i = 0; i< numlist; i++)
	{
		const StringT& side_ID = fside_set_IDs[i];
		
		/*get element geometry*/
		StringT elemID;   /*element group ID*/
		elemID = model.SideSetGroupID(side_ID);

		/*read and store sideset*/
		iArray2DT& side_set = fdata_sides[i];
		side_set = model.SideSet(side_ID);
		int num_sides = side_set.MajorDim();
		
		/* all facets in set must have the same number of nodes */
		if(i==0)
			fNumFacetNodes = num_facet_nodes[side_set(0,1)];
		for (int f = 0; f < num_sides; f++)
			if (num_facet_nodes[side_set(f,1)] != fNumFacetNodes)
				ExceptionT::BadInputValue(caller, "faces side set \"%s\" have different numbers of nodes",
					side_ID.Pointer());

		
		/*count number of nodes and facets*/
		fdata_global_nodes[i].Dimension(num_sides, fNumFacetNodes);
		fdata_local_nodes[i].Dimension(num_sides, fNumFacetNodes);
		fdata_eqnos[i].Dimension(num_sides, fNumFacetNodes*nsd);
		
		countfacets += num_sides;	
	}
	fNumFacets = countfacets;

	fsurf_nodes.Dimension(fNumFacetNodes);
	feqnos.Dimension(fNumFacetNodes*nsd);
	fQ.Dimension(nsd);
	fnormal.Dimension(nsd);
	fjacobian.Dimension(nsd, nsd-1);

	/* allocate */
	/* look for adjoint field*/
	const FieldT* primal = ElementSupport().Field("displacement");
	if (primal) {
					
		fLocSurfDisp.Dimension(fNumFacetNodes,NumDOF());
		fLocSurfCoords.Dimension(fNumFacetNodes,NumDOF());

		/* register */
		primal->RegisterLocal(fLocSurfDisp);
	}
	else {
		cout << "\nMust define a field for the primal (elasticity) problem. ";
		ExceptionT::GeneralFail(caller);
	}
	fLocData.Dimension(fNumFacetNodes,NumDOF());
	ElementSupport().RegisterCoordinates(fLocSurfCoords);
	
	fsurf_rhs.Dimension(fNumFacetNodes*NumDOF());
	fweight_cost.Dimension(NumSD());
}
	

/***********************************************************************
 * Protected
 ***********************************************************************/
/* calculate the cost function */
void FSFiber_OptSurf::Compute_Cost(void)
{
	const char caller[] = "FSFiber_SurfOpt::Compute_Cost";

	int numlist = fside_set_IDs.Length();
	int nsd = NumSD();
	int ndof = NumDOF();

	double tot = 0.0;
	for (int i = 0; i< numlist; i++)
	{
		const StringT& side_ID = fside_set_IDs[i];
		const iArray2DT& side_set = fdata_sides[i];
		const iArray2DT& global_nodes = fdata_global_nodes[i];    /* number facets x num_facet nodes (global numbering) */
		const iArray2DT& local_nodes = fdata_local_nodes[i];     /* number facets x num_facet nodes (global numbering) */
		
		int numfacets = side_set.MajorDim();
		for (int j = 0; j < numfacets; j++)
		{
			/*construct integration domain*/
			int elem = side_set(j,0);

			int facet = side_set(j,1);
			const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(facet);

			/*collect displacements in local ordering*/
			global_nodes.RowCopy(j,fsurf_nodes);
			fLocSurfDisp.SetLocal(fsurf_nodes);
			fLocSurfCoords.SetLocal(fsurf_nodes);

			/*collect data*/
			for (int k = 0; k < fNumFacetNodes; k++)
			{
				int a = fsurf_nodes[k];
				for(int m = 0; m < ndof && ElementSupport().StepNumber() > 0; m++)
				{
//					cout << "\ndiff: "<<fabs(fData_Coords(a,m)-fLocSurfCoords(k,m));
					if(fabs((fData_Coords(a,m) - fLocSurfCoords(k,m))/fData_Coords(a,m)) > tol)
					{
						cout <<"\nMismatch coords for facet "<<facet<< " of element "<<elem<<endl;
						ExceptionT::GeneralFail(caller);	
					}

					if(fData(a,m) == flag)
						ExceptionT::GeneralFail(caller, "No data provided for node %",a);
					fLocData(k,m) = fData(a,m);
				}
			}
						
			int facet_nip = surf_shape.NumIP();
			const double* w = surf_shape.Weight(); /*integration weights*/
			
//			fRHS = 0.0;
			for (int k = 0; k < facet_nip; k++)
			{
				/*calculate surface jacobian*/
				surf_shape.DomainJacobian(fLocSurfCoords, k, fjacobian);

				/*integration weight*/
				double detj = surf_shape.SurfaceJacobian(fjacobian,fQ);
				double scale = detj*w[k];

				surf_shape.Interpolate(fLocSurfDisp, fip_residual, k);
				surf_shape.Interpolate(fLocData, fip_data, k);
				fip_residual -= fip_data;
				
				for (int j = 0; j < nsd; j++)			
					tot += 0.5*fip_residual[j]*fip_residual[j]*scale*fweight_cost[j];

			} /*for nip*/

		} /*for nsides*/
	} /*for nlist*/
	
	fCostFunction += tot;
}


/* calculate the cost function */
void FSFiber_OptSurf::Compute_Gradients(void)
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
			int nblocks = fblock_IDs.Length();
			int nfp = fField_Params.Length();
			
			const StringT& blockID = ElementBlockID(CurrElementNumber());
			int blocknum = atoi(blockID)-1;
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
					int dex =  nparams +  fblock_index[blocknum]*(fparam_map[i]-1);
					fGradients[dex] += (product)*temp;
				}
				else
					fGradients[i] += (product)*temp;

			} /*for nparams*/
		}/*while ip*/
	}/*while elem*/
}


void FSFiber_OptSurf::RHSDriver(void)
{
	/* inherited */
	UpLagFiberCompT::RHSDriver();

	/* element contribution */
	ApplySurfaceForce();
}

/* calculate the internal force contribution ("-k*d") */
void FSFiber_OptSurf::FormKd(double constK)
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
	while ( fCurrShapes->NextIP() && ElementSupport().StepNumber() > 0)
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
	
}

/* compute contribution to RHS from traction BC's */
void FSFiber_OptSurf::ApplySurfaceForce(void)
{
	const char caller[] = "FSFiber_SurfOpt::ApplySurfaceForce";

	int numlist = fside_set_IDs.Length();
	int nsd = NumSD();
	int ndof = NumDOF();
	for (int i = 0; i< numlist; i++)
	{
		const StringT& side_ID = fside_set_IDs[i];
		const iArray2DT& side_set = fdata_sides[i];
		const iArray2DT& global_nodes = fdata_global_nodes[i];    /* number facets x num_facet nodes (global numbering) */
		const iArray2DT& local_nodes = fdata_local_nodes[i];     /* number facets x num_facet nodes (global numbering) */
		const iArray2DT& local_eqnos = fdata_eqnos[i];
		
		int numfacets = side_set.MajorDim();
		for (int j = 0; j < numfacets; j++)
		{
			/*construct integration domain*/
			int elem = side_set(j,0);

			int facet = side_set(j,1);
			const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(facet);

			/*collect displacements in local ordering*/
			local_eqnos.RowCopy(j, feqnos);
			global_nodes.RowCopy(j,fsurf_nodes);
			fLocSurfDisp.SetLocal(fsurf_nodes);
			fLocSurfCoords.SetLocal(fsurf_nodes);

			/*collect data*/
			for (int k = 0; k < fNumFacetNodes  && ElementSupport().StepNumber() > 0; k++)
			{
				int a = fsurf_nodes[k];
				for(int m = 0; m < ndof && ElementSupport().StepNumber() > 0; m++)
				{
					if(fabs((fData_Coords(a,m) - fLocSurfCoords(k,m))/fData_Coords(a,m)) > tol)
					{
						cout <<"\nMismatch coords for facet "<<facet<< " of element "<<elem<<endl;
						ExceptionT::GeneralFail(caller);	
					}

					if(fData(a,m) == flag)
						ExceptionT::GeneralFail(caller, "No data provided for node %",a);
					fLocData(k,m) = fData(a,m);
				}
			}
						
			int facet_nip = surf_shape.NumIP();
			const double* w = surf_shape.Weight(); /*integration weights*/
			
			fsurf_rhs = 0.0;
			for (int k = 0; k < facet_nip && ElementSupport().StepNumber() > 0; k++)
			{
				/*calculate surface jacobian*/
				surf_shape.DomainJacobian(fLocSurfCoords, k, fjacobian);

				/*integration weight*/
				double detj = surf_shape.SurfaceJacobian(fjacobian,fQ);
				double scale = -detj*w[k];

				surf_shape.Interpolate(fLocSurfDisp, fip_residual, k);
				surf_shape.Interpolate(fLocData, fip_data, k);
				fip_residual -= fip_data;
				
				const double* Na = surf_shape.Shape(k);
				for (int a = 0; a < fNumFacetNodes; a++)
				{
					for (int j = 0; j < nsd; j++)			
					{
						double q = a*nsd + j;
						fsurf_rhs[q] += scale*Na[a]*fip_residual[j]*fweight_cost[j];
					} /*for nsd*/
				} /*for nfn*/
			} /*for nip*/

			/* assemble */
			ElementSupport().AssembleRHS(Group(), fsurf_rhs, feqnos);

		} /*for nsides*/
	} /*for nlist*/
}


/* compute the measures of strain/deformation over the element */
/* called before internal force vector and stiffness matrix eval*/
void FSFiber_OptSurf::SetGlobalShape(void)
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
