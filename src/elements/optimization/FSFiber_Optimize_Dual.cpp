/* $Id: FSFiber_Optimize_Dual.cpp,v 1.1 2009-04-23 03:03:43 thao Exp $ */

#include "FSFiber_Optimize_Dual.h"
#include "ShapeFunctionT.h"
#include "SolidMaterialT.h"

/* materials lists */
#include "FSFiberOptimize_MatListT.h"
#include "FSFiberOptimize_MatT.h"

/*material support*/
#include "FSFiberMatSupportT.h"

#include "ParameterContainerT.h"
#include "ModelManagerT.h"

using namespace Tahoe;

/* constructor */
FSFiber_Optimize_Dual::FSFiber_Optimize_Dual(const ElementSupportT& support):
	UpLagFiberCompT(support),
	fPrimal_Element(NULL),
	fOptimize_Mat(NULL),
	fLocData(LocalArrayT::kUnspecified),
	fLocPrimalDisp(LocalArrayT::kDisp),
	fLocPrimalDisp_last(LocalArrayT::kLastDisp)
{
	SetName("uplag_fiber_optimize_dual");
	
}

/* destructor */
FSFiber_Optimize_Dual::~FSFiber_Optimize_Dual(void)
{
}

void FSFiber_Optimize_Dual::WriteOutput(void)
{
	UpLagFiberCompT::WriteOutput();
	
	if (foutput)
	{
		fOutput.open(fOutFile); 

		Write_Dakota_Output(fOutput);	

		fOutput.close();

		foutput_step++;
	}
}

void FSFiber_Optimize_Dual::InitialCondition(void)
{
	const char caller[] = "FSFiber_Optimize_Dual::InitialCondition";

	/*inherited*/
	UpLagFiberCompT::InitialCondition();
	
	const ModelManagerT& model = fFiberSupport->ModelManager();
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
	
	fGradients.Dimension(num_params);
	fConstraint_Grad.Dimension(num_params);
	fParamVals.Dimension(num_params);

	fData.Dimension(nnd, NumDOF());
	fData_Coords.Dimension(nnd,NumDOF());
	fData = 0.0;	
	
	/*initialize the cost function and gradients*/
	fCostFunction = 0.0;
	fconstraint = 0.0;
	fGradients = 0.0;
	
	
	foutput_step = 0;
}

void FSFiber_Optimize_Dual::InitStep(void)
{
	const char caller[] = "FSFiber_Optimize_Dual::InitStep";

	/*inherited*/
//	cout << caller << endl;
	
	UpLagFiberCompT::InitStep();
	
//	foutput = ElementSupport().WriteOutput();
	/*only read in data for time >0*/
	fDataFile = fDataFileRoot;
	fDataFile.Append(".");
	fDataFile.Append(foutput_step,4);
	
	fDataInput.open(fDataFile);
	
	int nnd, ndof;
	double time;
	
	fDataInput >> nnd;
	fDataInput >> ndof;
	fDataInput >> time;
	foutput = (fabs(ElementSupport().Time()-time) < 1.0e-4*time);
		
	if (foutput)
	{
		for (int ind = 0; ind < ndof; ind++) 
			fDataInput >> fweight_cost[ind];
		fData = -1;
	//		cout << "\nfoutput_step: "<<foutput_step;

/*		cout << "\nstep number: "<<ElementSupport().StepNumber();
		cout << "\nRoot: "<<fDataFileRoot;
		cout << "\nfDataInput: "<<fDataFile;
	   cout << "\ntime: "<<time;
	   cout << "\nanalysis time: "<< ElementSupport().Time();
*/
		const ModelManagerT& model = fFiberSupport->ModelManager();
		/*check time stamp*/
//		if (fabs(ElementSupport().Time()-time) > 1.0e-4*time)
//			ExceptionT::GeneralFail(caller, "time stamp of data %f and simulation %f do not match", time, ElementSupport().Time());

		if (nnd > model.NumNodes() || ndof > NumDOF())
			ExceptionT::SizeMismatch(caller, "number of node %d and dof %d of datafile lager than expected by model.", nnd, ndof);

		for (int n = 0 ; n < nnd; n++)
		{
			int node;
			fDataInput >> node;
			
			for (int m = 0; m < ndof; m++)
				fDataInput >> fData_Coords(n,m);
			for (int m = 0; m < ndof; m++)
				fDataInput >> fData(n,m);
		}
	}
}

void FSFiber_Optimize_Dual::CloseStep(void)
{
	UpLagFiberCompT::CloseStep();
}
/* implementation of the ParameterInterfaceT interface */
void FSFiber_Optimize_Dual::DefineParameters(ParameterListT& list) const
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
void FSFiber_Optimize_Dual::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidElementT::DefineSubs(sub_list);	
	
	/* element block/material specification */
	sub_list.AddSub("fiber_opto_element_block", ParameterListT::OnePlus);
	sub_list.AddSub("fiber_orientations", ParameterListT::OnePlus);
}


/* return the description of the given inline subordinate parameter list */
ParameterInterfaceT* FSFiber_Optimize_Dual::NewSub(const StringT& name) const
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
	else /* inherited */
		return UpLagFiberCompT::NewSub(name);
}

void FSFiber_Optimize_Dual::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FSFiber_Optimize_Dual::TakeParameterList";

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

	/*construct output file*/
	fOutFile = list.GetParameter("dakota_output_file");
	fDataFileRoot = list.GetParameter ("data_input_file_root");
	
	fmat.Dimension(NumSD());
	fmat2.Dimension(NumSD());
	fGradU.Dimension(NumSD());
	fweight_cost.Dimension(NumSD());
	foutput = false;
}
	
/* extract the list of material parameters */
void FSFiber_Optimize_Dual::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
	const char caller[] = "FSFiber_Optimize_Dual::CollectMaterialInfo";
	
	/* initialize */
	mat_params.Clear();
	mat_params.SetName("fiber_opto_material");
	
	/* collected material parameters */
	int num_blocks = all_params.NumLists("fiber_opto_element_block");
	for (int i = 0; i < num_blocks; i++) {

		/* block information */	
		const ParameterListT& block = all_params.GetList("fiber_opto_element_block", i);
		
		/* collect material parameters */
		const ParameterListT& mat_list = block.GetList(mat_params.Name());
		const ArrayT<ParameterListT>& mat = mat_list.Lists();
		mat_params.AddList(mat[0]);
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void FSFiber_Optimize_Dual::Write_Dakota_Output(ofstreamT& output)
{
	/*compute the objective function and gradients*/
	
	Compute_Cost();
	
	Compute_Gradients();
	
	output << setprecision(12)<< fCostFunction<<endl;

	cout << "\n"<<fDataFile<<"\tcost: "<<fCostFunction;
	output << "[ ";
	
	for (int i = 0; i< fGradients.Length(); i++)
			output << fGradients[i]<< " ";
		output << "]"<<endl;

}


/* return a pointer to a new material list */
MaterialListT* FSFiber_Optimize_Dual::NewMaterialList(const StringT& name, int size)
{
	/* resolve number of spatial dimensions */
	/* no match */
	if (name != "fiber_opto_material")
		return NULL;
	/* full list */
	if (size > 0)
	{
		/* material support */
		if (!fFiberSupport) {
			fFiberSupport = TB_DYNAMIC_CAST(FSFiberMatSupportT*, NewMaterialSupport());
			if (!fFiberSupport)
				ExceptionT::GeneralFail("UpLagFiberCompT::NewMaterialList");
		}
		
		fFiberSupport->SetElementCards(&fElementCards);
		
		return new FSFiberOptimize_MatListT(size, *fFiberSupport);
	}
	else
	{
		return new FSFiberOptimize_MatListT();
	}
	
	/* no match */
	return NULL;
}

/***********************************************************************
 * Protected
 ***********************************************************************/
/* calculate the cost function */
void FSFiber_Optimize_Dual::Compute_Cost(void)
{
	int nsd = NumSD();
	double tot = 0.0;
	/*accumulate the next time step*/
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
			fCurrShapes->InterpolateU(fLocPrimalDisp, fip_residual);
			fCurrShapes->InterpolateU(fLocData, fip_data);
			fip_residual -= fip_data;
			
			const dMatrixT& F = fF_List[CurrIP()];
			double J = F.Det();
			
			double temp = (*det++)*(*weight++)/J;
			for (int i = 0; i < nsd; i++)
				tot += 0.5*fip_residual[i]*fip_residual[i]*temp*fweight_cost[i];
		}/*while ip*/
	}/*while elem*/
	fCostFunction += tot;
}

/* calculate the cost function */
void FSFiber_Optimize_Dual::Compute_Gradients(void)
{
	const char caller[] = "FSFiber_Optimize_Dual::Compute_Gradients";
	int nsd = NumSD();
	
	/*accumulate the next time step*/
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
			
//			cout << "\nstress grad: "<<stress_grad;
//			cout << "\ndual_grad: "<<dual_grad;
			for (int i = 0; i < NumParams(); i++)
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
				fGradients[i] += (product)*temp;
			} /*for nparams*/
//			cout << "\ngradients: "<<fGradients;
		}/*while ip*/
	}/*while elem*/
}

/* calculate the internal force contribution ("-k*d") */
void FSFiber_Optimize_Dual::FormKd(double constK)
{
	const char caller[] = "FSFiber_Optimize_Dual::FormKd";
	
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
//		cout << "\ndual_strain: "<<dual_strain;
//		cout << "\nmodulus: "<<modulus;
		fmat.A_ijkl_B_kl(modulus, dual_strain);
//		cout << "\nfmat: "<<fmat;
		fB.MultTx(fmat, fNEEvec);
		fRHS.AddScaled(scale, fNEEvec);		
//		cout << "\nfRHS_dual: "<<fRHS;
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
	}
//	cout << "\nRHS: "<<fRHS;
		/* integration */
	const double* Det2    = fCurrShapes->IPDets();
	const double* Weight2 = fCurrShapes->IPWeights();

	/* initialize */
	fCurrShapes->TopIP();
	while ( fCurrShapes->NextIP() && foutput)
	{
		double scale = constK*(*Det2++)*(*Weight2++);

	/*COST FUNCTION*/
		/*interpolate displacement and target_data to calculate residual at integration point*/		
		fCurrShapes->InterpolateU(fLocPrimalDisp, fip_residual);
		fCurrShapes->InterpolateU(fLocData, fip_data);
		fip_residual -= fip_data;
		/* accumulate in element residual force vector */				
		const double* Na = fCurrShapes->IPShapeU();

		const dMatrixT& F = fF_List[CurrIP()];
		double J = F.Det();
		/*integration factor*/
		for (int a = 0; a < nen; a++)
		{
			for (int j = 0; j < nsd; j++)			
			{
				double q = a*nsd + j;
				fRHS[q] += scale*1.0/J*Na[a]*fip_residual[j]*fweight_cost[j];
			}
		}
	}/*while*/
		
}


/* initialize local arrays */
void FSFiber_Optimize_Dual::SetLocalArrays(void)
{
	const char caller[] = "FSFiber_Optimize_Dual::SetLocalArrays";

	/* inherited */
	UpLagFiberCompT::SetLocalArrays();

	/* allocate */
	int nen = NumElementNodes();
	/* look for adjoint field*/
	const FieldT* primal = ElementSupport().Field("displacement");
	if (primal) {
	
		fLocPrimalDisp.Dimension(nen, NumDOF());
		fLocPrimalDisp_last.Dimension(nen, NumDOF());
		
		/* register */
		primal->RegisterLocal(fLocPrimalDisp);
		primal->RegisterLocal(fLocPrimalDisp_last);
	}
	else {
		cout << "\nMust define a field for the primal (elasticity) problem. ";
		ExceptionT::GeneralFail(caller);
	}
	fLocData.Dimension(nen,NumDOF());
}

/* compute the measures of strain/deformation over the element */
/* called before internal force vector and stiffness matrix eval*/
void FSFiber_Optimize_Dual::SetGlobalShape(void)
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
			
//			cout << "\nW: "<<fLocDisp;
//			cout<<"\ngradW: "<<fDualGrad_List[ip];
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

	const iArrayT& nodes = CurrentElement().NodesX();
	for (int n = 0; n < nodes.Length(); n++)
	{
		for (int m = 0; m < NumDOF()  && foutput; m++)
		{
			double a = nodes[n];
			if (fabs(fData_Coords(a,m) - fLocInitCoords(n,m))>kSmall)
			{
				cout << "\nMismatch coords for elem: "<<CurrElementNumber()
					<< fData_Coords(a,m)<<fLocInitCoords(n,m)<<endl;
				ExceptionT::GeneralFail(caller);
			}
			fLocData(n, m) = fData(nodes[n], m);
		}
		
	}
	
}

/* current element operations */
bool FSFiber_Optimize_Dual::NextElement(void)
{
	const char caller[] = "FSFiber_Optimize_Dual::NextElement";
	/* inherited */
	bool result = SolidElementT::NextElement();
	
	/* get material pointer */
	if (result)
	{
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[CurrentElement().MaterialNumber()];
	
		/* cast is safe since class contructs materials list */
		fOptimize_Mat =  TB_DYNAMIC_CAST(FSFiberOptimize_MatT*, pcont_mat);
		
		if(!fOptimize_Mat)
			ExceptionT::GeneralFail(caller, "Unable to cast material");
	}
	
	return result;
}


