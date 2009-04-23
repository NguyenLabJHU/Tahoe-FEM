/* $Id: SS_Optimize_Dual.cpp,v 1.1 2009-04-23 03:03:43 thao Exp $ */
#include "SS_Optimize_Dual.h"
#include "ShapeFunctionT.h"
//#include "SSOptimize_MatT.h"
#include "SSSolidMatT.h"

/* materials lists */
#include "SSOptimize_MatList3DT.h"
/*material support*/
#include "SSMatSupportT.h"

#include "ParameterContainerT.h"
#include "ModelManagerT.h"

using namespace Tahoe;

/* constructor */
SS_Optimize_Dual::SS_Optimize_Dual(const ElementSupportT& support):
	SmallStrainT(support),
	fPrimal_Element(NULL),
	fOptimize_Mat(NULL),
	fLocData(LocalArrayT::kUnspecified),
	fLocPrimalDisp(LocalArrayT::kDisp),
	fLocPrimalDisp_last(LocalArrayT::kLastDisp)
{
	SetName("small_strain_optimize_dual");
}

/* destructor */
SS_Optimize_Dual::~SS_Optimize_Dual(void)
{
}

void SS_Optimize_Dual::WriteOutput(void)
{
	SmallStrainT::WriteOutput();
	
	fOutput.open(fOutFile); 

	Write_Dakota_Output(fOutput);	

	fOutput.close();
}

void SS_Optimize_Dual::InitialCondition(void)
{
	const char caller[] = "SS_Optimize_Dual::InitialCondition";

	/*inherited*/
	SmallStrainT::InitialCondition();
	
	const ModelManagerT& model = fSSMatSupport->ModelManager();
	const int nnd = model.NumNodes();

	/*assumes all elements in group have the same parameters*/
	ContinuumMaterialT* pcont_mat = (*fMaterialList)[0];
	/* cast is safe since class contructs materials list */
	fOptimize_Mat =  TB_DYNAMIC_CAST(SSOptimize_MatT*, pcont_mat);
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
}

void SS_Optimize_Dual::InitStep(void)
{
	const char caller[] = "SS_Optimize_Dual::InitStep";

	/*inherited*/
//	cout << caller << endl;
	
	SmallStrainT::InitStep();
	
	/*only read in data for time >0*/
	fDataFile = fDataFileRoot;
	if (ElementSupport().StepNumber() > 0)
	{
		const int time_step = ElementSupport().StepNumber();
		fDataFile.Append(".");
		fDataFile.Append(time_step-1,4);
	
		fDataInput.open(fDataFile);
	
		int nnd, ndof;
		double time;
	
		fDataInput >> nnd;
		fDataInput >> ndof;
		fDataInput >> time;
		fDataInput >> fweight_cost;
		fData = -1;

/*		cout << "\nstep number: "<<ElementSupport().StepNumber();
		cout << "\nRoot: "<<fDataFileRoot;
		cout << "\nfDataInput: "<<fDataFile;
	   cout << "\ntime: "<<time;
	   cout << "\nanalysis time: "<< ElementSupport().Time();
*/
		const ModelManagerT& model = fSSMatSupport->ModelManager();
		/*check time stamp*/
		if (fabs(ElementSupport().Time()-time) > 1.0e-4*time)
			ExceptionT::GeneralFail(caller, "time stamp of data %f and simulation %f do not match", time, ElementSupport().Time());

		if (nnd > model.NumNodes() || ndof > NumDOF())
			ExceptionT::SizeMismatch(caller, "number of node %d and dof %d of datafile lager than expected by model.", nnd, ndof);

		for (int n = 0 ; n < nnd; n++)
		{
			for (int m = 0; m < ndof; m++)
				fDataInput >> fData_Coords(n,m);
			for (int m = 0; m < ndof; m++)
				fDataInput >> fData(n,m);
		}
	}
}

/* implementation of the ParameterInterfaceT interface */
void SS_Optimize_Dual::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SmallStrainT::DefineParameters(list);

	list.SetDescription("Solves adjoint problem.  Must define SS_Optimize_Primal first.");
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
void SS_Optimize_Dual::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidElementT::DefineSubs(sub_list);	

	/* element block/material specification */
	sub_list.AddSub("ss_opto_element_block", ParameterListT::OnePlus);
}


/* return the description of the given inline subordinate parameter list */
ParameterInterfaceT* SS_Optimize_Dual::NewSub(const StringT& name) const
{
	if (name == "ss_opto_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("ss_opto_material", ParameterListT::Once);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	else /* inherited */
		return SolidElementT::NewSub(name);
}

void SS_Optimize_Dual::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SS_Optimize_Dual::TakeParameterList";

	/* resolve primal problem element group before calling inherited */
	int primal_element_group = list.GetParameter("primal_element_group");
	primal_element_group--;
	ElementBaseT& element = ElementSupport().ElementGroup(primal_element_group);	
	fPrimal_Element = TB_DYNAMIC_CAST(SmallStrainT*, &element);		

	if(!fPrimal_Element)
		ExceptionT::GeneralFail(caller, "Element group %d for primal problem not found", primal_element_group+1);


	/* inherited */
	SmallStrainT::TakeParameterList(list);
		
	/* what's needed */
	bool need_strain = false;
	bool need_strain_last = false;
	for (int i = 0; i < fMaterialNeeds.Length(); i++) {
		const ArrayT<bool>& needs = fMaterialNeeds[i];
		need_strain = need_strain || needs[fNeedsOffset + kstrain];
		need_strain_last = need_strain_last || needs[fNeedsOffset + kstrain_last];
	}

	/* allocate strain list */
	if (need_strain) {
		fDual_Strain_List.Dimension(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fDual_Strain_List[i].Dimension(NumSD());
	}	

	fip_residual.Dimension(NumSD());
	fip_data.Dimension(NumSD());


	/*construct output file*/
	fOutFile = list.GetParameter("dakota_output_file");
	fDataFileRoot = list.GetParameter ("data_input_file_root");
	
	fmat.Dimension(NumSD());
}

/*const int SS_Optimize_Dual::NumParams(void) const
{
	const char caller[] = "SS_Optimize_Dual::NumParams";
	if(fPrimal_Element)
		return(fPrimal_Element->NumParams());
	else
		ExceptionT::GeneralFail(caller, "No element group found for primal problems specified");	
}
*/
	
/* extract the list of material parameters */
void SS_Optimize_Dual::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
	const char caller[] = "SS_Optimize_Dual::CollectMaterialInfo";
	
	/* initialize */
	mat_params.Clear();
	mat_params.SetName("ss_opto_material");
	
	/* collected material parameters */
	int num_blocks = all_params.NumLists("ss_opto_element_block");
	for (int i = 0; i < num_blocks; i++) {

		/* block information */	
		const ParameterListT& block = all_params.GetList("ss_opto_element_block", i);
		
		/* collect material parameters */
		const ParameterListT& mat_list = block.GetList(mat_params.Name());
		const ArrayT<ParameterListT>& mat = mat_list.Lists();
		mat_params.AddList(mat[0]);
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void SS_Optimize_Dual::Write_Dakota_Output(ofstreamT& output)
{
	/*compute the objective function and gradients*/
	Compute_Cost();
	
	Compute_Gradients();
	
	output << fCostFunction<<endl;

	output << "[ ";
	for (int i = 0; i< fGradients.Length(); i++)
		output << fGradients[i]<< " ";
	output << "]"<<endl;

	output.close();
}


/* return a pointer to a new material list */
MaterialListT* SS_Optimize_Dual::NewMaterialList(const StringT& name, int size)
{
	/* resolve number of spatial dimensions */
	/* no match */
	if (name != "ss_opto_material")
		return NULL;
	/* full list */
	if (size > 0)
	{
		/* material support */
		if (!fSSMatSupport) {
			fSSMatSupport = TB_DYNAMIC_CAST(SSMatSupportT*, NewMaterialSupport());
			if (!fSSMatSupport)
				ExceptionT::GeneralFail("SmallStrainT::NewMaterialList");
		}
		
		fSSMatSupport->SetElementCards(&fElementCards);
		
		return new SSOptimize_MatList3DT(size, *fSSMatSupport);
	}
	else
	{
		return new SSOptimize_MatList3DT();
	}
	
	/* no match */
	return NULL;
}

/***********************************************************************
 * Protected
 ***********************************************************************/
/* calculate the cost function */
void SS_Optimize_Dual::Compute_Cost(void)
{
	int nsd = NumSD();
	double tot = 0.0;
	/*accumulate the next time step*/
	Top();
	while (NextElement())
	{
		/* set shape function derivatives */
		SetGlobalShape();
	
		const double* det    = fShapes->IPDets();
		const double* weight = fShapes->IPWeights();

		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/*interpolate displacement and target_data to calculate residual at integration point*/		
			fShapes->InterpolateU(fLocPrimalDisp, fip_residual);
			fShapes->InterpolateU(fLocData, fip_data);
			fip_residual -= fip_data;
			double temp = (*det++)*(*weight++);
			for (int i = 0; i < nsd; i++)
				tot += 0.5*fip_residual[i]*fip_residual[i]*temp;
		}/*while ip*/
	}/*while elem*/
	fCostFunction += tot*fweight_cost;
}

/* calculate the cost function */
void SS_Optimize_Dual::Compute_Gradients(void)
{
	const char caller[] = "SS_Optimize_Dual::Compute_Gradients";
	int nsd = NumSD();
	
	/*accumulate the next time step*/
	Top();
	while (NextElement())
	{
		/* set shape function derivatives */
		SetGlobalShape();
	
		const double* det    = fShapes->IPDets();
		const double* weight = fShapes->IPWeights();

		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/*interpolate displacement and target_data to calculate residual at integration point*/		
			const dSymMatrixT& dual_strain = LinearDualStrain();
			
			const dArray2DT& stress_grad = fOptimize_Mat->ds_ij_dlambda_q();
//			cout << "\nstress grad: "<<stress_grad<<endl;;
//			cout << "\ndual strain: "<<dual_strain<<endl;
			int dex = 0;
			double temp = (*det++)*(*weight++);
			double product=0;
			for (int i = 0; i < NumParams(); i++)
			{
				dex++;
				
				if (NumSD() == 2){
					product = dual_strain[0]*stress_grad(0,i) 
							+ dual_strain[1]*stress_grad(1,i)
							+ 2.0*dual_strain[2]*stress_grad(2,i);
				}
				else if (NumSD() == 3) {
					product = dual_strain[0]*stress_grad(0,i) 
							+ dual_strain[1]*stress_grad(1,i) 
							+ dual_strain[2]*stress_grad(2,i) 
							+ 2.0*dual_strain[3]*stress_grad(3,i) 
							+ 2.0*dual_strain[4]*stress_grad(4,i) 
							+ 2.0*dual_strain[5]*stress_grad(5,i);
				}
				else 
					ExceptionT::GeneralFail(caller, "nsd = %d not supported", NumSD());
				
				fGradients[i] += (product)*temp*fweight_cost;
			} /*for nparams*/
		}/*while ip*/
	}/*while elem*/
}

/* form the element stiffness matrix */
void SS_Optimize_Dual::FormStiffness(double constK)
{
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integrate element stiffness */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/********DEBUG*******/
	bool print = false; 
	int pos = fElementCards.Position(); 
	if (pos == 1&&0)  
	  print = true; 
	/*******************/
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{

		double scale = constK*(*Det++)*(*Weight++);
	
		/* strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
		if (print) cout << "\nmodulus: "<<fCurrMaterial->c_ijkl();
							
		/* multiply b(transpose) * db, taking account of symmetry, */
		/* and accumulate in elstif */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
}

/* calculate the internal force contribution ("-k*d") */
void SS_Optimize_Dual::FormKd(double constK)
{
	const char caller[] = "SS_Optimize_Dual::FormKd";
//	cout << caller <<endl;
	
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	int nsd = NumSD();
	int  nen = NumElementNodes();
	fRHS = 0.0;

	/*don't form residual when calculating inital conditions*/
	fShapes->TopIP();
	while (fShapes->NextIP() && ElementSupport().StepNumber() > 0) 
	{
		/*calculate stress and history variables*/
		const dSymMatrixT& stress = fCurrMaterial->s_ij();
		
		double scale = constK*(*Det++)*(*Weight++);	
		/* strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);

		/* get D matrix */
		const dMatrixT& modulus = fCurrMaterial->c_ijkl();

		const dSymMatrixT& dual_strain = fDual_Strain_List[CurrIP()];
		fmat.A_ijkl_B_kl(modulus, dual_strain);
		fB.MultTx(fmat, fNEEvec);
		fRHS.AddScaled(scale, fNEEvec);		
	
		/*interpolate displacement and target_data to calculate residual at integration point*/		
		fShapes->InterpolateU(fLocPrimalDisp, fip_residual);
		fShapes->InterpolateU(fLocData, fip_data);
		fip_residual -= fip_data;
		/* accumulate in element residual force vector */				
		const double* Na = fShapes->IPShapeU();

		/*integration factor*/
		for (int a = 0; a < nen; a++)
		{
			for (int j = 0; j < nsd; j++)			
			{
				int p = a*nsd + j;
				fRHS[p] += scale*Na[a]*fip_residual[j];
			}
		}/*for*/

	}/*while*/
//	cout << "\nRHS: "<<fRHS;
}

/* initialize local arrays */
void SS_Optimize_Dual::SetLocalArrays(void)
{
	const char caller[] = "SS_Optimize_Dual::SetLocalArrays";

	/* inherited */
	SmallStrainT::SetLocalArrays();

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
void SS_Optimize_Dual::SetGlobalShape(void)
{
	const char caller[] = "SS_Optimize_Dual::SetGlobalShape";
	/* inherited */
	SolidElementT::SetGlobalShape();
	
	/* material information */
	int material_number = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];

	/* material dependent local arrays.  Obtain element values (in local node ordering) from global fields*/
	if (needs[fNeedsOffset + kstrain]) SetLocalU(fLocPrimalDisp);	
	if (needs[fNeedsOffset + kstrain_last]) SetLocalU(fLocPrimalDisp_last);	


	/*SetLocal fLocData*/
	const iArrayT& nodes = CurrentElement().NodesX();
	for (int n = 0; n < nodes.Length(); n++)
		for (int m = 0; m < NumDOF()  && ElementSupport().StepNumber() > 0; m++)
		{
			double a = nodes[n];
			if (fabs(fData_Coords(a,m) - fLocInitCoords(n,m))>kSmall)
				cout << caller<<"\nMismatch coords for elem: "<<CurrElementNumber()
					<< fData_Coords(a,m)<<fLocInitCoords(n,m)<<endl;
			fLocData(n, m) = fData(nodes[n], m);
		}
	/* using B-bar */
	if (fStrainDispOpt == kMeanDilBbar)
	{
		/*B-bar was calculated in SmallStrainT*/
		/* loop over integration points */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			int ip = fShapes->CurrIP();
			Set_B_bar(fShapes->Derivatives_U(ip), fMeanGradient, fB);
	
			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			{
				/* transpose displacement array */
				fLocPrimalDisp.ReturnTranspose(fLocDispTranspose);
				/* compute strain using B-bar */
				dSymMatrixT& strain = fStrain_List[ip];
				fB.Multx(fLocDispTranspose, strain);
				strain.ScaleOffDiagonal(0.5);

				/* deformation gradient */
				/* transpose displacement array */
				fLocDisp.ReturnTranspose(fLocDispTranspose);
				/* compute strain using B-bar */
				dSymMatrixT& dual_strain = fDual_Strain_List[ip];
				fB.Multx(fLocDispTranspose, dual_strain);
				dual_strain.ScaleOffDiagonal(0.5);
			}

			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				/* transpose displacement array */
				fLocPrimalDisp_last.ReturnTranspose(fLocDispTranspose);

				/* compute strain using B-bar */
				dSymMatrixT& strain = fStrain_last_List[ip];
				fB.Multx(fLocDispTranspose, strain);
				strain.ScaleOffDiagonal(0.5);
			}
				
		}		
	}
	else
	{
		/* loop over integration points */
		for (int i = 0; i < NumIP(); i++)
		{
			if (needs[fNeedsOffset + kstrain])
			{
				/* deformation gradient */
				fShapes->GradU(fLocPrimalDisp, fGradU, i);

				/* symmetric part */
				fStrain_List[i].Symmetrize(fGradU);

				/* displacement gradient */
				fShapes->GradU(fLocDisp, fGradU, i);

				/* symmetric part */
				fDual_Strain_List[i].Symmetrize(fGradU);
			}
			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				/* displacement gradient */
				fShapes->GradU(fLocPrimalDisp_last, fGradU, i);

				/* symmetric part */
				 fStrain_last_List[i].Symmetrize(fGradU);
			}			
			
		}
	}
}

/* current element operations */
bool SS_Optimize_Dual::NextElement(void)
{
	const char caller[] = "SS_Optimize_Dual::NextElement";
	/* inherited */
	bool result = SolidElementT::NextElement();
	
	/* get material pointer */
	if (result)
	{
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[CurrentElement().MaterialNumber()];
	
		/* cast is safe since class contructs materials list */
		fOptimize_Mat =  TB_DYNAMIC_CAST(SSOptimize_MatT*, pcont_mat);
		
		if(!fOptimize_Mat)
			ExceptionT::GeneralFail(caller, "Unable to cast material");
	}
	
	return result;
}

