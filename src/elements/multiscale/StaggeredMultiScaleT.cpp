/* $Id: StaggeredMultiScaleT.cpp,v 1.39 2004-07-15 08:28:27 paklein Exp $ */
#include "StaggeredMultiScaleT.h"

#include "ShapeFunctionT.h"
#include "Traction_CardT.h"
#include "ifstreamT.h"

#include "VMF_Virtual_Work_EqT.h"
#include "VMS_BCJT.h"
//#include "VMS_EZ3T.h"
#include "E_Pr_MatlT.h"
#include "Iso_MatlT.h"
#include "BCJ_MatlT.h"
#include "OutputSetT.h"

#include "ofstreamT.h"

using namespace Tahoe;

/* parameters */
static int knum_d_state = 1; // double's needed per ip
static int knum_i_state = 1; // int's needed per ip

//---------------------------------------------------------------------

/* constructor */
StaggeredMultiScaleT::StaggeredMultiScaleT(const ElementSupportT& support, const FieldT& coarse, const FieldT& fine):
	ElementBaseT(support),
	ua(LocalArrayT::kDisp),
	ua_n(LocalArrayT::kLastDisp),
	ub(LocalArrayT::kDisp),
	DDub(LocalArrayT::kAcc),
	ub_n(LocalArrayT::kLastDisp),
	fInitCoords(LocalArrayT::kInitCoords),
	fCurrCoords(LocalArrayT::kCurrCoords),
	fBodySchedule(NULL),
	fBody(NumDOF()),
	fTractionBCSet(0),
	fCoarse(coarse),
	fFine(fine),
	fKa_1(ElementMatrixT::kNonSymmetric),
	fKb_1(ElementMatrixT::kNonSymmetric),
	fKa_2(ElementMatrixT::kNonSymmetric),
	fKb_2(ElementMatrixT::kNonSymmetric),
	fEquation_1(NULL),
	fEquation_2(NULL),
	render_settings_file_name(32),
	surface_file_name(32),
	write_file_name(32),
	iDesired_Force_Node_Num(0),
	iDesired_Force_Direction(0),
	bStep_Complete(0)
{
ExceptionT::GeneralFail("StaggeredMultiScaleT::StaggeredMultiScaleT", "out of date");
#if 0	
	int i;
	/* check - some code below assumes that both fields have the
	 * same dimension. TEMP? */
	    if (fCoarse.NumDOF() != fFine.NumDOF()) 
				ExceptionT::BadInputValue("StaggeredMultiScaleT::StaggeredMultiScaleT");

	/* read parameters from input */
	ifstreamT& in = ElementSupport().Input();
	in >> fGeometryCode; //TEMP - should actually come from the geometry database
	in >> fNumIP;

	bLogical_Switches.Dimension ( kNUM_LOGICAL_SWITCHES );
	iMaterial_Data.Dimension ( kNUM_IMAT_TERMS );
	fMaterial_Data.Dimension ( kNUM_FMAT_TERMS );

	in >> render_switch; 
	in >> render_settings_file_name; 
	in >> surface_file_name; 
	in >> render_time;
	in >> bLogical_Switches[k__Diagnosis_Variables];	
	in >> render_variable_group;	// {0,1} --> { Coarse, Fine }

	in >> num_tensors_to_render; 
	Render_Tensor_Names.Dimension ( num_tensors_to_render ); 
	for (i=0; i<num_tensors_to_render; i++) {
		Render_Tensor_Names[i].Dimension (32); 
		in >> Render_Tensor_Names[i];
	}
	
	in >> num_scalars_to_render; 
	Render_Scalar_Names.Dimension ( num_scalars_to_render ); 
	for (i=0; i<num_scalars_to_render; i++) {
		Render_Scalar_Names[i].Dimension (32); 
		in >> Render_Scalar_Names[i];
	}

	in >> write_file_name;
	in >> component_i; 							
	in >> component_j; 							
	in >> Elmt2Write; 							// ELmt for data dump 
	in >> ElmtIP2Write; 						// IP for data dump 
	in >> bLog_Strain;
	in >> cube_bottom_elmt; 
	in >> cube_bottom_elmt_bottom_local_node; 
	in >> cube_top_elmt; 
	in >> cube_top_elmt_top_local_node; 

	//-- Convert natural numbers to C style start at 0
	
	component_i--; 							
	component_j--; 							
	Elmt2Write--; 							
	ElmtIP2Write--; 						
	cube_top_elmt--; 
	cube_top_elmt_top_local_node--; 
	cube_bottom_elmt--; 
	cube_bottom_elmt_bottom_local_node--; 

	//in >> render_variable_order;  // {1,2,4} --> { Scalar (S[]), Matrix (A[]), 4th Order Tendor (T4[]) }

	render_data_stored = 0;  // obsolete very soon now that iter dumps fixed

	in >> iFineScaleModelType; 

	//-- General parameters
	in >> fMaterial_Data[k__E];
	in >> fMaterial_Data[k__Pr];
	in >> fMaterial_Data[k__Density];

	//-- Classic BCJ parameters
	in >> fMaterial_Data[k__f];
	in >> fMaterial_Data[k__V];
	in >> fMaterial_Data[k__Y];

	//-- CC^ee Parameters
	in >> bLogical_Switches[k__Control_Eb]; 
	in >> fMaterial_Data[k__Gamma_b];  // Eb strain rate at yeild for .001 true strain (Eb is usually close to E^true)
	in >> fMaterial_Data[k__Yield_Strain];
	in >> fMaterial_Data[k__Beta_tilde];
	in >> fMaterial_Data[k__Rho_tilde];
	in >> fMaterial_Data[k__AlphaY];
	in >> fMaterial_Data[k__E2];
	fMaterial_Data[k__Pi] 	= fMaterial_Data[k__Rho_tilde]  * fMaterial_Data[k__Yield_Strain];
	fMaterial_Data[k__Rho] 	= fMaterial_Data[k__Beta_tilde] / fMaterial_Data[k__Yield_Strain];
	fMaterial_Data[k__E1]		= fMaterial_Data[k__E];
	fMaterial_Data[k__Pr1]	= fMaterial_Data[k__Pr];
	fMaterial_Data[k__Pr2]	= fMaterial_Data[k__Pr];

	//-- Backstress Parameters
	in >> iMaterial_Data[k__BS_Type];
	in >> fMaterial_Data[k__c_zeta];
	in >> fMaterial_Data[k__l];

	//-- Isotropic Hardening Parameters
	in >> iMaterial_Data[k__IH_Type];
	in >> fMaterial_Data[k__K];
	in >> fMaterial_Data[k__H];

	//-- Extra Integer Parameters for Future Developments (no need to modify input decks)
	in >> num_extra_integer_vars; 
	Extra_Integer_Vars.Dimension (  num_extra_integer_vars ); 
	for (i=0; i<num_extra_integer_vars; i++) 
		in >> Extra_Integer_Vars[i];

	if 	(num_extra_integer_vars >= 1) 
		bLogical_Switches[k__Del_Curl_sE] = Extra_Integer_Vars[0];

	if 	(num_extra_integer_vars >= 3) {
		iDesired_Force_Node_Num 	= Extra_Integer_Vars[1];
		iDesired_Force_Direction 	= Extra_Integer_Vars[2];
		iDesired_Force_Node_Num--; 
		iDesired_Force_Direction--; 
	}

	//-- Extra Double Parameters for Future Developments (no need to modify input decks)
	in >> num_extra_double_vars; 
	Extra_Double_Vars.Dimension ( num_extra_double_vars ); 
	for (i=0; i<num_extra_double_vars; i++) 
		in >> Extra_Double_Vars[i];
	
	Echo_Input_Data();
	
	/* allocate the global stack object (once) */
	extern FEA_StackT* fStack;
	if (!fStack) fStack = new FEA_StackT;
#endif
}

//---------------------------------------------------------------------

/* destructor */
StaggeredMultiScaleT::~StaggeredMultiScaleT(void) 
{  
	delete fShapes;
	delete fEquation_1; 
	delete fEquation_2; 
	delete fCoarseMaterial; 
	delete fFineMaterial;

	/* free the global stack object (once) */
	extern FEA_StackT* fStack;
	if (fStack) {
		delete fStack;
		fStack = NULL;
	}

	var_plot_file.close(); 
}

//--------------------------------------------------------------------

void StaggeredMultiScaleT::Echo_Input_Data(void) {

	int i;
	cout << "#######################################################" << endl; 
	cout << "############### ECHO VMS DATA #########################" << endl; 
	cout << "#######################################################" << endl; 

	//################## rendering data ##################

	cout << "render_switch " 										<< render_switch 											<< endl; 
	cout << "render_settings_file_name "				<< render_settings_file_name 					<< endl; 
	cout << "surface_file_name "								<< surface_file_name 									<< endl; 
	cout << "render_time " 											<< render_time 												<< endl;
	cout << "bDiagnosis_variables " 						<< bLogical_Switches[k__Diagnosis_Variables] << endl;
	cout << "render_variable_group " 						<< render_variable_group 							<< endl;
	cout << "num_tensors_to_render "  					<< num_tensors_to_render							<< endl;
	for (i=0; i<num_tensors_to_render; i++) 
		cout << "Render_Tensor_Names[i] "					<< Render_Tensor_Names[i]							<< endl;
	cout << "num_scalar_to_render "  						<< num_scalars_to_render							<< endl;
	for (i=0; i<num_scalars_to_render; i++) 
		cout << "Render_Scalar_Names[i] "					<< Render_Scalar_Names[i]							<< endl;
	cout << "write_file_name "									<< write_file_name										<< endl;
	cout << "component_i "							 				<< component_i 												<< endl;
	cout << "component_j "							 				<< component_j 												<< endl;
	cout << "Elmt2Write "							 					<< Elmt2Write 												<< endl;
	cout << "ElmtIP2Write "						 					<< ElmtIP2Write 											<< endl;
	cout << "bLog_Strain "											<< bLog_Strain												<< endl;
	cout << "cube_bottom_elmt "									<< cube_bottom_elmt 									<< endl;
	cout << "cube_bottom_elmt_bottom_local_node " << cube_bottom_elmt_bottom_local_node << endl;
	cout << "cube_top_elmt" 										<< cube_top_elmt											<< endl;
	cout << "cube_top_elmt_top_local_node " 		<< cube_top_elmt_top_local_node 			<< endl;

	//################## material data ##################

	//-- General parameters
	cout << "iFineScaleModelType " 							<< iFineScaleModelType 								<< endl;  
	cout << "fMaterial_Data[k__E] "  						<< fMaterial_Data[k__E] 							<< endl;
	cout << "fMaterial_Data[k__Pr] " 						<< fMaterial_Data[k__Pr] 							<< endl;

	//-- Classic BCJ parameters
	cout << "fMaterial_Data[k__Density] " 			<< fMaterial_Data[k__Density] 				<< endl;
	cout << "fMaterial_Data[k__f] " 						<< fMaterial_Data[k__f] 							<< endl;
	cout << "fMaterial_Data[k__V] " 						<< fMaterial_Data[k__V] 							<< endl;
	cout << "fMaterial_Data[k__Y] " 						<< fMaterial_Data[k__Y] 							<< endl;

	//-- CC^ee Parameters
	cout << "bLogical_Switches[k__Control_Eb] " << bLogical_Switches[k__Control_Eb] 	<< endl; 
	cout << "fMaterial_Data[k__Gamma_b] "  			<< fMaterial_Data[k__Gamma_b] 				<< endl;  
	cout << "fMaterial_Data[k__Yield_Strain] " 	<< fMaterial_Data[k__Yield_Strain] 		<< endl;
	cout << "fMaterial_Data[k__Beta_tilde] " 		<< fMaterial_Data[k__Beta_tilde] 			<< endl;
	cout << "fMaterial_Data[k__Rho_tilde] " 		<< fMaterial_Data[k__Rho_tilde] 			<< endl;
	cout << "fMaterial_Data[k__AlphaY] " 				<< fMaterial_Data[k__AlphaY] 					<< endl;
	cout << "fMaterial_Data[k__E2] " 						<< fMaterial_Data[k__E2] 							<< endl;

	//-- Backstress Parameters
	cout << "iMaterial_Data[k__BS_Type] "				<< iMaterial_Data[k__BS_Type]					<< endl;
	cout << "fMaterial_Data[k__c_zeta] " 				<< fMaterial_Data[k__c_zeta]					<< endl;
	cout << "fMaterial_Data[k__l] " 						<< fMaterial_Data[k__l]								<< endl;

	//-- Isotropic Hardening Parameters
	cout << "iMaterial_Data[k__IH_Type] " 			<< iMaterial_Data[k__IH_Type]					<< endl;
	cout << "fMaterial_Data[k__K] " 						<< fMaterial_Data[k__K]								<< endl;
	cout << "fMaterial_Data[k__H] "							<< fMaterial_Data[k__H]								<< endl;

	for (i=0; i<num_extra_integer_vars; i++) 
		cout << "Extra Integer Variable Number "<<i<<" = "<< Extra_Integer_Vars[i] << endl;

	for (i=0; i<num_extra_double_vars; i++) 
		cout << "Extra Double Variable Number "<<i<<" = "<< Extra_Double_Vars[i] << endl;

}

void StaggeredMultiScaleT::Initialize(void)
{
ExceptionT::GeneralFail("StaggeredMultiScaleT::Initialize", "out of date");
#if 0
	/* inherited */
	ElementBaseT::Initialize();
	
	/* dimensions (notation as per Hughes' Book) */
	n_ip = fNumIP;
	n_sd = NumSD();
	n_df = NumDOF(); 
	n_en = NumElementNodes();
	n_np = ElementSupport().NumNodes();
	n_el = NumElements();

	n_en_x_n_df = n_en*n_df;

	/* set local arrays for coarse scale */
	ub.Dimension (n_en, n_df);
	DDub.Dimension (n_en, n_df);
	ub_n.Dimension (n_en, n_df);
	del_ub.Dimension (n_en, n_df);
	del_ub_vec.Dimension (n_en_x_n_df);
	fCoarse.RegisterLocal(ub);
	fCoarse.RegisterLocal(ub_n);

	/* set local arrays for fine scale */
	ua.Dimension (n_en, n_df);
	ua_n.Dimension (n_en, n_df);
	del_ua.Dimension (n_en, n_df);
	del_ua_vec.Dimension (n_en_x_n_df);
	fFine.RegisterLocal(ua);
	fFine.RegisterLocal(ua_n);

	/* set shape functions */
	fInitCoords.Dimension(n_en, n_sd);
	ElementSupport().RegisterCoordinates(fInitCoords);	
	fCurrCoords.Dimension(n_en, n_sd);
	fShapes = new ShapeFunctionT(fGeometryCode, fNumIP, fCurrCoords);
	fShapes->Initialize();
	
	/* allocate state variable storage */
	int num_ip = 1; //TEMP - need to decide where to set the number of integration
	fdState_new.Dimension(n_el, num_ip*knum_d_state);
	fdState.Dimension(n_el, num_ip*knum_d_state);
	fiState_new.Dimension(n_el, num_ip*knum_i_state);
	fiState.Dimension(n_el, num_ip*knum_i_state);
	
	/* storage for the fine scale equation numbers */
	fEqnos_fine.Dimension(n_el, n_en*n_df);
	fEqnos_fine = -1;

	/* initialize state variables */
	fdState = 0;
	fiState = 0;

	/* set cards to data in array - NOT NEEDED IF YOU'RE NOT
	 * GOING TO USE THE ElementCardT ARRAY? */
	for (int i= 0; i < fElementCards.Length(); i++)
		fElementCards[i].Set(fiState.MinorDim(), fiState(i), fdState.MinorDim(), fdState(i));
		                     
	/* construct the black boxs */  

	//Select_Equations ( CoarseScaleT::kVMF_Virtual_Work_Eq,	FineScaleT::kVMS_BCJ );
	Select_Equations ( CoarseScaleT::kVMF_Virtual_Work_Eq, iFineScaleModelType );
	fEquation_2 -> Initialize ( n_ip, n_sd, n_en, ElementSupport().StepNumber() );
	//step_number_last_iter = 0; 
	//step_number_last_iter = ElementSupport().StepNumber();  // This may crash or not work

	/* FEA Allocation */

	fGRAD_ua.FEA_Dimension 		( fNumIP, n_sd,n_sd );
	fGRAD_ub.FEA_Dimension 		( fNumIP, n_sd,n_sd );
	fGRAD_ua_n.FEA_Dimension 	( fNumIP, n_sd,n_sd );
	fGRAD_ub_n.FEA_Dimension 	( fNumIP, n_sd,n_sd );

	fKa_1.Dimension 		( n_en_x_n_df, n_en_x_n_df );
	fKb_1.Dimension 		( n_en_x_n_df, n_en_x_n_df );
	fKa_2.Dimension 		( n_en_x_n_df, n_en_x_n_df );
	fKb_2.Dimension 		( n_en_x_n_df, n_en_x_n_df );

	fFint_1.Dimension 	( n_en_x_n_df );
	fFext_1.Dimension 	( n_en_x_n_df );
	fR_2.Dimension 	( n_en_x_n_df );


	fFEA_Shapes.Construct	( fNumIP,n_sd,n_en );

	Render_Tensor.Dimension ( n_el );
	Render_Scalar.Dimension ( n_el );
	for (int e=0; e<n_el; e++) {
		Render_Tensor[e].Construct ( num_tensors_to_render, n_ip, n_sd, n_sd );	
		Render_Scalar[e].Construct ( num_scalars_to_render, n_ip );	
	}

  /* streams */
	ifstreamT& in  = ElementSupport().Input();
	ofstreamT& out = ElementSupport().Output();

	/* open plot file in the same directory as the input file */
	StringT path;
	path.FilePath(in.filename());
	StringT file_path = write_file_name;
	file_path.ToNativePathName();
	file_path.Prepend(path);
	var_plot_file.open(file_path);
	nodal_forces_file.open("nodal_forces.data");
	displacements_file.open("displacements.data");

	int n_wfld = (int) bLog_Strain + num_tensors_to_render + num_scalars_to_render;
	for (int v=0; v<n_wfld; v++)
		var_plot_file << 0.0 << " "; // Accounts for initial time
	var_plot_file << endl; 

	/* storage for integration point stresses */
	fIPVariable.Dimension (n_el, fNumIP*dSymMatrixT::NumValues(n_sd));
	fIPVariable = 0.0;

	/* allocate storage for nodal forces */
	fForces_at_Node.Dimension ( n_sd );

	/* render switch */
	if (render_switch)
		Init_Render();

	 /* body force specification */
	 fDOFvec.Dimension(n_df);
	 EchoBodyForce(in, out);

	 /* echo traction B.C.'s */
	 EchoTractionBC(in, out);
#endif
}

//---------------------------------------------------------------------

void StaggeredMultiScaleT::RHSDriver(void)	
{
	int curr_group = ElementSupport().CurrentGroup();

	/* traction boundary conditions acting on the coarse scale equations */
	if (curr_group == fCoarse.Group()) 
		ApplyTractionBC();

	/* choose solution method */
	if (fCoarse.Group() == fFine.Group())
	  RHSDriver_monolithic();
	 else
	  RHSDriver_staggered();
}
//---------------------------------------------------------------------

void StaggeredMultiScaleT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* doing monolithic solution */
	if (fCoarse.Group() == fFine.Group())
	{
		int ndof_fine = fFine.NumDOF();
		int ndof_coarse = fCoarse.NumDOF();
		int nen = NumElementNodes();
	
		/* loop over connectivity blocks */
		for (int i = 0; i < fEqnos.Length(); i++)
		{
			/* connectivities */
			const iArray2DT& connects = *(fConnectivities[i]);
			int nel = connects.MajorDim();
		
			/* dimension */
			fEqnos[i].Dimension(nel, nen*(ndof_coarse + ndof_fine));
			iArray2DT coarse_eq(nel, nen*ndof_fine);
			iArray2DT fine_eq(nel, nen*ndof_fine);
			
			/* get equation numbers */
			fCoarse.SetLocalEqnos(connects, coarse_eq);
			fFine.SetLocalEqnos(connects, fine_eq);
			
			/* write into one array */
			fEqnos[i].BlockColumnCopyAt(coarse_eq, 0);
			fEqnos[i].BlockColumnCopyAt(fine_eq, coarse_eq.MinorDim());

			/* add to list of equation numbers */
			eq_1.Append(&fEqnos[i]);
		}
	
		/* reset pointers to element cards */
		SetElementCards(fBlockData, fConnectivities, fEqnos, fElementCards);
	}
	else
	{
		/* ElementBaseT handles equation array for the coarse scale */
		if (ElementSupport().CurrentGroup() == fCoarse.Group())
			ElementBaseT::Equations(eq_1,eq_2);

		/* fine scale equations */
		if (ElementSupport().CurrentGroup() == fFine.Group())
		{
			/* collect local equation numbers */
			fFine.SetLocalEqnos(fConnectivities, fEqnos_fine);
		
			eq_1.Append(&fEqnos_fine);
		}
	}
}

#if 0
//---------------------------------------------------------------------

/* form group contribution to the stiffness matrix and RHS */
void StaggeredMultiScaleT::RHSDriver_staggered(void)	// LHS too!	(This was original RHSDriver()
{
 
	int curr_group = ElementSupport().CurrentGroup();

	/* stress output work space */
	int n_stress = dSymMatrixT::NumValues(NumSD());
	dArray2DT   out_variable_all;
	dSymMatrixT out_variable;

	/** Time Step Increment */
	double delta_t = ElementSupport().TimeStep();
						time = ElementSupport().Time();
						step_number = ElementSupport().StepNumber();

	iArrayT fine_eq;

	if ( curr_group == fCoarse.Group() )
		cout << "############### Coarse Group ###############\n";

	if ( curr_group == fFine.Group() )
		cout << "############### Fine Group ###############\n";

	//cout <<" s= " << render_switch <<"; t= "<< time << "; rt= " << render_time << "\n";
 
	/* loop over elements */
	int e,v,l;
	Top();
	while (NextElement())
	{
		e = CurrElementNumber();

		SetLocalU (ua);			 SetLocalU (ua_n);
		SetLocalU (ub);			 SetLocalU (ub_n);

		del_ua.DiffOf (ua, ua_n);
		del_ub.DiffOf (ub, ub_n);

	 	SetLocalX(fInitCoords); 
		fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, ua, 1.0, ub); 
		fShapes->SetDerivatives(); 

		if (bStep_Complete) { //-- Done iterating, get result data 
			if (bLog_Strain) {	//-- For calculation of Lf : epsilon = e^(Lf/Lo) 
				if ( e == cube_top_elmt ) 
					x_top = fCurrCoords ( cube_top_elmt_top_local_node, 1 ); // 1 is for the 2 direction
				if ( e == cube_bottom_elmt ) 
					x_bot = fCurrCoords ( cube_bottom_elmt_bottom_local_node, 1 ); // 1 is for the 2 direction
			}
		}
		
		/** repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradients 		( fShapes, 	ua, ua_n, fGRAD_ua, fGRAD_ua_n );
		Convert.Gradients 		( fShapes, 	ub, ub_n, fGRAD_ub, fGRAD_ub_n );
		Convert.Shapes				(	fShapes, 	fFEA_Shapes );
		Convert.Displacements	(	del_ua, 	del_ua_vec  );
		Convert.Displacements	(	del_ub, 	del_ub_vec  );

		Convert.Na						(	n_en, fShapes, 	fFEA_Shapes );

		/** Construct data used in BOTH FineScaleT and CoarseScaleT (F,Fa,Fb,grad_ua,...etc.)
		 * 	Presently, Tahoe cannot exploit this fact.  n and np1 are calculated for coarse field, then
		 * 	calculated again for fine field -- this is a waste and defeats the putpose of VMS_VariableT. 
		 *  Note: n is last time step (known data), no subscript,np1 or (n+1) is the 
		 *  next time step (what were solving for)   */

		VMS_VariableT np1(	fGRAD_ua, 	fGRAD_ub 	 ); // Many variables at time-step n+1
		VMS_VariableT   n(	fGRAD_ua_n, fGRAD_ub_n );	// Many variables at time-step n

#if 0
		cout <<" e = "<<e<<endl;
		if (e==0) {
			FEA_dMatrixT TEMP(n_ip,n_sd,n_sd);
			TEMP = np1.Get( VMS::kFb );
			TEMP.Print(" Fb TEMP ");
		}
#endif
		
		/* which field */
	  //SolverGroup 1 (gets field 2) <-- ub (obtained by a rearranged Equation I)
		if ( curr_group == fCoarse.Group() || (bStep_Complete && render_variable_group==0) )	
		{
			if (bStep_Complete) { //-- Done iterating, get result data from converged upon displacements 

				fEquation_1 -> Construct ( fFEA_Shapes, fCoarseMaterial, np1, n, step_number, delta_t );

				for (v=0; v<num_tensors_to_render; v++ ) 
					fEquation_1 -> Get ( Render_Tensor_Names[v], Render_Tensor[e][v] ); 
				for (v=0; v<num_scalars_to_render; v++ ) 
					fEquation_1 -> Get ( Render_Scalar_Names[v], Render_Scalar[e][v] ); 

				//-- Store/Register data in classic tahoe manner 
				out_variable_all.Set(fNumIP, n_stress, fIPVariable(CurrElementNumber()));
				for (l=0; l < fNumIP; l++) {
					out_variable.Set(NumSD(), out_variable_all(l));
					out_variable.FromMatrix(Render_Tensor[e][0][l]);
				} 
			}
			else { //-- Still Iterating

				/** Compute N-R matrix equations */
				fEquation_1 -> Construct ( fFEA_Shapes, fCoarseMaterial, np1, n, step_number, delta_t );
				fEquation_1 -> Form_LHS_Ka_Kb ( fKa_1, fKb_1 );
				fEquation_1 -> Form_RHS_F_int ( fFint_1 );

				/** Set coarse LHS */
				fLHS = fKb_1;

				/** Compute coarse RHS (or Fint_bar_2 in FAXed notes) */
				fKa_1.Multx ( del_ua_vec, fRHS );
				fRHS += fFint_1; 
				fRHS *= -1.0; 

				/** Compute Traction B.C. and Body Forces */
				Get_Fext_I ( fFext_1 );
				fRHS += fFext_1;
			
				/* add to global equations */
				ElementSupport().AssembleLHS	( fCoarse.Group(), fLHS, CurrentElement().Equations() );
				ElementSupport().AssembleRHS 	( fCoarse.Group(), fRHS, CurrentElement().Equations() );
			}
		}

		// SolverGroup 2 (gets field 1) <-- ua (obtained by a rearranged Equation 2)
		else if (curr_group == fFine.Group() || (bStep_Complete && render_variable_group==1) )	
		{

			if (bStep_Complete) { //-- Done iterating, get result data from converged upon displacements 

				fEquation_2 -> Construct ( fFEA_Shapes, fFineMaterial, np1, n, step_number, delta_t );

				for (v=0; v<num_tensors_to_render; v++ )  
					fEquation_2 -> Get ( Render_Tensor_Names[v], Render_Tensor[e][v] ); 
				for (v=0; v<num_scalars_to_render; v++ ) 
					fEquation_2 -> Get ( Render_Scalar_Names[v], Render_Scalar[e][v] ); 

				//-- Store/Register data in classic tahoe manner 
				out_variable_all.Set(fNumIP, n_stress, fIPVariable(CurrElementNumber()));
				for (l=0; l < fNumIP; l++) {
					out_variable.Set(NumSD(), out_variable_all(l));
					out_variable.FromMatrix(Render_Tensor[e][0][l]);
				} 
			}
			else { //-- Still Iterating

				/** Compute N-R matrix equations */
				fEquation_2 -> Construct ( fFEA_Shapes, fFineMaterial, np1, n, step_number, delta_t, FEA::kBackward_Euler );
				fEquation_2 -> Form_LHS_Ka_Kb ( fKa_2, 	fKb_2 );
				fEquation_2 -> Form_RHS_F_int ( fR_2 );

				//cout << "fKa_2 = "<< fKa_2 << endl;
				//cout << "fKb_2 = "<< fKb_2 << endl;

				/** Set LHS */
				fLHS = fKa_2;	
		
				/** Compute fine RHS (or Fint_bar_2 in FAXed notes)  */
				fKb_2.Multx ( del_ub_vec, fRHS );
				fRHS += fR_2; 
				fRHS *= -1.0; 
		
				/* fine scale equation numbers */
				fEqnos_fine.RowAlias ( CurrElementNumber(), fine_eq );

				/* add to global equations */
				ElementSupport().AssembleLHS ( fFine.Group(), fLHS, fine_eq );
				ElementSupport().AssembleRHS ( fFine.Group(), fRHS, fine_eq );
			}

		}
		else throw ExceptionT::kGeneralFail;
	}

}
#endif


//---------------------------------------------------------------------

void StaggeredMultiScaleT::LHSDriver(GlobalT::SystemTypeT)
{
  /** Everything done in RHSDriver for efficiency */
	//cout << "############### In LHS Driver ############### \n";

}

//---------------------------------------------------------------------

void StaggeredMultiScaleT::Select_Equations (const int &iCoarseScale,const int &iFineScale )
{
	/** Choices for Coarse-Scale Equation */

	switch ( iCoarseScale )	{

		case CoarseScaleT::kVMF_Virtual_Work_Eq :
			fEquation_1 		= new VMF_Virtual_Work_EqT;
			fCoarseMaterial = new Iso_MatlT;
			fCoarseMaterial -> Assign ( Iso_MatlT::kE, 				fMaterial_Data[k__E] 				); 
			fCoarseMaterial -> Assign ( Iso_MatlT::kPr, 			fMaterial_Data[k__Pr] 			); 
			fCoarseMaterial -> Assign ( Iso_MatlT::kE1, 			fMaterial_Data[k__E1] 			); 
			fCoarseMaterial -> Assign ( Iso_MatlT::kPr1, 			fMaterial_Data[k__Pr1] 			); 
			fCoarseMaterial -> Assign ( Iso_MatlT::kE2, 			fMaterial_Data[k__E2] 			); 
			fCoarseMaterial -> Assign ( Iso_MatlT::kPr2, 			fMaterial_Data[k__Pr2] 			); 
			fCoarseMaterial -> Assign ( Iso_MatlT::kPi, 			fMaterial_Data[k__Pi] 			); 
			fCoarseMaterial -> Assign ( Iso_MatlT::kRho, 			fMaterial_Data[k__Rho] 			); 
			fCoarseMaterial -> Assign ( Iso_MatlT::kGamma_b, 	fMaterial_Data[k__Gamma_b] 	); 
			fCoarseMaterial -> Assign ( Iso_MatlT::kAlphaY, 	fMaterial_Data[k__AlphaY] 	); 

			fCoarseMaterial -> E_Nu_2_Lamda_Mu	( Iso_MatlT::kE,			Iso_MatlT::kPr,	
																						Iso_MatlT::kLamda, 	Iso_MatlT::kMu 		);

			fCoarseMaterial -> E_Nu_2_Lamda_Mu	( Iso_MatlT::kE1,			Iso_MatlT::kPr1,	
																						Iso_MatlT::kLamda1, Iso_MatlT::kMu1 	);

			fCoarseMaterial -> E_Nu_2_Lamda_Mu	( Iso_MatlT::kE2,			Iso_MatlT::kPr2,	
																						Iso_MatlT::kLamda2, Iso_MatlT::kMu2 	);

			fEquation_1 -> bControl_Eb 						=  	bLogical_Switches [k__Control_Eb]; 	

			break;

		case CoarseScaleT::kLDV :
			//fEquation_1 		= new LDVT; 
			fCoarseMaterial = new Iso_MatlT;
			break;

		case CoarseScaleT::kStraight :
			//fEquation_1 		= new StraightT; 
			fCoarseMaterial = new E_Pr_MatlT;
			break;

		default :
			cout << " StaggeredMultiScaleT::Select_Equations() .. ERROR >> bad iCoarseScale \n";
			break;
	}

	/** Choices for Fine-Scale Equation */

	switch ( iFineScale )	{

		case FineScaleT::kVMS_BCJ :
			fEquation_2 	 = new VMS_BCJT;
			fFineMaterial  = new BCJ_MatlT;							
			fFineMaterial -> Assign (		BCJ_MatlT::kE, 				fMaterial_Data[k__E] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr, 			fMaterial_Data[k__Pr] 		); 	
			fFineMaterial -> Assign (		BCJ_MatlT::kE1, 			fMaterial_Data[k__E1] 		); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr1, 			fMaterial_Data[k__Pr1] 		); 	
			fFineMaterial -> Assign (		BCJ_MatlT::kE2, 			fMaterial_Data[k__E2] 		); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr2, 			fMaterial_Data[k__Pr2] 		); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 				fMaterial_Data[k__f] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kV, 				fMaterial_Data[k__V] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kY, 				fMaterial_Data[k__Y] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kl, 				fMaterial_Data[k__l] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kc_zeta, 	fMaterial_Data[k__c_zeta] ); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kH, 				fMaterial_Data[k__H] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPi, 			fMaterial_Data[k__Pi] 		); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kRho, 			fMaterial_Data[k__Rho] 		); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kGamma_b, 	fMaterial_Data[k__Gamma_b] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kAlphaY, 	fMaterial_Data[k__AlphaY] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPlastic_Modulus_K, 	fMaterial_Data[k__K] 	); 	

			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE,			BCJ_MatlT::kPr,	
																					BCJ_MatlT::kLamda, 	BCJ_MatlT::kMu 		);

			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE1,			BCJ_MatlT::kPr1,	
																					BCJ_MatlT::kLamda1, BCJ_MatlT::kMu1 	);

			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE2,			BCJ_MatlT::kPr2,	
																					BCJ_MatlT::kLamda2, BCJ_MatlT::kMu2 	);

			fEquation_2 -> bDiagnosis_Variables 	=  	bLogical_Switches [k__Diagnosis_Variables]; 	
			fEquation_2 -> bControl_Eb 					=  	bLogical_Switches [k__Control_Eb]; 	
			fEquation_2 -> Back_Stress_Type 			=  	iMaterial_Data[k__BS_Type]; 	
			fEquation_2 -> Iso_Hard_Type					= 	iMaterial_Data[k__IH_Type];

			//fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 			0.0000001 		); 	
			//fFineMaterial -> Assign ( 	BCJ_MatlT::kV, 			258870.0   		); 	
			break;


		case FineScaleT::kVMS_BCJ_X :
			fEquation_2 	 = new VMS_BCJ_XT;
			fFineMaterial  = new BCJ_MatlT;							
			fFineMaterial -> Assign (		BCJ_MatlT::kE, 				fMaterial_Data[k__E] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr, 			fMaterial_Data[k__Pr] 		); 	
			fFineMaterial -> Assign (		BCJ_MatlT::kE1, 			fMaterial_Data[k__E1] 		); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr1, 			fMaterial_Data[k__Pr1] 		); 	
			fFineMaterial -> Assign (		BCJ_MatlT::kE2, 			fMaterial_Data[k__E2] 		); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr2, 			fMaterial_Data[k__Pr2] 		); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 				fMaterial_Data[k__f] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kV, 				fMaterial_Data[k__V] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kY, 				fMaterial_Data[k__Y] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kl, 				fMaterial_Data[k__l] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kc_zeta, 	fMaterial_Data[k__c_zeta] ); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kH, 				fMaterial_Data[k__H] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPi, 			fMaterial_Data[k__Pi] 		); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kRho, 			fMaterial_Data[k__Rho] 		); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kGamma_b, 	fMaterial_Data[k__Gamma_b] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kAlphaY, 	fMaterial_Data[k__AlphaY] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPlastic_Modulus_K, 	fMaterial_Data[k__K] 	); 	

			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE,			BCJ_MatlT::kPr,	
																					BCJ_MatlT::kLamda, 	BCJ_MatlT::kMu 		);

			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE1,			BCJ_MatlT::kPr1,	
																					BCJ_MatlT::kLamda1, BCJ_MatlT::kMu1 	);

			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE2,			BCJ_MatlT::kPr2,	
																					BCJ_MatlT::kLamda2, BCJ_MatlT::kMu2 	);

			fEquation_2 -> bDiagnosis_Variables 	=  	bLogical_Switches [k__Diagnosis_Variables]; 	
			fEquation_2 -> bControl_Eb 					=  	bLogical_Switches [k__Control_Eb]; 	
			fEquation_2 -> Back_Stress_Type 			=  	iMaterial_Data[k__BS_Type]; 	
			fEquation_2 -> Iso_Hard_Type					= 	iMaterial_Data[k__IH_Type];

			//fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 			0.0000001 		); 	
			//fFineMaterial -> Assign ( 	BCJ_MatlT::kV, 			258870.0   		); 	
			break;


		case FineScaleT::kVMS_EZ : 
			fEquation_2 	= new VMS_EZT;
			fFineMaterial = new BCJ_MatlT; // Can use any material class really -- just to hold data 
			fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 				fMaterial_Data[k__f] 			); 	
			break;

#if 0
		case FineScaleT::kVMS_EZ2 : 
			fEquation_2 	= new VMS_EZ2T;
			fFineMaterial = new Iso_MatlT; // <-- not used
			break;
#endif

		case FineScaleT::kVMS_EZ3 : 
			fEquation_2 	= new VMS_EZ3T;
			fFineMaterial = new BCJ_MatlT; 
			fFineMaterial -> Assign (		BCJ_MatlT::kE, 		fMaterial_Data[k__E]	 		); 	// GPa
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr, 	fMaterial_Data[k__Pr] 		); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 		fMaterial_Data[k__f]  		); 	// 1.6e-5
			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE,			BCJ_MatlT::kPr,	
																					BCJ_MatlT::kLamda, 	BCJ_MatlT::kMu 	);
			break;
			
#if 0

		case FineScaleT::kVMS_EZ4 : 
			fEquation_2 	= new VMS_EZ4T;
			fFineMaterial = new BCJ_MatlT; 
			fFineMaterial -> Assign (		BCJ_MatlT::kE, 			168.0 	 		); 	// GPa
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr, 		0.34 				); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 			0.000016 		); 	// 1.6e-5
			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE,			BCJ_MatlT::kPr,	
																					BCJ_MatlT::kLamda, 	BCJ_MatlT::kMu 	);
			break;

		case FineScaleT::kVMS_EZ5 : 
			fEquation_2 	= new VMS_EZ5T;
			fFineMaterial = new BCJ_MatlT; 
			fFineMaterial -> Assign (		BCJ_MatlT::kE, 			168.0 	 		); 	// GPa
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr, 		0.34 				); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 			0.000016 		); 	// 1.6e-5
			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE,			BCJ_MatlT::kPr,	
																					BCJ_MatlT::kLamda, 	BCJ_MatlT::kMu 	);
			break;

		case FineScaleT::kPHEN :
			//fEquation_2 	= new PHENT;
			//fFineMaterial = new VMS_Phen_MaterialT;
			break;
#endif

		default :
			cout << " StaggeredMultiScaleT::Select_Equations() .. ERROR >> bad iFineScale \n";
			break;
	}

}

//---------------------------------------------------------------------

/* return true if the element contributes to the solution of the
 * given group. */
bool StaggeredMultiScaleT::InGroup(int group) const
{
	return group == fCoarse.Group() ||
	       group == fFine.Group();
}

//---------------------------------------------------------------------

/* close current time increment */
void StaggeredMultiScaleT::CloseStep(void)
{
	/* inherited */
	ElementBaseT::CloseStep();

	/* store more recently updated values */
	fdState = fdState_new;
	fiState = fiState_new;
}

//---------------------------------------------------------------------

void StaggeredMultiScaleT::SendOutput(int kincode)
{
#pragma unused(kincode)
	//not implemented
}
//---------------------------------------------------------------------

/* form of tangent matrix */
GlobalT::SystemTypeT StaggeredMultiScaleT::TangentType(void) const
{
	return GlobalT::kNonSymmetric; 
}

//#############################################################################
//#############################################################################
//#############################################################################
//#############################################################################
//############################### NODAL FORCE  ################################
//#############################################################################
//#############################################################################
//#############################################################################
//#############################################################################

/* accumulate the residual force on the specified node */
void StaggeredMultiScaleT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	const char caller[] = "StaggeredMultiScaleT::AddNodalForce";

	/* coarse, fine, or neither */
	bool is_coarse = false;
	dArrayT* element_force = NULL;
	int num_force = 0;
	if (field.Name() == fCoarse.Name()) {
		is_coarse = true;
		element_force = &fFint_1;
		num_force = fCoarse.NumDOF();
	}
	else if (field.Name() == fFine.Name()) {
		is_coarse = false;
		element_force = &fR_2;
		num_force = fFine.NumDOF();
	}
	else
		return;

	/* time Step Increment */
	double delta_t = ElementSupport().TimeStep();
	time = ElementSupport().Time();
	step_number = ElementSupport().StepNumber();

 	/* has (coarse scale) body forces */
	int formBody = 0;
	if (fBodySchedule && fBody.Magnitude() > kSmall)
		formBody = 1;

	/* temp for nodal force */
	dArrayT nodalforce;

	/* loop over elements */
	int e;
	Top();
	while (NextElement())
	{
		int nodeposition;
		const iArrayT& nodes_u = CurrentElement().NodesU();
		if (nodes_u.HasValue(node, nodeposition))
		{
			e = CurrElementNumber();

			SetLocalU (ua);			 SetLocalU (ua_n);
			SetLocalU (ub);			 SetLocalU (ub_n);

			del_ua.DiffOf (ua, ua_n);
			del_ub.DiffOf (ub, ub_n);

			SetLocalX(fInitCoords); 
			fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, ua, 1.0, ub); 
			fShapes->SetDerivatives(); 
		
			/* repackage data to forms compatible with FEA classes (very little cost in big picture) */
			Convert.Gradients 		( fShapes, 	ua, ua_n, fGRAD_ua, fGRAD_ua_n );
			Convert.Gradients 		( fShapes, 	ub, ub_n, fGRAD_ub, fGRAD_ub_n );
			Convert.Shapes				(	fShapes, 	fFEA_Shapes );
			Convert.Displacements	(	del_ua, 	del_ua_vec  );
			Convert.Displacements	(	del_ub, 	del_ub_vec  );
			Convert.Na(	n_en, fShapes, 	fFEA_Shapes );

			/* Construct data used in BOTH FineScaleT and CoarseScaleT (F,Fa,Fb,grad_ua,...etc.)
			 * 	Presently, Tahoe cannot exploit this fact.  n and np1 are calculated for coarse field, then
			 * 	calculated again for fine field -- this is a waste and defeats the putpose of VMS_VariableT. 
			 *  Note: n is last time step (known data), no subscript,np1 or (n+1) is the 
			 *  next time step (what were solving for)   */
			VMS_VariableT np1(	fGRAD_ua, 	fGRAD_ub 	 ); // Many variables at time-step n+1
			VMS_VariableT   n(	fGRAD_ua_n, fGRAD_ub_n );	// Many variables at time-step n

			/* calculate coarse scale nodal force */
			if (is_coarse)
			{
				/* residual and tangent for coarse scale */
				fEquation_1 -> Construct ( fFEA_Shapes, fCoarseMaterial, np1, n, step_number, delta_t );
				fEquation_1 -> Form_LHS_Ka_Kb ( fKa_1, fKb_1 );
				fEquation_1 -> Form_RHS_F_int ( fFint_1 );
				fFint_1 *= -1.0;  

				/* add body force */
				if (formBody) {
//					double density = fCoarseMaterial->Retrieve(Iso_MatlT::kDensity);
					double density = 1.0;
					DDub = 0.0;
					AddBodyForce(DDub);
				
					/* add body force to fRHS */
					fRHS = 0.0;
					FormMa(kConsistentMass, -density, &DDub, NULL);
					fFint_1 += fRHS;
				}
			}
			else /* fine scale nodal force */
			{
				/* residual and tangent for fine scale */
				fEquation_2 -> Construct ( fFEA_Shapes, fFineMaterial, np1, n, step_number, delta_t, FEA::kBackward_Euler );
				fEquation_2 -> Form_LHS_Ka_Kb ( fKa_2, 	fKb_2 );
				fEquation_2 -> Form_RHS_F_int ( fR_2 );
				fR_2 *= -1.0;
			}

			/* loop over nodes (double-noding OK) */
			int dex = 0;
			for (int i = 0; i < nodes_u.Length(); i++)
			{
				if (nodes_u[i] == node)
				{
					/* components for node */
					nodalforce.Set(num_force, element_force->Pointer(dex));
	
					/* accumulate */
					force += nodalforce;
				}
				dex += NumDOF();
			}			
		}
	}

	if (is_coarse && node==0) // Default is the 1st node #1 
		nodal_forces_file << force << endl;
}

//---------------------------------------------------------------------

double StaggeredMultiScaleT::InternalEnergy ( void )
{
	//not implemented
return 0.0;
}

//---------------------------------------------------------------------

/* write restart data to the output stream */
void StaggeredMultiScaleT::WriteRestart(ostream& out) const
{
	/* inherited */
	ElementBaseT::WriteRestart(out);

	/* write state variable data */
	out << fdState;
}

//---------------------------------------------------------------------

/* read restart data to the output stream */
void StaggeredMultiScaleT::ReadRestart(istream& in)
{
	/* inherited */
	ElementBaseT::ReadRestart(in);

	/* write state variable data */
	in >> fdState;
}

//---------------------------------------------------------------------

void StaggeredMultiScaleT::RegisterOutput(void)
{
	/* collect block ID's */
	ArrayT<StringT> block_ID(fBlockData.Length());
	for (int i = 0; i < block_ID.Length(); i++)
		block_ID[i] = fBlockData[i].ID();

	/* output per element - stresses at the integration points */
	int n_stress = dSymMatrixT::NumValues(NumSD());
	ArrayT<StringT> e_labels(fNumIP*n_stress);

	/* over integration points */
	const char* slabels2D[] = {"s11", "s22", "s12"};
	const char* slabels3D[] = {"s11", "s22", "s33", "s23", "s13", "s12"};
	const char** slabels = (NumSD() == 2) ? slabels2D : slabels3D;
	int count = 0;
	for (int j = 0; j < fNumIP; j++)
	{
		StringT ip_label;
		ip_label.Append("ip", j+1);
			
		/* over stress components */
		for (int i = 0; i < n_stress; i++)
		{
			e_labels[count].Clear();
			e_labels[count].Append(ip_label, ".", slabels[i]);
			count++;
		}
	}		

	/* output per node */
	int num_node_output = fCoarse.NumDOF() + fFine.NumDOF() + n_stress;
	ArrayT<StringT> n_labels(num_node_output);
	count = 0;

	/* labels from fine scale */
	const ArrayT<StringT>& fine_labels = fFine.Labels();
	for (int i = 0; i < fine_labels.Length(); i++)
		n_labels[count++] = fine_labels[i];

	/* labels from coarse scale */
	const ArrayT<StringT>& coarse_labels = fCoarse.Labels();
	for (int i = 0; i < coarse_labels.Length(); i++)
		n_labels[count++] = coarse_labels[i];

	/* labels from stresses at the nodes */
	for (int i = 0; i < n_stress; i++)
		n_labels[count++] = slabels[i];

	/* set output specifier */
	OutputSetT output_set(fGeometryCode, block_ID, fConnectivities, n_labels, e_labels, false);
		
	/* register and get output ID */
	fOutputID = ElementSupport().RegisterOutput(output_set);
}

//#############################################################################
//#############################################################################
//#############################################################################
//#############################################################################
//############################### WRITE OUTPUT ################################
//#############################################################################
//#############################################################################
//#############################################################################
//#############################################################################

void StaggeredMultiScaleT::WriteOutput(void)
{

	cout << "####################################################################" << endl; 
	cout << "####################################################################" << endl; 
	cout << "####################################################################" << endl; 
	cout << "########################## STEP COMPLETE ###########################" << endl; 
	cout << "####################################################################" << endl; 
	cout << "####################################################################" << endl; 
	cout << "####################################################################" << endl; 

	bStep_Complete=1;
	RHSDriver();
	bStep_Complete=0;

	/* my output set */
	const OutputSetT& output_set = ElementSupport().OutputSet(fOutputID);
	
	/* my nodes used */
	const iArrayT& nodes_used = output_set.NodesUsed();

	/* smooth stresses to nodes */
	int n_stress = dSymMatrixT::NumValues(NumSD());
	ElementSupport().ResetAverage(n_stress);
	dArray2DT out_variable_all;
	dSymMatrixT out_variable;
	dArray2DT nd_stress(NumElementNodes(), n_stress);
	Top();
	while (NextElement())
	{
		/* extrapolate */
		nd_stress = 0.0;
		out_variable_all.Set(fNumIP, n_stress, fIPVariable(CurrElementNumber()));
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			out_variable.Set(NumSD(), out_variable_all(fShapes->CurrIP()));
			fShapes->Extrapolate(out_variable, nd_stress);
		}
	
	/* accumulate - extrapolation done from ip's to corners => X nodes  */
		ElementSupport().AssembleAverage(CurrentElement().NodesX(), nd_stress);
	}

	/* get nodally averaged values */
	dArray2DT extrap_values;
	ElementSupport().OutputUsedAverage(extrap_values);

	/* temp space for group displacements */
	int num_node_output = fCoarse.NumDOF() + fFine.NumDOF() + n_stress;
	dArray2DT n_values(nodes_used.Length(), num_node_output);

	/* collect nodal values */
	const dArray2DT& fUa = fFine[0];
	const dArray2DT& fUb = fCoarse[0];
	for (int i = 0; i < nodes_used.Length(); i++)
	{
		int node = nodes_used[i];
		double* row = n_values(i);
		for (int j = 0; j < fUa.MinorDim(); j++) {
			*row++ = fUa(node,j);
			displacements_file <<  fUa(node,j) <<" "; 
		}
		displacements_file <<"   "; 

		for (int j = 0; j < fUb.MinorDim(); j++) {
			*row++ = fUb(node,j);
			displacements_file <<  fUb(node,j) <<" "; 
		}
		displacements_file <<"   "; 

		for (int j = 0; j < fUb.MinorDim(); j++) 
			displacements_file << fUa(node,j) + fUb(node,j) <<" "; 
		
		double* p_stress = extrap_values(i);
		for (int j = 0; j < n_stress; j++)
			*row++ = p_stress[j];

		displacements_file <<  endl;
	}


	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, fIPVariable);

	//--------------------- Write plot file

	int k;
	double epsilon,Lf,Lo=8.0; // Lo = 1.0;  Lo = 1.0e-6

	if (bLog_Strain) {
		Lf = x_top - x_bot;
		epsilon = log ( Lf/Lo ); // epsilon = e^(Lf/Lo) but Lo=1 for unit cube
		var_plot_file << epsilon << " ";
		cout << "x_top= "<<x_top<<" x_bot= "<<x_bot<<" Lf= "<<Lf<<" epsilon = Log(Lf/Lo) = "<<epsilon<< endl;
	}

	for (k=0; k<num_tensors_to_render; k++) 
		var_plot_file << Render_Tensor[Elmt2Write][k][ElmtIP2Write]( component_i, component_j ) << " ";
	
	for (k=0; k<num_scalars_to_render; k++) 
		var_plot_file << Render_Scalar[Elmt2Write][k][ElmtIP2Write] << " ";

#if 0
	cout << "FLAG 3 "<<endl;
	cout << "Node Forces "<< fForces_at_Node << endl;

	if  (iDesired_Force_Direction > -1 && 	iDesired_Force_Node_Num > -1)
		var_plot_file << fForces_at_Node[iDesired_Force_Direction] << " "; 
#endif

	var_plot_file <<endl; 

	//--------------------- Rendering access (here for now)

	double diff = render_time-time;
	diff *= diff; // square to get smaller and pos
	double tiny = 0.00000001;
	int sflag = (diff < tiny) ? 1 : 0;

	if (render_switch==1 && sflag==1 )  
		RenderOutput();

		//-------------------- End Rendering
	
}	

//---------------------------------------------------------------------

void StaggeredMultiScaleT::RenderOutput(void)
{
	/* my output set */
	const OutputSetT& output_set = ElementSupport().OutputSet(fOutputID);
	
	/* my nodes used */
	const iArrayT& nodes_used = output_set.NodesUsed();

	int e,a,i;

	//------- Gather Geometry Data -----------

	Top();
	while (NextElement()) {

		e = CurrElementNumber();
		const iArrayT& node = CurrentElement().NodesX(); // global node numbers start at 0 !
	 	SetLocalX(fInitCoords); 

  	for (a=0; a<n_en; a++) {
			Geometry.Element_Set[e_set].IEN ( e,a ) = node[a];
   		for (i=0; i<n_sd; i++)  
				Geometry.Xo ( node[a],i ) = fInitCoords ( a,i ); 
		}
	}

	//------- Smooth Stresses to Nodes -------
		/* smooth stresses to nodes */
	int n_stress = dSymMatrixT::NumValues(NumSD());
	ElementSupport().ResetAverage(n_stress);
	dArray2DT out_variable_all;
	dSymMatrixT out_variable;
	dArray2DT nd_stress(NumElementNodes(), n_stress);
	
	Top();
	while (NextElement())
	{
		/* extrapolate */
		nd_stress = 0.0;
		out_variable_all.Set(fNumIP, n_stress, fIPVariable(CurrElementNumber()));
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			out_variable.Set(NumSD(), out_variable_all(fShapes->CurrIP()));
			fShapes->Extrapolate(out_variable, nd_stress);
		}
	
		/* accumulate - extrapolation done from ip's to corners => X nodes  */
		ElementSupport().AssembleAverage(CurrentElement().NodesX(), nd_stress);
	}

	/* get nodally averaged values */
	dArray2DT extrap_values;
	ElementSupport().OutputUsedAverage(extrap_values);

#if RENDER

	Render_Boss.time = time;
	Render_Boss.Assemble_Facets 	( Geometry );
	Render_Boss.Color_to_T ( extrap_values );
	Render_Boss.Render ( );

#endif

}

//---------------------------------------------------------------------

void 	StaggeredMultiScaleT::Init_Render ( void )
{
	e_set=0; 
	int n_es=1; 

	Geometry.n_np = n_np;
	Geometry.n_sd = n_sd;
	Geometry.Xo.Dimension ( n_np,n_sd );
	Geometry.n_es = n_es; 
	Geometry.Element_Set.Dimension (n_es);
	Geometry.Element_Set[e_set].n_el = n_el; 
	Geometry.Element_Set[e_set].n_en = n_en; 
	Geometry.Element_Set[e_set].IEN.Dimension ( n_el, n_en ); 

	Geometry.Read_Surface_Data ( surface_file_name,e_set );

	if ( n_sd == 2 )  
		Geometry.Element_Set[e_set].element_type = ContinuumT::kQuad; 
	if ( n_sd == 3 ) 	 
		Geometry.Element_Set[e_set].element_type = ContinuumT::kHex; 
		
#if  RENDER 

	Render_Boss.Read_Render_Settings 	( render_settings_file_name );
	Render_Boss.active_field_component = 1; // Default is 22 component
	//Render_Boss.Print_Render_Settings ( );

	// The following fails during destruction because Map isa nMatrixT <int> that has
	// no inate destructor.  Answer, change it to an iArray2DT which does.
	//
	//FEA_Data_ProcessorT Data_ProX;
	//Data_ProX.Form_Order_Reduction_Map( n_sd );
	//Render_Boss.active_field_component = Data_ProX.Map ( component_i, component_j );
	//cout << "Flag 3 Map(i,j) = " << Data_ProX.Map ( component_i, component_j ) << "\n"; <-- works fine

#endif

}

//---------------------------------------------------------------------

void 	StaggeredMultiScaleT::Get_Fext_I ( dArrayT &fFext_1 )
{
	fFext_1 = 0.0;
}


//#############################################################################
//#############################################################################
//#############################################################################
//#############################################################################
//###################### Actual Solver Routines Below  ########################
//#############################################################################
//#############################################################################
//#############################################################################
//#############################################################################
	
/*************************************************************************
 * Private
 *************************************************************************/

/* form group contribution to the stiffness matrix and RHS */
void StaggeredMultiScaleT::RHSDriver_staggered(void)
{
	const char caller[] = "StaggeredMultiScaleT::RHSDriver_staggered";
	if (fCoarse.Group() == fFine.Group())
		ExceptionT::GeneralFail(caller, "coarse and fine group must be different: %d == %d",
			fCoarse.Group(), fFine.Group());

	int curr_group = ElementSupport().CurrentGroup();

	/* stress output work space */
	int n_stress = dSymMatrixT::NumValues(NumSD());
	dArray2DT   out_variable_all;
	dSymMatrixT out_variable;

	/* time Step Increment */
	double delta_t = ElementSupport().TimeStep();
	time = ElementSupport().Time();
	step_number = ElementSupport().StepNumber();
	iArrayT fine_eq;

	//cout <<" s= " << render_switch <<"; t= "<< time << "; rt= " << render_time << "\n";
 
 	/* has (coarse scale) body forces */
	int formBody = 0;
	if (fBodySchedule && fBody.Magnitude() > kSmall)
		formBody = 1;

	/* loop over elements */
	int e,v,l;
	Top();
	while (NextElement())
	{
		e = CurrElementNumber();

		SetLocalU (ua);			 SetLocalU (ua_n);
		SetLocalU (ub);			 SetLocalU (ub_n);

		del_ua.DiffOf (ua, ua_n);
		del_ub.DiffOf (ub, ub_n);

	 	SetLocalX(fInitCoords); 
		fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, ua, 1.0, ub); 
		fShapes->SetDerivatives(); 

		if (bStep_Complete) { //-- Done iterating, get result data 
			if (bLog_Strain) {	//-- For calculation of Lf : epsilon = e^(Lf/Lo) 
				if ( e == cube_top_elmt ) 
					x_top = fCurrCoords ( cube_top_elmt_top_local_node, 1 ); // 1 is for the 2 direction
				if ( e == cube_bottom_elmt ) 
					x_bot = fCurrCoords ( cube_bottom_elmt_bottom_local_node, 1 ); // 1 is for the 2 direction
			}
		}
		
		/** repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradients 		( fShapes, 	ua, ua_n, fGRAD_ua, fGRAD_ua_n );
		Convert.Gradients 		( fShapes, 	ub, ub_n, fGRAD_ub, fGRAD_ub_n );
		Convert.Shapes				(	fShapes, 	fFEA_Shapes );
		Convert.Displacements	(	del_ua, 	del_ua_vec  );
		Convert.Displacements	(	del_ub, 	del_ub_vec  );

		Convert.Na						(	n_en, fShapes, 	fFEA_Shapes );

		/* Construct data used in BOTH FineScaleT and CoarseScaleT (F,Fa,Fb,grad_ua,...etc.)
		 * 	Presently, Tahoe cannot exploit this fact.  n and np1 are calculated for coarse field, then
		 * 	calculated again for fine field -- this is a waste and defeats the putpose of VMS_VariableT. 
		 *  Note: n is last time step (known data), no subscript,np1 or (n+1) is the 
		 *  next time step (what were solving for)   */

		VMS_VariableT np1(	fGRAD_ua, 	fGRAD_ub 	 ); // Many variables at time-step n+1
		VMS_VariableT   n(	fGRAD_ua_n, fGRAD_ub_n );	// Many variables at time-step n
		
		/* which field */
	  //SolverGroup 1 (gets field 2) <-- ub (obtained by a rearranged Equation I)
		if ( curr_group == fCoarse.Group() || (bStep_Complete && render_variable_group==0) )	
		{

			if (bStep_Complete) { //-- Done iterating, get result data from converged upon displacements 

				fEquation_1 -> Construct ( fFEA_Shapes, fCoarseMaterial, np1, n, step_number, delta_t );

				for (v=0; v<num_tensors_to_render; v++ ) 
					fEquation_1 -> Get ( Render_Tensor_Names[v], Render_Tensor[e][v] ); 
				for (v=0; v<num_scalars_to_render; v++ ) 
					fEquation_1 -> Get ( Render_Scalar_Names[v], Render_Scalar[e][v] ); 

				//-- Store/Register data in classic tahoe manner 
				out_variable_all.Set(fNumIP, n_stress, fIPVariable(CurrElementNumber()));
				for (l=0; l < fNumIP; l++) {
					out_variable.Set(NumSD(), out_variable_all(l));
					out_variable.FromMatrix(Render_Tensor[e][0][l]);
				} 
			}
			else { //-- Still Iterating

				/** Compute N-R matrix equations */
				fEquation_1 -> Construct ( fFEA_Shapes, fCoarseMaterial, np1, n, step_number, delta_t );
				fEquation_1 -> Form_LHS_Ka_Kb ( fKa_1, fKb_1 );
				fEquation_1 -> Form_RHS_F_int ( fFint_1 );

				/** Set coarse LHS */
				fLHS = fKb_1;

				/** Compute coarse RHS (or Fint_bar_2 in FAXed notes) */
				fKa_1.Multx ( del_ua_vec, fRHS );
				fRHS += fFint_1; 
				fRHS *= -1.0; 

				/** Compute Traction B.C. and Body Forces */
				Get_Fext_I ( fFext_1 );
				fRHS += fFext_1;
				
				/* add body forces */
				if (formBody) {
//					double density = fCoarseMaterial->Retrieve(Iso_MatlT::kDensity);
					double density = 1.0;
					DDub = 0.0;
					AddBodyForce(DDub);
					FormMa(kConsistentMass, -density, &DDub, NULL);				
				}
			
				/* add to global equations */
				ElementSupport().AssembleLHS	( fCoarse.Group(), fLHS, CurrentElement().Equations() );
				ElementSupport().AssembleRHS 	( fCoarse.Group(), fRHS, CurrentElement().Equations() );
			}
		}

		// SolverGroup 2 (gets field 1) <-- ua (obtained by a rearranged Equation 2)
		else if (curr_group == fFine.Group() || (bStep_Complete && render_variable_group==1) )	
		{

			if (bStep_Complete) { //-- Done iterating, get result data from converged upon displacements 

				fEquation_2 -> Construct ( fFEA_Shapes, fFineMaterial, np1, n, step_number, delta_t ); 

				for (v=0; v<num_tensors_to_render; v++ )  
					fEquation_2 -> Get ( Render_Tensor_Names[v], Render_Tensor[e][v] ); 
				for (v=0; v<num_scalars_to_render; v++ ) 
					fEquation_2 -> Get ( Render_Scalar_Names[v], Render_Scalar[e][v] ); 

				//-- Store/Register data in classic tahoe manner 
				out_variable_all.Set(fNumIP, n_stress, fIPVariable(CurrElementNumber()));
				for (l=0; l < fNumIP; l++) {
					out_variable.Set(NumSD(), out_variable_all(l));
					out_variable.FromMatrix(Render_Tensor[e][0][l]);
				} 
			}
			else { //-- Still Iterating

				/** Compute N-R matrix equations */
				fEquation_2 -> Construct ( fFEA_Shapes, fFineMaterial, np1, n, step_number, delta_t, FEA::kBackward_Euler );
				fEquation_2 -> Form_LHS_Ka_Kb ( fKa_2, 	fKb_2 );
				fEquation_2 -> Form_RHS_F_int ( fR_2 );

				/** Set LHS */
				fLHS = fKa_2;	
		
				/** Compute fine RHS (or Fint_bar_2 in FAXed notes)  */
				fKb_2.Multx ( del_ub_vec, fRHS );
				fRHS += fR_2; 
				fRHS *= -1.0; 
		
				/* fine scale equation numbers */
				fEqnos_fine.RowAlias ( CurrElementNumber(), fine_eq );

				/* add to global equations */
				ElementSupport().AssembleLHS ( fFine.Group(), fLHS, fine_eq );
				ElementSupport().AssembleRHS ( fFine.Group(), fRHS, fine_eq );
			}

		}
		else ExceptionT::GeneralFail(caller);
	}
}

//---------------------------------------------------------------------
/* form group contribution to the stiffness matrix and RHS */
void StaggeredMultiScaleT::RHSDriver_monolithic(void)
{
	const char caller[] = "StaggeredMultiScaleT::RHSDriver_monolithic";
	if (fCoarse.Group() != fFine.Group())
		ExceptionT::GeneralFail(caller, "coarse and fine group must be the same: %d != %d",
			fCoarse.Group(), fFine.Group());

	int curr_group = ElementSupport().CurrentGroup();

	/* stress output work space */
	int n_stress = dSymMatrixT::NumValues(NumSD());
	dArray2DT   out_variable_all;
	dSymMatrixT out_variable;

	/* time Step Increment */
	double delta_t = ElementSupport().TimeStep();
	time = ElementSupport().Time();
	step_number = ElementSupport().StepNumber();
	iArrayT coarse_eq, fine_eq;

	//cout <<" s= " << render_switch <<"; t= "<< time << "; rt= " << render_time << "\n";
 
 	/* has (coarse scale) body forces */
	int formBody = 0;
	if (fBodySchedule && fBody.Magnitude() > kSmall)
		formBody = 1;

	/* loop over elements */
	int e,v,l;
	Top();
	while (NextElement())
	{
		e = CurrElementNumber();

		SetLocalU (ua);			 SetLocalU (ua_n);
		SetLocalU (ub);			 SetLocalU (ub_n);

		del_ua.DiffOf (ua, ua_n);
		del_ub.DiffOf (ub, ub_n);

	 	SetLocalX(fInitCoords); 
		fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, ua, 1.0, ub); 
		fShapes->SetDerivatives(); 

		if (bStep_Complete) { //-- Done iterating, get result data 
			if (bLog_Strain) {	//-- For calculation of Lf : epsilon = e^(Lf/Lo) 
				if ( e == cube_top_elmt ) 
					x_top = fCurrCoords ( cube_top_elmt_top_local_node, 1 ); // 1 is for the 2 direction
				if ( e == cube_bottom_elmt ) 
					x_bot = fCurrCoords ( cube_bottom_elmt_bottom_local_node, 1 ); // 1 is for the 2 direction
			}
		}
		
		/* repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradients 		( fShapes, 	ua, ua_n, fGRAD_ua, fGRAD_ua_n );
		Convert.Gradients 		( fShapes, 	ub, ub_n, fGRAD_ub, fGRAD_ub_n );
		Convert.Shapes				(	fShapes, 	fFEA_Shapes );
		Convert.Displacements	(	del_ua, 	del_ua_vec  );
		Convert.Displacements	(	del_ub, 	del_ub_vec  );
		Convert.Na(	n_en, fShapes, 	fFEA_Shapes );

		/* Construct data used in BOTH FineScaleT and CoarseScaleT (F,Fa,Fb,grad_ua,...etc.)
		 * 	Presently, Tahoe cannot exploit this fact.  n and np1 are calculated for coarse field, then
		 * 	calculated again for fine field -- this is a waste and defeats the putpose of VMS_VariableT. 
		 *  Note: n is last time step (known data), no subscript,np1 or (n+1) is the 
		 *  next time step (what were solving for)   */
		VMS_VariableT np1(	fGRAD_ua, 	fGRAD_ub 	 ); // Many variables at time-step n+1
		VMS_VariableT   n(	fGRAD_ua_n, fGRAD_ub_n );	// Many variables at time-step n

		if (bStep_Complete) { //-- Done iterating, get result data from converged upon displacements 

			if (render_variable_group == 0)
			{
				fEquation_1 -> Construct ( fFEA_Shapes, fCoarseMaterial, np1, n, step_number, delta_t );

				for (v=0; v<num_tensors_to_render; v++ ) 
					fEquation_1 -> Get ( Render_Tensor_Names[v], Render_Tensor[e][v] ); 
				for (v=0; v<num_scalars_to_render; v++ ) 
					fEquation_1 -> Get ( Render_Scalar_Names[v], Render_Scalar[e][v] ); 

				//-- Store/Register data in classic tahoe manner 
				out_variable_all.Set(fNumIP, n_stress, fIPVariable(CurrElementNumber()));
				for (l=0; l < fNumIP; l++) {
					out_variable.Set(NumSD(), out_variable_all(l));
					out_variable.FromMatrix(Render_Tensor[e][0][l]);
				}
			}
			else if (render_variable_group == 1)
			{
				fEquation_2 -> Construct ( fFEA_Shapes, fFineMaterial, np1, n, step_number, delta_t );

				for (v=0; v<num_tensors_to_render; v++ )  
					fEquation_2 -> Get ( Render_Tensor_Names[v], Render_Tensor[e][v] ); 
				for (v=0; v<num_scalars_to_render; v++ ) 
					fEquation_2 -> Get ( Render_Scalar_Names[v], Render_Scalar[e][v] ); 

				//-- Store/Register data in classic tahoe manner 
				out_variable_all.Set(fNumIP, n_stress, fIPVariable(CurrElementNumber()));
				for (l=0; l < fNumIP; l++) {
					out_variable.Set(NumSD(), out_variable_all(l));
					out_variable.FromMatrix(Render_Tensor[e][0][l]);
				}
			}
			else ExceptionT::GeneralFail(caller, "inrecognized render group %d", render_variable_group);
		}
		else { //-- Still Iterating

			/* residual and tangent for coarse scale */
			fEquation_1 -> Construct ( fFEA_Shapes, fCoarseMaterial, np1, n, step_number, delta_t );
			fEquation_1 -> Form_LHS_Ka_Kb ( fKa_1, fKb_1 );
			fEquation_1 -> Form_RHS_F_int ( fFint_1 );
			fFint_1 *= -1.0;

			/* add body force */
			if (formBody) {
//				double density = fCoarseMaterial->Retrieve(Iso_MatlT::kDensity);
				double density = 1.0;
				DDub = 0.0;
				AddBodyForce(DDub);
				
				/* add body force to fRHS */
				fRHS = 0.0;
				FormMa(kConsistentMass, -density, &DDub, NULL);
				fFint_1 += fRHS;
			}

			/* residual and tangent for fine scale */
			fEquation_2 -> Construct ( fFEA_Shapes, fFineMaterial, np1, n, step_number, delta_t, FEA::kBackward_Euler );
			fEquation_2 -> Form_LHS_Ka_Kb ( fKa_2, 	fKb_2 );
			fEquation_2 -> Form_RHS_F_int ( fR_2 );
			fR_2 *= -1.0;

			/* equations numbers */
			const iArrayT& all_eq = CurrentElement().Equations();
			coarse_eq.Alias(fFint_1.Length(), all_eq.Pointer());
			fine_eq.Alias(fR_2.Length(), all_eq.Pointer(fFint_1.Length()));

			/* assemble residuals */
			ElementSupport().AssembleRHS(curr_group, fFint_1, coarse_eq);
			ElementSupport().AssembleRHS(curr_group, fR_2, fine_eq);

			/* assemble components of the tangent */
			ElementSupport().AssembleLHS(curr_group, fKb_1, coarse_eq);
			ElementSupport().AssembleLHS(curr_group, fKa_2, fine_eq);
			ElementSupport().AssembleLHS(curr_group, fKa_1, coarse_eq, fine_eq);
			ElementSupport().AssembleLHS(curr_group, fKb_2, fine_eq, coarse_eq);
		}
	}
}

