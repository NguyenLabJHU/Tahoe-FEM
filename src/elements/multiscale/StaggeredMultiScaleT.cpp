/* $Id: StaggeredMultiScaleT.cpp,v 1.24 2003-03-07 22:23:57 creigh Exp $ */
//DEVELOPMENT
#include "StaggeredMultiScaleT.h"

#include "ShapeFunctionT.h"
#include "Traction_CardT.h"
#include "ifstreamT.h"

#include "VMF_Virtual_Work_EqT.h"
#include "VMS_BCJT.h"
#include "E_Pr_MatlT.h"
#include "Iso_MatlT.h"
#include "BCJ_MatlT.h"
#include "OutputSetT.h"

#include "ofstreamT.h"

using namespace Tahoe;

/* parameters */
int knum_d_state = 1; // double's needed per ip
int knum_i_state = 1; // int's needed per ip

//---------------------------------------------------------------------

/* constructor */
StaggeredMultiScaleT::StaggeredMultiScaleT(const ElementSupportT& support, const FieldT& coarse, const FieldT& fine):
	ElementBaseT(support, coarse), //pass the coarse scale field to the base class
	ua(LocalArrayT::kDisp),
	ua_n(LocalArrayT::kLastDisp),
	ub(LocalArrayT::kDisp),
	ub_n(LocalArrayT::kLastDisp),
	fInitCoords(LocalArrayT::kInitCoords),
	fCurrCoords(LocalArrayT::kCurrCoords),
	fBodySchedule(NULL),
	fBody(NumDOF()),
	fTractionBCSet(0),
	fCoarse(coarse),
	fFine(fine),
	fEquation_I(NULL),
	fEquation_II(NULL)
{
	/* check - some code below assumes that both fields have the
	 * same dimension. TEMP? */
	if (fCoarse.NumDOF() != fFine.NumDOF()) throw ExceptionT::kBadInputValue;

	/* read parameters from input */
	ifstreamT& in = ElementSupport().Input();
	in >> fGeometryCode; //TEMP - should actually come from the geometry database
	in >> fNumIP;

	in >> render_switch; 
	in >> render_settings_file_name; 
	in >> surface_file_name; 
	in >> render_time; 
	in >> render_variable; // code for variable to render (ex. Fa)
	in >> component_ij; // ex. {0} --> {11} component of render variable 
	in >> render_variable_group; 
	in >> render_variable_order;
	in >> render_displ;

	render_data_stored = 0; 

	fMaterial_Data.Dimension ( kNUM_FMAT_TERMS );
	iMaterial_Data.Dimension ( kNUM_IMAT_TERMS );

	in >> fMaterial_Data[k__E];
	in >> fMaterial_Data[k__Pr];
	in >> fMaterial_Data[k__f];
	in >> fMaterial_Data[k__V];
	in >> fMaterial_Data[k__Y];

	in >> iMaterial_Data[k__BS_Type];
	in >> fMaterial_Data[k__c_zeta];
	in >> fMaterial_Data[k__l];

	in >> iMaterial_Data[k__IH_Type];
	in >> fMaterial_Data[k__K];
	in >> fMaterial_Data[k__H];

	/* allocate the global stack object (once) */
	extern FEA_StackT* fStack;
	if (!fStack) fStack = new FEA_StackT;
}

//---------------------------------------------------------------------

/* destructor */
StaggeredMultiScaleT::~StaggeredMultiScaleT(void) 
{  
	delete fShapes;
	delete fEquation_I; 
	delete fEquation_II; 
	delete fCoarseMaterial; 
	delete fFineMaterial;

	/* free the global stack object (once) */
	extern FEA_StackT* fStack;
	if (fStack) {
		delete fStack;
		fStack = NULL;
	}
}

//---------------------------------------------------------------------

void StaggeredMultiScaleT::Initialize(void)
{
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
		fElementCards[i].Set(fiState.MinorDim(), fiState(i),
		                     fdState.MinorDim(), fdState(i));
		                     
	/* construct the black boxs */  

	Select_Equations ( CoarseScaleT::kVMF_Virtual_Work_Eq,	FineScaleT::kVMS_BCJ );
	fEquation_II -> Initialize ( n_ip, n_sd, n_en, ElementSupport().StepNumber() );
	//step_number_last_iter = 0; 
	//step_number_last_iter = ElementSupport().StepNumber();  // This may crash or not work

	/* FEA Allocation */

	fSigma.FEA_Dimension 			( fNumIP,n_sd,n_sd );

	fGRAD_ua.FEA_Dimension 		( fNumIP,n_sd,n_sd );
	fGRAD_ub.FEA_Dimension 		( fNumIP,n_sd,n_sd );
	fGRAD_ua_n.FEA_Dimension 	( fNumIP,n_sd,n_sd );
	fGRAD_ub_n.FEA_Dimension 	( fNumIP,n_sd,n_sd );

	fKa_I.Dimension 	( n_en_x_n_df, n_en_x_n_df );
	fKb_I.Dimension 	( n_en_x_n_df, n_en_x_n_df );
	fKa_II.Dimension 	( n_en_x_n_df, n_en_x_n_df );
	fKb_II.Dimension 	( n_en_x_n_df, n_en_x_n_df );

	fFint_I.Dimension 	( n_en_x_n_df );
	fFext_I.Dimension 	( n_en_x_n_df );
	fFint_II.Dimension 	( n_en_x_n_df );


	fFEA_Shapes.Construct	( fNumIP,n_sd,n_en 			);

	/* storage for integration point stresses */
	fIPStress.Dimension (n_el, fNumIP*dSymMatrixT::NumValues(n_sd));
	fIPStress = 0.0;

	if (render_switch)
		Init_Render();

}

//---------------------------------------------------------------------


/* form group contribution to the stiffness matrix and RHS */
void StaggeredMultiScaleT::RHSDriver(void)	// LHS too!	
{
 
	int curr_group = ElementSupport().CurrentGroup();

	/* stress output work space */
	int n_stress = dSymMatrixT::NumValues(NumSD());
	dArray2DT   ip_stress_all;
	dSymMatrixT ip_stress;

	/** Time Step Increment */
	double delta_t = ElementSupport().TimeStep();
						time = ElementSupport().Time();
						step_number = ElementSupport().StepNumber();

	iArrayT fine_eq;

	//cout <<" s= " << render_switch <<"; t= "<< time << "; rt= " << render_time << "\n";
 
	/* loop over elements */
	int e=0;
	Top();
	while (NextElement())
	{
		e++;
		SetLocalU (ua);			 SetLocalU (ua_n);
		SetLocalU (ub);			 SetLocalU (ub_n);

		del_ua.DiffOf (ua, ua_n);
		del_ub.DiffOf (ub, ub_n);

	 	SetLocalX(fInitCoords); // dNdX
		fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, ua, 1.0, ub); 
		fShapes->SetDerivatives(); 
	
		/** repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradiants 		( fShapes, 	ua, ua_n, fGRAD_ua, fGRAD_ua_n );
		Convert.Gradiants 		( fShapes, 	ub, ub_n, fGRAD_ub, fGRAD_ub_n );
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
		
		/* which field */
		if (curr_group == fCoarse.Group())  //SolverGroup 1 (gets field 2) <-- ub (obtained by a rearranged Equation I)
		{
			/** Compute N-R matrix equations */

			fEquation_I -> Construct ( fFEA_Shapes, fCoarseMaterial, np1, n, FEA::kBackward_Euler );
			fEquation_I -> Form_LHS_Ka_Kb ( fKa_I, fKb_I );
			fEquation_I -> Form_RHS_F_int ( fFint_I );

			/* copy/store stresses */
			ip_stress_all.Set(fNumIP, n_stress, fIPStress(CurrElementNumber()));
			fEquation_I -> Get (VMF_Virtual_Work_EqT::kSigma, fSigma, 2); // 2 because Sigma is 2nd order tensor
			for (int i = 0; i < fNumIP; i++)
			{
				ip_stress.Set(NumSD(), ip_stress_all(i));
				ip_stress.FromMatrix(fSigma[i]);
			} 

			/** Set coarse LHS */
			fLHS = fKb_I;

			/** Compute coarse RHS (or Fint_bar_II in FAXed notes) */
			fKa_I.Multx ( del_ua_vec, fRHS );
			fRHS += fFint_I; 
			fRHS *= -1.0; 

			/** Compute Traction B.C. and Body Forces */
			Get_Fext_I ( fFext_I );
			fRHS += fFext_I;
			
			/* add to global equations */
			ElementSupport().AssembleLHS	( fCoarse.Group(), fLHS, CurrentElement().Equations() );
			ElementSupport().AssembleRHS 	( fCoarse.Group(), fRHS, CurrentElement().Equations() );
		}
		else if (curr_group == fFine.Group())	// SolverGroup 2 (gets field 1) <-- ua (obtained by a rearranged Equation II)
		{

			/** Compute N-R matrix equations */
			fEquation_II -> Construct ( fFEA_Shapes, fFineMaterial, np1, n, step_number, delta_t, FEA::kBackward_Euler );
			fEquation_II -> Form_LHS_Ka_Kb ( fKa_II, 	fKb_II );
			fEquation_II -> Form_RHS_F_int ( fFint_II );

			/** Set LHS */
			fLHS = fKa_II;	
		
			/** Compute fine RHS (or Fint_bar_II in FAXed notes)  */
			fKb_II.Multx ( del_ub_vec, fRHS );
			fRHS += fFint_II; 
			fRHS *= -1.0; 
		
			/* fine scale equation numbers */
			fEqnos_fine.RowAlias ( CurrElementNumber(), fine_eq );

			/* add to global equations */
			ElementSupport().AssembleLHS ( fFine.Group(), fLHS, fine_eq );
			ElementSupport().AssembleRHS ( fFine.Group(), fRHS, fine_eq );

		}
		else throw ExceptionT::kGeneralFail;
	}

}

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
			fEquation_I 		= new VMF_Virtual_Work_EqT;
			fCoarseMaterial = new Iso_MatlT;
			//fCoarseMaterial -> Assign ( Iso_MatlT::kE, 	 	168.0 	); // 100.0
			//fCoarseMaterial -> Assign ( Iso_MatlT::kPr, 	0.34 		); // .25 
			fCoarseMaterial -> Assign ( Iso_MatlT::kE, 	 	236.4 	); // 100.0
			fCoarseMaterial -> Assign ( Iso_MatlT::kPr, 	0.249 	); // .25 
			fCoarseMaterial -> E_Nu_2_Lamda_Mu	( Iso_MatlT::kE,			Iso_MatlT::kPr,	
																						Iso_MatlT::kLamda, 	Iso_MatlT::kMu 	);
			break;

		case CoarseScaleT::kLDV :
			//fEquation_I 		= new LDVT; 
			fCoarseMaterial = new Iso_MatlT;
			break;

		case CoarseScaleT::kStraight :
			//fEquation_I 		= new StraightT; 
			fCoarseMaterial = new E_Pr_MatlT;
			break;

		default :
			cout << " StaggeredMultiScaleT::Select_Equations() .. ERROR >> bad iCoarseScale \n";
			break;
	}

	/** Choices for Fine-Scale Equation */

	switch ( iFineScale )	{

		case FineScaleT::kVMS_BCJ :
			fEquation_II 	 = new VMS_BCJT;
			fFineMaterial  = new BCJ_MatlT;							
			fFineMaterial -> Assign (		BCJ_MatlT::kE, 				fMaterial_Data[k__E] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr, 			fMaterial_Data[k__Pr] 		); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 				fMaterial_Data[k__f] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kV, 				fMaterial_Data[k__V] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kY, 				fMaterial_Data[k__Y] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kl, 				fMaterial_Data[k__l] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kc_zeta, 	fMaterial_Data[k__c_zeta] ); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kH, 				fMaterial_Data[k__H] 			); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPlastic_Modulus_K, 	fMaterial_Data[k__K] 	); 	

			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE,			BCJ_MatlT::kPr,	
																					BCJ_MatlT::kLamda, 	BCJ_MatlT::kMu 	);

			fEquation_II -> Back_Stress_Type 	=  	iMaterial_Data[k__BS_Type]; 	
			fEquation_II -> Iso_Hard_Type			= 	iMaterial_Data[k__IH_Type];

			//fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 			0.0000001 		); 	
			//fFineMaterial -> Assign ( 	BCJ_MatlT::kV, 			258870.0   		); 	
			break;

		case FineScaleT::kVMS_EZ : 
			fEquation_II 	= new VMS_EZT;
			fFineMaterial = new Iso_MatlT; // <-- not used
			break;

#if 0

		case FineScaleT::kVMS_BCJ :
			fEquation_II 	= new VMS_BCJT;
			fFineMaterial = new BCJ_MatlT;																	// Tantalum 
			fFineMaterial -> Assign (		BCJ_MatlT::kE, 			168.0		 		); 	// 1.68e11 (Pa) 
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr, 		0.34 				); 	// .34 
			fFineMaterial -> Assign ( 	BCJ_MatlT::kl, 			0.001 			); 
			fFineMaterial -> Assign ( 	BCJ_MatlT::kc_zeta, 0.001 			); 
			fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 			0.000016 		); 	// 1.6e-5
			fFineMaterial -> Assign ( 	BCJ_MatlT::kV, 			13.95    		); 	// 9.78e6 
			fFineMaterial -> Assign ( 	BCJ_MatlT::kY, 			1.0			 		); 	// 2.59e7
			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE,			BCJ_MatlT::kPr,	
																					BCJ_MatlT::kLamda, 	BCJ_MatlT::kMu 	);
			break;

		case FineScaleT::kVMS_EZ : 
			fEquation_II 	= new VMS_EZT;
			fFineMaterial = new Iso_MatlT; // <-- not used
			break;


		case FineScaleT::kVMS_EZ2 : 
			fEquation_II 	= new VMS_EZ2T;
			fFineMaterial = new Iso_MatlT; // <-- not used
			break;

		case FineScaleT::kVMS_EZ3 : 
			fEquation_II 	= new VMS_EZ3T;
			fFineMaterial = new BCJ_MatlT; 
			fFineMaterial -> Assign (		BCJ_MatlT::kE, 			168.0 	 		); 	// GPa
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr, 		0.34 				); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 			0.000016 		); 	// 1.6e-5
			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE,			BCJ_MatlT::kPr,	
																					BCJ_MatlT::kLamda, 	BCJ_MatlT::kMu 	);
			break;

		case FineScaleT::kVMS_EZ4 : 
			fEquation_II 	= new VMS_EZ4T;
			fFineMaterial = new BCJ_MatlT; 
			fFineMaterial -> Assign (		BCJ_MatlT::kE, 			168.0 	 		); 	// GPa
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr, 		0.34 				); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 			0.000016 		); 	// 1.6e-5
			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE,			BCJ_MatlT::kPr,	
																					BCJ_MatlT::kLamda, 	BCJ_MatlT::kMu 	);
			break;

		case FineScaleT::kVMS_EZ5 : 
			fEquation_II 	= new VMS_EZ5T;
			fFineMaterial = new BCJ_MatlT; 
			fFineMaterial -> Assign (		BCJ_MatlT::kE, 			168.0 	 		); 	// GPa
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr, 		0.34 				); 	
			fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 			0.000016 		); 	// 1.6e-5
			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE,			BCJ_MatlT::kPr,	
																					BCJ_MatlT::kLamda, 	BCJ_MatlT::kMu 	);
			break;

		case FineScaleT::kPHEN :
			//fEquation_II 	= new PHENT;
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

void StaggeredMultiScaleT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* ElementBaseT handles equation array for the coarse scale */
	if (ElementSupport().CurrentGroup() == fCoarse.Group())
		ElementBaseT::Equations(eq_1,eq_2);
	else if (ElementSupport().CurrentGroup() == fFine.Group())
	{
		/* collect local equation numbers */
		fFine.SetLocalEqnos(fConnectivities, fEqnos_fine);
	
		eq_1.Append(&fEqnos_fine);
	}
	else throw ExceptionT::kGeneralFail;
}

//---------------------------------------------------------------------

/* write element group parameters to out */
void StaggeredMultiScaleT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

	out << " Coarse scale field. . . . . . . . . . . . . . . = \"" << fCoarse.Name() << "\"\n";
	out << " Fine scale field. . . . . . . . . . . . . . . . = \"" << fFine.Name() << "\"\n";
	out << " Element geometry code . . . . . . . . . . . . . = " << fGeometryCode << '\n';
	out << "    eq." << GeometryT::kPoint         << ", point\n";
	out << "    eq." << GeometryT::kLine          << ", line\n";
	out << "    eq." << GeometryT::kQuadrilateral << ", quadrilateral\n";
	out << "    eq." << GeometryT::kTriangle	  << ", triangle\n";
	out << "    eq." << GeometryT::kHexahedron	  << ", hexahedron\n";
	out << "    eq." << GeometryT::kTetrahedron   << ", tetrahedron\n";
	out << " Number of integration points. . . . . . . . . . = " << fNumIP    << '\n';

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

//---------------------------------------------------------------------

/* accumulate the residual force on the specified node */
void StaggeredMultiScaleT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
//not implemented
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
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

//---------------------------------------------------------------------

void StaggeredMultiScaleT::WriteOutput(void)
{

	//--------------------- Rendering access (here for now)
	static int once_flag=0;
	if (render_switch==1 && render_time==time && once_flag==0) { 
		cout << "rt = "<<render_time<<": time = "<<time<<"\n";
		RenderOutput();
		once_flag++;
	}	//-------------------- End Rendering
	
	
	/* my output set */
	const OutputSetT& output_set = ElementSupport().OutputSet(fOutputID);
	
	/* my nodes used */
	const iArrayT& nodes_used = output_set.NodesUsed();

	/* smooth stresses to nodes */
	int n_stress = dSymMatrixT::NumValues(NumSD());
	ElementSupport().ResetAverage(n_stress);
	dArray2DT ip_stress_all;
	dSymMatrixT ip_stress;
	dArray2DT nd_stress(NumElementNodes(), n_stress);
	Top();
	while (NextElement())
	{
		/* extrapolate */
		nd_stress = 0.0;
		ip_stress_all.Set(fNumIP, n_stress, fIPStress(CurrElementNumber()));
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			ip_stress.Set(NumSD(), ip_stress_all(fShapes->CurrIP()));
			fShapes->Extrapolate(ip_stress, nd_stress);
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
	const dArray2DT& ua = fFine[0];
	const dArray2DT& ub = fCoarse[0];
	for (int i = 0; i < nodes_used.Length(); i++)
	{
		int node = nodes_used[i];
		double* row = n_values(i);
		for (int j = 0; j < ua.MinorDim(); j++)
			*row++ = ua(node,j);

		for (int j = 0; j < ub.MinorDim(); j++)
			*row++ = ub(node,j);

		double* p_stress = extrap_values(i);
		for (int j = 0; j < n_stress; j++)
			*row++ = p_stress[j];
	}

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, fIPStress);

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
	dArray2DT ip_stress_all;
	dSymMatrixT ip_stress;
	dArray2DT nd_stress(NumElementNodes(), n_stress);
	
	Top();
	while (NextElement())
	{
		/* extrapolate */
		nd_stress = 0.0;
		ip_stress_all.Set(fNumIP, n_stress, fIPStress(CurrElementNumber()));
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			ip_stress.Set(NumSD(), ip_stress_all(fShapes->CurrIP()));
			fShapes->Extrapolate(ip_stress, nd_stress);
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
	Render_Boss.active_field_component = component_ij;
	//Render_Boss.Print_Render_Settings ( );

#endif

}

//---------------------------------------------------------------------

void 	StaggeredMultiScaleT::Get_Fext_I ( dArrayT &fFext_I )
{
	fFext_I = 0.0;

	//#############################################################################
	//#############################################################################
	//### Code Goes Here to apply Traction b.c. and Body Forces ###################
	//#############################################################################
	//#############################################################################
	
}

