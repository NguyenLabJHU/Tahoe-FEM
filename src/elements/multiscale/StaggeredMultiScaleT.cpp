/* $Id: StaggeredMultiScaleT.cpp,v 1.15 2002-12-17 02:44:34 creigh Exp $ */
//DEVELOPMENT
#include "StaggeredMultiScaleT.h"

#include "eControllerT.h"
#include "ShapeFunctionT.h"
#include "ifstreamT.h"

#include "VMF_Virtual_Work_EqT.h"
#include "VMS_BCJT.h"
#include "E_Pr_MatlT.h"
#include "Iso_MatlT.h"
#include "BCJ_MatlT.h"

#include "ofstreamT.h"

using namespace Tahoe;

/* parameters */
int knum_d_state = 1; // double's needed per ip
int knum_i_state = 1; // int's needed per ip

//---------------------------------------------------------------------

/* constructor */
StaggeredMultiScaleT::StaggeredMultiScaleT(const ElementSupportT& support, 
	const FieldT& coarse, const FieldT& fine):
	ElementBaseT(support, coarse), //pass the coarse scale field to the base class
	ua(LocalArrayT::kDisp),
	ua_n(LocalArrayT::kLastDisp),
	ub(LocalArrayT::kDisp),
	ub_n(LocalArrayT::kLastDisp),
	fInitCoords(LocalArrayT::kInitCoords),
	fCurrCoords(LocalArrayT::kCurrCoords),
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
	
	/* dimensions */
	int n_sd = NumSD();
	int n_df = NumDOF(); 
	int n_en = NumElementNodes();
	int n_en_x_n_df = n_en*n_df;

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
	int n_ip = 1; //TEMP - need to decide where to set the number of integration
	int n_e = NumElements();	
	fdState_new.Dimension(n_e, n_ip*knum_d_state);
	fdState.Dimension(n_e, n_ip*knum_d_state);
	fiState_new.Dimension(n_e, n_ip*knum_i_state);
	fiState.Dimension(n_e, n_ip*knum_i_state);
	
	/* storage for the fine scale equation numbers */
	fEqnos_fine.Dimension(n_e, n_en*n_df);
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

	//Select_Equations ( CoarseScaleT::kVMF_Virtual_Work_Eq,	FineScaleT::kVMS_BCJ ); 
	Select_Equations ( CoarseScaleT::kVMF_Virtual_Work_Eq,	FineScaleT::kVMS_EZ ); 
	
	/* FEA Allocation */

	fSigma.FEA_Dimension 			(fNumIP,n_sd,n_sd);

	fGRAD_ua.FEA_Dimension 		(fNumIP,n_sd,n_sd);
	fGRAD_ub.FEA_Dimension 		(fNumIP,n_sd,n_sd);
	fGRAD_ua_n.FEA_Dimension 	(fNumIP,n_sd,n_sd);
	fGRAD_ub_n.FEA_Dimension 	(fNumIP,n_sd,n_sd);

	fKa_I.Dimension 	( n_en_x_n_df, n_en_x_n_df );
	fKb_I.Dimension 	( n_en_x_n_df, n_en_x_n_df );
	fKa_II.Dimension 	( n_en_x_n_df, n_en_x_n_df );
	fKb_II.Dimension 	( n_en_x_n_df, n_en_x_n_df );

	fFint_I.Dimension 	( n_en_x_n_df );
	fFint_II.Dimension 	( n_en_x_n_df );

	fFEA_Shapes.Construct( fNumIP,n_sd,n_en );

	cout << "############################# Initialize#: \n"; 
}

//---------------------------------------------------------------------


/* form group contribution to the stiffness matrix and RHS */
void StaggeredMultiScaleT::RHSDriver(void)	// LHS too!	
{
 
	int curr_group = ElementSupport().CurrentGroup();

#define OFF      0
#define DISPL    1 
#define COARSE   2 
#define FINE     3 
#define ALL      4 

#define DEBUG COARSE 

#if (DEBUG)
	int debug_iteration=41;
	static int loop_num=0;
	static int current_group=0;
	//static ofstreamT myout("matrix");

	if (curr_group!=current_group) { // new group
		//loop_num =0;
		current_group = curr_group;
	}

  loop_num++;
	//loop_num =  ElementSupportT().IterationNumber();
	
	cout << "############################# ITERATION#: " << loop_num << "\n";
#endif

	/** Time Step Increment */
	double delta_t = ElementSupport().TimeStep();
	iArrayT fine_eq;

	/* loop over elements */

	int e=0;
	Top();
	while (NextElement())
	{
		e++;
		SetLocalU (ua);			 SetLocalU (ua_n);
		SetLocalU (ub);			 SetLocalU (ub_n);

#if (DEBUG==DISPL || DEBUG==ALL) // Debugging Code
		if ( e==1 && loop_num==debug_iteration ) {
			cout << "|||||||||||||||||| DISPL ||||||||||||||||| Elmt number = "<<e<<"\n";
			cout << "ua = \n" << ua << "\n\n";
			cout << "ub = \n" << ub << "\n\n";
			cout << "ua_n = \n" << ua_n << "\n\n";
			cout << "ub_n = \n" << ub_n << "\n\n";
		}
#endif

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

		/** Construct data used in BOTH FineScaleT and CoarseScaleT (F,Fa,Fb,grad_ua,...etc.)
		 * 	Presently, Tahoe cannot exploit this fact.  n and np1 are calculated for coarse field, then
		 * 	calculated again for fine field -- this is a waste and defeats the putpose of VMS_VariableT. 
		 *  Note: n is last time step (known data), no subscript,np1 or (n+1) is the 
		 *  next time step (what were solving for)   */

		VMS_VariableT np1(	fGRAD_ua, 	fGRAD_ub 	 ); // Many variables at time-step n+1
		VMS_VariableT   n(	fGRAD_ua_n, fGRAD_ub_n );	// Many variables at time-step n
		
		/* which field */
		if (curr_group == fCoarse.Group())  // <-- ub (obtained by a rearranged Equation I)
		{
			/** Compute N-R matrix equations */
			fEquation_I -> Construct ( fFEA_Shapes, fCoarseMaterial, np1, n, FEA::kBackward_Euler );
			fEquation_I -> Form_LHS_Ka_Kb ( fKa_I, fKb_I );
			fEquation_I -> Form_RHS_F_int ( fFint_I );

			//fKa_I.Random(1);
			//fKb_I.Random(2);

#if (DEBUG==ALL || DEBUG==COARSE) // Debugging Code

			if ( e==1 && loop_num==debug_iteration ) {
				cout << "|||||||||||||||||| COARSE ||||||||||||||||| Elmt number = "<<e<<"\n";
				cout << "  fKa_I = \n" << fKa_I << "\n\n";
				cout << "  fKb_I = \n" << fKb_I << "\n\n";
				cout << "  fFint_I = \n" << fFint_I << "\n\n";

			  fEquation_I -> Sigma ( fSigma );
				fSigma.Print("Sigma");
			}

#endif

			/** Set coarse LHS */
			fLHS = fKb_I;

			/** Compute coarse RHS (or Fint_bar_II in FAXed notes) */
			fKa_I.Multx ( del_ua_vec, fRHS );
			fRHS += fFint_I; 
			fRHS *= -1.0; 

			/* add to global equations */
			ElementSupport().AssembleLHS	( fCoarse.Group(), fLHS, CurrentElement().Equations() );
			ElementSupport().AssembleRHS 	( fCoarse.Group(), fRHS, CurrentElement().Equations() );
		}
		else if (curr_group == fFine.Group())	// <-- ua (obtained by a rearranged Equation II)
		{

			/** Compute N-R matrix equations */
			fEquation_II -> Construct ( fFEA_Shapes, fFineMaterial, np1, n, delta_t, FEA::kBackward_Euler );
			fEquation_II -> Form_LHS_Ka_Kb ( fKa_II, 	fKb_II );
			fEquation_II -> Form_RHS_F_int ( fFint_II );

#if (DEBUG==ALL || DEBUG==FINE) // Debugging Code

			if (e==1 && loop_num==debug_iteration ) {

			cout << ".................. FINE ................. Elmt number = "<<e<<"\n";
			cout << ">>>>>>>>>>>>> e = "<<e<<"   loop_num =" <<loop_num<<"\n";
			cout << "fKa_II = \n" << fKa_II << "\n";
			cout << "fKb_II = \n" << fKb_II << "\n";

			/* // Use for printing on 8.5 x 11
				FEA_dMatrixT KAII(1, fKa_II.Rows(), fKa_II.Cols() );
				FEA_dMatrixT KBII(1, fKb_II.Rows(), fKb_II.Cols() );

      	KAII[0] = fKa_II;
				KAII.print("Ka_II");

      	KBII[0] = fKb_II;
				KBII.print("Kb_II");

			  myout << fKa_II;
				*/
			}
#endif

			/** Set LHS */
			fLHS = fKa_II;	
		
			/** Compute fine RHS (or Fint_bar_II in FAXed notes) */
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

//--------------------------------------------------------------------

void StaggeredMultiScaleT::Select_Equations (const int &iCoarseScale,const int &iFineScale )
{
	/** Choices for Coarse-Scale Equation */

	switch ( iCoarseScale )	{

		case CoarseScaleT::kVMF_Virtual_Work_Eq :
			fEquation_I 		= new VMF_Virtual_Work_EqT;
			fCoarseMaterial = new Iso_MatlT;
			fCoarseMaterial -> Assign ( Iso_MatlT::kE, 	29000000.0 		);
			fCoarseMaterial -> Assign ( Iso_MatlT::kPr, 	0.30 				); 
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
			fEquation_II 	= new VMS_BCJT;
			fFineMaterial = new BCJ_MatlT;																	// Tantalum 
			fFineMaterial -> Assign (		BCJ_MatlT::kE, 			1.68e+11 		); 	// 1.68e11 (Pa) 
			fFineMaterial -> Assign ( 	BCJ_MatlT::kPr, 		0.34 				); 	// .34 
			fFineMaterial -> Assign ( 	BCJ_MatlT::kl, 			0.001 			); 
			fFineMaterial -> Assign ( 	BCJ_MatlT::kc_zeta, 0.001 			); 
			fFineMaterial -> Assign ( 	BCJ_MatlT::kf, 			1.60e-05 		); 	// 1.6e-5
			fFineMaterial -> Assign ( 	BCJ_MatlT::kV, 			9.78e+06 		); 	// 9.78e6 
			fFineMaterial -> Assign ( 	BCJ_MatlT::kY, 			2.59e+07 		); 	// 2.59e7
			fFineMaterial -> E_Nu_2_Lamda_Mu	( BCJ_MatlT::kE,			BCJ_MatlT::kPr,	
																					BCJ_MatlT::kLamda, 	BCJ_MatlT::kMu 	);
			break;

		case FineScaleT::kVMS_EZ : 
			fEquation_II 	= new VMS_EZT;
			fFineMaterial = new Iso_MatlT; // <-- not used
			break;

		case FineScaleT::kPHEN :
			//fEquation_II 	= new PHENT;
			//fFineMaterial = new VMS_Phen_MaterialT;
			break;

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

void StaggeredMultiScaleT::RegisterOutput(void)
{
	//not implemented
}

//---------------------------------------------------------------------

void StaggeredMultiScaleT::WriteOutput(void)
{
	//not implemented
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


