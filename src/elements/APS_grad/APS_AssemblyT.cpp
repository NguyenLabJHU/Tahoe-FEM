/* $Id: APS_AssemblyT.cpp,v 1.1 2003-07-10 17:14:47 raregue Exp $ */
#include "APS_AssemblyT.h"

#include "ShapeFunctionT.h"
#include "Traction_CardT.h"
#include "ifstreamT.h"

#include "APS_Bal_EqT.h"
#include "APS_BCJT.h"
#include "Shear_MatlT.h"
#include "APS_MatlT.h"
#include "OutputSetT.h"

#include "ofstreamT.h"

using namespace Tahoe;

/* parameters */
int knum_d_state = 1; // double's needed per ip
int knum_i_state = 1; // int's needed per ip

//---------------------------------------------------------------------

/* constructor */
APS_AssemblyT::APS_AssemblyT(const ElementSupportT& support, const FieldT& displ, const FieldT& gammap):
	ElementBaseT(support, displ), //pass the displacement field to the base class
	u(LocalArrayT::kDisp),
	u_n(LocalArrayT::kLastDisp),
	gamma_p(LocalArrayT::kDisp),
	gamma_p_n(LocalArrayT::kLastDisp),
	fInitCoords(LocalArrayT::kInitCoords),
	fCurrCoords(LocalArrayT::kCurrCoords),
	fBodySchedule(NULL),
	fBody(NumDOF()),
	fTractionBCSet(0),
	fDispl(displ),
	fPlast(gammap),
	fKd_I(ElementMatrixT::kNonSymmetric),
	fKeps_I(ElementMatrixT::kNonSymmetric),
	fKd_II(ElementMatrixT::kNonSymmetric),
	fKeps_II(ElementMatrixT::kNonSymmetric),
	fEquation_I(NULL),
	fEquation_II(NULL)
{
	int i;
	/* check - some code below assumes that both fields have the
	 * same dimension. TEMP?  get rid of this!!!!!! */ 
	    if (fDispl.NumDOF() != fPlast.NumDOF()) 
				ExceptionT::BadInputValue("APS_AssemblyT::APS_AssemblyT");

	/* read parameters from input */
	ifstreamT& in = ElementSupport().Input();
	in >> fGeometryCode; //TEMP - should actually come from the geometry database
	in >> fNumIP;

	fMaterial_Data.Dimension ( kNUM_FMAT_TERMS );

	in >> iPlastModelType; 

	//-- Elasticity parameters
	in >> fMaterial_Data[k__mu];

	//-- Plasticity parameters
	in >> fMaterial_Data[k__m_rate];
	in >> fMaterial_Data[k__gamma0_dot];
	in >> fMaterial_Data[k__m1];
	in >> fMaterial_Data[k__m2];

	//-- Backstress Parameters
	in >> fMaterial_Data[k__l];

	//-- Isotropic Hardening Parameters
	in >> fMaterial_Data[k__H];
	
	Echo_Input_Data();
	
	/* allocate the global stack object (once) */
	extern FEA_StackT* fStack;
	if (!fStack) fStack = new FEA_StackT;
}

//---------------------------------------------------------------------

/* destructor */
APS_AssemblyT::~APS_AssemblyT(void) 
{  
	delete fShapes;
	delete fEquation_I; 
	delete fEquation_II; 
	delete fBalLinMomMaterial; 
	delete fPlastMaterial;

	/* free the global stack object (once) */
	extern FEA_StackT* fStack;
	if (fStack) {
		delete fStack;
		fStack = NULL;
	}

	var_plot_file.close(); 
}

//--------------------------------------------------------------------

void APS_AssemblyT::Echo_Input_Data(void) {

	int i;
	cout << "#######################################################" << endl; 
	cout << "############### ECHO APS DATA #########################" << endl; 
	cout << "#######################################################" << endl; 

	//################## material data ##################

	cout << "iPlastModelType " 						<< iPlastModelType 						<< endl; 
	
	//-- Elasticity parameters 
	cout << "fMaterial_Data[k__mu] "  				<< fMaterial_Data[k__mu] 					<< endl;

	//-- Plasticity parameters
	cout << "fMaterial_Data[k__m_rate] " 			<< fMaterial_Data[k__m_rate] 				<< endl;
	cout << "fMaterial_Data[k__gamma0_dot] " 		<< fMaterial_Data[k__gamma0_dot] 			<< endl;
	cout << "fMaterial_Data[k__m1] " 				<< fMaterial_Data[k__m1] 					<< endl;
	cout << "fMaterial_Data[k__m2] " 				<< fMaterial_Data[k__m2] 					<< endl;

	//-- Backstress Parameters
	cout << "fMaterial_Data[k__l] " 				<< fMaterial_Data[k__l]						<< endl;

	//-- Isotropic Hardening Parameters
	cout << "fMaterial_Data[k__H] "					<< fMaterial_Data[k__H]						<< endl;
	
}

void APS_AssemblyT::Initialize(void)
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
	u.Dimension (n_en, n_df);
	u_n.Dimension (n_en, n_df);
	del_u.Dimension (n_en, n_df);
	del_u_vec.Dimension (n_en_x_n_df);
	fDispl.RegisterLocal(u);
	fDispl.RegisterLocal(u_n);

	/* set local arrays for fine scale */
	gamma_p.Dimension (n_en, 2);
	gamma_p_n.Dimension (n_en, 2);
	del_gamma_p.Dimension (n_en, 2);
	int dum = n_en*2;
	del_gamma_p_vec.Dimension (dum);
	fPlast.RegisterLocal(gamma_p);
	fPlast.RegisterLocal(gamma_p_n);

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
	fEqnos_plast.Dimension(n_el, n_en*2);
	fEqnos_plast = -1;

	/* initialize state variables */
	fdState = 0;
	fiState = 0;

	/* set cards to data in array - NOT NEEDED IF YOU'RE NOT
	 * GOING TO USE THE ElementCardT ARRAY? */
	for (int i= 0; i < fElementCards.Length(); i++)
		fElementCards[i].Set(fiState.MinorDim(), fiState(i), fdState.MinorDim(), fdState(i));
		                     
	/* construct the black boxs */  

	Select_Equations ( BalLinMomT::kAPS_Bal_Eq, iPlastModelType );
	fEquation_II -> Initialize ( n_ip, n_sd, n_en, ElementSupport().StepNumber() );
	//step_number_last_iter = 0; 
	//step_number_last_iter = ElementSupport().StepNumber();  // This may crash or not work

	/* FEA Allocation */

	// these dimensions should be different since want to use quadratic interp for u
	// and linear interp for gamma_p
	fKd_I.Dimension 		( n_en_x_n_df, n_en_x_n_df );
	fKeps_I.Dimension 		( n_en_x_n_df, n_en_x_n_df );
	fKd_II.Dimension 		( n_en_x_n_df, n_en_x_n_df );
	fKeps_II.Dimension 		( n_en_x_n_df, n_en_x_n_df );

	fFint_I.Dimension 	( n_en_x_n_df );
	fFext_I.Dimension 	( n_en_x_n_df );
	fFint_II.Dimension 	( n_en_x_n_df );
	fFext_II.Dimension 	( n_en_x_n_df );

	fFEA_Shapes.Construct	( fNumIP,n_sd,n_en );

	/* streams */
	ifstreamT& in  = ElementSupport().Input();
	ofstreamT& out = ElementSupport().Output();

	/* storage for integration point stresses */
	fIPVariable.Dimension (n_el, fNumIP*dSymMatrixT::NumValues(n_sd));
	fIPVariable = 0.0;

	/* allocate storage for nodal forces */
	fForces_at_Node.Dimension ( n_sd );

	/* body force specification */
	fDOFvec.Dimension(n_df);
	EchoBodyForce(in, out);

	/* echo traction B.C.'s */
	EchoTractionBC(in, out);
}

//---------------------------------------------------------------------

void APS_AssemblyT::RHSDriver(void)	
{
	int curr_group = ElementSupport().CurrentGroup();

	/* traction boundary conditions acting on the coarse scale equations */
	if (curr_group == fDispl.Group()) 
		ApplyTractionBC();

	/* choose solution method */
	if (fDispl.Group() == fPlast.Group())
	  RHSDriver_monolithic();
	 else
	  RHSDriver_staggered();
}
//---------------------------------------------------------------------

void APS_AssemblyT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* doing monolithic solution */
	if (fDispl.Group() == fPlast.Group())
	{
		int ndof_plast = fPlast.NumDOF();
		int ndof_displ = fDispl.NumDOF();
		int nen = NumElementNodes();
	
		/* loop over connectivity blocks */
		for (int i = 0; i < fEqnos.Length(); i++)
		{
			/* connectivities */
			const iArray2DT& connects = *(fConnectivities[i]);
			int nel = connects.MajorDim();
		
			/* dimension */
			fEqnos[i].Dimension(nel, nen*(ndof_displ + ndof_plast));
			iArray2DT displ_eq(nel, nen*ndof_plast);
			iArray2DT plast_eq(nel, nen*ndof_plast);
			
			/* get equation numbers */
			fDispl.SetLocalEqnos(connects, displ_eq);
			fPlast.SetLocalEqnos(connects, plast_eq);
			
			/* write into one array */
			fEqnos[i].BlockColumnCopyAt(displ_eq, 0);
			fEqnos[i].BlockColumnCopyAt(plast_eq, displ_eq.MinorDim());

			/* add to list of equation numbers */
			eq_1.Append(&fEqnos[i]);
		}
	
		/* reset pointers to element cards */
		SetElementCards();	
	}
	else
	{
		/* ElementBaseT handles equation array for the coarse scale */
		if (ElementSupport().CurrentGroup() == fDispl.Group())
			ElementBaseT::Equations(eq_1,eq_2);

		/* fine scale equations */
		if (ElementSupport().CurrentGroup() == fPlast.Group())
		{
			/* collect local equation numbers */
			fPlast.SetLocalEqnos(fConnectivities, fEqnos_plast);
		
			eq_1.Append(&fEqnos_plast);
		}
	}
}

#if 0
//---------------------------------------------------------------------

/* form group contribution to the stiffness matrix and RHS */
void APS_AssemblyT::RHSDriver_staggered(void)	// LHS too!	(This was original RHSDriver()
{
 
	int curr_group = ElementSupport().CurrentGroup();

	/* stress output work space */
	int n_stress = dSymMatrixT::NumValues(NumSD());
	dArray2DT   out_variable_all;
	dSymMatrixT out_variable;

	/** Time Step Increment */
	double delta_t = ElementSupport().TimeStep();
	double time = ElementSupport().Time();
	int step_number = ElementSupport().StepNumber();

	iArrayT plast_eq;

	if ( curr_group == fDispl.Group() )
		cout << "############### Displacement Group ###############\n";

	if ( curr_group == fPlast.Group() )
		cout << "############### Plasticity Group ###############\n";
 
	/* loop over elements */
	int e,v,l;
	Top();
	while (NextElement())
	{
		e = CurrElementNumber();

		SetLocalU (u);
		SetLocalU (u_n);
		SetLocalU (gamma_p);
		SetLocalU (gamma_p_n);

		del_u.DiffOf (u, u_n);
		del_gamma_p.DiffOf (gamma_p, gamma_p_n);

	 	SetLocalX(fInitCoords); 
		fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, ua, 1.0, ub); 
		fShapes->SetDerivatives(); 
		
		// ?????
		/** repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradients 		( fShapes, 	u, u_n );
		Convert.Gradients 		( fShapes, 	gamma_p, gamma_p_n );
		Convert.Shapes			(	fShapes, 	fFEA_Shapes );
		Convert.Displacements	(	del_u, 	del_u_vec  );
		Convert.Displacements	(	del_gamma_p, 	del_gamma_p_vec  );
		Convert.Na				(	n_en, fShapes, 	fFEA_Shapes );

		/** Construct data used in BOTH PlastT and BalLinMomT (grad_u and gammap)
		 * 	Presently, Tahoe cannot exploit this fact.  n and np1 are calculated for coarse field, then
		 * 	calculated again for fine field -- this is a waste and defeats the purpose of APS_VariableT. 
		 *  Note: n is last time step (known data), no subscript, np1 or (n+1) is the 
		 *  next time step (what were solving for)   */
		
		/* which field */
	  //SolverGroup 1 (gets field 2) <-- gamma_p (obtained by a rearranged Equation I)
		if ( curr_group == fDispl.Group() )	
		{
			if (bStep_Complete) { 
			// do nothing
			}
			else { //-- Still Iterating

				/** Compute N-R matrix equations */
				fEquation_I -> Construct ( fFEA_Shapes, fBalLinMomMaterial, np1, n, step_number, delta_t );
				fEquation_I -> Form_LHS_Keps_Kd ( fKeps_I, fKd_I );
				fEquation_I -> Form_RHS_F_int ( fFint_I );

				/** Set coarse LHS */
				fLHS = fKd_I;

				/** Compute coarse RHS */
				fKeps_I.Multx ( del_gamma_p_vec, fRHS );
				fRHS += fFint_I; 
				fRHS *= -1.0; 

				/** Compute Traction B.C. and Body Forces */
				Get_Fext_I ( fFext_I );
				fRHS += fFext_I;
			
				/* add to global equations */
				ElementSupport().AssembleLHS	( fDispl.Group(), fLHS, CurrentElement().Equations() );
				ElementSupport().AssembleRHS 	( fDispl.Group(), fRHS, CurrentElement().Equations() );
			}
		}

		// SolverGroup 2 (gets field 1) <-- ua (obtained by a rearranged Equation II)
		else if (curr_group == fPlast.Group() )	
		{

			if (bStep_Complete) { 
			//do nothing
			}
			else { //-- Still Iterating

				/** Compute N-R matrix equations */
				fEquation_II -> Construct ( fFEA_Shapes, fPlastMaterial, np1, n, step_number, delta_t, FEA::kBackward_Euler );
				fEquation_II -> Form_LHS_Keps_Kd ( fKeps_II, 	fKd_II );
				fEquation_II -> Form_RHS_F_int ( fFint_II );

				/** Set LHS */
				fLHS = fKeps_II;	
		
				/** Compute fine RHS (or Fint_bar_II in FAXed notes)  */
				fKd_II.Multx ( del_u_vec, fRHS );
				fRHS += fFint_II; 
				fRHS *= -1.0; 
		
				/* fine scale equation numbers */
				fEqnos_plast.RowAlias ( CurrElementNumber(), plast_eq );

				/* add to global equations */
				ElementSupport().AssembleLHS ( fPlast.Group(), fLHS, plast_eq );
				ElementSupport().AssembleRHS ( fPlast.Group(), fRHS, plast_eq );
			}

		}
		else throw ExceptionT::kGeneralFail;
	}

}
#endif


//---------------------------------------------------------------------

void APS_AssemblyT::LHSDriver(GlobalT::SystemTypeT)
{
  /** Everything done in RHSDriver for efficiency */
	//cout << "############### In LHS Driver ############### \n";

}

//---------------------------------------------------------------------

void APS_AssemblyT::Select_Equations (const int &iBalScale,const int &iPlastScale )
{
	/** Choices for Coarse-Scale Equation */

	switch ( iBalScale )	{

		case BalLinMomT::kAPS_Bal_Eq :
			fEquation_I 	= new APS_Bal_EqT;
			fBalLinMomMaterial = new Shear_MatlT;
			fBalLinMomMaterial -> Assign ( Shear_MatlT::kmu, fMaterial_Data[k__mu] );
			break;

		default :
			cout << " APS_AssemblyT::Select_Equations() .. ERROR >> bad iBalScale \n";
			break;
	}

	/** Choices for Fine-Scale Equation */

	switch ( iPlastScale )	{

		case PlastT::kAPS_BCJ :
			fEquation_II	= new APS_BCJT;
			fPlastMaterial	= new APS_MatlT;		
			fPlastMaterial -> Assign (	APS_MatlT::kMu, 			fMaterial_Data[k__mu] 			); 	
			fPlastMaterial -> Assign ( 	APS_MatlT::km_rate, 		fMaterial_Data[k__m_rate] 		); 	
			fPlastMaterial -> Assign ( 	APS_MatlT::kgamma0_dot, 	fMaterial_Data[k__gamma0_dot] 	); 	
			fPlastMaterial -> Assign ( 	APS_MatlT::km1, 			fMaterial_Data[k__m1] 			); 
			fPlastMaterial -> Assign ( 	APS_MatlT::km2, 			fMaterial_Data[k__m2] 			); 	
			fPlastMaterial -> Assign ( 	APS_MatlT::kl, 				fMaterial_Data[k__l] 			); 	
			fPlastMaterial -> Assign ( 	APS_MatlT::kH, 				fMaterial_Data[k__H] 			); 	
			break;

#endif

		default :
			cout << " APS_AssemblyT::Select_Equations() .. ERROR >> bad iPlastScale \n";
			break;
	}

}

//---------------------------------------------------------------------

/* return true if the element contributes to the solution of the
 * given group. */
bool APS_AssemblyT::InGroup(int group) const
{
	return group == fDispl.Group() ||
	       group == fPlast.Group();
}

//---------------------------------------------------------------------

/* close current time increment */
void APS_AssemblyT::CloseStep(void)
{
	/* inherited */
	ElementBaseT::CloseStep();

	/* store more recently updated values */
	fdState = fdState_new;
	fiState = fiState_new;
}

//---------------------------------------------------------------------

/* write element group parameters to out */
void APS_AssemblyT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

	out << " Displacement field. . . . . . . . . . . . . . . = \"" << fDispl.Name() << "\"\n";
	out << " Plastic gradient field. . . . . . . . . . . . . . . . = \"" << fPlast.Name() << "\"\n";
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

void APS_AssemblyT::SendOutput(int kincode)
{
#pragma unused(kincode)
	//not implemented
}
//---------------------------------------------------------------------

/* form of tangent matrix */
GlobalT::SystemTypeT APS_AssemblyT::TangentType(void) const
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
void APS_AssemblyT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	const char caller[] = "APS_AssemblyT::AddNodalForce";

	/* coarse, fine, or neither */
	bool is_coarse = false;
	dArrayT* element_force = NULL;
	int num_force = 0;
	if (field.Name() == fDispl.Name()) {
		is_coarse = true;
		element_force = &fFint_I;
		num_force = fDispl.NumDOF();
	}
	else if (field.Name() == fPlast.Name()) {
		is_coarse = false;
		element_force = &fFint_II;
		num_force = fPlast.NumDOF();
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
	int e,v,l;
	Top();
	while (NextElement())
	{
		int nodeposition;
		const iArrayT& nodes_u = CurrentElement().NodesU();
		if (nodes_u.HasValue(node, nodeposition))
		{
			e = CurrElementNumber();

		SetLocalU (u);
		SetLocalU (u_n);
		SetLocalU (gamma_p);
		SetLocalU (gamma_p_n);

		del_u.DiffOf (u, u_n);
		del_gamma_p.DiffOf (gamma_p, gamma_p_n);

	 	SetLocalX(fInitCoords); 
	 	// coordinates do not change for anti-plane shear
		fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, u); 
		fShapes->SetDerivatives(); 
		
		// ?????
		/** repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradients 		( fShapes, 	u, u_n );
		Convert.Gradients 		( fShapes, 	gamma_p, gamma_p_n );
		Convert.Shapes			(	fShapes, 	fFEA_Shapes );
		Convert.Displacements	(	del_u, 	del_u_vec  );
		Convert.Displacements	(	del_gamma_p, 	del_gamma_p_vec  );
		Convert.Na				(	n_en, fShapes, 	fFEA_Shapes );

		/** Construct data used in BOTH PlastT and BalLinMomT (F,Fa,Fb,grad_ua,...etc.)
		 * 	Presently, Tahoe cannot exploit this fact.  n and np1 are calculated for coarse field, then
		 * 	calculated again for fine field -- this is a waste and defeats the purpose of VMS_VariableT. 
		 *  Note: n is last time step (known data), no subscript, np1 or (n+1) is the 
		 *  next time step (what were solving for)   */

			/* calculate coarse scale nodal force */
			if (is_coarse)
			{
				/* residual and tangent for coarse scale */
				fEquation_I -> Construct ( fFEA_Shapes, fBalLinMomMaterial, np1, n, step_number, delta_t );
				fEquation_I -> Form_LHS_Keps_Kd ( fKeps_I, fKd_I );
				fEquation_I -> Form_RHS_F_int ( fFint_I );
				fFint_I *= -1.0;  

				/* add body force */
				if (formBody) {
//					double density = fBalLinMomMaterial->Retrieve(Iso_MatlT::kDensity);
					double density = 1.0;
					DDub = 0.0;
					AddBodyForce(DDub);
				
					/* add body force to fRHS */
					fRHS = 0.0;
					FormMa(kConsistentMass, -density, &DDub, NULL);
					fFint_I += fRHS;
				}
			}
			else /* fine scale nodal force */
			{
				/* residual and tangent for fine scale */
				fEquation_II -> Construct ( fFEA_Shapes, fPlastMaterial, np1, n, step_number, delta_t, FEA::kBackward_Euler );
				fEquation_II -> Form_LHS_Keps_Kd ( fKeps_II, 	fKd_II );
				fEquation_II -> Form_RHS_F_int ( fFint_II );
				fFint_II *= -1.0;
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
//	cout << "F_int = \n" << fFint_I << endl;
}

//---------------------------------------------------------------------

double APS_AssemblyT::InternalEnergy ( void )
{
	//not implemented
return 0.0;
}

//---------------------------------------------------------------------

/* write restart data to the output stream */
void APS_AssemblyT::WriteRestart(ostream& out) const
{
	/* inherited */
	ElementBaseT::WriteRestart(out);

	/* write state variable data */
	out << fdState;
}

//---------------------------------------------------------------------

/* read restart data to the output stream */
void APS_AssemblyT::ReadRestart(istream& in)
{
	/* inherited */
	ElementBaseT::ReadRestart(in);

	/* write state variable data */
	in >> fdState;
}

//---------------------------------------------------------------------

void APS_AssemblyT::RegisterOutput(void)
{
	/* collect block ID's */
	ArrayT<StringT> block_ID(fBlockData.Length());
	for (int i = 0; i < block_ID.Length(); i++)
		block_ID[i] = fBlockData[i].ID();

	/* output per element - stresses at the integration points */
	int n_stress = dSymMatrixT::NumValues(NumSD());
	ArrayT<StringT> e_labels(fNumIP*n_stress);

	/* over integration points */
	// what stress output?????
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
	int num_node_output = fDispl.NumDOF() + fPlast.NumDOF() + n_stress;
	ArrayT<StringT> n_labels(num_node_output);
	count = 0;

	/* labels from fine scale */
	const ArrayT<StringT>& fine_labels = fPlast.Labels();
	for (int i = 0; i < fine_labels.Length(); i++)
		n_labels[count++] = fine_labels[i];

	/* labels from coarse scale */
	const ArrayT<StringT>& coarse_labels = fDispl.Labels();
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

void APS_AssemblyT::WriteOutput(void)
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
	int num_node_output = fDispl.NumDOF() + fPlast.NumDOF() + n_stress;
	dArray2DT n_values(nodes_used.Length(), num_node_output);

	/* collect nodal values */
	const dArray2DT& fgamma_p = fPlast[0];
	const dArray2DT& fU = fDispl[0];
	for (int i = 0; i < nodes_used.Length(); i++)
	{
		int node = nodes_used[i];
		double* row = n_values(i);
		for (int j = 0; j < fgamma_p.MinorDim(); j++)
			*row++ = fgamma_p(node,j);

		for (int j = 0; j < fU.MinorDim(); j++)
			*row++ = fU(node,j);

		double* p_stress = extrap_values(i);
		for (int j = 0; j < n_stress; j++)
			*row++ = p_stress[j];
	}

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, fIPVariable);


}	



//---------------------------------------------------------------------

void 	APS_AssemblyT::Get_Fext_I ( dArrayT &fFext_I )
{
	fFext_I = 0.0;
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
void APS_AssemblyT::RHSDriver_staggered(void)
{
	const char caller[] = "APS_AssemblyT::RHSDriver_staggered";
	if (fDispl.Group() == fPlast.Group())
		ExceptionT::GeneralFail(caller, "coarse and fine group must be different: %d == %d",
			fDispl.Group(), fPlast.Group());

	int curr_group = ElementSupport().CurrentGroup();

	/* stress output work space */
	int n_stress = dSymMatrixT::NumValues(NumSD());
	dArray2DT   out_variable_all;
	dSymMatrixT out_variable;

	/* time Step Increment */
	double delta_t = ElementSupport().TimeStep();
	time = ElementSupport().Time();
	step_number = ElementSupport().StepNumber();
	iArrayT plast_eq;

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

		SetLocalU (u);			 
		SetLocalU (u_n);
		SetLocalU (gamma_p);			 
		SetLocalU (gamma_p_n);

		del_u.DiffOf (u, u_n);
		del_gamma_p.DiffOf (gamma_p, gamma_p_n);

	 	SetLocalX(fInitCoords); 
	 	// for anti-plane shear, coordinates do not change
		fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, u); 
		fShapes->SetDerivatives(); 
		
		// ?????
		/** repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradients 		( fShapes, 	u, u_n );
		Convert.Gradients 		( fShapes, 	gamma_p, gamma_p_n );
		Convert.Shapes			(	fShapes, 	fFEA_Shapes );
		Convert.Displacements	(	del_u, 	del_u_vec  );
		Convert.Displacements	(	del_gamma_p, 	del_gamma_p_vec  );
		Convert.Na				(	n_en, fShapes, 	fFEA_Shapes );
		
		/* which field */
	  //SolverGroup 1 (gets field 2) <-- u (obtained by a rearranged Equation I)
		if ( curr_group == fDispl.Group()  )	
		{

			if (bStep_Complete) {
			// do nothing
			}
			else { //-- Still Iterating

				/** Compute N-R matrix equations */
				fEquation_I -> Construct ( fFEA_Shapes, fBalLinMomMaterial, np1, n, step_number, delta_t );
				fEquation_I -> Form_LHS_Keps_Kd ( fKeps_I, fKd_I );
				fEquation_I -> Form_RHS_F_int ( fFint_I );

				/** Set coarse LHS */
				fLHS = fKd_I;

				/** Compute coarse RHS */
				fKeps_I.Multx ( del_gamma_p_vec, fRHS );
				fRHS += fFint_I; 
				fRHS *= -1.0; 

				/** Compute Traction B.C. and Body Forces */
				Get_Fext_I ( fFext_I );
				fRHS += fFext_I;
				
				/* add body forces */
				if (formBody) {
//					double density = fBalLinMomMaterial->Retrieve(Iso_MatlT::kDensity);
					double density = 1.0;
					DDub = 0.0;
					AddBodyForce(DDub);
					FormMa(kConsistentMass, -density, &DDub, NULL);				
				}
			
				/* add to global equations */
				ElementSupport().AssembleLHS	( fDispl.Group(), fLHS, CurrentElement().Equations() );
				ElementSupport().AssembleRHS 	( fDispl.Group(), fRHS, CurrentElement().Equations() );
			}
		}

		// SolverGroup 2 (gets field 1) <-- ua (obtained by a rearranged Equation II)
		else if (curr_group == fPlast.Group() || (bStep_Complete && render_variable_group==1) )	
		{

			if (bStep_Complete) { 
			// do nothing
			}
			else { //-- Still Iterating

				/** Compute N-R matrix equations */
				fEquation_II -> Construct ( fFEA_Shapes, fPlastMaterial, np1, n, step_number, delta_t, FEA::kBackward_Euler );
				fEquation_II -> Form_LHS_Keps_Kd ( fKeps_II, 	fKd_II );
				fEquation_II -> Form_RHS_F_int ( fFint_II );

				/** Set LHS */
				fLHS = fKeps_II;	
		
				/** Compute fine RHS (or Fint_bar_II in FAXed notes)  */
				fKd_II.Multx ( del_u_vec, fRHS );
				fRHS += fFint_II; 
				fRHS *= -1.0; 
		
				/* fine scale equation numbers */
				fEqnos_plast.RowAlias ( CurrElementNumber(), plast_eq );

				/* add to global equations */
				ElementSupport().AssembleLHS ( fPlast.Group(), fLHS, plast_eq );
				ElementSupport().AssembleRHS ( fPlast.Group(), fRHS, plast_eq );
			}

		}
		else ExceptionT::GeneralFail(caller);
	}
}

//---------------------------------------------------------------------
/* form group contribution to the stiffness matrix and RHS */
void APS_AssemblyT::RHSDriver_monolithic(void)
{
	const char caller[] = "APS_AssemblyT::RHSDriver_monolithic";
	if (fDispl.Group() != fPlast.Group())
		ExceptionT::GeneralFail(caller, "coarse and fine group must be the same: %d != %d",
			fDispl.Group(), fPlast.Group());

	int curr_group = ElementSupport().CurrentGroup();

	/* stress output work space */
	int n_stress = dSymMatrixT::NumValues(NumSD());
	dArray2DT   out_variable_all;
	dSymMatrixT out_variable;

	/* time Step Increment */
	double delta_t = ElementSupport().TimeStep();
	time = ElementSupport().Time();
	step_number = ElementSupport().StepNumber();
	iArrayT displ_eq, plast_eq;

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

		SetLocalU (u);			 
		SetLocalU (u_n);
		SetLocalU (gamma_p);			 
		SetLocalU (gamma_p_n);

		del_u.DiffOf (u, u_n);
		del_gamma_p.DiffOf (gamma_p, gamma_p_n);

	 	SetLocalX(fInitCoords); 
	 	// not change in coordinates
		fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, u); 
		fShapes->SetDerivatives(); 
		
		/* repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradients 		( fShapes, 	u, u_n );
		Convert.Gradients 		( fShapes, 	gamma_p, gamma_p_n );
		Convert.Shapes			(	fShapes, 	fFEA_Shapes );
		Convert.Displacements	(	del_u, 	del_u_vec  );
		Convert.Displacements	(	del_gamma_p, 	del_gamma_p_vec  );
		Convert.Na				(	n_en, fShapes, 	fFEA_Shapes );

		if (bStep_Complete) { 
		//nothing done
		}
		else { //-- Still Iterating

			/* residual and tangent for coarse scale */
			fEquation_I -> Construct ( fFEA_Shapes, fBalLinMomMaterial, np1, n, step_number, delta_t );
			fEquation_I -> Form_LHS_Keps_Kd ( fKeps_I, fKd_I );
			fEquation_I -> Form_RHS_F_int ( fFint_I );
			fFint_I *= -1.0;

			/* add body force */
			if (formBody) {
//				double density = fBalLinMomMaterial->Retrieve(Iso_MatlT::kDensity);
				double density = 1.0;
				DDub = 0.0;
				AddBodyForce(DDub);
				
				/* add body force to fRHS */
				fRHS = 0.0;
				FormMa(kConsistentMass, -density, &DDub, NULL);
				fFint_I += fRHS;
			}

			/* residual and tangent for fine scale */
			fEquation_II -> Construct ( fFEA_Shapes, fPlastMaterial, np1, n, step_number, delta_t, FEA::kBackward_Euler );
			fEquation_II -> Form_LHS_Keps_Kd ( fKeps_II, 	fKd_II );
			fEquation_II -> Form_RHS_F_int ( fFint_II );
			fFint_II *= -1.0;

			/* equations numbers */
			const iArrayT& all_eq = CurrentElement().Equations();
			displ_eq.Set(fFint_I.Length(), all_eq.Pointer());
			plast_eq.Set(fFint_II.Length(), all_eq.Pointer(fFint_I.Length()));

			/* assemble residuals */
			ElementSupport().AssembleRHS(curr_group, fFint_I, displ_eq);
			ElementSupport().AssembleRHS(curr_group, fFint_II, plast_eq);

			/* assemble components of the tangent */
			ElementSupport().AssembleLHS(curr_group, fKd_I, displ_eq);
			ElementSupport().AssembleLHS(curr_group, fKeps_II, plast_eq);
			ElementSupport().AssembleLHS(curr_group, fKeps_I, displ_eq, plast_eq);
			ElementSupport().AssembleLHS(curr_group, fKd_II, plast_eq, displ_eq);
		}
	}
}

