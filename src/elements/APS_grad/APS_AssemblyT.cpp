/* $Id: APS_AssemblyT.cpp,v 1.39 2003-10-12 02:51:19 raregue Exp $ */
#include "APS_AssemblyT.h"

#include "ShapeFunctionT.h"
#include "Traction_CardT.h"

#include "ifstreamT.h"
#include "ofstreamT.h"

#include "APS_MatlT.h"
#include "Shear_MatlT.h"

#include "APS_Bal_EqT.h"
#include "APS_BCJT.h"

#include "OutputSetT.h"

using namespace Tahoe;

/* parameters */
//static int knum_d_state = 2; // double's needed per ip
//static int knum_i_state = 0; // int's needed per ip

//---------------------------------------------------------------------

/* constructor */
APS_AssemblyT::APS_AssemblyT(const ElementSupportT& support, const FieldT& displ, 
							const FieldT& gammap):
	ElementBaseT(support, displ), //pass the displacement field to the base class
	u(LocalArrayT::kDisp),
	u_n(LocalArrayT::kLastDisp),
	DDu(LocalArrayT::kAcc),
	gamma_p(LocalArrayT::kDisp),
	gamma_p_n(LocalArrayT::kLastDisp),
	fInitCoords(LocalArrayT::kInitCoords),
	fCurrCoords(LocalArrayT::kCurrCoords),
	fBodySchedule(NULL),
	fBody(NumDOF()),
	fTractionBCSet(0),
	fDispl(displ),
	fPlast(gammap),
	fKdd(ElementMatrixT::kNonSymmetric),
	fKdd_face(ElementMatrixT::kNonSymmetric),
	fKdeps(ElementMatrixT::kNonSymmetric),
	fKepsd(ElementMatrixT::kNonSymmetric),
	fKepseps(ElementMatrixT::kNonSymmetric),
	fEquation_d(NULL),
	fEquation_eps(NULL),
	bStep_Complete(0)
{
	
	knum_d_state = 6; // double's needed per ip
	knum_i_state = 0; // int's needed per ip
	
	knumstrain = 3; 
	knumstress = 2; 
	
	output = "out";

	/* check - some code below assumes that both fields have the
	 * same dimension. TEMP?  get rid of this!!!!!! */ 
	/*    if (fDispl.NumDOF() != fPlast.NumDOF()) 
				ExceptionT::BadInputValue("APS_AssemblyT::APS_AssemblyT"); */

	/* read parameters from input */
	ifstreamT& in = ElementSupport().Input();
	in >> fGeometryCode; //TEMP - should actually come from the geometry database
	in >> fNumIP;
	in >> fGeometryCodeSurf; //TEMP - should actually come from the geometry database
	in >> fNumIPSurf;

	fMaterial_Data.Dimension ( kNUM_FMAT_TERMS );

	in >> iPlastModelType; 

	//-- Elasticity parameters
	in >> fMaterial_Data[kMu];

	//-- Plasticity parameters
	in >> fMaterial_Data[km_rate];
	in >> fMaterial_Data[kgamma0_dot_1];
	in >> fMaterial_Data[kgamma0_dot_2];
	in >> fMaterial_Data[km1_x];
	in >> fMaterial_Data[km1_y];
	in >> fMaterial_Data[km2_x];
	in >> fMaterial_Data[km2_y];

	//-- Backstress Parameters
	in >> fMaterial_Data[kl];

	//-- Isotropic Hardening Parameters
	in >> fMaterial_Data[kH];
	
	//-- Initial value of state variable
	in >> fMaterial_Data[kkappa0_1];
	in >> fMaterial_Data[kkappa0_2];
	
	/* allocate the global stack object (once) */
	extern FEA_StackT* fStack;
	if (!fStack) fStack = new FEA_StackT;

	/* prescribed plastic gradient at surface */
	in >> num_sidesets;
	fSideSetID.Dimension(num_sidesets);
	fSideSetElements.Dimension(num_sidesets);
	fSideSetFaces.Dimension(num_sidesets);
	fPlasticGradientWght.Dimension(num_sidesets);
	fPlasticGradientFaces.Dimension(num_sidesets);
	fPlasticGradientFaceEqnos.Dimension(num_sidesets);
	
	//enable the model manager
	ModelManagerT& model = ElementSupport().Model();

	ArrayT<GeometryT::CodeT> facet_geom;
	iArrayT facet_nodes;
	for (int i = 0; i < num_sidesets; i++)
	{
		in >> fSideSetID[i];
		in >> fPlasticGradientWght[i];
	
		// get nodes-on-faces
		model.SideSet(fSideSetID[i], facet_geom, facet_nodes, fPlasticGradientFaces[i]);
		
		// get the side set information {element, face number} for each
		// face in the set
		const iArray2DT& side_set = model.SideSet(fSideSetID[i]);
		
		// get side set elements - element numbers in zeroth column
		fSideSetElements[i].Dimension(side_set.MajorDim());
		fSideSetFaces[i].Dimension(side_set.MajorDim());
		side_set.ColumnCopy(0, fSideSetElements[i]);
		side_set.ColumnCopy(1, fSideSetFaces[i]);
	}
	
	Echo_Input_Data();
}

//---------------------------------------------------------------------

/* destructor */
APS_AssemblyT::~APS_AssemblyT(void) 
{  
	delete fShapes;
	delete fEquation_d; 
	delete fEquation_eps; 
	delete fBalLinMomMaterial; 
	delete fPlastMaterial;

	/* free the global stack object (once) */
	extern FEA_StackT* fStack;
	if (fStack) {
		delete fStack;
		fStack = NULL;
	}

}

//--------------------------------------------------------------------

void APS_AssemblyT::Echo_Input_Data(void) {

	cout << "#######################################################" << endl; 
	cout << "############### ECHO APS DATA #########################" << endl; 
	cout << "#######################################################" << endl; 

	//################## material data ##################

	cout << "iPlastModelType " 						<< iPlastModelType 				<< endl; 
	
	//-- Elasticity parameters 
	cout << "fMaterial_Data[kMu] "  				<< fMaterial_Data[kMu] 		<< endl;

	//-- Plasticity parameters
	cout << "fMaterial_Data[km_rate] " 				<< fMaterial_Data[km_rate] 	<< endl;
	cout << "fMaterial_Data[kgamma0_dot_1] " 		<< fMaterial_Data[kgamma0_dot_1] << endl;
	cout << "fMaterial_Data[kgamma0_dot_2] " 		<< fMaterial_Data[kgamma0_dot_2] << endl;
	cout << "fMaterial_Data[km1_x] " 				<< fMaterial_Data[km1_x] 		<< endl;
	cout << "fMaterial_Data[km1_y] " 				<< fMaterial_Data[km1_y] 		<< endl;
	cout << "fMaterial_Data[km2_x] " 				<< fMaterial_Data[km2_x] 		<< endl;
	cout << "fMaterial_Data[km2_y] " 				<< fMaterial_Data[km2_y] 		<< endl;

	//-- Backstress Parameters
	cout << "fMaterial_Data[kl] " 					<< fMaterial_Data[kl]			<< endl;

	//-- Isotropic Hardening Parameters
	cout << "fMaterial_Data[kH] "					<< fMaterial_Data[kH]			<< endl;
	
	//-- Initial state variable
	cout << "fMaterial_Data[kkappa0_1] "			<< fMaterial_Data[kkappa0_1]	<< endl;
	cout << "fMaterial_Data[kkappa0_2] "			<< fMaterial_Data[kkappa0_2]	<< endl;
	
}

void APS_AssemblyT::Initialize(void)
{
	/* inherited */
	ElementBaseT::Initialize();
	
	/* dimensions (notation as per Hughes' Book) */
	n_ip = fNumIP;
	n_sd = NumSD();
	//n_df = NumDOF(); 
	n_df = 1+n_sd; 
	n_en = NumElementNodes();
	n_np = ElementSupport().NumNodes();
	n_el = NumElements();
	
	n_sd_surf = n_sd;
	n_en_surf = 2;

	n_en_x_n_df = n_en*n_df;
	n_en_x_n_sd = n_en*n_sd;

	/* set local arrays for coarse scale */
	int dum=1;
	u.Dimension (n_en, dum);
	u_n.Dimension (n_en, dum);
	DDu.Dimension (n_en, dum);
	del_u.Dimension (n_en, dum);
	del_u_vec.Dimension (n_en);
	fDispl.RegisterLocal(u);
	fDispl.RegisterLocal(u_n);

	/* set local arrays for fine scale */
	gamma_p.Dimension (n_en, n_sd);
	gamma_p_n.Dimension (n_en, n_sd);
	del_gamma_p.Dimension (n_en, n_sd);
	dum = n_en*n_sd;
	del_gamma_p_vec.Dimension (dum);
	fPlast.RegisterLocal(gamma_p);
	fPlast.RegisterLocal(gamma_p_n);

	/* set shape functions */
	fInitCoords.Dimension(n_en, n_sd);
	ElementSupport().RegisterCoordinates(fInitCoords);	
	fCurrCoords.Dimension(n_en, n_sd);
	fShapes = new ShapeFunctionT(fGeometryCode, fNumIP, fCurrCoords);
	fShapes->Initialize();
	
	//fNormal.Dimension ( n_sd );
	
	/* allocate state variable storage */
	int num_ip = fNumIP;
	fdState_new.Dimension(n_el, num_ip*knum_d_state);
	fdState.Dimension(n_el, num_ip*knum_d_state);
	fiState_new.Dimension(n_el, num_ip*knum_i_state);
	fiState.Dimension(n_el, num_ip*knum_i_state);
	
	/* storage for the fine scale equation numbers */
	fEqnos_plast.Dimension(n_el, n_en*n_sd);
	fEqnos_plast = -1;

	/* initialize state variables */
	fdState = 0;
	fdState_new = 0;
	fiState = 0;
	fiState_new = 0;
	
	/* set cards to data in array - NOT NEEDED IF YOU'RE NOT
	 * GOING TO USE THE ElementCardT ARRAY? */
	for (int i= 0; i < fElementCards.Length(); i++)
		fElementCards[i].Set(fiState.MinorDim(), fiState(i), fdState.MinorDim(), fdState(i));
		                     
	/* construct the black boxs */  

	Select_Equations ( BalLinMomT::kAPS_Bal_Eq, iPlastModelType );
	fEquation_eps -> Initialize ( n_ip, n_sd, n_en, knum_d_state, ElementSupport().StepNumber() );
	//step_number_last_iter = 0; 
	//step_number_last_iter = ElementSupport().StepNumber();  // This may crash or not work

	/* FEA Allocation */

	// these dimensions should be different since want to use quadratic interp for u
	// and linear interp for gamma_p
	
	//fgrad_u.FEA_Dimension 			( fNumIP, n_sd );
	dum=1;
	fgrad_u.FEA_Dimension 			( fNumIP, dum, n_sd );
	fgrad_u_surf.FEA_Dimension 		( fNumIPSurf, dum, n_sd );
	fgamma_p.FEA_Dimension 			( fNumIP, n_sd );
	fgamma_p_surf.FEA_Dimension 	( fNumIPSurf, n_sd );
	fgrad_gamma_p.FEA_Dimension 	( fNumIP, n_sd,n_sd );
	//fgrad_u_n.FEA_Dimension 		( fNumIP, n_sd );
	fgrad_u_n.FEA_Dimension 		( fNumIP, dum, n_sd );
	fgrad_u_surf_n.FEA_Dimension 	( fNumIPSurf, dum, n_sd );
	fgamma_p_n.FEA_Dimension 		( fNumIP, n_sd );
	fgamma_p_surf_n.FEA_Dimension 	( fNumIPSurf, n_sd );
	fgrad_gamma_p_n.FEA_Dimension 	( fNumIP, n_sd,n_sd );
	
	fstate.FEA_Dimension 			( fNumIP, knum_d_state );
	fstate_n.FEA_Dimension 			( fNumIP, knum_d_state );

	//check these dims
	fKdd.Dimension 			( n_en, n_en );
	fKdeps.Dimension 		( n_en, n_en_x_n_sd );
	fKepsd.Dimension 		( n_en_x_n_sd, n_en );
	fKepseps.Dimension 		( n_en_x_n_sd, n_en_x_n_sd );

	fFd_int.Dimension 		( n_en );
	fFd_ext.Dimension 		( n_en );
	fFeps_int.Dimension 	( n_en_x_n_sd );
	fFeps_ext.Dimension 	( n_en_x_n_sd );
	
	fKdd_face.Dimension 		( n_en_surf, n_en_surf );
	fFd_int_face.Dimension 		( n_en_surf );

	fFEA_Shapes.Construct	( fNumIP,n_sd,n_en );
	fFEA_SurfShapes.Construct	( fNumIPSurf,n_sd_surf,n_en_surf );
	
	Render_Vector.Dimension ( n_el );
	for (int e=0; e<n_el; e++) {
		Render_Vector[e].Construct ( 1, n_ip, knumstrain+knumstress+knum_d_state );	
	}


	/* streams */
	ifstreamT& in  = ElementSupport().Input();
	ofstreamT& out = ElementSupport().Output();

	/* storage for integration point strain, stress, and ISVs*/
	fIPVariable.Dimension (n_el, fNumIP*(knumstrain+knumstress+knum_d_state));
	fIPVariable = 0.0;

	/* allocate storage for nodal forces */
	//fForces_at_Node.Dimension ( n_sd );

	/* body force specification */
	#pragma message("APS_AssemblyT::Initialize: careful, no body force for gammap ")
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

void APS_AssemblyT::Equations(AutoArrayT<const iArray2DT*>& eq_d,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_eps)
{
	/* doing monolithic solution */
	if (fDispl.Group() == fPlast.Group())
	{
		#pragma message("APS_AssemblyT::Equations: how determine NumDOF for each field?")
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
			iArray2DT displ_eq(nel, nen*ndof_displ);
			iArray2DT plast_eq(nel, nen*ndof_plast);
			
			/* get equation numbers */
			fDispl.SetLocalEqnos(connects, displ_eq);
			fPlast.SetLocalEqnos(connects, plast_eq);
			
			/* write into one array */
			fEqnos[i].BlockColumnCopyAt(displ_eq, 0);
			fEqnos[i].BlockColumnCopyAt(plast_eq, displ_eq.MinorDim());

			/* add to list of equation numbers */
			eq_d.Append(&fEqnos[i]);
		}
	
		/* reset pointers to element cards */
		SetElementCards();	
	}
	else
	{
		/* ElementBaseT handles equation array for the coarse scale */
		if (ElementSupport().CurrentGroup() == fDispl.Group())
			ElementBaseT::Equations(eq_d,eq_eps);

		/* fine scale equations */
		if (ElementSupport().CurrentGroup() == fPlast.Group())
		{
			/* collect local equation numbers */
			fPlast.SetLocalEqnos(fConnectivities, fEqnos_plast);
		
			eq_d.Append(&fEqnos_plast);
		}
	}
	
	/* get the equation number for the nodes on the faces */
	for (int i = 0; i < fPlasticGradientFaceEqnos.Length(); i++)
	{
		iArray2DT& faces = fPlasticGradientFaces[i];
		iArray2DT& eqnos = fPlasticGradientFaceEqnos[i];
		eqnos.Dimension(faces.MajorDim(), faces.MajorDim()*fDispl.NumDOF());
	
		fDispl.SetLocalEqnos(faces, eqnos);
	}
}


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
			fEquation_d = new APS_Bal_EqT;
			fBalLinMomMaterial = new Shear_MatlT;
			fBalLinMomMaterial -> Assign ( Shear_MatlT::kMu, fMaterial_Data[kMu] );
			break;

		default :
			cout << "APS_AssemblyT::Select_Equations() .. ERROR >> bad iBalScale \n";
			break;
	}

	/** Choices for Fine-Scale Equation */

	switch ( iPlastScale )	{

		case PlastT::kAPS_BCJ :
			fEquation_eps	= new APS_BCJT;
			fPlastMaterial	= new APS_MatlT;		
			fPlastMaterial -> Assign (	APS_MatlT::kMu, 		fMaterial_Data[kMu] 		); 	
			fPlastMaterial -> Assign ( 	APS_MatlT::km_rate, 	fMaterial_Data[km_rate] 	); 	
			fPlastMaterial -> Assign ( 	APS_MatlT::kgamma0_dot_1, fMaterial_Data[kgamma0_dot_1]); 
			fPlastMaterial -> Assign ( 	APS_MatlT::kgamma0_dot_2, fMaterial_Data[kgamma0_dot_2]); 	
			fPlastMaterial -> Assign ( 	APS_MatlT::km1_x, 		fMaterial_Data[km1_x] 		); 
			fPlastMaterial -> Assign ( 	APS_MatlT::km1_y, 		fMaterial_Data[km1_y] 		); 	
			fPlastMaterial -> Assign ( 	APS_MatlT::km2_x, 		fMaterial_Data[km2_x] 		); 
			fPlastMaterial -> Assign ( 	APS_MatlT::km2_y, 		fMaterial_Data[km2_y] 		); 	
			fPlastMaterial -> Assign ( 	APS_MatlT::kl, 			fMaterial_Data[kl] 		); 	
			fPlastMaterial -> Assign ( 	APS_MatlT::kH, 			fMaterial_Data[kH] 		);
			fPlastMaterial -> Assign ( 	APS_MatlT::kkappa0_1, 	fMaterial_Data[kkappa0_1] 	); 	
			fPlastMaterial -> Assign ( 	APS_MatlT::kkappa0_2, 	fMaterial_Data[kkappa0_2] 	); 	
			break;
			
		default :
			cout << "APS_AssemblyT::Select_Equations() .. ERROR >> bad iPlastScale \n";
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
	out << " Plastic gradient field. . . . . . . . . . . . . = \"" << fPlast.Name() << "\"\n";
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
		element_force = &fFd_int;
		num_force = fDispl.NumDOF();
		}
	else if (field.Name() == fPlast.Name()) {
		is_coarse = false;
		element_force = &fFeps_int;
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
	
	dArray2DT fdstatenew_all, fdstate_all;

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

		SetLocalU (u);
		SetLocalU (u_n);
		SetLocalU (gamma_p);
		SetLocalU (gamma_p_n);

		del_u.DiffOf (u, u_n);
		del_gamma_p.DiffOf (gamma_p, gamma_p_n);

	 	SetLocalX(fInitCoords); 
	 	// coordinates do not change for anti-plane shear
		//fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, u); 
		fCurrCoords = fInitCoords;
		fShapes->SetDerivatives(); 
		
		//update state variables
		fdstatenew_all.Alias(fNumIP, knum_d_state, fdState_new(CurrElementNumber()));
		fdstate_all.Alias(fNumIP, knum_d_state, fdState(CurrElementNumber()));
		
		/** repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradients 		( fShapes, 	u, u_n, fgrad_u, fgrad_u_n );
		Convert.Gradients 		( fShapes, 	gamma_p, gamma_p_n, fgrad_gamma_p, fgrad_gamma_p_n );
		Convert.Interpolate 	( fShapes, 	gamma_p, gamma_p_n, fgamma_p, fgamma_p_n );
		Convert.Shapes			(	fShapes, fFEA_Shapes );
		Convert.Displacements	(	del_u, 	del_u_vec  );
		Convert.Displacements	(	del_gamma_p, 	del_gamma_p_vec  );
		Convert.Na				(	n_en, fShapes, 	fFEA_Shapes );
		Convert.Copy			(	fNumIP, knum_d_state, fdstatenew_all, fstate );
		Convert.Copy			(	fNumIP, knum_d_state, fdstate_all, fstate_n );
		
		APS_VariableT np1(	fgrad_u, fgrad_u_surf, fgamma_p, fgamma_p_surf, fgrad_gamma_p, fstate ); // Many variables at time-step n+1
		APS_VariableT   n(	fgrad_u_n, fgrad_u_surf_n,fgamma_p_n, fgamma_p_surf_n, fgrad_gamma_p_n, fstate_n );	// Many variables at time-step n	 

			/* calculate coarse scale nodal force */
			if (is_coarse)
			{
				/* residual and tangent for coarse scale */
				fEquation_d -> Construct ( fFEA_Shapes, fBalLinMomMaterial, fPlastMaterial, np1, n, 
											step_number, delta_t );
				fEquation_d -> Form_LHS_Keps_Kd ( fKdeps, fKdd );
				fEquation_d -> Form_RHS_F_int ( fFd_int, np1 );
				fFd_int *= -1.0;  

				/* add body force */
				if (formBody) {
//					double density = fBalLinMomMaterial->Retrieve(Iso_MatlT::kDensity);
					double density = 1.0;
					DDu = 0.0;
					AddBodyForce(DDu);
				
					/* add body force to fRHS */
					fRHS = 0.0;
					FormMa(kConsistentMass, -density, &DDu, NULL);
					fFd_int += fRHS;
				}
			}
			else /* fine scale nodal force */
			{
				/* residual and tangent for fine scale */
				fEquation_eps -> Construct ( fFEA_Shapes, fPlastMaterial, np1, n, 
											step_number, delta_t, FEA::kBackward_Euler );
				fEquation_eps -> Form_LHS_Keps_Kd ( fKepseps, 	fKepsd );
				fEquation_eps -> Form_RHS_F_int ( fFeps_int );
				fFeps_int *= -1.0;
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
//	cout << "F_int = \n" << fFd_int << endl;
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

	/* output per element - strain, stress, and ISVs at the integration points */
	ArrayT<StringT> e_labels(fNumIP*(knumstrain+knumstress+knum_d_state));

	/* over integration points */
	const char* slabels2D[] = {"gamma_x", "gamma_y", "gammap_curl", "s_xz", "s_yz"};
	const char* svlabels2D[] = {"xi_1", "kappa_1", "gamma_dot_1", "xi_2", "kappa_2", "gamma_dot_2"};
	int count = 0;
	for (int j = 0; j < fNumIP; j++)
	{
		StringT ip_label;
		ip_label.Append("ip", j+1);
			
		/* over strain and stress components */
		for (int i = 0; i < knumstrain+knumstress; i++)
		{
			e_labels[count].Clear();
			e_labels[count].Append(ip_label, ".", slabels2D[i]);
			count++;
		}
		
		/* over state variables */
		for (int i = 0; i < knum_d_state; i++)
		{
			e_labels[count].Clear();
			e_labels[count].Append(ip_label, ".", svlabels2D[i]);
			count++;
		}
	}		

	/* output per node */
	int num_node_output = fDispl.NumDOF() + fPlast.NumDOF() + knumstrain + knumstress + knum_d_state;
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

	/* labels from strains and stresses at the nodes */
	for (int i = 0; i < knumstrain+knumstress; i++)
		n_labels[count++] = slabels2D[i];
		
	/* labels from state variables at the nodes */
	for (int i = 0; i < knum_d_state; i++)
		n_labels[count++] = svlabels2D[i];

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
	bStep_Complete=1;
	RHSDriver();
	bStep_Complete=0;

	/* my output set */
	const OutputSetT& output_set = ElementSupport().OutputSet(fOutputID);
	
	/* my nodes used */
	const iArrayT& nodes_used = output_set.NodesUsed();

	/* smooth stresses to nodes */
	ElementSupport().ResetAverage(knumstrain+knumstress+knum_d_state);
	dArray2DT out_variable_all;
	dArrayT out_variable;
	dArray2DT nd_var(NumElementNodes(), knumstrain+knumstress+knum_d_state);
	Top();
	while (NextElement())
	{
		/* extrapolate */
		nd_var = 0.0;
		out_variable_all.Alias(fNumIP, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(fShapes->CurrIP()));
			fShapes->Extrapolate(out_variable, nd_var);
		}
	
		/* accumulate - extrapolation done from ip's to corners => X nodes  */
		ElementSupport().AssembleAverage(CurrentElement().NodesX(), nd_var);
	}

	/* get nodally averaged values */
	dArray2DT extrap_values;
	ElementSupport().OutputUsedAverage(extrap_values);

	/* temp space for group displacements */
	int num_node_output = fDispl.NumDOF() + fPlast.NumDOF() + knumstrain + knumstress + knum_d_state;
	dArray2DT n_values(nodes_used.Length(), num_node_output);

	/* collect nodal values */
	const dArray2DT& fGamma_p = fPlast[0];
	const dArray2DT& fU = fDispl[0];
	for (int i = 0; i < nodes_used.Length(); i++)
	{
		int node = nodes_used[i];
		double* row = n_values(i);
		for (int j = 0; j < fGamma_p.MinorDim(); j++)
			*row++ = fGamma_p(node,j);

		for (int j = 0; j < fU.MinorDim(); j++)
			*row++ = fU(node,j);

		double* p_stress = extrap_values(i);
		for (int j = 0; j < (knumstrain+knumstress+knum_d_state); j++)
			*row++ = p_stress[j];
	}

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, fIPVariable);

}	



//---------------------------------------------------------------------

void 	APS_AssemblyT::Get_Fd_ext ( dArrayT &fFd_ext )
{
	fFd_ext = 0.0;
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
	dArray2DT out_variable_all, fdstatenew_all, fdstate_all;
	dArrayT out_variable;

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
	int e,l;
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
		//fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, u);
		fCurrCoords = fInitCoords; 
		fShapes->SetDerivatives(); 
				
		//-- Alias to update state variables and pass past ones to model subroutine 
		fdstatenew_all.Alias(fNumIP, knum_d_state, fdState_new(CurrElementNumber()));
		fdstate_all.Alias(fNumIP, knum_d_state, fdState(CurrElementNumber()));	
		
		//repackage data to forms compatible with FEA classes (very little cost in big picture)
		Convert.Gradients 		( fShapes, 	u, u_n, fgrad_u, fgrad_u_n );
		Convert.Gradients 		( fShapes, 	gamma_p, gamma_p_n, fgrad_gamma_p, fgrad_gamma_p_n );
		Convert.Interpolate 	( fShapes, 	gamma_p, gamma_p_n, fgamma_p, fgamma_p_n );
		Convert.Shapes			(	fShapes, 	fFEA_Shapes );
		Convert.Displacements	(	del_u, 	del_u_vec  );
		Convert.Displacements	(	del_gamma_p, 	del_gamma_p_vec  );
		Convert.Na				(	n_en, fShapes, 	fFEA_Shapes );
		Convert.Copy			(	fNumIP, knum_d_state, fdstatenew_all, fstate );
		Convert.Copy			(	fNumIP, knum_d_state, fdstate_all, fstate_n );
		
		APS_VariableT np1(	fgrad_u, fgrad_u_surf, fgamma_p, fgamma_p_surf, fgrad_gamma_p, fstate ); // Many variables at time-step n+1
		APS_VariableT   n(	fgrad_u_n, fgrad_u_surf_n, fgamma_p_n, fgamma_p_surf_n, fgrad_gamma_p_n, fstate_n );	// Many variables at time-step n
		
		/* which field */
	  	//SolverGroup 1 (gets field 1) <-- u (obtained by a rearranged Equation_d)
		if ( curr_group == fDispl.Group()  )	
		{

			if (bStep_Complete) {
			//do nothing here
			}
			else { //-- Still Iterating

				/** Compute N-R matrix equations */
				fEquation_d -> Construct ( fFEA_Shapes, fBalLinMomMaterial, fPlastMaterial, np1, n, 
											step_number, delta_t );
				fEquation_d -> Form_LHS_Keps_Kd ( fKdeps, fKdd );
				fEquation_d -> Form_RHS_F_int ( fFd_int, np1 );

				/** Set coarse LHS */
				fLHS = fKdd;

				/** Compute coarse RHS */
				fKdeps.Multx ( del_gamma_p_vec, fRHS );
				fRHS += fFd_int; 
				fRHS *= -1.0; 

				/** Compute Traction B.C. and Body Forces */
				Get_Fd_ext ( fFd_ext );
				fRHS += fFd_ext;
				
				/* add body forces */
				if (formBody) {
//					double density = fBalLinMomMaterial->Retrieve(Iso_MatlT::kDensity);
					double density = 1.0;
					DDu = 0.0;
					AddBodyForce(DDu);
					FormMa(kConsistentMass, -density, &DDu, NULL);				
				}
			
				/* add to global equations */
				ElementSupport().AssembleLHS ( fDispl.Group(), fLHS, CurrentElement().Equations() );
				ElementSupport().AssembleRHS ( fDispl.Group(), fRHS, CurrentElement().Equations() );
			}
		}

		// SolverGroup 2 (gets field 2) <-- gamma_p (obtained by a rearranged Equation_eps)
		else if (curr_group == fPlast.Group() )	
		{

			if (bStep_Complete) { 
			
				fEquation_eps -> Construct ( fFEA_Shapes, fPlastMaterial, np1, n, step_number, delta_t );
				fEquation_eps -> Get ( output, Render_Vector[e][0] );
			
				//-- Store/Register data in classic tahoe manner 
				out_variable_all.Alias(fNumIP, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
				for (l=0; l < fNumIP; l++) {
					out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(l));
					out_variable=Render_Vector[e][0][l];
					} 
			
			}
			else { //-- Still Iterating

				/** Compute N-R matrix equations */
				fEquation_eps -> Construct ( fFEA_Shapes, fPlastMaterial, np1, n, 
											step_number, delta_t, FEA::kBackward_Euler );
				fEquation_eps -> Form_LHS_Keps_Kd ( fKepseps, 	fKepsd );
				fEquation_eps -> Form_RHS_F_int ( fFeps_int );
				
				// update state variables
				np1.Update( APS::kstate, fstate );		
				Convert.Copy ( fNumIP, knum_d_state, fstate, fdstatenew_all );

				/** Set LHS */
				fLHS = fKepseps;	
		
				/** Compute fine RHS (or Fint_bar_II in FAXed notes)  */
				fKepsd.Multx ( del_u_vec, fRHS );
				fRHS += fFeps_int; 
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
	dArray2DT	out_variable_all, fdstatenew_all, fdstate_all;
	dArrayT		out_variable;

	/* time Step Increment */
	double delta_t = ElementSupport().TimeStep();
	time = ElementSupport().Time();
	step_number = ElementSupport().StepNumber();
	iArrayT displ_eq, plast_eq;

	/* work space for integration over faces */
	LocalArrayT face_coords(LocalArrayT::kInitCoords, n_en_surf, NumSD());
	ElementSupport().RegisterCoordinates(face_coords);
	iArrayT face_nodes, face_equations;
	dMatrixT face_jacobian(NumSD(), NumSD()-1);
	dMatrixT face_Q(NumSD());
	LocalArrayT face_gamma_p(LocalArrayT::kDisp, n_en_surf, NumSD());
	fPlast.RegisterLocal(face_gamma_p);

 	/* has (coarse scale) body forces */
	int formBody = 0;
	if (fBodySchedule && fBody.Magnitude() > kSmall)
		formBody = 1;

	/* loop over elements */
	int e,l;
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
	 	// no change in coordinates
		//fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, u); 
		fCurrCoords = fInitCoords;
		fShapes->SetDerivatives(); 
		
		//update state variables
		fdstatenew_all.Alias(fNumIP, knum_d_state, fdState_new(CurrElementNumber()));
		fdstate_all.Alias(fNumIP, knum_d_state, fdState(CurrElementNumber()));
		
		/* repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradients 		( fShapes, 	u, u_n, fgrad_u, fgrad_u_n );
		Convert.Gradients 		( fShapes, 	gamma_p, gamma_p_n, fgrad_gamma_p, fgrad_gamma_p_n );
		Convert.Interpolate 	( fShapes, 	gamma_p, gamma_p_n, fgamma_p, fgamma_p_n );
		Convert.Shapes			(	fShapes, 	fFEA_Shapes );
		Convert.Displacements	(	del_u, 	del_u_vec  );
		Convert.Displacements	(	del_gamma_p, 	del_gamma_p_vec  );
		Convert.Na				(	n_en, fShapes, 	fFEA_Shapes );
		Convert.Copy			(	fNumIP, knum_d_state, fdstatenew_all, fstate );
		Convert.Copy			(	fNumIP, knum_d_state, fdstate_all, fstate_n );
		
		APS_VariableT np1(	fgrad_u, fgrad_u_surf, fgamma_p, fgamma_p_surf, fgrad_gamma_p, fstate ); // Many variables at time-step n+1
		APS_VariableT   n(	fgrad_u_n, fgrad_u_surf_n, fgamma_p_n, fgamma_p_surf_n, fgrad_gamma_p_n, fstate_n );	// Many variables at time-step n
		
		if (bStep_Complete) { 
		
			fEquation_eps -> Construct ( fFEA_Shapes, fPlastMaterial, np1, n, step_number, delta_t );
			fEquation_eps -> Get ( output, Render_Vector[e][0] );
			
			//-- Store/Register data in classic tahoe manner 
			out_variable_all.Alias(fNumIP, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
			for (l=0; l < fNumIP; l++) {
				out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(l));
				out_variable=Render_Vector[e][0][l];
			} 
	
		}
		else { //-- Still Iterating
		
			/* residual and tangent for coarse scale */
			fEquation_d -> Construct (	fFEA_Shapes, fBalLinMomMaterial, fPlastMaterial, np1, n, 
										step_number, delta_t );
			fEquation_d -> Form_LHS_Keps_Kd ( fKdeps, fKdd );
			fEquation_d -> Form_RHS_F_int ( fFd_int, np1 );
			fFd_int *= -1.0;
			
			// add contribution from sidesets								
			for (int i = 0; i < num_sidesets; i++)
			{
				for (int j = 0; j < fSideSetElements[i].Length(); j++)
				{
					if (e == fSideSetElements[i][j])
					{
						/* collect coordinates over the face */
						fPlasticGradientFaces[i].RowAlias(j, face_nodes);
						face_coords.SetLocal(face_nodes);
						
						/* collect plastic strain over the face */
						face_gamma_p.SetLocal(face_nodes);
						
						/* shape functions over the given face */

						//int face = fSideSetElements[i][j];
						int face = fSideSetFaces[i][j];
						const ParentDomainT& surf_shape = fShapes->FacetShapeFunction(face);
						const ParentDomainT& parent = ShapeFunction().ParentDomain();
						iArrayT face_local_nodes(2);
						parent.NodesOnFacet(face, face_local_nodes);
						
						/* equations for the nodes on the face */
						fPlasticGradientFaceEqnos[i].RowAlias(j, face_equations);
						
						Convert.SurfShapeGradient	( n_en_surf, surf_shape, fFEA_SurfShapes, face_coords,
													parent, fInitCoords, *fShapes, u, u_n, fgrad_u_surf, fgrad_u_surf_n,
													face_gamma_p, fgamma_p_surf, face_local_nodes );
						APS_VariableT np1(	fgrad_u, fgrad_u_surf, fgamma_p, fgamma_p_surf, fgrad_gamma_p, fstate ); 
						fEquation_d -> Form_LHS_Kd_Surf ( fKdd_face, fFEA_SurfShapes );
						fEquation_d -> Form_RHS_F_int_Surf ( fFd_int_face, np1, fPlasticGradientWght[i] );

						ElementSupport().AssembleRHS(curr_group, fFd_int_face, face_equations);
						ElementSupport().AssembleLHS(curr_group, fKdd_face, face_equations);
					}
				}
			}
			
			
			/* add body force */
			if (formBody) {
//				double density = fBalLinMomMaterial->Retrieve(Iso_MatlT::kDensity);
				double density = 1.0;
				DDu = 0.0;
				AddBodyForce(DDu);
				
				/* add body force to fRHS */
				fRHS = 0.0;
				FormMa(kConsistentMass, -density, &DDu, NULL);
				fFd_int += fRHS;
			}

			/* residual and tangent for fine scale */
			fEquation_eps -> Construct ( fFEA_Shapes, fPlastMaterial, np1, n, step_number, 
										delta_t, FEA::kBackward_Euler );
			fEquation_eps -> Form_LHS_Keps_Kd ( fKepseps, 	fKepsd );
			fEquation_eps -> Form_RHS_F_int ( fFeps_int );
			fFeps_int *= -1.0;
			
			// update state variables
			np1.Update( APS::kstate, fstate );
			Convert.Copy ( fNumIP, knum_d_state, fstate, fdstatenew_all );

			/* equations numbers */
			const iArrayT& all_eq = CurrentElement().Equations();
			displ_eq.Set(fFd_int.Length(), all_eq.Pointer());
			plast_eq.Set(fFeps_int.Length(), all_eq.Pointer(fFd_int.Length()));

			/* assemble residuals */
			ElementSupport().AssembleRHS(curr_group, fFd_int, displ_eq);
			ElementSupport().AssembleRHS(curr_group, fFeps_int, plast_eq);

			/* assemble components of the tangent */
			ElementSupport().AssembleLHS(curr_group, fKdd, displ_eq);
			ElementSupport().AssembleLHS(curr_group, fKepseps, plast_eq);
			ElementSupport().AssembleLHS(curr_group, fKdeps, displ_eq, plast_eq);
			ElementSupport().AssembleLHS(curr_group, fKepsd, plast_eq, displ_eq);
		}
	}	
}

