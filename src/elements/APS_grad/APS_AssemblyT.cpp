/* $Id: APS_AssemblyT.cpp,v 1.60 2004-08-02 23:54:07 raregue Exp $ */
#include "APS_AssemblyT.h"

#include "APS_MatlT.h"
#include "Shear_MatlT.h"

#include "OutputSetT.h"
#include "ParameterContainerT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

//---------------------------------------------------------------------

/* constructor */
APS_AssemblyT::APS_AssemblyT(const ElementSupportT& support):
	ElementBaseT(support), //pass the displacement field to the base class
	u(LocalArrayT::kDisp),
	u_n(LocalArrayT::kLastDisp),
	gamma_p(LocalArrayT::kDisp),
	gamma_p_n(LocalArrayT::kLastDisp),
	fInitCoords_displ(LocalArrayT::kInitCoords),
	fCurrCoords_displ(LocalArrayT::kCurrCoords),
	fInitCoords_plast(LocalArrayT::kInitCoords),
	fCurrCoords_plast(LocalArrayT::kCurrCoords),
	fTractionBCSet(0),
	fDispl(NULL),
	fPlast(NULL),
	fShapes_displ(NULL),
	fShapes_plast(NULL),
	fBalLinMomMaterial(NULL),
	fPlastMaterial(NULL),
	fKdd(ElementMatrixT::kNonSymmetric),
	fKdd_face(ElementMatrixT::kNonSymmetric),
	fKdeps(ElementMatrixT::kNonSymmetric),
	fKepsd(ElementMatrixT::kNonSymmetric),
	fKepseps(ElementMatrixT::kNonSymmetric),
	fEquation_d(NULL),
	fEquation_eps(NULL),
	bStep_Complete(0)
{
	SetName("antiplane_shear_grad_plast");
}

//---------------------------------------------------------------------

/* destructor */
APS_AssemblyT::~APS_AssemblyT(void) 
{  
	delete fDispl;
	delete fPlast;
	delete fEquation_d; 
	delete fEquation_eps; 
	delete fShapes_displ;
	delete fShapes_plast; 
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

	cout << "iPlastModelType " 						<< iPlastModelType 			<< endl; 
	
	//-- Elasticity parameters 
	cout << "fMaterial_Data[kMu] "  				<< fMaterial_Data[kMu] 		<< endl;

	//-- Plasticity parameters
	cout << "fMaterial_Data[km_rate] " 				<< fMaterial_Data[km_rate] 	<< endl;
	cout << "fMaterial_Data[kgamma0_dot_1] " 		<< fMaterial_Data[kgamma0_dot_1] << endl;
	cout << "fMaterial_Data[kgamma0_dot_2] " 		<< fMaterial_Data[kgamma0_dot_2] << endl;
	cout << "fMaterial_Data[kgamma0_dot_3] " 		<< fMaterial_Data[kgamma0_dot_3] << endl;
	cout << "fMaterial_Data[km1_x] " 				<< fMaterial_Data[km1_x] 		<< endl;
	cout << "fMaterial_Data[km1_y] " 				<< fMaterial_Data[km1_y] 		<< endl;
	cout << "fMaterial_Data[km2_x] " 				<< fMaterial_Data[km2_x] 		<< endl;
	cout << "fMaterial_Data[km2_y] " 				<< fMaterial_Data[km2_y] 		<< endl;
	cout << "fMaterial_Data[km3_x] " 				<< fMaterial_Data[km3_x] 		<< endl;
	cout << "fMaterial_Data[km3_y] " 				<< fMaterial_Data[km3_y] 		<< endl;

	//-- Backstress Parameter
	cout << "fMaterial_Data[kl] " 					<< fMaterial_Data[kl]			<< endl;

	//-- Isotropic Hardening Parameter
	cout << "fMaterial_Data[kH] "					<< fMaterial_Data[kH]			<< endl;
	
	//-- Initial state variables
	cout << "fMaterial_Data[kkappa0_1] "			<< fMaterial_Data[kkappa0_1]	<< endl;
	cout << "fMaterial_Data[kkappa0_2] "			<< fMaterial_Data[kkappa0_2]	<< endl;
	cout << "fMaterial_Data[kkappa0_3] "			<< fMaterial_Data[kkappa0_3]	<< endl;
	
}


//---------------------------------------------------------------------

void APS_AssemblyT::RHSDriver(void)	
{
	int curr_group = ElementSupport().CurrentGroup();

	/* traction boundary conditions acting on displacement equations */
	if (curr_group == fDispl->Group()) 
		ApplyTractionBC();

	/* choose solution method */
	if (fDispl->Group() == fPlast->Group())
	  RHSDriver_monolithic();
	 else
	  RHSDriver_staggered();
}
//---------------------------------------------------------------------

void APS_AssemblyT::Equations(AutoArrayT<const iArray2DT*>& eq_d,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_eps)
{
	/* doing monolithic solution */
	if (fDispl->Group() == fPlast->Group())
	{
		int ndof_plast = fPlast->NumDOF();
		int ndof_displ = fDispl->NumDOF();
	
		/* loop over connectivity blocks */
		fEqnos_displ.Dimension(fEqnos.Length());
		fEqnos_plast.Dimension(fEqnos.Length());
		for (int i = 0; i < fEqnos.Length(); i++)
		{
			/* connectivities */
			const iArray2DT& connects_displ = *(fConnectivities_displ[i]);
			const iArray2DT& connects_plast = *(fConnectivities_plast[i]);
			int nel = connects_displ.MajorDim();
		
			/* dimension */ 
			fEqnos[i].Dimension(nel, n_en_displ*ndof_displ + n_en_plast*ndof_plast);
			iArray2DT& displ_eq = fEqnos_displ[i];
			iArray2DT& plast_eq = fEqnos_plast[i];
			displ_eq.Dimension(nel, n_en_displ*ndof_displ);
			plast_eq.Dimension(nel, n_en_plast*ndof_plast);
			
			/* get equation numbers */
			fDispl->SetLocalEqnos(connects_displ, displ_eq);
			fPlast->SetLocalEqnos(connects_plast, plast_eq);
			
			/* write into one array */
			fEqnos[i].BlockColumnCopyAt(displ_eq, 0);
			fEqnos[i].BlockColumnCopyAt(plast_eq, displ_eq.MinorDim());

			/* add to list of equation numbers */
			eq_d.Append(&fEqnos[i]);
		}
	
		/* reset pointers to element cards */
		SetElementCards(fBlockData, fConnectivities_displ, fEqnos_displ, fElementCards_displ);
		SetElementCards(fBlockData, fConnectivities_plast, fEqnos_plast, fElementCards_plast);
	}
	else
	{
#pragma message("correct initialization for staggered solution")
	
		/* ElementBaseT handles equation array for displacements */
		if (ElementSupport().CurrentGroup() == fDispl->Group())
			ElementBaseT::Equations(eq_d, eq_eps);

		/* plasticity equation */
		if (ElementSupport().CurrentGroup() == fPlast->Group())
		{
			/* collect local equation numbers */
			//fPlast.SetLocalEqnos(fConnectivities_plast, fEqnos_plast);
		
			//eq_d.Append(&fEqnos_plast);
		}
	}
	
	/* get the equation number for the nodes on the faces */
	for (int i = 0; i < fPlasticGradientFaceEqnos.Length(); i++)
	{
		iArray2DT& faces = fPlasticGradientFaces[i];
		iArray2DT& eqnos = fPlasticGradientFaceEqnos[i];
		eqnos.Dimension(faces.MajorDim(), faces.MajorDim()*fDispl->NumDOF());
	
		fDispl->SetLocalEqnos(faces, eqnos);
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
			fBalLinMomMaterial->Assign ( Shear_MatlT::kMu, fMaterial_Data[kMu] );
			break;

		default :
			cout << "APS_AssemblyT::Select_Equations() .. ERROR >> bad iBalScale \n";
			break;
	}

	/** Choices for Fine-Scale Equation */

	switch ( iPlastScale )	{

		case PlastT::kAPS_kappa_alpha :
			fEquation_eps	= new APS_kappa_alphaT;
			fPlastMaterial	= new APS_MatlT;		
			fPlastMaterial->Assign (	APS_MatlT::kMu, 		fMaterial_Data[kMu] 		); 	
			fPlastMaterial->Assign ( 	APS_MatlT::km_rate, 	fMaterial_Data[km_rate] 	); 	
			fPlastMaterial->Assign ( 	APS_MatlT::kgamma0_dot_1, fMaterial_Data[kgamma0_dot_1]); 
			fPlastMaterial->Assign ( 	APS_MatlT::kgamma0_dot_2, fMaterial_Data[kgamma0_dot_2]);
			fPlastMaterial->Assign ( 	APS_MatlT::kgamma0_dot_3, fMaterial_Data[kgamma0_dot_3]); 	
			fPlastMaterial->Assign ( 	APS_MatlT::km1_x, 		fMaterial_Data[km1_x] 		); 
			fPlastMaterial->Assign ( 	APS_MatlT::km1_y, 		fMaterial_Data[km1_y] 		); 	
			fPlastMaterial->Assign ( 	APS_MatlT::km2_x, 		fMaterial_Data[km2_x] 		); 
			fPlastMaterial->Assign ( 	APS_MatlT::km2_y, 		fMaterial_Data[km2_y] 		); 
			fPlastMaterial->Assign ( 	APS_MatlT::km3_x, 		fMaterial_Data[km3_x] 		); 
			fPlastMaterial->Assign ( 	APS_MatlT::km3_y, 		fMaterial_Data[km3_y] 		); 	
			fPlastMaterial->Assign ( 	APS_MatlT::kl, 			fMaterial_Data[kl] 			); 	
			fPlastMaterial->Assign ( 	APS_MatlT::kH, 			fMaterial_Data[kH] 			);
			fPlastMaterial->Assign ( 	APS_MatlT::kkappa0_1, 	fMaterial_Data[kkappa0_1] 	); 	
			fPlastMaterial->Assign ( 	APS_MatlT::kkappa0_2, 	fMaterial_Data[kkappa0_2] 	); 	
			fPlastMaterial->Assign ( 	APS_MatlT::kkappa0_3, 	fMaterial_Data[kkappa0_3] 	); 	
			break;

		case PlastT::kAPS_kappa_gp :
			fEquation_eps	= new APS_kappa_gpT;
			fPlastMaterial	= new APS_MatlT;		
			fPlastMaterial->Assign (	APS_MatlT::kMu, 		fMaterial_Data[kMu] 		); 	
			fPlastMaterial->Assign ( 	APS_MatlT::km_rate, 	fMaterial_Data[km_rate] 	); 	
			fPlastMaterial->Assign ( 	APS_MatlT::kgamma0_dot_1, fMaterial_Data[kgamma0_dot_1]); 
			fPlastMaterial->Assign ( 	APS_MatlT::kgamma0_dot_2, fMaterial_Data[kgamma0_dot_2]);
			fPlastMaterial->Assign ( 	APS_MatlT::kgamma0_dot_3, fMaterial_Data[kgamma0_dot_3]); 	
			fPlastMaterial->Assign ( 	APS_MatlT::km1_x, 		fMaterial_Data[km1_x] 		); 
			fPlastMaterial->Assign ( 	APS_MatlT::km1_y, 		fMaterial_Data[km1_y] 		); 	
			fPlastMaterial->Assign ( 	APS_MatlT::km2_x, 		fMaterial_Data[km2_x] 		); 
			fPlastMaterial->Assign ( 	APS_MatlT::km2_y, 		fMaterial_Data[km2_y] 		); 
			fPlastMaterial->Assign ( 	APS_MatlT::km3_x, 		fMaterial_Data[km3_x] 		); 
			fPlastMaterial->Assign ( 	APS_MatlT::km3_y, 		fMaterial_Data[km3_y] 		); 	
			fPlastMaterial->Assign ( 	APS_MatlT::kl, 			fMaterial_Data[kl] 			); 	
			fPlastMaterial->Assign ( 	APS_MatlT::kH, 			fMaterial_Data[kH] 			);
			fPlastMaterial->Assign ( 	APS_MatlT::kkappa0_1, 	fMaterial_Data[kkappa0_1] 	); 	
			fPlastMaterial->Assign ( 	APS_MatlT::kkappa0_2, 	fMaterial_Data[kkappa0_2] 	); 	
			fPlastMaterial->Assign ( 	APS_MatlT::kkappa0_3, 	fMaterial_Data[kkappa0_3] 	); 	
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
	return group == fDispl->Group() ||
	       group == fPlast->Group();
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


void APS_AssemblyT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented
}


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

	/* displ, plast, or neither */
	bool is_displ = false;
	dArrayT* element_force = NULL;
	int num_force = 0;
	if (field.Name() == fDispl->Name()) {
		is_displ = true;
		element_force = &fFd_int;
		num_force = fDispl->NumDOF();
		}
	else if (field.Name() == fPlast->Name()) {
		is_displ = false;
		element_force = &fFeps_int;
		num_force = fPlast->NumDOF();
		}
	else
		return;

	/* time Step Increment */
	double delta_t = ElementSupport().TimeStep();
	time = ElementSupport().Time();
	step_number = ElementSupport().StepNumber();

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

	 	// coordinates do not change for anti-plane shear
		//fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, u); 
	 	SetLocalX(fInitCoords_displ); 
		fCurrCoords_displ = fInitCoords_displ;
		fShapes_displ->SetDerivatives(); 
		//
		SetLocalX(fInitCoords_plast); 
		fCurrCoords_plast = fInitCoords_plast;
		fShapes_plast->SetDerivatives(); 
		
		//update state variables
		fdstatenew_all.Alias(fNumIP_plast, knum_d_state, fdState_new(CurrElementNumber()));
		fdstate_all.Alias(fNumIP_plast, knum_d_state, fdState(CurrElementNumber()));
		
		/** repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradients 		( fShapes_displ, 	u, u_n, fgrad_u, fgrad_u_n );
		Convert.Gradients 		( fShapes_plast, 	gamma_p, gamma_p_n, fgrad_gamma_p, fgrad_gamma_p_n );
		Convert.Interpolate 	( fShapes_plast, 	gamma_p, gamma_p_n, fgamma_p, fgamma_p_n );
		Convert.Shapes			(	fShapes_displ, fFEA_Shapes_displ );
		Convert.Shapes			(	fShapes_plast, fFEA_Shapes_plast );
		Convert.Displacements	(	del_u, 	del_u_vec  );
		Convert.Displacements	(	del_gamma_p, 	del_gamma_p_vec  );
		Convert.Na				(	n_en_displ, fShapes_displ, 	fFEA_Shapes_displ );
		Convert.Na				(	n_en_plast, fShapes_plast, 	fFEA_Shapes_plast );
		Convert.Copy			(	fNumIP_plast, knum_d_state, fdstatenew_all, fstate );
		Convert.Copy			(	fNumIP_plast, knum_d_state, fdstate_all, fstate_n );
		
		APS_VariableT np1(	fgrad_u, fgrad_u_surf, fgamma_p, fgamma_p_surf, fgrad_gamma_p, fstate ); // variables at time-step n+1
		APS_VariableT   n(	fgrad_u_n, fgrad_u_surf_n,fgamma_p_n, fgamma_p_surf_n, fgrad_gamma_p_n, fstate_n );	// variables at time-step n	 

			/* calculate displacement nodal force */
			if (is_displ)
			{
				/* residual and tangent for displacement field */
				fEquation_d->Construct ( fNumIPSurf_displ, n_en_surf, fFEA_Shapes_displ, fFEA_Shapes_plast, fBalLinMomMaterial, 
										fPlastMaterial, np1, n, step_number, delta_t );
				fEquation_d->Form_LHS_Keps_Kd ( fKdeps, fKdd );
				fEquation_d->Form_RHS_F_int ( fFd_int, np1 );
				fFd_int *= -1.0;  
			}
			else /* plasticity nodal force */
			{
				/* residual and tangent for plastic gradient field */
				fEquation_eps->Construct ( fFEA_Shapes_displ, fFEA_Shapes_plast, fPlastMaterial, np1, n, 
											step_number, delta_t, FEA::kBackward_Euler );
				fEquation_eps->Form_LHS_Keps_Kd ( fKepseps, fKepsd );
				fEquation_eps->Form_RHS_F_int ( fFeps_int );
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
	ArrayT<StringT> e_labels(fNumIP_plast*(knumstrain+knumstress+knum_d_state));

	/* over integration points */
	const char* slabels2D[] = {"gamma_x", "gamma_y", "gammap_curl", "effstr", "s_xz", "s_yz", "J2"};
	const char* svlabels2D[] = {"xi_1", "kappa_1", "gamma_dot_1", "xi_2", "kappa_2", "gamma_dot_2", "xi_3", "kappa_3", "gamma_dot_3"};
	int count = 0;
	for (int j = 0; j < fNumIP_plast; j++)
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
	int num_node_output = fDispl->NumDOF() + fPlast->NumDOF() + knumstrain + knumstress + knum_d_state;
	ArrayT<StringT> n_labels(num_node_output);
	count = 0;

	/* labels from plastic gradient */
	const ArrayT<StringT>& plast_labels = fPlast->Labels();
	for (int i = 0; i < plast_labels.Length(); i++)
		n_labels[count++] = plast_labels[i];

	/* labels from displacement */
	const ArrayT<StringT>& displ_labels = fDispl->Labels();
	for (int i = 0; i < displ_labels.Length(); i++)
		n_labels[count++] = displ_labels[i];

	/* labels from strains and stresses at the nodes */
	for (int i = 0; i < knumstrain+knumstress; i++)
		n_labels[count++] = slabels2D[i];
		
	/* labels from state variables at the nodes */
	for (int i = 0; i < knum_d_state; i++)
		n_labels[count++] = svlabels2D[i];

	/* set output specifier */
	#pragma message("APS_AssemblyT::RegisterOutput: is this right? ")
	OutputSetT output_set(fGeometryCode_displ, block_ID, fConnectivities, n_labels, e_labels, false);
		
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
		out_variable_all.Alias(fNumIP_plast, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
		fShapes_displ->TopIP();
		while (fShapes_displ->NextIP())
		{
			out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(fShapes_displ->CurrIP()));
			fShapes_displ->Extrapolate(out_variable, nd_var);
		}
	
		/* accumulate - extrapolation done from ip's to corners => X nodes  */
		ElementSupport().AssembleAverage(CurrentElement().NodesX(), nd_var);
	}

	/* get nodally averaged values */
	dArray2DT extrap_values;
	ElementSupport().OutputUsedAverage(extrap_values);

	/* temp space for group displacements */
	int num_node_output = fDispl->NumDOF() + fPlast->NumDOF() + knumstrain + knumstress + knum_d_state;
	dArray2DT n_values(nodes_used.Length(), num_node_output);

	/* collect nodal values */
	const dArray2DT& fGamma_p = (*fPlast)[0];
	const dArray2DT& fU = (*fDispl)[0];
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
	if (fDispl->Group() == fPlast->Group())
		ExceptionT::GeneralFail(caller, "displacement and plastic gradient groups must be different: %d == %d",
			fDispl->Group(), fPlast->Group());

	int curr_group = ElementSupport().CurrentGroup();

	/* stress output work space */
	dArray2DT out_variable_all, fdstatenew_all, fdstate_all;
	dArrayT out_variable;

	/* time Step Increment */
	double delta_t = ElementSupport().TimeStep();
	time = ElementSupport().Time();
	step_number = ElementSupport().StepNumber();
	
	/* loop over elements */
	int e,l;
	Top();
	while (NextElement())
	{
		e = CurrElementNumber();
		const iArrayT& nodes_displ = fElementCards_displ[e].NodesU();
		const iArrayT& nodes_plast = fElementCards_plast[e].NodesU();

		u.SetLocal(nodes_displ);
		u_n.SetLocal(nodes_displ);
		gamma_p.SetLocal(nodes_plast);
		gamma_p_n.SetLocal(nodes_plast);

		del_u.DiffOf (u, u_n);
		del_gamma_p.DiffOf (gamma_p, gamma_p_n);
		
		// coordinates do not change for anti-plane shear
		//fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, u); 
	 	fInitCoords_displ.SetLocal(nodes_displ);
		fCurrCoords_displ = fInitCoords_displ;
		fShapes_displ->SetDerivatives(); 
		//
		fInitCoords_plast.SetLocal(nodes_plast); 
		fCurrCoords_plast = fInitCoords_plast;
		fShapes_plast->SetDerivatives(); 
		
		//update state variables
		fdstatenew_all.Alias(fNumIP_plast, knum_d_state, fdState_new(CurrElementNumber()));
		fdstate_all.Alias(fNumIP_plast, knum_d_state, fdState(CurrElementNumber()));
		
		/** repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradients 		( fShapes_displ, 	u, u_n, fgrad_u, fgrad_u_n );
		Convert.Gradients 		( fShapes_plast, 	gamma_p, gamma_p_n, fgrad_gamma_p, fgrad_gamma_p_n );
		Convert.Interpolate 	( fShapes_plast, 	gamma_p, gamma_p_n, fgamma_p, fgamma_p_n );
		Convert.Shapes			(	fShapes_displ, fFEA_Shapes_displ );
		Convert.Shapes			(	fShapes_plast, fFEA_Shapes_plast );
		Convert.Displacements	(	del_u, 	del_u_vec  );
		Convert.Displacements	(	del_gamma_p, 	del_gamma_p_vec  );
		Convert.Na				(	n_en_displ, fShapes_displ, 	fFEA_Shapes_displ );
		Convert.Na				(	n_en_plast, fShapes_plast, 	fFEA_Shapes_plast );
		Convert.Copy			(	fNumIP_plast, knum_d_state, fdstatenew_all, fstate );
		Convert.Copy			(	fNumIP_plast, knum_d_state, fdstate_all, fstate_n );
		
		// variables at time-step n+1
		APS_VariableT np1(	fgrad_u, fgrad_u_surf, fgamma_p, fgamma_p_surf, fgrad_gamma_p, fstate ); 
		// variables at time-step n
		APS_VariableT   n(	fgrad_u_n, fgrad_u_surf_n,fgamma_p_n, fgamma_p_surf_n, fgrad_gamma_p_n, fstate_n );		 

		/* which field */
	  	//SolverGroup 1 (gets field 1) <-- u (obtained by a rearranged Equation_d)
		if ( curr_group == fDispl->Group()  )	
		{

			if (bStep_Complete) {
			//do nothing here
			}
			else { //-- Still Iterating

				/** Compute N-R matrix equations */
				fEquation_d->Construct ( fNumIPSurf_displ, n_en_surf, fFEA_Shapes_displ, fFEA_Shapes_plast, fBalLinMomMaterial, 
											fPlastMaterial, np1, n, step_number, delta_t );
				fEquation_d->Form_LHS_Keps_Kd ( fKdeps, fKdd );
				fEquation_d->Form_RHS_F_int ( fFd_int, np1 );

				/** Set displacement LHS */
				fLHS = fKdd;

				/** Compute displacement RHS */
				fKdeps.Multx ( del_gamma_p_vec, fRHS );
				fRHS += fFd_int; 
				fRHS *= -1.0; 

				/** Compute Traction B.C. */
				Get_Fd_ext ( fFd_ext );
				fRHS += fFd_ext;
			
				/* add to global equations */
				ElementSupport().AssembleLHS ( fDispl->Group(), fLHS, CurrentElement().Equations() );
				ElementSupport().AssembleRHS ( fDispl->Group(), fRHS, CurrentElement().Equations() );
			}
		}

		// SolverGroup 2 (gets field 2) <-- gamma_p (obtained by a rearranged Equation_eps)
		else if (curr_group == fPlast->Group() )	
		{

			if (bStep_Complete) { 
			
				fEquation_eps->Construct ( fFEA_Shapes_displ, fFEA_Shapes_plast, fPlastMaterial, np1, n, step_number, delta_t );
				fEquation_eps->Get ( output, Render_Vector[e][0] );
			
				//-- Store/Register data in classic tahoe manner 
				out_variable_all.Alias(fNumIP_plast, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
				for (l=0; l < fNumIP_plast; l++) {
					out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(l));
					out_variable=Render_Vector[e][0][l];
					} 
			
			}
			else { //-- Still Iterating

				/** Compute N-R matrix equations */
				fEquation_eps->Construct ( fFEA_Shapes_displ, fFEA_Shapes_plast, fPlastMaterial, np1, n, 
											step_number, delta_t, FEA::kBackward_Euler );
				fEquation_eps->Form_LHS_Keps_Kd ( fKepseps, fKepsd );
				fEquation_eps->Form_RHS_F_int ( fFeps_int );
				
				// update state variables
				np1.Update( APS::kstate, fstate );		
				Convert.Copy ( fNumIP_plast, knum_d_state, fstate, fdstatenew_all );

				/** Set LHS */
				fLHS = fKepseps;	
		
				/** Compute plasticity RHS  */
				fKepsd.Multx ( del_u_vec, fRHS );
				fRHS += fFeps_int; 
				fRHS *= -1.0; 
		
				/* plastic gradient equation numbers */
				const iArrayT& plast_eq = fElementCards_plast[e].Equations();

				/* add to global equations */
				ElementSupport().AssembleLHS ( fPlast->Group(), fLHS, plast_eq );
				ElementSupport().AssembleRHS ( fPlast->Group(), fRHS, plast_eq );
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
	if (fDispl->Group() != fPlast->Group())
		ExceptionT::GeneralFail(caller, "displacement and plastic gradient groups must be the same: %d != %d",
			fDispl->Group(), fPlast->Group());

	int curr_group = ElementSupport().CurrentGroup();

	/* stress output work space */
	dArray2DT	out_variable_all, fdstatenew_all, fdstate_all;
	dArrayT		out_variable;

	/* time Step Increment */
	double delta_t = ElementSupport().TimeStep();
	time = ElementSupport().Time();
	step_number = ElementSupport().StepNumber();

	/* work space for integration over faces */
	LocalArrayT face_coords(LocalArrayT::kInitCoords, n_en_surf, NumSD());
	ElementSupport().RegisterCoordinates(face_coords);
	iArrayT face_nodes, face_equations;
	dMatrixT face_jacobian(NumSD(), NumSD()-1);
	dMatrixT face_Q(NumSD());
	LocalArrayT face_gamma_p(LocalArrayT::kDisp, n_en_surf, NumSD());
	fPlast->RegisterLocal(face_gamma_p);

	/* loop over elements */
	int e,l;
	Top();
	while (NextElement())
	{
		e = CurrElementNumber();
		const iArrayT& nodes_displ = fElementCards_displ[e].NodesU();
		const iArrayT& nodes_plast = fElementCards_plast[e].NodesU();

		u.SetLocal(nodes_displ);
		u_n.SetLocal(nodes_displ);
		gamma_p.SetLocal(nodes_plast);
		gamma_p_n.SetLocal(nodes_plast);

		del_u.DiffOf (u, u_n);
		del_gamma_p.DiffOf (gamma_p, gamma_p_n);
		
		// coordinates do not change for anti-plane shear
		//fCurrCoords.SetToCombination (1.0, fInitCoords, 1.0, u); 
	 	fInitCoords_displ.SetLocal(nodes_displ); 
		fCurrCoords_displ = fInitCoords_displ;
		fShapes_displ->SetDerivatives(); 
		//
		fInitCoords_plast.SetLocal(nodes_plast); 
		fCurrCoords_plast = fInitCoords_plast;
		fShapes_plast->SetDerivatives(); 
		
		//update state variables
		fdstatenew_all.Alias(fNumIP_plast, knum_d_state, fdState_new(CurrElementNumber()));
		fdstate_all.Alias(fNumIP_plast, knum_d_state, fdState(CurrElementNumber()));
		
		/** repackage data to forms compatible with FEA classes (very little cost in big picture) */
		Convert.Gradients 		( fShapes_displ, u, u_n, fgrad_u, fgrad_u_n );
		Convert.Gradients 		( fShapes_plast, gamma_p, gamma_p_n, fgrad_gamma_p, fgrad_gamma_p_n );
		Convert.Interpolate 	( fShapes_plast, gamma_p, gamma_p_n, fgamma_p, fgamma_p_n );
		Convert.Shapes			( fShapes_displ, fFEA_Shapes_displ );
		Convert.Shapes			( fShapes_plast, fFEA_Shapes_plast );
		Convert.Displacements	( del_u, del_u_vec  );
		Convert.Displacements	( del_gamma_p, del_gamma_p_vec  );
		Convert.Na				( n_en_displ, fShapes_displ, fFEA_Shapes_displ );
		Convert.Na				( n_en_plast, fShapes_plast, fFEA_Shapes_plast );
		Convert.Copy			( fNumIP_plast, knum_d_state, fdstatenew_all, fstate );
		Convert.Copy			( fNumIP_plast, knum_d_state, fdstate_all, fstate_n );
		
		// variables at time-step n+1
		APS_VariableT np1(	fgrad_u, fgrad_u_surf, fgamma_p, fgamma_p_surf, fgrad_gamma_p, fstate ); 
		// variables at time-step n
		APS_VariableT   n(	fgrad_u_n, fgrad_u_surf_n,fgamma_p_n, fgamma_p_surf_n, fgrad_gamma_p_n, fstate_n );		 

				
		if (bStep_Complete) { 
		
			fEquation_eps->Construct ( fFEA_Shapes_displ, fFEA_Shapes_plast, fPlastMaterial, np1, n, step_number, delta_t );
			fEquation_eps->Get ( output, Render_Vector[e][0] );
			
			//-- Store/Register data in classic tahoe manner 
			out_variable_all.Alias(fNumIP_plast, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
			for (l=0; l < fNumIP_plast; l++) {
				out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(l));
				out_variable=Render_Vector[e][0][l];
			} 
	
		}
		else { //-- Still Iterating
		
			/* residual and tangent for displacements */
			fEquation_d->Construct ( fNumIPSurf_displ, n_en_surf, fFEA_Shapes_displ, fFEA_Shapes_plast, fBalLinMomMaterial, 
									fPlastMaterial, np1, n, step_number, delta_t );
			fEquation_d->Form_LHS_Keps_Kd ( fKdeps, fKdd );
			fEquation_d->Form_RHS_F_int ( fFd_int, np1 );
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
						int face = fSideSetFaces[i][j];
						const ParentDomainT& surf_shape = fShapes_displ->FacetShapeFunction(face);
						const ParentDomainT& parent = ShapeFunction().ParentDomain();
						iArrayT face_local_nodes(2);
						parent.NodesOnFacet(face, face_local_nodes);
						
						/* equations for the nodes on the face */
						fPlasticGradientFaceEqnos[i].RowAlias(j, face_equations);
						
						Convert.SurfShapeGradient ( n_en_surf, surf_shape, fFEA_SurfShapes, face_coords,
													parent, fInitCoords_displ, *fShapes_displ, u, u_n, fgrad_u_surf, 
													fgrad_u_surf_n, face_gamma_p, fgamma_p_surf, face_local_nodes );
						APS_VariableT np1_surf(	fgrad_u, fgrad_u_surf, fgamma_p, fgamma_p_surf, fgrad_gamma_p, fstate ); 
						fEquation_d->Form_LHS_Kd_Surf ( fKdd_face, fFEA_SurfShapes );
						double wght = fPlasticGradientWght[i];
						fEquation_d->Form_RHS_F_int_Surf ( fFd_int_face, np1_surf, wght );
						
						ElementSupport().AssembleRHS(curr_group, fFd_int_face, face_equations);
						ElementSupport().AssembleLHS(curr_group, fKdd_face, face_equations);
					}
				}
			}

			/* residual and tangent for plasticity */
			fEquation_eps->Construct ( fFEA_Shapes_displ, fFEA_Shapes_plast, fPlastMaterial, np1, n, step_number, 
										delta_t, FEA::kBackward_Euler );
			fEquation_eps->Form_LHS_Keps_Kd ( fKepseps, fKepsd );
			fEquation_eps->Form_RHS_F_int ( fFeps_int );
			fFeps_int *= -1.0;
			
			// update state variables
			np1.Update( APS::kstate, fstate );
			Convert.Copy ( fNumIP_plast, knum_d_state, fstate, fdstatenew_all );

			/* equations numbers */
			const iArrayT& displ_eq = fElementCards_displ[e].Equations();
			const iArrayT& plast_eq = fElementCards_plast[e].Equations();

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



/* form global shape function derivatives */
void APS_AssemblyT::SetGlobalShape(void)
{
	/* fetch (initial) coordinates */
	SetLocalX(fLocInitCoords);
	
	/* compute shape function derivatives */
	fShapes_displ->SetDerivatives();
	fShapes_plast->SetDerivatives();
}



/* describe the parameters needed by the interface */
void APS_AssemblyT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

	/* displacement field */
	//already done in ElementBaseT
	//list.AddParameter(ParameterT::Word, "displ_field_name");
	
	/* plastic gradient field */
	list.AddParameter(ParameterT::Word, "plastic_grad_field_name");
		
	list.AddParameter(fGeometryCode_displ_int, "GeometryCode_displ");
	list.AddParameter(fNumIP_displ, "NumIP_displ");
	list.AddParameter(fGeometryCodeSurf_displ_int, "GeometryCodeSurf_displ");
	list.AddParameter(fNumIPSurf_displ, "NumIPSurf_displ");
	list.AddParameter(n_en_displ, "n_en_displ");
	list.AddParameter(n_en_plast, "n_en_plast");
	
	list.AddParameter(iPlastModelType, "plast_mod_type");
			
	/*		
	// elasticity
	list.AddParameter(fMaterial_Data[kMu], "shear_modulus");
	
	// plasticity
	list.AddParameter(fMaterial_Data[km_rate], "rate_param");
	list.AddParameter(fMaterial_Data[kgamma0_dot_1], "gamma0_dot_1");
	list.AddParameter(fMaterial_Data[kgamma0_dot_2], "gamma0_dot_2");
	list.AddParameter(fMaterial_Data[kgamma0_dot_3], "gamma0_dot_3");
	list.AddParameter(fMaterial_Data[km1_x], "m1_x");
	list.AddParameter(fMaterial_Data[km1_y], "m1_y");
	list.AddParameter(fMaterial_Data[km2_x], "m2_x");
	list.AddParameter(fMaterial_Data[km2_y], "m2_y");
	list.AddParameter(fMaterial_Data[km3_x], "m3_x");
	list.AddParameter(fMaterial_Data[km3_y], "m3_y");

	//-- length scale, Backstress Parameter
	list.AddParameter(fMaterial_Data[kl], "length_scale");

	//-- Isotropic Hardening Parameter
	list.AddParameter(fMaterial_Data[kH], "H_param");
	
	//-- Initial values of state variable along various slip systems
	list.AddParameter(fMaterial_Data[kkappa0_1], "kappa0_1");
	list.AddParameter(fMaterial_Data[kkappa0_2], "kappa0_2");
	list.AddParameter(fMaterial_Data[kkappa0_3], "kappa0_3");
	*/
	
	double shearMu, m_rate, gamma0_dot_1, gamma0_dot_2, gamma0_dot_3,
			m1_x,m1_y,m2_x,m2_y,m3_x,m3_y, length_scale, H_param, kappa0_1,
			kappa0_2,kappa0_3;
			
	// elasticity
	list.AddParameter(shearMu, "shear_modulus");
	
	// plasticity
	list.AddParameter(m_rate, "rate_param");
	list.AddParameter(gamma0_dot_1, "gamma0_dot_1");
	list.AddParameter(gamma0_dot_2, "gamma0_dot_2");
	list.AddParameter(gamma0_dot_3, "gamma0_dot_3");
	list.AddParameter(m1_x, "m1_x");
	list.AddParameter(m1_y, "m1_y");
	list.AddParameter(m2_x, "m2_x");
	list.AddParameter(m2_y, "m2_y");
	list.AddParameter(m3_x, "m3_x");
	list.AddParameter(m3_y, "m3_y");

	//-- length scale, Backstress Parameter
	list.AddParameter(length_scale, "length_scale");

	//-- Isotropic Hardening Parameter
	list.AddParameter(H_param, "H_param");
	
	//-- Initial values of state variable along various slip systems
	list.AddParameter(kappa0_1, "kappa0_1");
	list.AddParameter(kappa0_2, "kappa0_2");
	list.AddParameter(kappa0_3, "kappa0_3");
	
	// number of side sets
	list.AddParameter(num_sidesets, "num_sidesets");

	//doesn't appear to work
	/*
	int num_sidesets_tmp=5; //dummy memory??
	fPlasticGradientWght.Dimension(num_sidesets_tmp);
	
 	list.AddParameter(ParameterT::Word, "fSideSetID[0]");
 	list.AddParameter(fPlasticGradientWght[0], "fPlasticGradientWght[0]");
 	
 	list.AddParameter(ParameterT::Word, "fSideSetID[1]");
 	list.AddParameter(fPlasticGradientWght[1], "fPlasticGradientWght[1]");
 	
 	list.AddParameter(ParameterT::Word, "fSideSetID[2]");
 	list.AddParameter(fPlasticGradientWght[2], "fPlasticGradientWght[2]");
 	
 	list.AddParameter(ParameterT::Word, "fSideSetID[3]");
 	list.AddParameter(fPlasticGradientWght[3], "fPlasticGradientWght[3]");
 	
 	list.AddParameter(ParameterT::Word, "fSideSetID[4]");
 	list.AddParameter(fPlasticGradientWght[4], "fPlasticGradientWght[4]");
 	*/

}


/* accept parameter list */
void APS_AssemblyT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "APS_AssemblyT::TakeParameterList";
	
	/* inherited */
	ElementBaseT::TakeParameterList(list);

	/* get form of tangent */
	GlobalT::SystemTypeT type = TangentType();
	
	/* set form of element stiffness matrix */
	if (type == GlobalT::kSymmetric)
		fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
	else if (type == GlobalT::kNonSymmetric)
		fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
	else if (type == GlobalT::kDiagonal)
		fLHS.SetFormat(ElementMatrixT::kDiagonal);
	
	/* get displacement field */
	/*
	const StringT& displ_field_name = list.GetParameter("displ_field_name");
	fDispl = ElementSupport().Field(displ_field_name);
	if (!fDispl)
		ExceptionT::GeneralFail(caller, "could not resolve \"%s\" displ_field", 
		displ_field_name.Pointer());
		*/
	const StringT& displ_field_name = list.GetParameter("field_name");
	fDispl = ElementSupport().Field(displ_field_name);
	if (!fDispl)
		ExceptionT::GeneralFail(caller, "could not resolve \"%s\" displ_field", 
		displ_field_name.Pointer());	

	/* get plastic gradient field */
	const StringT& plastic_grad_field_name = list.GetParameter("plastic_grad_field_name");
	fPlast = ElementSupport().Field(plastic_grad_field_name);
	if (!fPlast)
		ExceptionT::GeneralFail(caller, "could not resolve \"%s\" plastic_grad_field", 
		plastic_grad_field_name.Pointer());

	fGeometryCode_displ_int = list.GetParameter("GeometryCode_displ");
	fGeometryCode_displ = GeometryT::int2CodeT(fGeometryCode_displ_int);
	fNumIP_displ = list.GetParameter("NumIP_displ");
	fGeometryCodeSurf_displ_int = list.GetParameter("GeometryCodeSurf_displ");
	fGeometryCodeSurf_displ = GeometryT::int2CodeT(fGeometryCodeSurf_displ_int);
	fNumIPSurf_displ = list.GetParameter("NumIPSurf_displ");
	n_en_displ = list.GetParameter("n_en_displ");
	n_en_plast = list.GetParameter("n_en_plast");
	
	fGeometryCode_plast = fGeometryCode_displ; 
	fNumIP_plast = fNumIP_displ;
	fGeometryCodeSurf_plast = fGeometryCodeSurf_displ;
	fNumIPSurf_plast = fNumIPSurf_displ;
	
	iPlastModelType = list.GetParameter("plast_mod_type");
	
	fMaterial_Data.Dimension ( kNUM_FMAT_TERMS );
	
	fMaterial_Data[kMu] = list.GetParameter("shear_modulus");
	
	fMaterial_Data[km_rate] = list.GetParameter("rate_param");
	fMaterial_Data[kgamma0_dot_1] = list.GetParameter("gamma0_dot_1");
	fMaterial_Data[kgamma0_dot_2] = list.GetParameter("gamma0_dot_2");
	fMaterial_Data[kgamma0_dot_3] = list.GetParameter("gamma0_dot_3");
	fMaterial_Data[km1_x] = list.GetParameter("m1_x");
	fMaterial_Data[km1_y] = list.GetParameter("m1_y");
	fMaterial_Data[km2_x] = list.GetParameter("m2_x");
	fMaterial_Data[km2_y] = list.GetParameter("m2_y");
	fMaterial_Data[km3_x] = list.GetParameter("m3_x");
	fMaterial_Data[km3_y] = list.GetParameter("m3_y");
	
	fMaterial_Data[kl] = list.GetParameter("length_scale");
	
	fMaterial_Data[kH] = list.GetParameter("H_param");
	
	fMaterial_Data[kkappa0_1] = list.GetParameter("kappa0_1");
	fMaterial_Data[kkappa0_2] = list.GetParameter("kappa0_2");
	fMaterial_Data[kkappa0_3] = list.GetParameter("kappa0_3");
	
	num_sidesets = list.GetParameter("num_sidesets");
	num_sidesets = 0;
	
	/* prescribed plastic gradient at surface */
	fSideSetID.Dimension(num_sidesets);
	fSideSetElements.Dimension(num_sidesets);
	fSideSetFaces.Dimension(num_sidesets);
	fPlasticGradientWght.Dimension(num_sidesets);
	fPlasticGradientFaces.Dimension(num_sidesets);
	fPlasticGradientFaceEqnos.Dimension(num_sidesets);
	
	// enable the model manager
	ModelManagerT& model = ElementSupport().ModelManager();

	// get sideset info for plastic gradient field
	/*
	fSideSetID[0] = list.GetParameter("fSideSetID[0]");
	fPlasticGradientWght[0] = list.GetParameter("fPlasticGradientWght[0]");
	fSideSetID[1] = list.GetParameter("fSideSetID[1]");
	fPlasticGradientWght[1] = list.GetParameter("fPlasticGradientWght[1]");
	fSideSetID[2] = list.GetParameter("fSideSetID[2]");
	fPlasticGradientWght[2] = list.GetParameter("fPlasticGradientWght[2]");
	fSideSetID[3] = list.GetParameter("fSideSetID[3]");
	fPlasticGradientWght[3] = list.GetParameter("fPlasticGradientWght[3]");
	fSideSetID[4] = list.GetParameter("fSideSetID[4]");
	fPlasticGradientWght[4] = list.GetParameter("fPlasticGradientWght[4]");
	*/

	int n_en_surf=0;
	ArrayT<GeometryT::CodeT> facet_geom;
	iArrayT facet_nodes;
	for (int i = 0; i < num_sidesets; i++)
	{
		// get nodes-on-faces
		model.SideSet(fSideSetID[i], facet_geom, facet_nodes, fPlasticGradientFaces[i]);
		// get number of facet nodes, assuming each facet has an equal number
		n_en_surf=facet_nodes[0];

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
	
	knum_d_state = 9; // double's needed per ip, state variables
	knum_i_state = 0; // int's needed per ip, state variables
	
	knumstrain = 4; // number of strain outputs
	knumstress = 3; // number of stress outputs
	
	output = "out";
	
	/* dimensions (notation as per Hughes' Book) */
	int& n_ip_displ = fNumIP_displ;
	int& n_ip_plast = fNumIP_plast;
	n_sd = NumSD();
	//n_df = NumDOF(); 
	//n_df = 1+n_sd; 		
	int nen = NumElementNodes(); /* number of nodes/element in the mesh */
	//n_en_plast = n_en_displ = nen;
	//n_np = ElementSupport().NumNodes();

	/* initialize connectivities */
	fConnectivities_displ.Alias(fConnectivities);
	fConnectivities_plast.Alias(fConnectivities);

	/* pick element interpolations based on available number of element nodes
	 * and the specified number of integration points */
	// only implemented for 2D, quadratic quads
	//if (n_sd == 2 && nen == 8 && fGeometryCode_displ == GeometryT::kQuadrilateral) 
	if (n_sd == 2 && n_en_plast != n_en_displ && fGeometryCode_displ == GeometryT::kQuadrilateral) 
	{
		// don't expect reduced integration for both fields 
		// if (n_ip_displ == 4 && n_ip_plast == 4)
		//	ExceptionT::GeneralFail(caller, "not expecting 4 ips for both fields");
		//else if (n_ip_displ == 4 || n_ip_plast == 4) // create reduced connectivities
		//{ 
			// reduce the number of element nodes based on the number ip's
			int& nen_red = (n_ip_displ == 4) ? n_en_displ : n_en_plast;
			nen_red = 4;
			ArrayT<const iArray2DT*>& connects_red = (n_ip_displ == 4) ? 
				fConnectivities_displ : 
				fConnectivities_plast;
		
			//create reduced connectivities
			connects_red.Dimension(0);
			connects_red.Dimension(fConnectivities.Length());
			fConnectivities_reduced.Dimension(fConnectivities.Length());
			for (int i = 0; i < fConnectivities_reduced.Length(); i++) {
				iArray2DT& connects_red_store = fConnectivities_reduced[i];
				const iArray2DT& connects = *(fConnectivities[i]);
				connects_red_store.Dimension(connects.MajorDim(), nen_red);				
				connects_red[i] = &connects_red_store;
				
				//take 1st four element nodes (columns)
				for (int j = 0; j < nen_red; j++)
					connects_red_store.ColumnCopy(j, connects, j);
			}
		//}
	}
	

	n_el = NumElements();	
	n_sd_surf = n_sd;
	
	/* set shape functions */
	// u
	fInitCoords_displ.Dimension(n_en_displ, n_sd);
	ElementSupport().RegisterCoordinates(fInitCoords_displ);	
	fCurrCoords_displ.Dimension(n_en_displ, n_sd);
	fShapes_displ = new ShapeFunctionT(fGeometryCode_displ, fNumIP_displ, fCurrCoords_displ);
	//fShapes_displ->Initialize();
	// gamma_p
	fInitCoords_plast.Dimension(n_en_plast, n_sd);
	ElementSupport().RegisterCoordinates(fInitCoords_plast);	
	fCurrCoords_plast.Dimension(n_en_plast, n_sd);
	fShapes_plast = new ShapeFunctionT(fGeometryCode_plast, fNumIP_plast, fCurrCoords_plast);
	//fShapes_plast = new ShapeFunctionT(fGeometryCode_plast, fNumIP_plast, fCurrCoords_displ);
	//fShapes_plast->Initialize();

	/* set local arrays for displacement field */
	int dum=1;
	u.Dimension (n_en_displ, dum);
	u_n.Dimension (n_en_displ, dum);
	del_u.Dimension (n_en_displ, dum);
	del_u_vec.Dimension (n_en_displ);
	//ElementSupport().RegisterCoordinates(fInitCoords_displ);
	fDispl->RegisterLocal(u);
	fDispl->RegisterLocal(u_n);

	/* set local arrays for plastic gradient field */
	gamma_p.Dimension (n_en_plast, n_sd);
	gamma_p_n.Dimension (n_en_plast, n_sd);
	del_gamma_p.Dimension (n_en_plast, n_sd);
	n_en_plast_x_n_sd = n_en_plast*n_sd;
	del_gamma_p_vec.Dimension (n_en_plast_x_n_sd);
	//ElementSupport().RegisterCoordinates(fInitCoords_plast);
	fPlast->RegisterLocal(gamma_p);
	fPlast->RegisterLocal(gamma_p_n);
	
	/* allocate state variable storage */
	// state variables are calculated at IPs for gamma_p field
	int num_ip = fNumIP_plast;
	fdState_new.Dimension(n_el, num_ip*knum_d_state);
	fdState.Dimension(n_el, num_ip*knum_d_state);
	fiState_new.Dimension(n_el, num_ip*knum_i_state);
	fiState.Dimension(n_el, num_ip*knum_i_state);
	
	/* initialize equations */
	fEqnos_displ.Alias(fEqnos_displ);
	fEqnos_plast.Dimension(fConnectivities_plast.Length());

	/* initialize state variables */
	fdState = 0;
	fdState_new = 0;
	fiState = 0;
	fiState_new = 0;

	/* initialize element cards */
	fElementCards_displ.Alias(fElementCards);
	fElementCards_plast.Dimension(fElementCards.Length());
	
	/* set cards to data in array - NOT NEEDED IF YOU'RE NOT
	 * GOING TO USE THE ElementCardT ARRAY? */
	for (int i= 0; i < fElementCards.Length(); i++)
		fElementCards[i].Set(fiState.MinorDim(), fiState(i), fdState.MinorDim(), fdState(i));
		                     
	/* allocate the global stack object (once) */
	extern FEA_StackT* fStack;
	if (!fStack) fStack = new FEA_StackT;

	Select_Equations ( BalLinMomT::kAPS_Bal_Eq, iPlastModelType );
	dum=knumstrain+knumstress;
	fEquation_eps->Initialize ( n_ip_plast, n_sd, n_en_displ, n_en_plast, 
					knum_d_state, dum, ElementSupport().StepNumber() );

	/* FEA Allocation */
	dum=1;
	fgrad_u.FEA_Dimension 			( fNumIP_displ, dum, n_sd );
	fgrad_u_surf.FEA_Dimension 		( fNumIPSurf_displ, dum, n_sd );
	fgamma_p.FEA_Dimension 			( fNumIP_plast, n_sd );
	//need gamma_p at the surface of the displ eqs
	fgamma_p_surf.FEA_Dimension 	( fNumIPSurf_displ, n_sd );
	fgrad_gamma_p.FEA_Dimension 	( fNumIP_plast, n_sd, n_sd );
	fgrad_u_n.FEA_Dimension 		( fNumIP_displ, dum, n_sd );
	fgrad_u_surf_n.FEA_Dimension 	( fNumIPSurf_displ, dum, n_sd );
	fgamma_p_n.FEA_Dimension 		( fNumIP_plast, n_sd );
	fgamma_p_surf_n.FEA_Dimension 	( fNumIPSurf_displ, n_sd );
	fgrad_gamma_p_n.FEA_Dimension 	( fNumIP_plast, n_sd,n_sd );
	
	fstate.FEA_Dimension 			( fNumIP_plast, knum_d_state );
	fstate_n.FEA_Dimension 			( fNumIP_plast, knum_d_state );

	fKdd.Dimension 			( n_en_displ, n_en_displ );
	fKdeps.Dimension 		( n_en_displ, n_en_plast_x_n_sd );
	fKepsd.Dimension 		( n_en_plast_x_n_sd, n_en_displ );
	fKepseps.Dimension 		( n_en_plast_x_n_sd, n_en_plast_x_n_sd );

	fFd_int.Dimension 		( n_en_displ );
	fFd_ext.Dimension 		( n_en_displ );
	fFeps_int.Dimension 	( n_en_plast_x_n_sd );
	fFeps_ext.Dimension 	( n_en_plast_x_n_sd );
	
	/* only allow this dimensioning if there are sidesets */
	if (n_en_surf > 0) {
		fKdd_face.Dimension 		( n_en_surf, n_en_surf );
		fFd_int_face.Dimension 		( n_en_surf );
		fFEA_SurfShapes.Construct	( fNumIPSurf_displ,n_sd_surf,n_en_surf );
	}

	fFEA_Shapes_displ.Construct	( fNumIP_displ,n_sd,n_en_displ );
	fFEA_Shapes_plast.Construct	( fNumIP_plast,n_sd,n_en_plast );
	
	Render_Vector.Dimension ( n_el );
	for (int e=0; e<n_el; e++) {
		Render_Vector[e].Construct ( 1, n_ip_plast, knumstrain+knumstress+knum_d_state );	
	}

	/* streams */
	ofstreamT& out = ElementSupport().Output();

	/* storage for integration point strain, stress, and ISVs*/
	fIPVariable.Dimension (n_el, fNumIP_plast*(knumstrain+knumstress+knum_d_state));
	fIPVariable = 0.0;

	/* allocate storage for nodal forces */
	//fForces_at_Node.Dimension ( n_sd );
	
	/* extract natural boundary conditions */
	TakeNaturalBC(list);
}



/* information about subordinate parameter lists */
void APS_AssemblyT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);
	
	/* tractions */
	sub_list.AddSub("aps_natural_bc", ParameterListT::Any);
}

/* return the description of the given inline subordinate parameter list */
void APS_AssemblyT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	ElementBaseT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* APS_AssemblyT::NewSub(const StringT& name) const
{
	/* create non-const this */
	APS_AssemblyT* non_const_this = const_cast<APS_AssemblyT*>(this);

	if (name == "aps_natural_bc") /* traction bc */
	{
		ParameterContainerT* natural_bc = new ParameterContainerT(name);

		natural_bc->AddParameter(ParameterT::Word, "side_set_ID");
		natural_bc->AddParameter(ParameterT::Integer, "schedule");

		ParameterT coord_sys(ParameterT::Enumeration, "coordinate_system");
		coord_sys.AddEnumeration("global", Traction_CardT::kCartesian);
		coord_sys.AddEnumeration( "local", Traction_CardT::kLocal);
		coord_sys.SetDefault(Traction_CardT::kCartesian);
		natural_bc->AddParameter(coord_sys);

		natural_bc->AddSub("DoubleList", ParameterListT::OnePlus); 		
		
		return natural_bc;
	}
	else /* inherited */
		return ElementBaseT::NewSub(name);
}



//##################################################################################
//###### Traction B.C. Methods (Cut and Paste from ContinuumElementT) ##############
//##################################################################################

//---------------------------------------------------------------------

//---------------------------------------------------------------------

/* update traction BC data */
void APS_AssemblyT::SetTractionBC(void)
{
//NOTE: With the possibility of variable global node numbers and
//		and equations, we assume as little as possible here with
//      regard to the validity of the node/equation numbers, requiring
//      only that NodesX in the element cards has the correct global
//      node numbers.

	/* dimensions */
	int ndof = NumDOF();

	/* echo values */
	iArray2DT nd_tmp, eq_tmp;
	for (int i = 0; i < fTractionList.Length(); i++)
	{
		Traction_CardT& BC_card = fTractionList[i];
			
		/* traction element/facet */
		int elem, facet;
		BC_card.Destination(elem, facet);

		/* set global node numbers */
		const iArrayT& loc_nodes = BC_card.LocalNodeNumbers();
		int nnd = loc_nodes.Length();
		
		iArrayT& nodes = BC_card.Nodes();
		nodes.Dimension(nnd);
		nodes.Collect(loc_nodes, fElementCards[elem].NodesX());
		
		/* set global equation numbers */
		iArrayT& eqnos = BC_card.Eqnos();
		eqnos.Dimension(ndof*nnd);
		
		/* get from node manager */
		nd_tmp.Set(1, nnd, nodes.Pointer());
		eq_tmp.Set(1, ndof*nnd, eqnos.Pointer());
		fDispl->SetLocalEqnos(nd_tmp, eq_tmp);
	}

	/* set flag */
	fTractionBCSet = 1;
}



/* extract natural boundary condition information */
void APS_AssemblyT::TakeNaturalBC(const ParameterListT& list)
{
	const char caller[] = "APS_AssemblyT::TakeTractionBC";

	int num_natural_bc = list.NumLists("natural_bc");
	if (num_natural_bc > 0)
	{
		/* model manager */
		ModelManagerT& model = ElementSupport().ModelManager();
	
		/* temp space */
		ArrayT<StringT> block_ID(num_natural_bc);
	    ArrayT<iArray2DT> localsides(num_natural_bc);
	    iArrayT LTf(num_natural_bc);
	    ArrayT<Traction_CardT::CoordSystemT> coord_sys(num_natural_bc);
	    ArrayT<dArray2DT> values(num_natural_bc);

	    /* nodes on element facets */
	    iArrayT num_facet_nodes;
	    fShapes_displ->NumNodesOnFacets(num_facet_nodes);
	    
	    /* loop over natural BC's */
	    int tot_num_sides = 0;
	    for (int i = 0; i < num_natural_bc; i++) 
	   	{
	    	const ParameterListT& natural_bc = list.GetList("natural_bc", i);
	    
	    	/* side set */
	    	const StringT& ss_ID = natural_bc.GetParameter("side_set_ID");
			localsides[i] = model.SideSet(ss_ID);
			int num_sides = localsides[i].MajorDim();
			tot_num_sides += num_sides;
			if (num_sides > 0)
			{
				block_ID[i] = model.SideSetGroupID(ss_ID);
				LTf[i] = natural_bc.GetParameter("schedule");
				coord_sys[i] = Traction_CardT::int2CoordSystemT(natural_bc.GetParameter("coordinate_system"));

				/* switch to elements numbering within the group */
				iArray2DT& side_set = localsides[i];
				iArrayT elems(num_sides);
				side_set.ColumnCopy(0, elems);
				BlockToGroupElementNumbers(elems, block_ID[i]);
				side_set.SetColumn(0, elems);

				/* all facets in set must have the same number of nodes */
				int num_nodes = num_facet_nodes[side_set(0,1)];
				for (int f = 0; f < num_sides; f++)
					if (num_facet_nodes[side_set(f,1)] != num_nodes)
						ExceptionT::BadInputValue(caller, "faces side set \"%s\" have different numbers of nodes",
							ss_ID.Pointer());

				/* read traction nodal values */
				dArray2DT& nodal_values = values[i];
				nodal_values.Dimension(num_nodes, NumDOF());
				int num_traction_vectors = natural_bc.NumLists("DoubleList");
				if (num_traction_vectors != 1 && num_traction_vectors != num_nodes)
					ExceptionT::GeneralFail(caller, "expecting 1 or %d vectors not %d",
						num_nodes, num_traction_vectors);
						
				/* constant over the face */
				if (num_traction_vectors == 1) {
					const ParameterListT& traction_vector = natural_bc.GetList("DoubleList");
					int dim = traction_vector.NumLists("Double");
					if (dim != NumDOF())
						ExceptionT::GeneralFail(caller, "expecting traction vector length %d not %d",
							NumDOF(), dim);

					/* same for all face nodes */
					for (int f = 0; f < NumDOF(); f++) {
						double t = traction_vector.GetList("Double", f).GetParameter("value");
						nodal_values.SetColumn(f, t);
					}
				}
				else
				{
					/* read separate vector for each face node */
					dArrayT t;
					for (int f = 0; f < num_nodes; f++) {
						const ParameterListT& traction_vector = natural_bc.GetList("DoubleList", f);
					int dim = traction_vector.NumLists("Double");
						if (dim != NumDOF())
							ExceptionT::GeneralFail(caller, "expecting traction vector length %d not %d",
								NumDOF(), dim);

						nodal_values.RowAlias(f, t);
						for (int j = 0; j < NumDOF(); j++)
							t[j] = traction_vector.GetList("Double", j).GetParameter("value");
					}
				}
			}
	    }
#pragma message("OK with empty side sets?")

		/* allocate all traction BC cards */
	    fTractionList.Dimension(tot_num_sides);

	    /* correct numbering offset */
	    LTf--;

		/* define traction cards */
		if (tot_num_sides > 0)
		{
			iArrayT loc_node_nums;
			int dex = 0;
			for (int i = 0; i < num_natural_bc; i++)
			{
				/* set traction BC cards */
				iArray2DT& side_set = localsides[i];
				int num_sides = side_set.MajorDim();
				for (int j = 0; j < num_sides; j++)
				{					
					/* get facet local node numbers */
					fShapes_displ->NodesOnFacet(side_set(j, 1), loc_node_nums);
					
					/* set and echo */
					fTractionList[dex++].SetValues(ElementSupport(), side_set(j,0), side_set (j,1), LTf[i],
						 coord_sys[i], loc_node_nums, values[i]);
				}
			}
		}

		/* check coordinate system specifications */
		if (NumSD() != NumDOF())
			for (int i = 0; i < fTractionList.Length(); i++)
				if (fTractionList[i].CoordSystem() != Traction_CardT::kCartesian)
					ExceptionT::BadInputValue(caller, "coordinate system must be Cartesian if (nsd != ndof) for card %d", i+1);
	}
}


//---------------------------------------------------------------------

/* compute contribution to RHS from traction BC's */
void APS_AssemblyT::ApplyTractionBC(void)
{
	if (fTractionList.Length() > 0)
	{
		/* dimensions */
		int nsd = NumSD();
		int ndof = NumDOF();
	
		/* update equation numbers */
		if (!fTractionBCSet) SetTractionBC();
	
		/* force vector */
		dArrayT rhs;
		VariArrayT<double> rhs_man(25, rhs);
		
		/* local coordinates */
		LocalArrayT coords(LocalArrayT::kInitCoords);
		VariLocalArrayT coord_man(25, coords, nsd);
		ElementSupport().RegisterCoordinates(coords);
		
		/* nodal tractions */
		LocalArrayT tract(LocalArrayT::kUnspecified);
		VariLocalArrayT tract_man(25, tract, ndof);

		/* integration point tractions */
		dArray2DT ip_tract;
		nVariArray2DT<double> ip_tract_man(25, ip_tract, ndof);
		dArrayT tract_loc, tract_glb(ndof);
		dMatrixT Q(ndof);
		
		/* Jacobian of the surface mapping */
		dMatrixT jacobian(nsd, nsd-1);
		
		for (int i = 0; i < fTractionList.Length(); i++)
		{
			const Traction_CardT& BC_card = fTractionList[i];

			/* dimension */
			const iArrayT& nodes = BC_card.Nodes();
			int nnd = nodes.Length();
			rhs_man.SetLength(nnd*ndof, false);
			coord_man.SetNumberOfNodes(nnd);
			tract_man.SetNumberOfNodes(nnd);
			
			/* local coordinates */
			coords.SetLocal(nodes);

			/* nodal traction vectors: (ndof x nnd) */
			BC_card.CurrentValue(tract);
			
			/* BC destination */
			int elem, facet;
			BC_card.Destination(elem, facet);
			
			/* default thickness */
			double thick = 1.0;
			
			/* boundary shape functions */
			const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(facet);
			int nip = surf_shape.NumIP();
			
			/* all ip tractions: (nip x ndof) */
			ip_tract_man.SetMajorDimension(nip, false);
			surf_shape.Interpolate(tract, ip_tract);

			/* traction vector coordinate system */
			if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
			{
				/* integrate */			
				rhs = 0.0;
				const double* w = surf_shape.Weight();
				for (int j = 0; j < nip; j++)
				{
					/* coordinate mapping */
					surf_shape.DomainJacobian(coords, j, jacobian);
					double detj = surf_shape.SurfaceJacobian(jacobian);
	
					/* ip weight */
					double jwt = detj*w[j]*thick;
					
					/* ip traction */
					const double* tj = ip_tract(j);
					
					/* accumulate */
					for (int l = 0; l < ndof; l++)
					{
						/* nodal shape function */
						const double* Na = surf_shape.Shape(j);
					
						double* prhs = rhs.Pointer(l);
						double  fact = jwt*(*tj++);
						for (int k = 0; k < nnd; k++)
						{
							*prhs += fact*(*Na++);
							prhs += ndof;
						}
					}				
				}
			}
			else if (BC_card.CoordSystem() == Traction_CardT::kLocal)
			{
				/* integrate */			
				rhs = 0.0;
				const double* w = surf_shape.Weight();
				for (int j = 0; j < nip; j++)
				{
					/* coordinate mapping */
					surf_shape.DomainJacobian(coords, j, jacobian);
					double detj = surf_shape.SurfaceJacobian(jacobian, Q);
	
					/* ip weight */
					double jwt = detj*w[j]*thick;
					
					/* transform ip traction out of local frame */
					ip_tract.RowAlias(j, tract_loc);
					Q.Multx(tract_loc, tract_glb);

					/* ip traction */
					const double* tj = tract_glb.Pointer();
					
					/* accumulate */
					for (int l = 0; l < ndof; l++)
					{
						/* nodal shape function */
						const double* Na = surf_shape.Shape(j);
					
						double* prhs = rhs.Pointer(l);
						double  fact = jwt*(*tj++);
						for (int k = 0; k < nnd; k++)
						{
							*prhs += fact*(*Na++);
							prhs += ndof;
						}
					}				
				}
			}
			else
				throw ExceptionT::kGeneralFail;

			/* assemble into displacement equations */
			ElementSupport().AssembleRHS(fDispl->Group(), rhs, BC_card.Eqnos());
		}
	}
}

