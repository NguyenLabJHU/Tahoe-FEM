/* $Id: FSSolidFluidMixT.cpp,v 1.3 2006-09-23 17:23:58 regueiro Exp $ */
#include "FSSolidFluidMixT.h"

#include "OutputSetT.h"
#include "ParameterContainerT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* constructor */
FSSolidFluidMixT::FSSolidFluidMixT(const ElementSupportT& support):
	ElementBaseT(support), //pass the solid displacement field to the base class
	u(LocalArrayT::kDisp),
	u_n(LocalArrayT::kLastDisp),
	press(LocalArrayT::kDisp),
	press_n(LocalArrayT::kLastDisp),
	fInitCoords_displ(LocalArrayT::kInitCoords),
	fCurrCoords_displ(LocalArrayT::kCurrCoords),
	fInitCoords_press(LocalArrayT::kInitCoords),
	fCurrCoords_press(LocalArrayT::kCurrCoords),
	fTractionBCSet(0),
	fDispl(NULL),
	fPress(NULL),
	fShapes_displ(NULL),
	fShapes_press(NULL),
	fKdd(ElementMatrixT::kNonSymmetric),
	fKdd_face(ElementMatrixT::kNonSymmetric),
	fKdtheta(ElementMatrixT::kNonSymmetric),
	fKthetad(ElementMatrixT::kNonSymmetric),
	fKthetatheta(ElementMatrixT::kNonSymmetric),
	bStep_Complete(0)
{
	SetName("total_lagrangian_solid_fluid_mix");
}

/* destructor */
FSSolidFluidMixT::~FSSolidFluidMixT(void) 
{  
	delete fShapes_displ;
	delete fShapes_press;
}


void FSSolidFluidMixT::Echo_Input_Data(void) {

	cout << "#######################################################" << endl; 
	cout << "############### ECHO FSSolidFluidMix DATA #########################" << endl; 
	cout << "#######################################################" << endl; 

	//################## material parameters ##################
	cout << "iConstitutiveModelType " 				<< iConstitutiveModelType 	<< endl; 
	
	//-- Elasticity parameters for solid
	cout << "fMaterial_Params[kMu] "  				<< fMaterial_Params[kMu] 	 << endl;
	cout << "fMaterial_Params[kLambda] "  			<< fMaterial_Params[kLambda] << endl;
	
	//-- Viscoelasticity parameters for solid
	cout << "fMaterial_Params[kAlpha] "  			<< fMaterial_Params[kAlpha] << endl;
	
	//-- Elasticity parameters for fluid
	cout << "fMaterial_Params[kKf] "  				<< fMaterial_Params[kKf] 	<< endl;

	//-- Initial real (intrinsic) mass densities
	cout << "fMaterial_Params[kRho_sR0] " 			<< fMaterial_Params[kRho_sR0] << endl;
	cout << "fMaterial_Params[kRho_fR0] " 			<< fMaterial_Params[kRho_fR0] << endl;
	
	//-- Initial volume fractions
	cout << "fMaterial_Params[kPhi_s0] " 			<< fMaterial_Params[kPhi_s0] << endl;
	cout << "fMaterial_Params[kPhi_f0] " 			<< fMaterial_Params[kPhi_f0] << endl;
	
	//################## Newmark time integration parameters ##################
	cout << "fIntegration_Params[kBeta] " 		<< fIntegration_Params[kBeta] 	<< endl;
	cout << "fIntegration_Params[kGamma] " 		<< fIntegration_Params[kGamma] 	<< endl;
}


//---------------------------------------------------------------------

void FSSolidFluidMixT::RHSDriver(void)	
{
	int curr_group = ElementSupport().CurrentGroup();

	/* traction boundary conditions acting on displacement equations */
	if (curr_group == fDispl->Group()) 
		ApplyTractionBC();

	/* choose solution method */
	if (fDispl->Group() == fPress->Group())
	  RHSDriver_monolithic();
	 else
	  RHSDriver_staggered();
}
//---------------------------------------------------------------------

void FSSolidFluidMixT::Equations(AutoArrayT<const iArray2DT*>& eq_d,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_theta)
{
	/* doing monolithic solution */
	if (fDispl->Group() == fPress->Group())
	{
		int ndof_press = fPress->NumDOF();
		int ndof_displ = fDispl->NumDOF();
	
		/* loop over connectivity blocks */
		fEqnos_displ.Dimension(fEqnos.Length());
		fEqnos_press.Dimension(fEqnos.Length());
		for (int i = 0; i < fEqnos.Length(); i++)
		{
			/* connectivities */
			const iArray2DT& connects_displ = *(fConnectivities_displ[i]);
			const iArray2DT& connects_press = *(fConnectivities_press[i]);
			int nel = connects_displ.MajorDim();
		
			/* dimension */ 
			fEqnos[i].Dimension(nel, n_en_displ*ndof_displ + n_en_press*ndof_press);
			iArray2DT& displ_eq = fEqnos_displ[i];
			iArray2DT& press_eq = fEqnos_press[i];
			displ_eq.Dimension(nel, n_en_displ*ndof_displ);
			press_eq.Dimension(nel, n_en_press*ndof_press);
			
			/* get equation numbers */
			fDispl->SetLocalEqnos(connects_displ, displ_eq);
			fPress->SetLocalEqnos(connects_press, press_eq);
			
			/* write into one array */
			fEqnos[i].BlockColumnCopyAt(displ_eq, 0);
			fEqnos[i].BlockColumnCopyAt(press_eq, displ_eq.MinorDim());

			/* add to list of equation numbers */
			eq_d.Append(&fEqnos[i]);
		}
	
		/* reset pointers to element cards */
		SetElementCards(fBlockData, fConnectivities_displ, fEqnos_displ, fElementCards_displ);
		SetElementCards(fBlockData, fConnectivities_press, fEqnos_press, fElementCards_press);
	}
	else
	/* doing staggered */
	{
#pragma message("initialization for staggered solution needs to be corrected")
	
		/* ElementBaseT handles equation array for displacements */
		if (ElementSupport().CurrentGroup() == fDispl->Group())
			ElementBaseT::Equations(eq_d, eq_theta);

		/* pore pressure equation */
		if (ElementSupport().CurrentGroup() == fPress->Group())
		{
			/* collect local equation numbers */
			//fPress.SetLocalEqnos(fConnectivities_press, fEqnos_press);
		
			//eq_d.Append(&fEqnos_press);
		}
	}
	
	/* get the equation number for the nodes on the faces */
	for (int i = 0; i < fPorePressureFaceEqnos.Length(); i++)
	{
		iArray2DT& faces = fPorePressureFaces[i];
		iArray2DT& eqnos = fPorePressureFaceEqnos[i];
		eqnos.Dimension(faces.MajorDim(), faces.MajorDim()*fDispl->NumDOF());
	
		fDispl->SetLocalEqnos(faces, eqnos);
	}
}


//---------------------------------------------------------------------

void FSSolidFluidMixT::LHSDriver(GlobalT::SystemTypeT)
{
/** Everything done in RHSDriver for efficiency */
//cout << "############### In LHS Driver ############### \n";
}

//---------------------------------------------------------------------

void FSSolidFluidMixT::Select_Equations (const int &iBalLinChoice, const int &iMassBalChoice )
{
	/** Choices for Linear Momentum Balance Equation */

	switch ( iBalLinChoice )	{

		default :
			cout << "FSSolidFluidMixT::Select_Equations() .. currently only one linear momentum balance for mixture \n";
			break;
	}

	/** Choices for Mass Balance Equation */

	switch ( iMassBalChoice )	{

		default :
			cout << "FSSolidFluidMixT::Select_Equations() .. currently only one mass balance equation for mixture \n";
			break;
	}

}

//---------------------------------------------------------------------

/* return true if the element contributes to the solution of the
 * given group. */
bool FSSolidFluidMixT::InGroup(int group) const
{
	return group == fDispl->Group() ||
	       group == fPress->Group();
}

//---------------------------------------------------------------------

/* close current time increment */
void FSSolidFluidMixT::CloseStep(void)
{
	/* inherited */
	ElementBaseT::CloseStep();

	/* store more recently updated values */
	fdState = fdState_new;
	fiState = fiState_new;
}


void FSSolidFluidMixT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented
}


/* form of tangent matrix */
GlobalT::SystemTypeT FSSolidFluidMixT::TangentType(void) const
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
void FSSolidFluidMixT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	const char caller[] = "FSSolidFluidMixT::AddNodalForce";

	/* displ, press, or neither */
	bool is_displ = false;
	dArrayT* element_force = NULL;
	int num_force = 0;
	if (field.FieldName() == fDispl->FieldName()) {
		is_displ = true;
		element_force = &fFd_int;
		num_force = fDispl->NumDOF();
		}
	else if (field.FieldName() == fPress->FieldName()) {
		is_displ = false;
		element_force = &fFtheta_int;
		num_force = fPress->NumDOF();
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
		const iArrayT& nodes_displ = fElementCards_displ[e].NodesU();
		const iArrayT& nodes_press = fElementCards_press[e].NodesU();

		u.SetLocal(nodes_displ);
		u_n.SetLocal(nodes_displ);
		press.SetLocal(nodes_press);
		press_n.SetLocal(nodes_press);

		del_u.DiffOf (u, u_n);
		del_press.DiffOf (press, press_n);

	 	// calculate derivatives based on reference coordinates
		fInitCoords_displ.SetLocal(fElementCards_displ[e].NodesX());
		//fCurrCoords_displ.SetToCombination (1.0, fInitCoords_displ, 1.0, u); 
		fCurrCoords_displ=fInitCoords_displ;
		fShapes_displ->SetDerivatives(); 
		//
		fInitCoords_press.SetLocal(fElementCards_press[e].NodesX());
		fCurrCoords_press=fInitCoords_press;
		//fCurrCoords_press.SetToCombination (1.0, fInitCoords_press, 1.0, u); 
		fShapes_press->SetDerivatives(); 
		
		//update state variables
		fdstatenew_all.Alias(fNumIP_press, knum_d_state, fdState_new(CurrElementNumber()));
		fdstate_all.Alias(fNumIP_press, knum_d_state, fdState(CurrElementNumber()));

		const double* Det    = fShapes_displ->IPDets();
		const double* Weight = fShapes_displ->IPWeights();
		/* calculate displacement nodal force */
		if (is_displ)
		{
			/* residual for displacement field */
			//generate this vector fFd_int 
			fShapes_displ->TopIP();
			while (fShapes_displ->NextIP())
			{
				//nothing right now
				fFd_int=0.0;
			}
		}
		else /* pressure nodal force */
		{
			/* residual for pore pressure field */ 
			// generate this vector fFtheta_int
			fShapes_displ->TopIP();
			while (fShapes_displ->NextIP())
			{
				//nothing right now
				fFtheta_int=0.0;
			}
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

double FSSolidFluidMixT::InternalEnergy ( void )
{
//not implemented
return 0.0;
}

//---------------------------------------------------------------------

/* write restart data to the output stream */
void FSSolidFluidMixT::WriteRestart(ostream& out) const
{
	/* inherited */
	ElementBaseT::WriteRestart(out);

	/* write state variable data */
	out << fdState;
}

//---------------------------------------------------------------------

/* read restart data to the output stream */
void FSSolidFluidMixT::ReadRestart(istream& in)
{
	/* inherited */
	ElementBaseT::ReadRestart(in);

	/* write state variable data */
	in >> fdState;
}

//---------------------------------------------------------------------

void FSSolidFluidMixT::RegisterOutput(void)
{
	/* collect block ID's */
	ArrayT<StringT> block_ID(fBlockData.Length());
	for (int i = 0; i < block_ID.Length(); i++)
		block_ID[i] = fBlockData[i].ID();

	/* output per element - strain, stress, and ISVs at the integration points */
	ArrayT<StringT> e_labels(fNumIP_press*(knumstrain+knumstress+knum_d_state));

	/* over integration points */
	//enter what values you need at integration points
	// stress and strain
	const char* slabels3D[] = {"s11", "s22", "s33"};
	// state variables; kappa is fictitious right now
	const char* svlabels3D[] = {"kappa"};
	int count = 0;
	for (int j = 0; j < fNumIP_press; j++)
	{
		StringT ip_label;
		ip_label.Append("ip", j+1);
			
		/* over strain and stress components */
		for (int i = 0; i < knumstrain+knumstress; i++)
		{
			e_labels[count].Clear();
			e_labels[count].Append(ip_label, ".", slabels3D[i]);
			count++;
		}
		
		/* over state variables */
		for (int i = 0; i < knum_d_state; i++)
		{
			e_labels[count].Clear();
			e_labels[count].Append(ip_label, ".", svlabels3D[i]);
			count++;
		}
	}		

	/* output per node */
	int num_node_output = fDispl->NumDOF() + fPress->NumDOF() + knumstrain + knumstress + knum_d_state;
	ArrayT<StringT> n_labels(num_node_output);
	count = 0;

	/* labels from pressic gradient */
	const ArrayT<StringT>& press_labels = fPress->Labels();
	for (int i = 0; i < press_labels.Length(); i++)
		n_labels[count++] = press_labels[i];

	/* labels from displacement */
	const ArrayT<StringT>& displ_labels = fDispl->Labels();
	for (int i = 0; i < displ_labels.Length(); i++)
		n_labels[count++] = displ_labels[i];

	/* labels from strains and stresses at the nodes */
	for (int i = 0; i < knumstrain+knumstress; i++)
		n_labels[count++] = slabels3D[i];
		
	/* labels from state variables at the nodes */
	for (int i = 0; i < knum_d_state; i++)
		n_labels[count++] = svlabels3D[i];

	/* set output specifier */
	#pragma message("FSSolidFluidMixT::RegisterOutput: is this right? ")
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

void FSSolidFluidMixT::WriteOutput(void)
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
		out_variable_all.Alias(fNumIP_press, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
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
	int num_node_output = fDispl->NumDOF() + fPress->NumDOF() + knumstrain + knumstress + knum_d_state;
	dArray2DT n_values(nodes_used.Length(), num_node_output);

	/* collect nodal values */
	const dArray2DT& fPressure = (*fPress)[0];
	const dArray2DT& fU = (*fDispl)[0];
	for (int i = 0; i < nodes_used.Length(); i++)
	{
		int node = nodes_used[i];
		double* row = n_values(i);
		for (int j = 0; j < fPressure.MinorDim(); j++)
			*row++ = fPressure(node,j);

		for (int j = 0; j < fU.MinorDim(); j++)
			*row++ = fU(node,j);

		double* p_stress = extrap_values(i);
		for (int j = 0; j < (knumstrain+knumstress+knum_d_state); j++)
			*row++ = p_stress[j];
	}

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, fIPVariable);

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
void FSSolidFluidMixT::RHSDriver_staggered(void)
{
	const char caller[] = "FSSolidFluidMixT::RHSDriver_staggered";
#pragma message("staggered solution not implemented")
}

//---------------------------------------------------------------------
/* form group contribution to the stiffness matrix and RHS */
void FSSolidFluidMixT::RHSDriver_monolithic(void)
{
	const char caller[] = "FSSolidFluidMixT::RHSDriver_monolithic";
	if (fDispl->Group() != fPress->Group())
		ExceptionT::GeneralFail(caller, "displacement and pore pressure groups must be the same: %d != %d",
			fDispl->Group(), fPress->Group());

	int curr_group = ElementSupport().CurrentGroup();

	/* stress output work space */
	dArray2DT	out_variable_all, fdstatenew_all, fdstate_all;
	dArrayT		out_variable;

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
		const iArrayT& nodes_press = fElementCards_press[e].NodesU();

		u.SetLocal(nodes_displ);
		u_n.SetLocal(nodes_displ);
		press.SetLocal(nodes_press);
		press_n.SetLocal(nodes_press);

		del_u.DiffOf (u, u_n);
		del_press.DiffOf (press, press_n);
		
		// calculate derivatives based on reference coordinates
		fInitCoords_displ.SetLocal(fElementCards_displ[e].NodesX());
		fCurrCoords_displ=fInitCoords_displ;
		//fCurrCoords_displ.SetToCombination (1.0, fInitCoords_displ, 1.0, u); 
		fShapes_displ->SetDerivatives(); 
		//
		fInitCoords_press.SetLocal(fElementCards_press[e].NodesX());
		fCurrCoords_press=fInitCoords_press;
		//fCurrCoords_press.SetToCombination (1.0, fInitCoords_press, 1.0, u); 
		fShapes_press->SetDerivatives(); 
		
		//update state variables
		fdstatenew_all.Alias(fNumIP_press, knum_d_state, fdState_new(CurrElementNumber()));
		fdstate_all.Alias(fNumIP_press, knum_d_state, fdState(CurrElementNumber()));
			
		if (bStep_Complete) { 
			
			//-- Store/Register data in classic tahoe manner 
			out_variable_all.Alias(fNumIP_press, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
			for (l=0; l < fNumIP_press; l++) 
			{
				out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(l));
				//out_variable=??;
			} 
	
		}
		else { //-- Still Iterating
		
			/* residual and tangent for displacements */
			const double* Det    = fShapes_displ->IPDets();
			const double* Weight = fShapes_displ->IPWeights();
			fShapes_displ->TopIP();
			while (fShapes_displ->NextIP())
			{
				//nothing right now for fKdtheta, fKdd, fFd_int
				fKdd=0.0;
				fKdtheta=0.0;
				fFd_int=0.0;
			}
			
			fShapes_displ->TopIP();
			while (fShapes_displ->NextIP())
			{
				//nothing right now for fKthetatheta, fKthetad, fFtheta_int
				fKthetad=0.0;
				fKthetatheta=0.0;
				fFtheta_int=0.0;
			}

			/* equations numbers */
			const iArrayT& displ_eq = fElementCards_displ[e].Equations();
			const iArrayT& press_eq = fElementCards_press[e].Equations();

			/* assemble residuals */
			ElementSupport().AssembleRHS(curr_group, fFd_int, displ_eq);
			ElementSupport().AssembleRHS(curr_group, fFtheta_int, press_eq);

			/* assemble components of the tangent */
			ElementSupport().AssembleLHS(curr_group, fKdd, displ_eq);
			ElementSupport().AssembleLHS(curr_group, fKthetatheta, press_eq);
			ElementSupport().AssembleLHS(curr_group, fKdtheta, displ_eq, press_eq);
			ElementSupport().AssembleLHS(curr_group, fKthetad, press_eq, displ_eq);
		}
	}	
}



/* form global shape function derivatives */
void FSSolidFluidMixT::SetGlobalShape(void)
{
	/* fetch (initial) coordinates */
	SetLocalX(fLocInitCoords);
	
	/* compute shape function derivatives */
	fShapes_displ->SetDerivatives();
	fShapes_press->SetDerivatives();
}



/* describe the parameters needed by the interface */
void FSSolidFluidMixT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

	/* displacement field */
	//already done in ElementBaseT
	//list.AddParameter(ParameterT::Word, "displ_field_name");
	
	/* pore pressure field */
	list.AddParameter(ParameterT::Word, "pore_pressure_field_name");
		
	list.AddParameter(fGeometryCode_displ_int, "GeometryCode_displ");
	list.AddParameter(fNumIP_displ, "NumIP_displ");
	list.AddParameter(fGeometryCodeSurf_displ_int, "GeometryCodeSurf_displ");
	list.AddParameter(fNumIPSurf_displ, "NumIPSurf_displ");
	list.AddParameter(n_en_displ, "n_en_displ");
	list.AddParameter(n_en_press, "n_en_press");
	
	list.AddParameter(iConstitutiveModelType, "constitutive_mod_type");
	
	double shearMu, sLambda, sAlpha, Rho_sR0, Rho_fR0, Phi_s0, Phi_f0, bulkK, 
			newBeta, newGamma;
			
	// solid elasticity
	list.AddParameter(shearMu, "mu");
	list.AddParameter(sLambda, "lambda");
	
	// fluid elasticity
	list.AddParameter(bulkK, "Kf");
	
	// viscoelasticity
	list.AddParameter(sAlpha, "alpha");
	
	// initial real mass density
	list.AddParameter(Rho_sR0, "rho_sR0");
	list.AddParameter(Rho_fR0, "rho_fR0");
	
	// initial volume fractions
	list.AddParameter(Phi_s0, "phi_s0");
	list.AddParameter(Phi_f0, "phi_f0");
	
	// Newmark time integration parameters
	list.AddParameter(newBeta, "beta");
	list.AddParameter(newGamma, "gamma");

}


/* accept parameter list */
void FSSolidFluidMixT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FSSolidFluidMixT::TakeParameterList";
	
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

	/* get pore pressure field */
	const StringT& pore_pressure_field_name = list.GetParameter("pore_pressure_field_name");
	fPress = ElementSupport().Field(pore_pressure_field_name);
	if (!fPress)
		ExceptionT::GeneralFail(caller, "could not resolve \"%s\" pore_pressure_field", 
		pore_pressure_field_name.Pointer());

	fGeometryCode_displ_int = list.GetParameter("GeometryCode_displ");
	fGeometryCode_displ = GeometryT::int2CodeT(fGeometryCode_displ_int);
	fNumIP_displ = list.GetParameter("NumIP_displ");
	fGeometryCodeSurf_displ_int = list.GetParameter("GeometryCodeSurf_displ");
	fGeometryCodeSurf_displ = GeometryT::int2CodeT(fGeometryCodeSurf_displ_int);
	fNumIPSurf_displ = list.GetParameter("NumIPSurf_displ");
	n_en_displ = list.GetParameter("n_en_displ");
	n_en_press = list.GetParameter("n_en_press");
	
	fGeometryCode_press = fGeometryCode_displ; 
	fNumIP_press = fNumIP_displ;
	fGeometryCodeSurf_press = fGeometryCodeSurf_displ;
	fNumIPSurf_press = fNumIPSurf_displ;
	
	iConstitutiveModelType = list.GetParameter("constitutive_mod_type");
	
	fMaterial_Params.Dimension ( kNUM_FMATERIAL_TERMS );
	fIntegration_Params.Dimension ( kNUM_FINTEGRATE_TERMS );
	
	fMaterial_Params[kMu] = list.GetParameter("mu");
	fMaterial_Params[kLambda] = list.GetParameter("lambda");
	fMaterial_Params[kKf] = list.GetParameter("Kf");
	fMaterial_Params[kAlpha] = list.GetParameter("alpha");
	fMaterial_Params[kRho_sR0] = list.GetParameter("rho_sR0");
	fMaterial_Params[kRho_fR0] = list.GetParameter("rho_fR0");
	fMaterial_Params[kPhi_s0] = list.GetParameter("phi_s0");
	fMaterial_Params[kPhi_f0] = list.GetParameter("phi_f0");
	
	fIntegration_Params[kBeta] = list.GetParameter("beta");
	fIntegration_Params[kGamma] = list.GetParameter("gamma");
	
	Echo_Input_Data();
	
	knum_d_state = 1; // double's needed per ip, state variables
	knum_i_state = 0; // int's needed per ip, state variables
	
	knumstrain = 0; // number of strain outputs
	knumstress = 3; // number of stress outputs
	
	output = "out";
	
	/* dimensions (notation as per Hughes' Book) */
	int& n_ip_displ = fNumIP_displ;
	int& n_ip_press = fNumIP_press;
	n_sd = NumSD();
	//n_df = NumDOF(); 
	//n_df = 1+n_sd; 		
	int nen = NumElementNodes(); /* number of nodes/element in the mesh */
	//n_en_press = n_en_displ = nen;
	//n_np = ElementSupport().NumNodes();

	/* initialize connectivities */
	fConnectivities_displ.Alias(fConnectivities);
	fConnectivities_press.Alias(fConnectivities);

	/* pick element interpolations based on available number of element nodes
	 * and the specified number of integration points */
	// only implemented for 3D, quadratic hexs
	//if (n_sd == 2 && nen == 8 && fGeometryCode_displ == GeometryT::kQuadrilateral) 
	if (n_sd == 3 && n_en_press != n_en_displ && fGeometryCode_displ == GeometryT::kHexahedron) 
	{
		// don't expect reduced integration for both fields 
		// if (n_ip_displ == 4 && n_ip_press == 4)
		//	ExceptionT::GeneralFail(caller, "not expecting 4 ips for both fields");
		//else if (n_ip_displ == 4 || n_ip_press == 4) // create reduced connectivities
		//{ 
			// reduce the number of element nodes based on the number ip's
			int& nen_red = (n_ip_displ == 8) ? n_en_displ : n_en_press;
			nen_red = 8;
			ArrayT<const iArray2DT*>& connects_red = (n_ip_displ == 8) ? 
				fConnectivities_displ : 
				fConnectivities_press;
		
			//create reduced connectivities
			connects_red.Dimension(0);
			connects_red.Dimension(fConnectivities.Length());
			fConnectivities_reduced.Dimension(fConnectivities.Length());
			for (int i = 0; i < fConnectivities_reduced.Length(); i++) {
				iArray2DT& connects_red_store = fConnectivities_reduced[i];
				const iArray2DT& connects = *(fConnectivities[i]);
				connects_red_store.Dimension(connects.MajorDim(), nen_red);				
				connects_red[i] = &connects_red_store;
				
				//take 1st eight element nodes (columns)
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
	// press
	fInitCoords_press.Dimension(n_en_press, n_sd);
	ElementSupport().RegisterCoordinates(fInitCoords_press);	
	fCurrCoords_press.Dimension(n_en_press, n_sd);
	fShapes_press = new ShapeFunctionT(fGeometryCode_press, fNumIP_press, fCurrCoords_press);
	//fShapes_press = new ShapeFunctionT(fGeometryCode_press, fNumIP_press, fCurrCoords_displ);
	//fShapes_press->Initialize();

	/* set local arrays for displacement field */
	u.Dimension (n_en_displ, n_sd);
	u_n.Dimension (n_en_displ, n_sd);
	del_u.Dimension (n_en_displ, n_sd);
	n_en_displ_x_n_sd = n_en_displ*n_sd;
	del_u_vec.Dimension (n_en_displ_x_n_sd);
	//ElementSupport().RegisterCoordinates(fInitCoords_displ);
	fDispl->RegisterLocal(u);
	fDispl->RegisterLocal(u_n);

	/* set local arrays for pore pressure field */
	int dum=1;
	press.Dimension (n_en_press, dum);
	press_n.Dimension (n_en_press, dum);
	del_press.Dimension (n_en_press, dum);
	del_press_vec.Dimension (n_en_press);
	//ElementSupport().RegisterCoordinates(fInitCoords_press);
	fPress->RegisterLocal(press);
	fPress->RegisterLocal(press_n);
	
	/* allocate state variable storage */
	// state variables are calculated at IPs for press field
	int num_ip = fNumIP_press;
	fdState_new.Dimension(n_el, num_ip*knum_d_state);
	fdState.Dimension(n_el, num_ip*knum_d_state);
	fiState_new.Dimension(n_el, num_ip*knum_i_state);
	fiState.Dimension(n_el, num_ip*knum_i_state);
	
	/* initialize equations */
	fEqnos_displ.Alias(fEqnos_displ);
	fEqnos_press.Dimension(fConnectivities_press.Length());

	/* initialize state variables */
	fdState = 0;
	fdState_new = 0;
	fiState = 0;
	fiState_new = 0;

	/* initialize element cards */
	fElementCards_displ.Alias(fElementCards);
	fElementCards_press.Dimension(fElementCards.Length());
	
	/* set cards to data in array - NOT NEEDED IF YOU'RE NOT
	 * GOING TO USE THE ElementCardT ARRAY? */
	for (int i= 0; i < fElementCards.Length(); i++)
		fElementCards[i].Set(fiState.MinorDim(), fiState(i), fdState.MinorDim(), fdState(i));

	fKdd.Dimension 			( n_en_displ_x_n_sd, n_en_displ_x_n_sd );
	fKdtheta.Dimension 		( n_en_displ_x_n_sd, n_en_press );
	fKthetad.Dimension 		( n_en_press, n_en_displ_x_n_sd );
	fKthetatheta.Dimension 	( n_en_press, n_en_press );

	fFd_int.Dimension 		( n_en_displ_x_n_sd );
	fFd_ext.Dimension 		( n_en_displ_x_n_sd );
	fFtheta_int.Dimension 	( n_en_press );
	fFtheta_ext.Dimension 	( n_en_press );

	/* streams */
	ofstreamT& out = ElementSupport().Output();

	/* storage for integration point strain, stress, and ISVs*/
	fIPVariable.Dimension (n_el, fNumIP_press*(knumstrain+knumstress+knum_d_state));
	fIPVariable = 0.0;

	/* allocate storage for nodal forces */
	//fForces_at_Node.Dimension ( n_sd );
	
	/* extract natural boundary conditions */
	TakeNaturalBC(list);
}



/* information about subordinate parameter lists */
void FSSolidFluidMixT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

	/* element blocks */
	sub_list.AddSub("total_lagrangian_solid_fluid_mix_element_block");
	
	/* tractions */
	sub_list.AddSub("total_lagrangian_solid_fluid_mix_natural_bc", ParameterListT::Any);
}



/* return the description of the given inline subordinate parameter list */
void FSSolidFluidMixT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	ElementBaseT::DefineInlineSub(name, order, sub_lists);
}



/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSSolidFluidMixT::NewSub(const StringT& name) const
{
	/* create non-const this */
	FSSolidFluidMixT* non_const_this = const_cast<FSSolidFluidMixT*>(this);

	if (name == "total_lagrangian_solid_fluid_mix_natural_bc") /* traction bc */
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
	else if (name == "total_lagrangian_solid_fluid_mix_element_block")
	{
		ParameterContainerT* element_block = new ParameterContainerT(name);
		element_block->AddSub("block_ID_list");
		return element_block;
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
void FSSolidFluidMixT::SetTractionBC(void)
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
void FSSolidFluidMixT::TakeNaturalBC(const ParameterListT& list)
{
	const char caller[] = "FSSolidFluidMixT::TakeTractionBC";

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
void FSSolidFluidMixT::ApplyTractionBC(void)
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

