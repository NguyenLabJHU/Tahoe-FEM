/* Id: FSSolidFluidMixT.cpp,v 1.6 2006/10/10 19:55:23 regueiro Exp $ */
#include "FSSolidFluidMixT.h"

#include "OutputSetT.h"
#include "ParameterContainerT.h"
#include "CommunicatorT.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
FSSolidFluidMixT::FSSolidFluidMixT(const ElementSupportT& support):
    ElementBaseT(support), //pass the solid displacement field to the base class
    u(LocalArrayT::kDisp),
    u_dot(LocalArrayT::kVel),
    u_dot_n(LocalArrayT::kLastVel),
    u_dotdot(LocalArrayT::kAcc),
    u_dotdot_n(LocalArrayT::kLastAcc),
    u_n(LocalArrayT::kLastDisp),
    press(LocalArrayT::kDisp),
    press_dot(LocalArrayT::kVel),
    press_dot_n(LocalArrayT::kLastVel),
    press_dotdot(LocalArrayT::kAcc),
    press_dotdot_n(LocalArrayT::kLastAcc),
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


    //-- Type of analysis
    cout << "kAnalysisType"  				<< kAnalysisType	 << endl;
	
    //-- Type of initial condition
    cout << "kInitialConditionType"  				<< kInitialConditionType	 << endl;

    //-- Elasticity parameters for solid
    cout << "fMaterial_Params[kMu] "  				<< fMaterial_Params[kMu] 	 << endl;
    cout << "fMaterial_Params[kLambda] "  			<< fMaterial_Params[kLambda] << endl;
	
    //-- Viscoelasticity parameters for solid
    cout << "fMaterial_Params[kAlpha] "  			<< fMaterial_Params[kAlpha] << endl;
	
    //-- Elasticity parameters for fluid
    cout << "fMaterial_Params[kKf] "  				<< fMaterial_Params[kKf] 	<< endl;

    //-- Hydraulic Conductivity
    cout << "fMaterial_Params[kK] "  				<< fMaterial_Params[kK] 	<< endl;

    //-- gravity
    cout << "fMaterial_Params[kg] "  				<< fMaterial_Params[kg] 	<< endl;

    //-- gravity in each direction (depends on the coordinate system which we have chosen for the problem)
    cout << "fMaterial_Params[kg1] "  				<< fMaterial_Params[kg1] 	<< endl;
    cout << "fMaterial_Params[kg2] "  				<< fMaterial_Params[kg2] 	<< endl;
    cout << "fMaterial_Params[kg3] "  				<< fMaterial_Params[kg3] 	<< endl;

    //-- Initial real (intrinsic) mass densities
    cout << "fMaterial_Params[kRho_sR0] " 			<< fMaterial_Params[kRho_sR0] << endl;
    cout << "fMaterial_Params[kRho_fR0] " 			<< fMaterial_Params[kRho_fR0] << endl;
	
    //-- Initial volume fractions
    cout << "fMaterial_Params[kPhi_s0] " 			<< fMaterial_Params[kPhi_s0] << endl;
    cout << "fMaterial_Params[kPhi_f0] " 			<< fMaterial_Params[kPhi_f0] << endl;
	
    //################## Newmark time integration parameters ##################
//    cout << "fIntegration_Params[kBeta] " 		<< fIntegration_Params[kBeta] 	<< endl;
//    cout << "fIntegration_Params[kGamma] " 		<< fIntegration_Params[kGamma] 	<< endl;
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
    /*
    for (int i = 0; i < fPorePressureFaceEqnos.Length(); i++)
    {
		iArray2DT& faces = fPorePressureFaces[i];
		iArray2DT& eqnos = fPorePressureFaceEqnos[i];
		eqnos.Dimension(faces.MajorDim(), faces.MajorDim()*fDispl->NumDOF());
		
		fDispl->SetLocalEqnos(faces, eqnos);
    }
    */
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
    return group == fDispl->Group() || group == fPress->Group();
}

//---------------------------------------------------------------------


/* initialize/finalize step */
void FSSolidFluidMixT::InitStep(void)
{
	/* inherited */
	ElementBaseT::InitStep();
}


/* close current time increment */
void FSSolidFluidMixT::CloseStep(void)
{
    /* inherited */
    ElementBaseT::CloseStep();
    

    //-- Store/Register initial values in classic tahoe manner 
    if ( ElementSupport().Time()==0 )
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

    /* zero first derivative of fields which are created at time=0 during calculating geostatic equilibrium(Terapezoidal rule) */    
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==1)
    {
	FieldT* fpress = const_cast <FieldT*> (fPress);
	FieldT* fdispl = const_cast <FieldT*> (fDispl);

	(*fdispl)[1] = 0;
	(*fpress)[1] = 0;
    }   

    /* zero second derivative of fields which are created at time=0 during calculating geostatic equilibrium(Newmark method) */
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==2)
    {
	FieldT* fpress = const_cast <FieldT*> (fPress);
	FieldT* fdispl = const_cast <FieldT*> (fDispl);
	(*fdispl)[2] = 0;
	(*fpress)[2] = 0;
    }

    /* reassign initial 2nd time derivative of pressure to 1st derivative */
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==2)
    {
	FieldT* fpress = const_cast <FieldT*> (fPress);
	(*fpress)[1] = (*fpress)[2];
	(*fpress)[2] = 0;
    }
   

    /* store more recently updated values */
    fdState = fdState_new;
    fiState = fiState_new;
	
//    fs_mix_out	<< endl 
//		<< setw(outputFileWidth) << "time_step"
//		<< endl;
    step_number = ElementSupport().StepNumber();
//    fs_mix_out	<< setw(outputFileWidth) << step_number
//		<< endl;
//    fs_mix_out	<< endl << "**********************************************************************************************";
//    fs_mix_out	<< endl << "**********************************************************************************************" << endl;
}


/* resets to the last converged solution */
/*
GlobalT::RelaxCodeT FSSolidFluidMixT::ResetStep(void)
{
	const char caller[] = "FSSolidFluidMixT::ResetStep";
	
	// inherited 
	GlobalT::RelaxCodeT relax = ElementBaseT::ResetStep();

	// update material internal variables 
	//needs to be implemented
#pragma message("reseting internal variables not implemented")	
	//ExceptionT::GeneralFail(caller, "reseting internal variables not implemented");

	return relax;
}
*/

/* element level reconfiguration for the current time increment */
/*
GlobalT::RelaxCodeT FSSolidFluidMixT::RelaxSystem(void)
{
	const char caller[] = "FSSolidFluidMixT::RelaxSystem";
	
	// inherited 
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	// loop over materials 
	//needs to be implemented
#pragma message("relax step for materials not implemented")	
	//ExceptionT::GeneralFail(caller, "relax step for materials not implemented");

	return relax;
}
*/


void FSSolidFluidMixT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented
}


/* return geometry and number of nodes on each facet */
void FSSolidFluidMixT::FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry, 
	iArrayT& num_facet_nodes) const
{
	/* from integration domain */
	ShapeFunctionDispl().FacetGeometry(facet_geometry, num_facet_nodes);
}


/* form of tangent matrix */
GlobalT::SystemTypeT FSSolidFluidMixT::TangentType(void) const
{
    return GlobalT::kNonSymmetric; 
}

/*
void FSSolidFluidMixT::SetStatus(const ArrayT<ElementCardT::StatusT>& status)
{
	// loop over elements and initial state variables
	int elem_num = 0;
	Top();
	while (NextElement())
	{
		// current element
		ElementCardT::StatusT& flag = CurrentElement().Flag();
		flag = status[elem_num++];

		if (flag == ElementCardT::kMarkON)
			flag = ElementCardT::kON;
		else if (flag == ElementCardT::kMarkOFF)
			flag = ElementCardT::kOFF;
	}
}
*/

/* initial condition/restart functions (per time sequence) */
void FSSolidFluidMixT::InitialCondition(void)
{
	/* inherited */
	ElementBaseT::InitialCondition();
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
    if (field.FieldName() == fDispl->FieldName()) 
    {
		is_displ = true;
		element_force = &fFd_int;
		num_force = fDispl->NumDOF();
    }
    else if (field.FieldName() == fPress->FieldName()) 
    {
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
		    fShapes_displ->SetDerivatives_DN_DDN(); 

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
		    if (is_displ)
		    {
			    for (int i = 0; i < nodes_displ.Length(); i++)
			    {
					if (nodes_displ[i] == node)
					{
					    /* components for node */
					    nodalforce.Set(num_force, element_force->Pointer(dex));

					    /* accumulate */
					    force += nodalforce;
					}
					dex += fDispl->NumDOF();
			    }		
		    }
		    else /* pressure nodal dof */
			{
			    for (int i = 0; i < nodes_press.Length(); i++)
			    {
					if (nodes_press[i] == node)
					{
					    /* components for node */
					    nodalforce.Set(num_force, element_force->Pointer(dex));

					    /* accumulate */
					    force += nodalforce;
					}
					dex += fPress->NumDOF();
			    }		
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
    // enter what values you need at integration points
    // stress and strain
    const char* slabels3D[] = {"s11", "s22", "s33","s23","s13","s12","p","e11","e22","e33","e23","e13","e12"};
    // state variables; phi_s and phi_f
    const char* svlabels3D[] = {"phi_s","phi_f","J"};
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

    /* print time */
//     fs_mix_out	<<"delta_t "<<delta_t << endl ;
//     fs_mix_out	<<"time "<<time << endl ;



    /* loop over elements */
    int e,l;
    Top();

    /* {fGravity_vector} will be formed */
    fGravity_vector[0]= fMaterial_Params[kg1];
    fGravity_vector[1]= fMaterial_Params[kg2];
    fGravity_vector[2]= fMaterial_Params[kg3];


    /* [fGravity_column_matrix] will be formed */
    for (int i=0; i<n_sd; i++)
	fGravity_column_matrix(i,0)=fGravity_vector[i];


//   fs_mix_out	<<"kInitialConditionType "<<kInitialConditionType << endl ;
//   fs_mix_out	<<"kAnalysisType "<<kAnalysisType << endl ;

    /* at time=0 when geostatic initial condition is calculated, trapezoidal integrator will calculate first time derivative of fields which by setting alpha_delta_t = 1 will be changed to 
       displacement and pressure which should be assigned to them, note that at time=0, delta_t=0 and Trapezoidal scheme which is embeded in the integrator will do nothing by itself(in changing previous iteration values)*/
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==1)
    {
	FieldT* fpress = const_cast <FieldT*> (fPress);
	FieldT* fdispl = const_cast <FieldT*> (fDispl);
	(*fpress)[0] = (*fpress)[1];
	(*fdispl)[0] = (*fdispl)[1];

    }  

    /* at time=0 when geostatic initial condition is calculated, dynamic Newmark integrator will calculate second time derivative of fields which by setting beta_delta_t2 = 1 will be changed to 
       displacement and pressure which should be assigned to them, note that at time=0, delta_t=0 and Newmark scheme which is embeded in dynamic integrator will do nothing by itself(in changing previous iteration value)*/
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==2)
    {
	FieldT* fpress = const_cast <FieldT*> (fPress);
	FieldT* fdispl = const_cast <FieldT*> (fDispl);
	(*fpress)[0] = (*fpress)[2];
	(*fdispl)[0] = (*fdispl)[2];

    } 
    while (NextElement())
    {
	fFd_int_N1_vector = 0.0;
	fFd_int_N2_vector = 0.0;
	fFtheta_int_N1_vector = 0.0;
	fFtheta_int_N2_vector = 0.0;
	fK_dd_G3_1_matrix = 0.0;
	fK_dd_G3_2_matrix = 0.0;
	fK_dd_G3_3_matrix = 0.0;
	fK_dd_G3_4_matrix = 0.0;
	fK_dd_G3_5_matrix = 0.0;
	fK_dtheta_G3_matrix = 0.0;
	fK_thetad_H3_1_matrix =0.0;
	fK_thetad_H3_2_matrix = 0.0;
	fK_thetad_H3_3_matrix = 0.0;
	fK_thetad_H3_4_matrix = 0.0;
	fK_thetatheta_H3_1_matrix = 0.0;
	fK_thetatheta_H3_2_matrix = 0.0;
	fK_dd_BTDB_matrix = 0.0;
	fFd_int_smallstrain_vector = 0.0;
	fM_dd_matrix = 0.0;
	fC_dd_matrix = 0.0;
	fFd_int_G4_vector = 0.0;
	fFd_int_M_vector = 0.0;
	fFd_int_C_vector = 0.0;
	fM_thetad_matrix = 0.0;
	fC_thetatheta_matrix = 0.0;
	fC_thetad_matrix = 0.0;
	fFtheta_int_H4_vector = 0.0;
	fFtheta_int_M_vector = 0.0;
	fFtheta_int_C1_vector = 0.0;
	fFtheta_int_C2_vector = 0.0;
	fK_dd_G1_1_matrix = 0.0;
	fK_dd_G1_2_matrix = 0.0;
	fK_dtheta_G1_matrix = 0.0;
	fK_dd_G2_1_matrix = 0.0;
	fK_dd_G2_2_matrix = 0.0;
	fK_dd_G2_3_matrix = 0.0;
	fK_dd_G2_4_matrix = 0.0;
	fK_dd_G4_matrix = 0.0;	
	fK_dtheta_G4_matrix = 0.0;
	fK_thetad_H1_1_matrix = 0.0; 
	fK_thetad_H1_2_matrix = 0.0;
	fK_thetad_H1_3_matrix = 0.0; 
	fK_thetad_H1_4_matrix = 0.0;
	fK_thetatheta_H1_matrix = 0.0;
	fK_thetad_H2_1_matrix = 0.0;
	fK_thetad_H2_2_matrix = 0.0;
	fK_thetad_H2_3_matrix = 0.0;
	fK_thetad_H2_4_matrix = 0.0;
	fK_thetad_H2_5_matrix = 0.0;
	fK_thetatheta_H2_1_matrix = 0.0;
	fK_thetatheta_H2_2_matrix = 0.0;
	fK_thetatheta_H2_3_matrix = 0.0;
	fK_thetad_H4_1_matrix = 0.0;
	fK_thetad_H4_2_matrix = 0.0;
	fK_thetad_H4_3_matrix = 0.0;
	fK_thetatheta_H4_matrix = 0.0;
	
		
	e = CurrElementNumber();
	const iArrayT& nodes_displ = fElementCards_displ[e].NodesU();
	const iArrayT& nodes_press = fElementCards_press[e].NodesU();

	u.SetLocal(nodes_displ);
	u_n.SetLocal(nodes_displ);
	if (u_dot.IsRegistered()) u_dot.SetLocal(nodes_displ);
	if (u_dot_n.IsRegistered()) u_dot_n.SetLocal(nodes_displ);
	if (u_dotdot.IsRegistered()) u_dotdot.SetLocal(nodes_displ);
	if (u_dotdot_n.IsRegistered())u_dotdot_n.SetLocal(nodes_displ);

	press.SetLocal(nodes_press);
	press_n.SetLocal(nodes_press);
	if (press_dot.IsRegistered()) press_dot.SetLocal(nodes_press);
	if (press_dot_n.IsRegistered()) press_dot_n.SetLocal(nodes_press);
	if (press_dotdot.IsRegistered()) press_dotdot.SetLocal(nodes_press);
	if (press_dotdot_n.IsRegistered()) press_dotdot_n.SetLocal(nodes_press);

	/* print solid displacement from previous step (u)*/
/*	fs_mix_out	<<"nodal solid displacement from previous step(u)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    fs_mix_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		fs_mix_out << u(i,j) << "\t";
	    fs_mix_out	<< endl ;
	}
*/


	/* print solid displacement from previous step (u_n)*/
/*	fs_mix_out	<<"nodal solid displacement from previous step(u_n)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    fs_mix_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		fs_mix_out << u_n(i,j) << "\t";
	    fs_mix_out	<< endl ;
	}
*/

	
	/* print solid velocity from previous step (u_dot)*/
/*	fs_mix_out	<<"nodal solid velocity from previous step(u_dot)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    fs_mix_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		fs_mix_out << u_dot(i,j) << "\t";
	    fs_mix_out	<< endl ;
	}
*/


	/* print solid velocity from previous step (u_dot_n)*/
/*	fs_mix_out	<<"nodal solid velocity from previous step(u_dot_n)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    fs_mix_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		fs_mix_out << u_dot_n(i,j) << "\t";
	    fs_mix_out	<< endl ;
	}
*/



	/* print solid acceleration from previous step (u_dotdot)*/
/*	fs_mix_out	<<"nodal solid velocity from previous step(u_dotdot)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    fs_mix_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		fs_mix_out << u_dotdot(i,j) << "\t";
	    fs_mix_out	<< endl ;
	}
*/


	/* print solid acceleration from previous step (u_dotdot_n)*/
/*	fs_mix_out	<<"nodal solid velocity from previous step(u_dotdot_n)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    fs_mix_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		fs_mix_out << u_dotdot_n(i,j) << "\t";
	    fs_mix_out	<< endl ;
	}
*/



		 
	/* print fluid pressure from previous step (press)*/
/*	fs_mix_out	<<"nodal fluid pressure from previous step(press)"<< endl ;
	for (int i=0; i<n_en_press; i++)
	{
	    fs_mix_out	<< "node number " << i+1 <<" :  " ;	
	    fs_mix_out	<< press(i,0) << endl;
	}
*/

	/* print fluid pressure from previous step (press_n)*/
/*	fs_mix_out	<<"nodal fluid pressure from previous step(press_n)"<< endl ;
	for (int i=0; i<n_en_press; i++)
	{
	    fs_mix_out	<< "node number " << i+1 <<" :  " ;	
	    fs_mix_out	<< press_n(i,0) << endl;
	}
*/

	
	/* print first derivative of pressure from previous step (press_dot)*/
/*	fs_mix_out	<<"first derivative of nodal fluid pressure from previous step(press_dot)"<< endl ;
	for (int i=0; i<n_en_press; i++)
	{
	    fs_mix_out	<< "node number " << i+1 <<" :  " ;	
	    fs_mix_out	<< press_dot(i,0) << endl;
	}
*/


	/* print first derivative of pressure from previous step (press_dot_n)*/
/*	fs_mix_out	<<"first derivative of nodal fluid pressure from previous step(press_dot_n)"<< endl ;
	for (int i=0; i<n_en_press; i++)
	{
	    fs_mix_out	<< "node number " << i+1 <<" :  " ;	
	    fs_mix_out	<< press_dot_n(i,0) << endl;
	}
*/



	/* print second derivative of pressure from previous step (press_dotdot)*/
/*	fs_mix_out	<<"second derivative of nodal fluid pressure from previous step(press_dotdot)"<< endl ;
	for (int i=0; i<n_en_press; i++)
	{
	    fs_mix_out	<< "node number " << i+1 <<" :  " ;	
	    fs_mix_out	<< press_dotdot(i,0) << endl;
	}
*/

	/* print second derivative of pressure from previous step (press_dotdot_n)*/
/*	fs_mix_out	<<"second derivative of nodal fluid pressure from previous step(press_dotdot_n)"<< endl ;
	for (int i=0; i<n_en_press; i++)
	{
	    fs_mix_out	<< "node number " << i+1 <<" :  " ;	
	    fs_mix_out	<< press_dotdot_n(i,0) << endl;
	}
*/

	
	
	/* populate solid displacement,solid velocity and 
	   solid accelration in vector form*/
	int index = 0;
	for (int i=0; i<n_en_displ; i++)
	{
	    for (int j=0; j<n_sd; j++)
	    {
		u_vec[index] = u(i,j);
		u_dot_vec[index] = u_dot(i,j);
		u_dotdot_vec[index] = u_dotdot(i,j);
		index += 1;
	    }
	}

        /* [u_dot_column_matrix] will be formed */
	for (int i=0; i<n_en_displ_x_n_sd; i++)
	    u_dot_column_matrix(i,0) = u_dot_vec[i];

        /* [u_dot_column_matrix_Transpose] will be formed */
	u_dot_column_matrix_Transpose.Transpose(u_dot_column_matrix);

	/* [u_dotdot_column_matrix] will be formed */
	for (int i=0; i<n_en_displ_x_n_sd; i++)
	    u_dotdot_column_matrix(i,0) = u_dotdot_vec[i];
	
	/* populate fluid pressure, first and second derivatives of fluid pressure 
	   in vector form*/
	for (int i=0; i<n_en_press; i++) 
	{
	    press_vec[i] = press(i,0);
	    press_dot_vec[i] = press_dot(i,0);
	    press_dotdot_vec[i] = press_dotdot(i,0);
	}

        /* [press_dot_column_matrix] will be formed */	
	for (int i=0; i<n_en_press; i++)
	    press_dot_column_matrix(i,0) = press_dot_vec[i];



	del_u.DiffOf (u, u_n);
	del_press.DiffOf (press, press_n);
	
	// calculate derivatives based on reference coordinates
	fInitCoords_displ.SetLocal(fElementCards_displ[e].NodesX());
	fCurrCoords_displ=fInitCoords_displ;
	//fCurrCoords_displ.SetToCombination (1.0, fInitCoords_displ, 1.0, u); 
	fShapes_displ->SetDerivatives_DN_DDN(); 
	//
	fInitCoords_press.SetLocal(fElementCards_press[e].NodesX());
	fCurrCoords_press=fInitCoords_press;
	//fCurrCoords_press.SetToCombination (1.0, fInitCoords_press, 1.0, u); 
	fShapes_press->SetDerivatives(); 
	
	//update state variables
	fdstatenew_all.Alias(fNumIP_press, knum_d_state, fdState_new(CurrElementNumber()));
	fdstate_all.Alias(fNumIP_press, knum_d_state, fdState(CurrElementNumber()));
	
	if (bStep_Complete) 
	{ 
	    //-- Store/Register data in classic tahoe manner 
	    out_variable_all.Alias(fNumIP_press, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
	    for (l=0; l < fNumIP_press; l++) 
	    {
		out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(l));
		Put_values_In_dArrayT_vector(fCauchy_effective_stress_Elements_IPs, e,l,fTemp_six_values);
		out_variable.CopyIn(0,fTemp_six_values);
		out_variable[6]=fPhysical_pore_water_pressure_Elements_IPs(e,l);
		Put_values_In_dArrayT_vector(fEulerian_effective_strain_Elements_IPs, e,l,fTemp_six_values);
		out_variable.CopyIn(7,fTemp_six_values);
		out_variable[13]=fState_variables_Elements_IPs(e,l*3+0);
		out_variable[14]=fState_variables_Elements_IPs(e,l*3+1);
		out_variable[15]=fState_variables_Elements_IPs(e,l*3+2);
	    } 
	}
	else 
	{ //-- Still Iterating

		if (time == 0 &&  kInitialConditionType==1)
		{
		    /* residual and tangent for displacements */
                    /* at time=0 program will solve for geostatic initial equilibrium by satisfying static equilibrium([M]=[C]=0)
		       it will create initial pore water pressure for fully saturated condition, displacement of solid phase and 
		       initial effective stresses due to geostatic equilibrium */
		    const double* Det    = fShapes_displ->IPDets();
		    const double* Weight = fShapes_displ->IPWeights();
		    fShapes_displ->TopIP();
		    fShapes_press->TopIP();
		    while (fShapes_displ->NextIP() && fShapes_press->NextIP())
		    {
			double scale_const = (*Weight++)*(*Det++);
			double beta_delta_t,beta_delta_t2,gamma_delta_t;
			beta_delta_t2 = 1;
			beta_delta_t = 1;
			gamma_delta_t = 1;

		
			const int IP = fShapes_displ->CurrIP();	
			dArrayT SolidIPCoordinate(n_sd),FluidIPCoordinate(n_sd);
			fShapes_displ->IPCoords(SolidIPCoordinate);
			fShapes_press->IPCoords(FluidIPCoordinate);


			const double* shapes_displ_X = fShapes_displ->IPShapeX();
			/* [fShapeSolid]will be formed */
			Form_solid_shape_functions(shapes_displ_X);
			
			
			fShapes_displ->GradNa(fShapeSolidGrad_temp);
			/* [fShapeSolidGrad] will be formed */
			Form_Gradient_of_solid_shape_functions(fShapeSolidGrad_temp);
			
			
			/* [fShapeSolidGrad_t] and [fShapeSolidGrad_t_Transpose] will be formed */
			Form_Gradient_t_of_solid_shape_functions(fShapeSolidGrad_temp);
			fShapeSolidGrad_t_Transpose.Transpose(fShapeSolidGrad_t);
		
		
			const double* shapes_press_X = fShapes_press->IPShapeX();
			/* {fShapeFluid} will be formed */
			Form_fluid_shape_functions(shapes_press_X);
			
			/* [fShapeFluid_row_matrix] will be formed */				
			for (int i=0; i<n_en_press ; i++)
			    fShapeFluid_row_matrix(0,i) = fShapeFluid[i];
			
			
			/* [fShapeFluidGrad] will be formed */
			fShapes_press->GradNa(fShapeFluidGrad);
			
			
			
			/* [fDeformation_Gradient] will be formed */
			Form_deformation_gradient_tensor();

                        /* [fDefGradT_9x9_matrix] will be formed */
			Form_fDefGradT_9x9_matrix();
			
			
			/* [fIdentity_matrix] will be formed */
			fIdentity_matrix = 0.0;			
			for (int i=0; i<n_sd ; i++)
			    fIdentity_matrix(i,i) =1.0;
		
		
			/* [fDeformation_Gradient_Inverse] and [fDeformation_Gradient_Transpose] and [fDeformation_Gradient_Inverse_Transpose] will be formed */
			if (fDeformation_Gradient.Det()==0)
			    fDeformation_Gradient = fIdentity_matrix; 
			fDeformation_Gradient_Inverse.Inverse(fDeformation_Gradient);
			fDeformation_Gradient_Inverse_Transpose.Transpose(fDeformation_Gradient_Inverse);
			fDeformation_Gradient_Transpose.Transpose(fDeformation_Gradient);
			
			
			/* {fDefGradInv_vector} will be formed */
			Form_deformation_gradient_inv_vector();
			
			/* [fDefGradInv_column_matrix] will be formed */
			for (int i=0; i<n_sd_x_n_sd; i++)
			    fDefGradInv_column_matrix(i,0)=fDefGradInv_vector[i];
			
			/* [fDefGradInv_column_matrix_Transpose] will be formed */		
			fDefGradInv_column_matrix_Transpose.Transpose(fDefGradInv_column_matrix);
			
			/* [fDefGradInv_Grad_grad] will be formed */
			Form_Grad_grad_transformation_matrix();
			
			
			/* [fDefGradInv_Grad_grad_Transpose] will be formed */
			fDefGradInv_Grad_grad_Transpose.Transpose(fDefGradInv_Grad_grad);
			
			/* Calculating theta */
			theta = fShapeFluid[0]*press_vec[0];
			for (int i=1; i<8; i++)
			    theta += fShapeFluid[i]*press_vec[i];
			
			
			/* Calculating Jacobian */
			double J = fDeformation_Gradient.Det();


			/* Jacobian for the current IP will be saved */
			fState_variables_IPs(IP,2)=J;

			
			/* Calculating fP_f */
			double fP_f=theta/J;

                        /* Physical pore water pressure for the current IP will be saved */
			fPhysical_pore_water_pressure_IPs(IP,0)=fP_f;

			
			/* Calculating fRho_f */
			fRho_f = fMaterial_Params[kRho_fR0]*exp((fP_f-fPf_0_matrix(CurrElementNumber(),IP))/
								fMaterial_Params[kKf]);
			
			
			
			/* Calculating phi_s and phi_f, volume fractions */ 
			phi_s = fMaterial_Params[kPhi_s0]/J;
			phi_f = 1.0 - phi_s;
			
			
			/*  Calculating fRho */
			fRho = phi_f*fRho_f+ phi_s*fMaterial_Params[kRho_sR0];
			
			/* Calculating fRho_0 */
			fRho_0 = J*fRho;
			
			
			/* [fRight_Cauchy_Green_tensor] will be formed */
			fRight_Cauchy_Green_tensor.MultATB(fDeformation_Gradient, fDeformation_Gradient);
		
		
			/* [fRight_Cauchy_Green_tensor_Inverse] will be formed */
			if (fRight_Cauchy_Green_tensor.Det()==0)
			    fRight_Cauchy_Green_tensor = fIdentity_matrix;
			fRight_Cauchy_Green_tensor_Inverse.Inverse(fRight_Cauchy_Green_tensor);
			
			
			/* [fLeft_Cauchy_Green_tensor] will be formed */
			fLeft_Cauchy_Green_tensor.MultABT(fDeformation_Gradient, fDeformation_Gradient);
			/* [fLeft_Cauchy_Green_tensor_Inverse] will be formed */
			if (fLeft_Cauchy_Green_tensor.Det()==0)
			    fLeft_Cauchy_Green_tensor = fIdentity_matrix;
			fLeft_Cauchy_Green_tensor_Inverse.Inverse(fLeft_Cauchy_Green_tensor);
			
			
			/* [fEulerian_effective_strain_tensor_current_IP] will be formed */
			fEulerian_effective_strain_tensor_current_IP = fLeft_Cauchy_Green_tensor_Inverse;
			fEulerian_effective_strain_tensor_current_IP *= -1; 
			fEulerian_effective_strain_tensor_current_IP += fIdentity_matrix;
			fEulerian_effective_strain_tensor_current_IP *= 0.5;
			
			
			/* extract six values of strain from symmetric eulerian strain tensor */
			Extract_six_values_from_symmetric_tensor(fEulerian_effective_strain_tensor_current_IP,fTemp_six_values);
			
			/* Save Euilerian effective strain tensor of the current IP */ 
			fEulerian_effective_strain_IPs.SetRow(IP,fTemp_six_values);
			
			
			/* Calculating J_Prim */
			if (fRight_Cauchy_Green_tensor.Det()==0)
			    fRight_Cauchy_Green_tensor = fIdentity_matrix; 
			double TempJ_Prim=fRight_Cauchy_Green_tensor.Det();
			double J_Prim=sqrt(fabs(TempJ_Prim));
			
			/* [fEffective_Second_Piola_tensor_inviscid] will be formed */
			fEffective_Second_Piola_tensor_inviscid.SetToScaled(fMaterial_Params[kLambda]*log(J_Prim)-fMaterial_Params[kMu],fRight_Cauchy_Green_tensor_Inverse); 
			fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fIdentity_matrix);
			fEffective_Second_Piola_tensor_inviscid += fTemp_matrix_nsd_x_nsd;
			
			
			/* [fEffective_Kirchhoff_tensor_inviscid] will be formed */
			fEffective_Kirchhoff_tensor_inviscid.MultABCT(fDeformation_Gradient,fEffective_Second_Piola_tensor_inviscid,fDeformation_Gradient);
			
			
			/* {fEffective_Kirchhoff_vector_inviscid} will be formed */
			Form_effective_kirchhoff_stress_vector_inviscid();
			
			
			/* [fIota_temp_matrix] will be formed */
			fIota_temp_matrix.MultATB(fShapeSolidGrad,fDefGradInv_Grad_grad);
			
			
			/* second derivatives of solid shape functions, [fShapeSolidGradGrad] will be formed */
			fShapes_displ->Grad_GradNa(fShapeSolidGradGrad);
			
			
			
			/* [fVarpi_temp_matrix] will be formed */
			Form_Varpi_temp_matrix();
			
			
			/* hydraulic conductivity matrix in the current coordinate, [k] will be formed */
			fK_hydraulic_conductivity_matrix.SetToScaled(fMaterial_Params[kK],fIdentity_matrix); 
			fk_hydraulic_conductivity_matrix.SetToScaled(1/J,fK_hydraulic_conductivity_matrix); 
			fTemp_matrix_nsd_x_nsd.MultABCT(fDeformation_Gradient,fk_hydraulic_conductivity_matrix,fDeformation_Gradient);
			fk_hydraulic_conductivity_matrix = fTemp_matrix_nsd_x_nsd;

//			fk_hydraulic_conductivity_matrix = fK_hydraulic_conductivity_matrix;
		
			
			/* [fLambda_temp_matrix] will be formed */
			fLambda_temp_matrix.MultATBC(fShapeFluidGrad,fDeformation_Gradient_Inverse,fk_hydraulic_conductivity_matrix);
			
			
			/* {fChi_temp_vector} will be formed */
			fVarpi_temp_matrix.Multx(u_vec,fChi_temp_vector);
			
			
			/* [fChi_temp_column_matrix] will be formed */
			for (int i=0; i<3 ; i++)
			    fChi_temp_column_matrix(i,0)= fChi_temp_vector[i];
			
			
			/* {fFd_int_N1_vector} will be formed */
			double scale = scale_const;
			fIota_temp_matrix.Multx(fEffective_Kirchhoff_vector_inviscid,fTemp_vector_ndof_se,scale);
			/* fFd_int_N1_vector for the current IP */
			/* accumulate */
			fFd_int_N1_vector += fTemp_vector_ndof_se;
			
			/* {fFd_int_N2_vector} will be formed */
			scale = -1.0*theta*scale_const;
			fShapeSolidGrad.MultTx(fDefGradInv_vector,fTemp_vector_ndof_se);
			fTemp_vector_ndof_se *= scale;
			/* accumulate */
			fFd_int_N2_vector += fTemp_vector_ndof_se; 
			
			
			/* state vaiables(volume fractions) for the current IP will be saved */
			fState_variables_IPs(IP,0)=phi_s;
			fState_variables_IPs(IP,1)=phi_f;
			
			
			/* {fFtheta_int_N1_vector} will be formed */
			double const1 = fMaterial_Params[kg]* phi_f;
			if (fabs(const1)>1.0e-16)
			    scale = theta*(1.0/fMaterial_Params[kg]) * ((phi_s/phi_f)-1) * scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
			fTemp_matrix_nen_press_x_nsd.Multx(fChi_temp_vector, fTemp_vector_nen_press,scale);
			/* accumulate */
			fFtheta_int_N1_vector += fTemp_vector_nen_press;
			
			
			/* {fFtheta_int_N2_vector} will be formed */
			fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_nsd,fShapeFluidGrad);
			scale = 1.0/fMaterial_Params[kg]*scale_const; 
			fTemp_matrix_nen_press_x_nen_press.Multx(press_vec, fTemp_vector_nen_press,scale);
			/* accumulate */
			fFtheta_int_N2_vector += fTemp_vector_nen_press;
			
			
			/* [fIm_temp_matrix] will be formed */
			Form_Im_temp_matrix();
			
			
			/* [fHbar_temp_matrix] will be formed */
			Form_Hbar_temp_matrix();
			
			
			/* [fEll_temp_matrix] will be formed */
			Form_Ell_temp_matrix();
			
			
			/* {fPi_temp_transpose_vector} will be formed */
			fShapeSolidGrad.MultTx(fDefGradInv_vector,fPi_temp_transpose_vector);
		
		
			/* [fPi_temp_row_matrix] will be formed */
			for (int i=0; i<n_en_displ_x_n_sd; i++)
			    fPi_temp_row_matrix(0,i) = fPi_temp_transpose_vector[i];
		
			/* [fK_dd_G3_1_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fIm_temp_matrix,fIota_temp_matrix);
			scale = -1*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_1_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fI_ij_column_matrix] will be formed */
			fI_ij_column_matrix = 0.0;
			fI_ij_column_matrix(0,0) = 1.0;
			fI_ij_column_matrix(4,0) = 1.0;
			fI_ij_column_matrix(8,0) = 1.0;
			
			
			/* [fK_dd_G3_2_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fHbar_temp_matrix,fIota_temp_matrix);
			scale = fMaterial_Params[kMu] * beta_delta_t2 * scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_2_matrix += fTemp_matrix_ndof_se_x_ndof_se;
		
		
			/* [fK_dd_G3_3_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fEll_temp_matrix,fIota_temp_matrix);
			scale = fMaterial_Params[kMu] * beta_delta_t2 * scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_3_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dd_G3_4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABC(fIota_temp_matrix,fI_ij_column_matrix,fPi_temp_row_matrix);
			scale = fMaterial_Params[kLambda] * beta_delta_t2 * scale_const; 
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_4_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dd_G3_5_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fShapeSolidGrad_t_Transpose,fDefGradInv_Grad_grad_Transpose,fIota_temp_matrix);
			scale = theta * beta_delta_t2 * scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_5_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dtheta_G3_matrix] will be formed */
			fTemp_matrix_ndof_se_x_nen_press.MultATB(fPi_temp_row_matrix,fShapeFluid_row_matrix); 
			scale = -1*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_nen_press *= scale;
			/* accumulate */
			fK_dtheta_G3_matrix += fTemp_matrix_ndof_se_x_nen_press;
			
			/* {fGrad_1_J_vector} will be filled */
			fVarpi_temp_matrix.Multx(u_vec,fGrad_1_J_vector, -1.0/J);
			
			
			/* {fGrad_theta_vector} will be filled */
			fShapeFluidGrad.Multx(press_vec, fGrad_theta_vector);
		
		
			/* {fGrad_phi_f_vector} will be filled */
			fGrad_phi_f_vector.SetToScaled(-1* fMaterial_Params[kPhi_s0],fGrad_1_J_vector);
		
			
			/* {fGrad_Omega_vector} will be filled */
			fTemp_nsd_vector.SetToScaled(theta/J,fGrad_phi_f_vector) ; 
			fGrad_Omega_vector = fTemp_nsd_vector;
			fTemp_nsd_vector.SetToScaled(phi_f/J,fGrad_theta_vector) ;
			fGrad_Omega_vector += fTemp_nsd_vector;
			fTemp_nsd_vector.SetToScaled(phi_f * theta,fGrad_1_J_vector) ;
			fGrad_Omega_vector += fTemp_nsd_vector;
			
			
			/* {fgrad_Omega_vector} will be formed */
			fDeformation_Gradient_Inverse_Transpose.Multx(fGrad_Omega_vector,fgrad_Omega_vector);
			
			
			/* [fJmath_temp_matrix] will be formed */
			Form_Jmath_temp_matrix(); 
			
			/* [fWp_temp_matrix] will be formed */
			Form_Wp_temp_matrix(); 
			
			
			/* [fK_thetad_H3_1_matrix] will be formed */
			const1 = fMaterial_Params[kg]*phi_f;
			if (fabs(const1) > 1e-16) 
			    scale = (-1*theta*phi_s/(fMaterial_Params[kg]*phi_f))*
				(phi_s/phi_f +1)*beta_delta_t2*scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
			fTemp_matrix_nen_press_x_ndof_se.MultABC(fTemp_matrix_nen_press_x_nsd,fChi_temp_column_matrix,fPi_temp_row_matrix);
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H3_1_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H3_2_matrix] will be formed */
			const1 = fMaterial_Params[kg]*phi_f ;
			if (fabs(const1) > 1e-16)
			    scale = 1*(theta/fMaterial_Params[kg])*
				(phi_s/phi_f -1)*beta_delta_t2*scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_ndof_se.MultAB(fTemp_matrix_nen_press_x_nsd,fVarpi_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H3_2_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H3_3_matrix] will be formed */
			const1 = fMaterial_Params[kg]*phi_f;
			if (fabs(const1) > 1e-16)
			    scale = -1*J/(fMaterial_Params[kg]*phi_f)*beta_delta_t2*scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nsd.MultATB(fShapeFluidGrad,fDeformation_Gradient_Inverse);
			fTemp_matrix_nen_press_x_ndof_se.MultABCT(fTemp_matrix_nen_press_x_nsd,fJmath_temp_matrix,fIota_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H3_3_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H3_4_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultABCT(fTemp_matrix_nen_press_x_nsd,fWp_temp_matrix,fIota_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H3_4_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetatheta_H3_1_matrix] will be formed */
			scale = (1/fMaterial_Params[kg])*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
			fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_nsd,fShapeFluidGrad);
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H3_1_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			/* [fK_thetatheta_H3_2_matrix] will be formed */
			const1 = fMaterial_Params[kg] * phi_f ;
			if (fabs(const1)> 1e-16)
			    scale =  (1/fMaterial_Params[kg])*
				(phi_s/phi_f -1)*beta_delta_t2*scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nen_press.MultABC(fTemp_matrix_nen_press_x_nsd,fChi_temp_column_matrix,fShapeFluid_row_matrix);
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H3_2_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			
			/* Creating Second tangential elasticity tensor in the Ref. coordinate [fC_matrix] */
			Form_C_matrix(J_Prim);

                        /* [fEffective_Second_Piola_tensor_viscous] will be formed */
			double fAlpha = fMaterial_Params[kAlpha];
			Form_fEffective_Second_Piola_tensor_viscous(fAlpha);

			/* [fEffective_Kirchhoff_tensor_viscous] will be formed */
			fEffective_Kirchhoff_tensor_viscous.MultABCT(fDeformation_Gradient,fEffective_Second_Piola_tensor_viscous,fDeformation_Gradient);

                        /* [fEffective_Kirchhoff_tensor] will be formed */
			fEffective_Kirchhoff_tensor = fEffective_Kirchhoff_tensor_inviscid;
			fEffective_Kirchhoff_tensor += fEffective_Kirchhoff_tensor_viscous;

			
			/* [fCauchy_effective_stress_tensor_current_IP] will be formed */
			fCauchy_effective_stress_tensor_current_IP = fEffective_Kirchhoff_tensor;
			fCauchy_effective_stress_tensor_current_IP *= 1/J;
			
			
			/* extract six values of stress from symmetric cauchy stress tensor */
			Extract_six_values_from_symmetric_tensor(fCauchy_effective_stress_tensor_current_IP,fTemp_six_values);
			
			
			
			/* Save Cauchy effective stress tensor of the current IP */ 
			fCauchy_effective_stress_IPs.SetRow(IP,fTemp_six_values); 
			

			/* Creating Second tangential elasticity tensor in the Current coordinate [fc_matrix]*/
			Form_c_matrix();
		
		
			/* [fIm_Prim_temp_matrix] will be formed */
			Form_Im_Prim_temp_matrix();
			
			/* [fUpsilon_temp_matrix] will be formed */
			fUpsilon_temp_matrix.MultATB(fShapeSolid,fShapeSolid);
			
			
			/* {fFd_int_G4_vector} will be formed */
			fShapeSolid.MultTx(fGravity_vector,fTemp_vector_ndof_se);
			scale = -1*fRho_0*scale_const;
			fTemp_vector_ndof_se *= scale; 
			/* accumulate */
			fFd_int_G4_vector += fTemp_vector_ndof_se;
			
			
			/* {fFtheta_int_H4_vector} will be formed */
			fLambda_temp_matrix.Multx(fGravity_vector,fTemp_vector_nen_press);
			scale = -1/(fMaterial_Params[kg])*J*fRho_f*scale_const;
			fTemp_vector_nen_press *= scale ;
			/* accumulate */
			fFtheta_int_H4_vector += fTemp_vector_nen_press;
			
			
			/* fC1, fC2 and fC3 will be formed */
			fC1 = phi_f*fRho_f/(fMaterial_Params[kKf]*J);
			fC2 = 1/J*(fRho_f*fMaterial_Params[kPhi_s0]-
				   phi_f*fRho_f*theta/fMaterial_Params[kKf]);
			fC3 = fC2 - fMaterial_Params[kRho_sR0]*fMaterial_Params[kPhi_s0]/J;
			
			
			/* {fgradv_vector} will be formed */
			Form_gradv_vector();
			
			
			/* [fXi_temp_matrix] will be formed */
			Form_Xi_temp_matrix();
			
			
			/* [fVarsigma_temp_matrix] will be formed */
			Form_Varsigma_temp_matrix();
			
			
			/* [fI_ijkl_matrix] will be formed */
			Form_I_ijkl_matrix();
			
			
			/* [fK_dd_G4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultATBC(fShapeSolid,fGravity_column_matrix,fPi_temp_row_matrix);
			scale = -1*J*(fC3 + fRho)*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G4_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dtheta_G4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_nen_press.MultATBC(fShapeSolid,fGravity_column_matrix,fShapeFluid_row_matrix);
			scale = -1*J*fC1*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_nen_press *= scale;
			/* accumulate */
			fK_dtheta_G4_matrix += fTemp_matrix_ndof_se_x_nen_press;
			

			/* [fAleph_temp_matrix] will be formed */
			Form_Aleph_temp_matrix(IP);
			

			/* [fImath_temp_matrix] will be formed */
			Form_Imath_temp_matrix();
			
			
			/* [fK_thetad_H4_1_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultABC(fLambda_temp_matrix,fGravity_column_matrix,fPi_temp_row_matrix);
			scale = 1/(fMaterial_Params[kg]*fMaterial_Params[kKf])*
			    fRho_f*theta*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H4_1_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H4_2_matrix] will be formed */
			fTemp_matrix_nsd_x_ndof_se.MultABCT(fDeformation_Gradient_Inverse,fImath_temp_matrix,fIota_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se.MultATB(fShapeFluidGrad,fTemp_matrix_nsd_x_ndof_se);
			scale = 1/(fMaterial_Params[kg])*J*fRho_f*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H4_2_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H4_3_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultABC(fLambda_temp_matrix,fGravity_column_matrix,fPi_temp_row_matrix);
			scale = -1/(fMaterial_Params[kg])*J*fRho_f*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H4_3_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetatheta_H4_matrix] will be formed */
			fTemp_matrix_nen_press_x_nen_press.MultABC(fLambda_temp_matrix,fGravity_column_matrix,fShapeFluid_row_matrix);
			scale = -1/(fMaterial_Params[kg]*fMaterial_Params[kKf])*
			    fRho_f*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H4_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			
		    }		
	    
		    /* saving eulerian effective strain for each IPs of the current element */
		    fEulerian_effective_strain_Elements_IPs.SetRow(e,fEulerian_effective_strain_IPs);
	    
	    
		    /* saving cauchy effective stress for each IPs of the current element */
		    fCauchy_effective_stress_Elements_IPs.SetRow(e,fCauchy_effective_stress_IPs);

                    /* saving physical pore water pressure for each IPs of the current element */
		    fPhysical_pore_water_pressure_Elements_IPs.SetRow(e,fPhysical_pore_water_pressure_IPs);
	    
	    
		    /* saving state variables for each IPs of the current element */
		    fState_variables_Elements_IPs.SetRow(e,fState_variables_IPs);

	    
		    /* {fFd_int} will be formed */
		    fFd_int = fFd_int_N1_vector;
		    fFd_int += fFd_int_N2_vector; 
		    fFd_int += fFd_int_G4_vector;
		    fFd_int *= -1;
	    

		    /* [fKdd] will be formed */
		    fKdd = fK_dd_G3_1_matrix;
		    fKdd += fK_dd_G3_2_matrix;
		    fKdd += fK_dd_G3_3_matrix;
		    fKdd += fK_dd_G3_4_matrix;
		    fKdd += fK_dd_G3_5_matrix; 
		    fKdd += fK_dd_G4_matrix; 



		    /* [fKdtheta] will be formed */
		    fKdtheta = fK_dtheta_G3_matrix;
		    fKdtheta += fK_dtheta_G4_matrix;
	    
	    
		    /* [fKthetad] will be formed */
		    fKthetad = fK_thetad_H3_1_matrix;
		    fKthetad += fK_thetad_H3_2_matrix;
		    fKthetad += fK_thetad_H3_3_matrix;
		    fKthetad += fK_thetad_H3_4_matrix;
		    fKthetad += fK_thetad_H4_1_matrix;
		    fKthetad += fK_thetad_H4_2_matrix;
		    fKthetad += fK_thetad_H4_3_matrix;


		    /* [fKthetatheta] will be formed */
		    fKthetatheta = fK_thetatheta_H3_1_matrix;
		    fKthetatheta += fK_thetatheta_H3_2_matrix;
		    fKthetatheta += fK_thetatheta_H4_matrix;	


		    /* {fFtheta_int} will be formed */
		    fFtheta_int = fFtheta_int_N1_vector;
		    fFtheta_int += fFtheta_int_N2_vector;
		    fFtheta_int += fFtheta_int_H4_vector;
		    fFtheta_int *= -1;
		    

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
		else if (time == 0 &&  kInitialConditionType==2)		
		{

                    /* at time=0 program will solve for second time derivative of solid phase and first time derivative of fluid phase based on initial
		       displacement and velocity for the solid phase and initial pore pressure for the fluid phase by satisfying dynamic equilibrium equations at t==0 
		       which contains [M] and [C]*/
		    /* residual and tangent for displacements */
		    const double* Det    = fShapes_displ->IPDets();
		    const double* Weight = fShapes_displ->IPWeights();
		    fShapes_displ->TopIP();
		    fShapes_press->TopIP();

		    while (fShapes_displ->NextIP() && fShapes_press->NextIP())
		    {
			const int IP = fShapes_displ->CurrIP();	
			/* Interpolating initial nodal values of pressures to integration points and saving them in a nel*num_ip matrix */
			fShapes_press->InterpolateU(press_n,fP0_temp_value);
			fPf_0_matrix(CurrElementNumber(),IP)=fP0_temp_value[0];

			double scale_const = (*Weight++)*(*Det++);
			

			dArrayT SolidIPCoordinate(n_sd),FluidIPCoordinate(n_sd);
			fShapes_displ->IPCoords(SolidIPCoordinate);
			fShapes_press->IPCoords(FluidIPCoordinate);
			
			
			const double* shapes_displ_X = fShapes_displ->IPShapeX();
			/* [fShapeSolid]will be formed */
			Form_solid_shape_functions(shapes_displ_X);
			
			
			fShapes_displ->GradNa(fShapeSolidGrad_temp);
			/* [fShapeSolidGrad] will be formed */
			Form_Gradient_of_solid_shape_functions(fShapeSolidGrad_temp);
			
			
			/* [fShapeSolidGrad_t] and [fShapeSolidGrad_t_Transpose] will be formed */
			Form_Gradient_t_of_solid_shape_functions(fShapeSolidGrad_temp);
			fShapeSolidGrad_t_Transpose.Transpose(fShapeSolidGrad_t);
			
			
			const double* shapes_press_X = fShapes_press->IPShapeX();
			/* {fShapeFluid} will be formed */
			Form_fluid_shape_functions(shapes_press_X);
			
			/* [fShapeFluid_row_matrix] will be formed */				
			for (int i=0; i<n_en_press ; i++)
			    fShapeFluid_row_matrix(0,i) = fShapeFluid[i];
			
			
			/* [fShapeFluidGrad] will be formed */
			fShapes_press->GradNa(fShapeFluidGrad);
			
			
			
			/* [fDeformation_Gradient] will be formed */
			Form_deformation_gradient_tensor();

                       /* [fDefGradT_9x9_matrix] will be formed */
			Form_fDefGradT_9x9_matrix();
			
			
			/* [fIdentity_matrix] will be formed */
			fIdentity_matrix = 0.0;			
			for (int i=0; i<n_sd ; i++)
			    fIdentity_matrix(i,i) =1.0;
			
			
			/* [fDeformation_Gradient_Inverse] and [fDeformation_Gradient_Transpose] and [fDeformation_Gradient_Inverse_Transpose] will be formed */
			if (fDeformation_Gradient.Det()==0)
			    fDeformation_Gradient = fIdentity_matrix; 
			fDeformation_Gradient_Inverse.Inverse(fDeformation_Gradient);
			fDeformation_Gradient_Inverse_Transpose.Transpose(fDeformation_Gradient_Inverse);
			fDeformation_Gradient_Transpose.Transpose(fDeformation_Gradient);
			
			
			/* {fDefGradInv_vector} will be formed */
			Form_deformation_gradient_inv_vector();
				
			/* [fDefGradInv_column_matrix] will be formed */
			for (int i=0; i<n_sd_x_n_sd; i++)
			    fDefGradInv_column_matrix(i,0)=fDefGradInv_vector[i];
			
			/* [fDefGradInv_column_matrix_Transpose] will be formed */		
			fDefGradInv_column_matrix_Transpose.Transpose(fDefGradInv_column_matrix);
			
			/* [fDefGradInv_Grad_grad] will be formed */
			Form_Grad_grad_transformation_matrix();
			
			
			/* [fDefGradInv_Grad_grad_Transpose] will be formed */
			fDefGradInv_Grad_grad_Transpose.Transpose(fDefGradInv_Grad_grad);
			
			/* Calculating theta */
			theta = fShapeFluid[0]*press_vec[0];
			for (int i=1; i<8; i++)
			    theta += fShapeFluid[i]*press_vec[i];
			
			
			/* Calculating Jacobian */
			double J = fDeformation_Gradient.Det();

			
			/* Jacobian for the current IP will be saved */
			fState_variables_IPs(IP,2)=J;

			
			/* Calculating fP_f */
			double fP_f=theta/J;


                        /* Physical pore water pressure for the current IP will be saved */
			fPhysical_pore_water_pressure_IPs(IP,0)=fP_f;
			
			/* Calculating fRho_f */
			fRho_f = fMaterial_Params[kRho_fR0]*exp((fP_f-fPf_0_matrix(CurrElementNumber(),IP))/
								fMaterial_Params[kKf]);
			
			
			
			/* Calculating phi_s and phi_f, volume fractions */ 
			phi_s = fMaterial_Params[kPhi_s0]/J;
			phi_f = 1.0 - phi_s;
			
			
			/*  Calculating fRho */
			fRho = phi_f*fRho_f+ phi_s*fMaterial_Params[kRho_sR0];
			
			/* Calculating fRho_0 */
			fRho_0 = J*fRho;
			
			
			/* [fRight_Cauchy_Green_tensor] will be formed */
			fRight_Cauchy_Green_tensor.MultATB(fDeformation_Gradient, fDeformation_Gradient);
			
			
			/* [fRight_Cauchy_Green_tensor_Inverse] will be formed */
			if (fRight_Cauchy_Green_tensor.Det()==0)
			    fRight_Cauchy_Green_tensor = fIdentity_matrix;
			fRight_Cauchy_Green_tensor_Inverse.Inverse(fRight_Cauchy_Green_tensor);
			
			
			/* [fLeft_Cauchy_Green_tensor] will be formed */
			fLeft_Cauchy_Green_tensor.MultABT(fDeformation_Gradient, fDeformation_Gradient);
			/* [fLeft_Cauchy_Green_tensor_Inverse] will be formed */
			if (fLeft_Cauchy_Green_tensor.Det()==0)
			    fLeft_Cauchy_Green_tensor = fIdentity_matrix;
			fLeft_Cauchy_Green_tensor_Inverse.Inverse(fLeft_Cauchy_Green_tensor);
			
			
			/* [fEulerian_effective_strain_tensor_current_IP] will be formed */
			fEulerian_effective_strain_tensor_current_IP = fLeft_Cauchy_Green_tensor_Inverse;
			fEulerian_effective_strain_tensor_current_IP *= -1; 
			fEulerian_effective_strain_tensor_current_IP += fIdentity_matrix;
			fEulerian_effective_strain_tensor_current_IP *= 0.5;
			
			
			/* extract six values of strain from symmetric eulerian strain tensor */
			Extract_six_values_from_symmetric_tensor(fEulerian_effective_strain_tensor_current_IP,fTemp_six_values);
			
			/* Save Euilerian effective strain tensor of the current IP */ 
			fEulerian_effective_strain_IPs.SetRow(IP,fTemp_six_values);
			
			
			/* Calculating J_Prim */
			if (fRight_Cauchy_Green_tensor.Det()==0)
			    fRight_Cauchy_Green_tensor = fIdentity_matrix; 
			double TempJ_Prim=fRight_Cauchy_Green_tensor.Det();
			double J_Prim=sqrt(fabs(TempJ_Prim));
			
			/* [fEffective_Second_Piola_tensor_inviscid] will be formed */
			fEffective_Second_Piola_tensor_inviscid.SetToScaled(fMaterial_Params[kLambda]*log(J_Prim)-fMaterial_Params[kMu],fRight_Cauchy_Green_tensor_Inverse); 
			fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fIdentity_matrix);
			fEffective_Second_Piola_tensor_inviscid += fTemp_matrix_nsd_x_nsd;
			
			
			/* [fEffective_Kirchhoff_tensor_inviscid] will be formed */
			fEffective_Kirchhoff_tensor_inviscid.MultABCT(fDeformation_Gradient,fEffective_Second_Piola_tensor_inviscid,fDeformation_Gradient);
			
			
			
			/* {fEffective_Kirchhoff_vector_inviscid} will be formed */
			Form_effective_kirchhoff_stress_vector_inviscid();
			
			
			/* [fIota_temp_matrix] will be formed */
			fIota_temp_matrix.MultATB(fShapeSolidGrad,fDefGradInv_Grad_grad);
			
			
			/* second derivatives of solid shape functions, [fShapeSolidGradGrad] will be formed */
			fShapes_displ->Grad_GradNa(fShapeSolidGradGrad);
			
			
			
			/* [fVarpi_temp_matrix] will be formed */
			Form_Varpi_temp_matrix();
			
			
			/* hydraulic conductivity matrix in the current coordinate, [k] will be formed */
			fK_hydraulic_conductivity_matrix.SetToScaled(fMaterial_Params[kK],fIdentity_matrix); 
			fk_hydraulic_conductivity_matrix.SetToScaled(1/J,fK_hydraulic_conductivity_matrix); 
			fTemp_matrix_nsd_x_nsd.MultABCT(fDeformation_Gradient,fk_hydraulic_conductivity_matrix,fDeformation_Gradient);
			fk_hydraulic_conductivity_matrix = fTemp_matrix_nsd_x_nsd;

//			fk_hydraulic_conductivity_matrix = fK_hydraulic_conductivity_matrix;
			
			/* [fLambda_temp_matrix] will be formed */
			fLambda_temp_matrix.MultATBC(fShapeFluidGrad,fDeformation_Gradient_Inverse,fk_hydraulic_conductivity_matrix);
			
			
			/* {fChi_temp_vector} will be formed */
			fVarpi_temp_matrix.Multx(u_vec,fChi_temp_vector);
			
			
			/* [fChi_temp_column_matrix] will be formed */
			for (int i=0; i<3 ; i++)
			    fChi_temp_column_matrix(i,0)= fChi_temp_vector[i];
			
			
			/* {fFd_int_N1_vector} will be formed */
			double scale = scale_const;
			fIota_temp_matrix.Multx(fEffective_Kirchhoff_vector_inviscid,fTemp_vector_ndof_se,scale);
			/* fFd_int_N1_vector for the current IP */
			/* accumulate */
			fFd_int_N1_vector += fTemp_vector_ndof_se;
			
			/* {fFd_int_N2_vector} will be formed */
			scale = -1.0*theta*scale_const;
			fShapeSolidGrad.MultTx(fDefGradInv_vector,fTemp_vector_ndof_se);
			fTemp_vector_ndof_se *= scale;
			/* accumulate */
			fFd_int_N2_vector += fTemp_vector_ndof_se; 
			
			
			/* state vaiables(volume fractions) for the current IP will be saved */
			fState_variables_IPs(IP,0)=phi_s;
			fState_variables_IPs(IP,1)=phi_f;
			
			
			/* {fFtheta_int_N1_vector} will be formed */
			double const1 = fMaterial_Params[kg]* phi_f;
			if (fabs(const1)>1.0e-16)
			    scale = theta*(1.0/fMaterial_Params[kg]) * ((phi_s/phi_f)-1) * scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
			fTemp_matrix_nen_press_x_nsd.Multx(fChi_temp_vector, fTemp_vector_nen_press,scale);
			/* accumulate */
			fFtheta_int_N1_vector += fTemp_vector_nen_press;
			
			
			/* {fFtheta_int_N2_vector} will be formed */
			fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_nsd,fShapeFluidGrad);
			scale = 1.0/fMaterial_Params[kg]*scale_const; 
			fTemp_matrix_nen_press_x_nen_press.Multx(press_vec, fTemp_vector_nen_press,scale);
			/* accumulate */
			fFtheta_int_N2_vector += fTemp_vector_nen_press;
			
			
			/* [fIm_temp_matrix] will be formed */
			Form_Im_temp_matrix();
			
			
			/* [fHbar_temp_matrix] will be formed */
			Form_Hbar_temp_matrix();
			
			
			/* [fEll_temp_matrix] will be formed */
			Form_Ell_temp_matrix();
			
			
			/* {fPi_temp_transpose_vector} will be formed */
			fShapeSolidGrad.MultTx(fDefGradInv_vector,fPi_temp_transpose_vector);
			
			
			/* [fPi_temp_row_matrix] will be formed */
			for (int i=0; i<n_en_displ_x_n_sd; i++)
			    fPi_temp_row_matrix(0,i) = fPi_temp_transpose_vector[i];
			
			/* [fI_ij_column_matrix] will be formed */
			fI_ij_column_matrix = 0.0;
			fI_ij_column_matrix(0,0) = 1.0;
			fI_ij_column_matrix(4,0) = 1.0;
			fI_ij_column_matrix(8,0) = 1.0;
			
			
			/* {fGrad_1_J_vector} will be filled */
			fVarpi_temp_matrix.Multx(u_vec,fGrad_1_J_vector, -1.0/J);
			
			/* {fGrad_theta_vector} will be filled */
			fShapeFluidGrad.Multx(press_vec, fGrad_theta_vector);
			
			
			/* {fGrad_phi_f_vector} will be filled */
			fGrad_phi_f_vector.SetToScaled(-1* fMaterial_Params[kPhi_s0],fGrad_1_J_vector);
			
			
			/* {fGrad_Omega_vector} will be filled */
			fTemp_nsd_vector.SetToScaled(theta/J,fGrad_phi_f_vector) ; 
			fGrad_Omega_vector = fTemp_nsd_vector;
			fTemp_nsd_vector.SetToScaled(phi_f/J,fGrad_theta_vector) ;
			fGrad_Omega_vector += fTemp_nsd_vector;
			fTemp_nsd_vector.SetToScaled(phi_f * theta,fGrad_1_J_vector) ;
			fGrad_Omega_vector += fTemp_nsd_vector;
			
			
			/* {fgrad_Omega_vector} will be formed */
			fDeformation_Gradient_Inverse_Transpose.Multx(fGrad_Omega_vector,fgrad_Omega_vector);
			
			
			/* [fJmath_temp_matrix] will be formed */
			Form_Jmath_temp_matrix(); 
			
			/* [fWp_temp_matrix] will be formed */
			Form_Wp_temp_matrix(); 
			
			
			/* Creating Second tangential elasticity tensor in the Ref. coordinate [fC_matrix] */
			Form_C_matrix(J_Prim);

                        /* [fEffective_Second_Piola_tensor_viscous] will be formed */
			double fAlpha = fMaterial_Params[kAlpha];
			Form_fEffective_Second_Piola_tensor_viscous(fAlpha);

			/* [fEffective_Kirchhoff_tensor_viscous] will be formed */
			fEffective_Kirchhoff_tensor_viscous.MultABCT(fDeformation_Gradient,fEffective_Second_Piola_tensor_viscous,fDeformation_Gradient);

                         /* [fEffective_Kirchhoff_tensor] will be formed */
			fEffective_Kirchhoff_tensor = fEffective_Kirchhoff_tensor_inviscid;
			fEffective_Kirchhoff_tensor += fEffective_Kirchhoff_tensor_viscous;


			/* [fCauchy_effective_stress_tensor_current_IP] will be formed */
			fCauchy_effective_stress_tensor_current_IP = fEffective_Kirchhoff_tensor;
			fCauchy_effective_stress_tensor_current_IP *= 1/J;
			
			
			/* extract six values of stress from symmetric cauchy stress tensor */
			Extract_six_values_from_symmetric_tensor(fCauchy_effective_stress_tensor_current_IP,fTemp_six_values);
			
			
			
			/* Save Cauchy effective stress tensor of the current IP */ 
			fCauchy_effective_stress_IPs.SetRow(IP,fTemp_six_values); 
			

			/* Creating Second tangential elasticity tensor in the Current coordinate [fc_matrix]*/
			Form_c_matrix();
			
			
			/* [fIm_Prim_temp_matrix] will be formed */
			Form_Im_Prim_temp_matrix();
			
			/* [fUpsilon_temp_matrix] will be formed */
			fUpsilon_temp_matrix.MultATB(fShapeSolid,fShapeSolid);
			
			/* [fM_dd_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se = fUpsilon_temp_matrix;
			scale = fRho_0*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fM_dd_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fC_dd_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fc_matrix,fIota_temp_matrix);
			scale = fMaterial_Params[kAlpha]*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fC_dd_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* {fFd_int_G4_vector} will be formed */
			fShapeSolid.MultTx(fGravity_vector,fTemp_vector_ndof_se);
			scale = -1*fRho_0*scale_const;
			fTemp_vector_ndof_se *= scale; 
			/* accumulate */
			fFd_int_G4_vector += fTemp_vector_ndof_se;
			
			
			/* [fM_thetad_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultAB(fLambda_temp_matrix,fShapeSolid);
			scale = 1/(fMaterial_Params[kg])*J*fRho_f*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fM_thetad_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fC_thetatheta_matrix] will be formed */
			fTemp_matrix_nen_press_x_nen_press.MultATB(fShapeFluid_row_matrix,fShapeFluid_row_matrix);
			scale = phi_f*fRho_f/(fMaterial_Params[kKf])*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fC_thetatheta_matrix += fTemp_matrix_nen_press_x_nen_press;

			
			/* [fC_thetad_matrix] will be formed */		
			fTemp_matrix_nen_press_x_ndof_se.MultATBC(fShapeFluid_row_matrix,fDefGradInv_column_matrix_Transpose,fShapeSolidGrad);
			scale = J*fRho_f*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fC_thetad_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* {fFtheta_int_H4_vector} will be formed */
			fLambda_temp_matrix.Multx(fGravity_vector,fTemp_vector_nen_press);
			scale = -1/(fMaterial_Params[kg])*J*fRho_f*scale_const;
			fTemp_vector_nen_press *= scale ;
			/* accumulate */
			fFtheta_int_H4_vector += fTemp_vector_nen_press;
			
			
			/* fC1, fC2 and fC3 will be formed */
			fC1 = phi_f*fRho_f/(fMaterial_Params[kKf]*J);
			fC2 = 1/J*(fRho_f*fMaterial_Params[kPhi_s0]-
				   phi_f*fRho_f*theta/fMaterial_Params[kKf]);
			fC3 = fC2 - fMaterial_Params[kRho_sR0]*fMaterial_Params[kPhi_s0]/J;
			
			
			/* {fgradv_vector} will be formed */
			Form_gradv_vector();
			
			
			/* [fXi_temp_matrix] will be formed */
			Form_Xi_temp_matrix();
			
			
			/* [fVarsigma_temp_matrix] will be formed */
			Form_Varsigma_temp_matrix();
			
			
			/* [fI_ijkl_matrix] will be formed */
			Form_I_ijkl_matrix();
			
			
			/* [fAleph_temp_matrix] will be formed */
			Form_Aleph_temp_matrix(IP);
			
			
		    }
		    
		    
		    /* saving eulerian effective strain for each IPs of the current element */
		    fEulerian_effective_strain_Elements_IPs.SetRow(e,fEulerian_effective_strain_IPs);
		    
		    
		    /* saving cauchy effective stress for each IPs of the current element */
		    fCauchy_effective_stress_Elements_IPs.SetRow(e,fCauchy_effective_stress_IPs);

                    /* saving physical pore water pressure for each IPs of the current element */
		    fPhysical_pore_water_pressure_Elements_IPs.SetRow(e,fPhysical_pore_water_pressure_IPs);
		    
		    /* saving state variables for each IPs of the current element */
		    fState_variables_Elements_IPs.SetRow(e,fState_variables_IPs);
		    
		    
		    /* {fFd_int_M_vector} will be formed */	    
		    fM_dd_matrix.Multx(u_dotdot_vec,fFd_int_M_vector);
			
			
		    /* {fFd_int_C_vector} will be formed */
		    fC_dd_matrix.Multx(u_dot_vec,fFd_int_C_vector);
			
		    /* {fFd_int} will be formed */
		    fFd_int = fFd_int_M_vector;
		    fFd_int += fFd_int_C_vector;
		    fFd_int += fFd_int_N1_vector;
		    fFd_int += fFd_int_N2_vector; 
		    fFd_int += fFd_int_G4_vector;
		    fFd_int *=-1;
		    
		    /* {fFtheta_int_M_vector} will be formed */
		    fM_thetad_matrix.Multx(u_dotdot_vec,fFtheta_int_M_vector);
		    
		    
		    /* {fFtheta_int_C2_vector} will be formed */
		    fC_thetad_matrix.Multx(u_dot_vec,fFtheta_int_C2_vector);
		    
		    /* {fFtheta_int_C1_vector} will be formed */
		    fC_thetatheta_matrix.Multx(press_dotdot_vec,fFtheta_int_C1_vector);

		    /* {fFtheta_int} will be formed */
		    fFtheta_int = fFtheta_int_C1_vector;
		    fFtheta_int += fFtheta_int_M_vector; 
		    fFtheta_int += fFtheta_int_C2_vector; 
		    fFtheta_int += fFtheta_int_N1_vector;
		    fFtheta_int += fFtheta_int_N2_vector;
		    fFtheta_int += fFtheta_int_H4_vector;
		    fFtheta_int *=-1;	
		    

		    fKdd = fM_dd_matrix;
		    
		    fKthetatheta = fC_thetatheta_matrix;

		    fKdtheta = 0.0;

		    fKthetad = fM_thetad_matrix;

		    /* equations numbers */
		    const iArrayT& displ_eq = fElementCards_displ[e].Equations();
		    const iArrayT& press_eq = fElementCards_press[e].Equations();

		    
		    /* assemble residuals */
		    ElementSupport().AssembleRHS(curr_group, fFd_int, displ_eq);
		    ElementSupport().AssembleRHS(curr_group, fFtheta_int, press_eq);
		    
		    /* assemble components of the tangent */
		    ElementSupport().AssembleLHS(curr_group, fKdd, displ_eq);
		    ElementSupport().AssembleLHS(curr_group, fKthetatheta, press_eq);
		    ElementSupport().AssembleLHS(curr_group, fKthetad, press_eq, displ_eq);
		    ElementSupport().AssembleLHS(curr_group, fKdtheta, displ_eq, press_eq);

		    
		}
		else if(kAnalysisType==1)
		{
		    /* residual and tangent for displacements */
                    /* solving consolidation problem */
		    const double* Det    = fShapes_displ->IPDets();
		    const double* Weight = fShapes_displ->IPWeights();
		    fShapes_displ->TopIP();
		    fShapes_press->TopIP();
		    while (fShapes_displ->NextIP() && fShapes_press->NextIP())
		    {
			double scale_const = (*Weight++)*(*Det++);
			double beta_delta_t,beta_delta_t2,gamma_delta_t,alpha_delta_t;
                        // alpha in Trapezoidal rule = 0.5
			alpha_delta_t = 0.5*delta_t;
                        /* Trapezoidal integrator will solve for velocity */
			gamma_delta_t = 1;

		
			const int IP = fShapes_displ->CurrIP();	
			dArrayT SolidIPCoordinate(n_sd),FluidIPCoordinate(n_sd);
			fShapes_displ->IPCoords(SolidIPCoordinate);
			fShapes_press->IPCoords(FluidIPCoordinate);


			const double* shapes_displ_X = fShapes_displ->IPShapeX();
			/* [fShapeSolid]will be formed */
			Form_solid_shape_functions(shapes_displ_X);
			
			
			fShapes_displ->GradNa(fShapeSolidGrad_temp);
			/* [fShapeSolidGrad] will be formed */
			Form_Gradient_of_solid_shape_functions(fShapeSolidGrad_temp);
			
			
			/* [fShapeSolidGrad_t] and [fShapeSolidGrad_t_Transpose] will be formed */
			Form_Gradient_t_of_solid_shape_functions(fShapeSolidGrad_temp);
			fShapeSolidGrad_t_Transpose.Transpose(fShapeSolidGrad_t);
		
		
			const double* shapes_press_X = fShapes_press->IPShapeX();
			/* {fShapeFluid} will be formed */
			Form_fluid_shape_functions(shapes_press_X);
			
			/* [fShapeFluid_row_matrix] will be formed */				
			for (int i=0; i<n_en_press ; i++)
			    fShapeFluid_row_matrix(0,i) = fShapeFluid[i];
			
			
			/* [fShapeFluidGrad] will be formed */
			fShapes_press->GradNa(fShapeFluidGrad);
			
			
			
			/* [fDeformation_Gradient] will be formed */
			Form_deformation_gradient_tensor();

                       /* [fDefGradT_9x9_matrix] will be formed */
			Form_fDefGradT_9x9_matrix();
			
			
			/* [fIdentity_matrix] will be formed */
			fIdentity_matrix = 0.0;			
			for (int i=0; i<n_sd ; i++)
			    fIdentity_matrix(i,i) =1.0;
		
		
			/* [fDeformation_Gradient_Inverse] and [fDeformation_Gradient_Transpose] and [fDeformation_Gradient_Inverse_Transpose] will be formed */
			if (fDeformation_Gradient.Det()==0)
			    fDeformation_Gradient = fIdentity_matrix; 
			fDeformation_Gradient_Inverse.Inverse(fDeformation_Gradient);
			fDeformation_Gradient_Inverse_Transpose.Transpose(fDeformation_Gradient_Inverse);
			fDeformation_Gradient_Transpose.Transpose(fDeformation_Gradient);
			
			
			/* {fDefGradInv_vector} will be formed */
			Form_deformation_gradient_inv_vector();
			
			/* [fDefGradInv_column_matrix] will be formed */
			for (int i=0; i<n_sd_x_n_sd; i++)
			    fDefGradInv_column_matrix(i,0)=fDefGradInv_vector[i];
			
			/* [fDefGradInv_column_matrix_Transpose] will be formed */		
			fDefGradInv_column_matrix_Transpose.Transpose(fDefGradInv_column_matrix);
			
			/* [fDefGradInv_Grad_grad] will be formed */
			Form_Grad_grad_transformation_matrix();
			
			
			/* [fDefGradInv_Grad_grad_Transpose] will be formed */
			fDefGradInv_Grad_grad_Transpose.Transpose(fDefGradInv_Grad_grad);
			
			/* Calculating theta */
			theta = fShapeFluid[0]*press_vec[0];
			for (int i=1; i<8; i++)
			    theta += fShapeFluid[i]*press_vec[i];
			
			
			/* Calculating Jacobian */
			double J = fDeformation_Gradient.Det();


			/* Jacobian for the current IP will be saved */
			fState_variables_IPs(IP,2)=J;

			
			/* Calculating fP_f */
			double fP_f=theta/J;

                        /* Physical pore water pressure for the current IP will be saved */
			fPhysical_pore_water_pressure_IPs(IP,0)=fP_f;

			
			/* Calculating fRho_f */
			fRho_f = fMaterial_Params[kRho_fR0]*exp((fP_f-fPf_0_matrix(CurrElementNumber(),IP))/
								fMaterial_Params[kKf]);
			
			
			
			/* Calculating phi_s and phi_f, volume fractions */ 
			phi_s = fMaterial_Params[kPhi_s0]/J;
			phi_f = 1.0 - phi_s;
			
			
			/*  Calculating fRho */
			fRho = phi_f*fRho_f+ phi_s*fMaterial_Params[kRho_sR0];
			
			/* Calculating fRho_0 */
			fRho_0 = J*fRho;
			
			
			/* [fRight_Cauchy_Green_tensor] will be formed */
			fRight_Cauchy_Green_tensor.MultATB(fDeformation_Gradient, fDeformation_Gradient);
		
		
			/* [fRight_Cauchy_Green_tensor_Inverse] will be formed */
			if (fRight_Cauchy_Green_tensor.Det()==0)
			    fRight_Cauchy_Green_tensor = fIdentity_matrix;
			fRight_Cauchy_Green_tensor_Inverse.Inverse(fRight_Cauchy_Green_tensor);
			
			
			/* [fLeft_Cauchy_Green_tensor] will be formed */
			fLeft_Cauchy_Green_tensor.MultABT(fDeformation_Gradient, fDeformation_Gradient);
			/* [fLeft_Cauchy_Green_tensor_Inverse] will be formed */
			if (fLeft_Cauchy_Green_tensor.Det()==0)
			    fLeft_Cauchy_Green_tensor = fIdentity_matrix;
			fLeft_Cauchy_Green_tensor_Inverse.Inverse(fLeft_Cauchy_Green_tensor);
			
			
			/* [fEulerian_effective_strain_tensor_current_IP] will be formed */
			fEulerian_effective_strain_tensor_current_IP = fLeft_Cauchy_Green_tensor_Inverse;
			fEulerian_effective_strain_tensor_current_IP *= -1; 
			fEulerian_effective_strain_tensor_current_IP += fIdentity_matrix;
			fEulerian_effective_strain_tensor_current_IP *= 0.5;
			
			
			/* extract six values of strain from symmetric eulerian strain tensor */
			Extract_six_values_from_symmetric_tensor(fEulerian_effective_strain_tensor_current_IP,fTemp_six_values);
			
			/* Save Euilerian effective strain tensor of the current IP */ 
			fEulerian_effective_strain_IPs.SetRow(IP,fTemp_six_values);
			
			
			/* Calculating J_Prim */
			if (fRight_Cauchy_Green_tensor.Det()==0)
			    fRight_Cauchy_Green_tensor = fIdentity_matrix; 
			double TempJ_Prim=fRight_Cauchy_Green_tensor.Det();
			double J_Prim=sqrt(fabs(TempJ_Prim));
			
			/* [fEffective_Second_Piola_tensor_inviscid] will be formed */
			fEffective_Second_Piola_tensor_inviscid.SetToScaled(fMaterial_Params[kLambda]*log(J_Prim)-fMaterial_Params[kMu],fRight_Cauchy_Green_tensor_Inverse); 
			fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fIdentity_matrix);
			fEffective_Second_Piola_tensor_inviscid += fTemp_matrix_nsd_x_nsd;
			
			
			/* [fEffective_Kirchhoff_tensor_inviscid] will be formed */
			fEffective_Kirchhoff_tensor_inviscid.MultABCT(fDeformation_Gradient,fEffective_Second_Piola_tensor_inviscid,fDeformation_Gradient);
			
			
			/* {fEffective_Kirchhoff_vector_inviscid} will be formed */
			Form_effective_kirchhoff_stress_vector_inviscid();
			
			
			/* [fIota_temp_matrix] will be formed */
			fIota_temp_matrix.MultATB(fShapeSolidGrad,fDefGradInv_Grad_grad);
			
			
			/* second derivatives of solid shape functions, [fShapeSolidGradGrad] will be formed */
			fShapes_displ->Grad_GradNa(fShapeSolidGradGrad);
			
			
			
			/* [fVarpi_temp_matrix] will be formed */
			Form_Varpi_temp_matrix();
			
			
			/* hydraulic conductivity matrix in the current coordinate, [k] will be formed */
			fK_hydraulic_conductivity_matrix.SetToScaled(fMaterial_Params[kK],fIdentity_matrix); 
			fk_hydraulic_conductivity_matrix.SetToScaled(1/J,fK_hydraulic_conductivity_matrix); 
			fTemp_matrix_nsd_x_nsd.MultABCT(fDeformation_Gradient,fk_hydraulic_conductivity_matrix,fDeformation_Gradient);
			fk_hydraulic_conductivity_matrix = fTemp_matrix_nsd_x_nsd;

//			fk_hydraulic_conductivity_matrix = fK_hydraulic_conductivity_matrix;
			
			
			/* [fLambda_temp_matrix] will be formed */
			fLambda_temp_matrix.MultATBC(fShapeFluidGrad,fDeformation_Gradient_Inverse,fk_hydraulic_conductivity_matrix);
			
			
			/* {fChi_temp_vector} will be formed */
			fVarpi_temp_matrix.Multx(u_vec,fChi_temp_vector);
			
			
			/* [fChi_temp_column_matrix] will be formed */
			for (int i=0; i<3 ; i++)
			    fChi_temp_column_matrix(i,0)= fChi_temp_vector[i];
			
			
			/* {fFd_int_N1_vector} will be formed */
			double scale = scale_const;
			fIota_temp_matrix.Multx(fEffective_Kirchhoff_vector_inviscid,fTemp_vector_ndof_se,scale);
			/* fFd_int_N1_vector for the current IP */
			/* accumulate */
			fFd_int_N1_vector += fTemp_vector_ndof_se;
			
			/* {fFd_int_N2_vector} will be formed */
			scale = -1.0*theta*scale_const;
			fShapeSolidGrad.MultTx(fDefGradInv_vector,fTemp_vector_ndof_se);
			fTemp_vector_ndof_se *= scale;
			/* accumulate */
			fFd_int_N2_vector += fTemp_vector_ndof_se; 
			
			
			/* state vaiables(volume fractions) for the current IP will be saved */
			fState_variables_IPs(IP,0)=phi_s;
			fState_variables_IPs(IP,1)=phi_f;
			
			
			/* {fFtheta_int_N1_vector} will be formed */
			double const1 = fMaterial_Params[kg]* phi_f;
			if (fabs(const1)>1.0e-16)
			    scale = theta*(1.0/fMaterial_Params[kg]) * ((phi_s/phi_f)-1) * scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
			fTemp_matrix_nen_press_x_nsd.Multx(fChi_temp_vector, fTemp_vector_nen_press,scale);
			/* accumulate */
			fFtheta_int_N1_vector += fTemp_vector_nen_press;
			
			
			/* {fFtheta_int_N2_vector} will be formed */
			fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_nsd,fShapeFluidGrad);
			scale = 1.0/fMaterial_Params[kg]*scale_const; 
			fTemp_matrix_nen_press_x_nen_press.Multx(press_vec, fTemp_vector_nen_press,scale);
			/* accumulate */
			fFtheta_int_N2_vector += fTemp_vector_nen_press;
			
			
			/* [fIm_temp_matrix] will be formed */
			Form_Im_temp_matrix();
			
			
			/* [fHbar_temp_matrix] will be formed */
			Form_Hbar_temp_matrix();
			
			
			/* [fEll_temp_matrix] will be formed */
			Form_Ell_temp_matrix();
			
			
			/* {fPi_temp_transpose_vector} will be formed */
			fShapeSolidGrad.MultTx(fDefGradInv_vector,fPi_temp_transpose_vector);
		
		
			/* [fPi_temp_row_matrix] will be formed */
			for (int i=0; i<n_en_displ_x_n_sd; i++)
			    fPi_temp_row_matrix(0,i) = fPi_temp_transpose_vector[i];
		
			/* [fK_dd_G3_1_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fIm_temp_matrix,fIota_temp_matrix);
			scale = -1*alpha_delta_t*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_1_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fI_ij_column_matrix] will be formed */
			fI_ij_column_matrix = 0.0;
			fI_ij_column_matrix(0,0) = 1.0;
			fI_ij_column_matrix(4,0) = 1.0;
			fI_ij_column_matrix(8,0) = 1.0;
			
			
			/* [fK_dd_G3_2_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fHbar_temp_matrix,fIota_temp_matrix);
			scale = fMaterial_Params[kMu] * alpha_delta_t * scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_2_matrix += fTemp_matrix_ndof_se_x_ndof_se;
		
		
			/* [fK_dd_G3_3_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fEll_temp_matrix,fIota_temp_matrix);
			scale = fMaterial_Params[kMu] * alpha_delta_t * scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_3_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dd_G3_4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABC(fIota_temp_matrix,fI_ij_column_matrix,fPi_temp_row_matrix);
			scale = fMaterial_Params[kLambda] * alpha_delta_t * scale_const; 
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_4_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dd_G3_5_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fShapeSolidGrad_t_Transpose,fDefGradInv_Grad_grad_Transpose,fIota_temp_matrix);
			scale = theta * alpha_delta_t * scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_5_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dtheta_G3_matrix] will be formed */
			fTemp_matrix_ndof_se_x_nen_press.MultATB(fPi_temp_row_matrix,fShapeFluid_row_matrix); 
			scale = -1*alpha_delta_t*scale_const;
			fTemp_matrix_ndof_se_x_nen_press *= scale;
			/* accumulate */
			fK_dtheta_G3_matrix += fTemp_matrix_ndof_se_x_nen_press;
			
			/* {fGrad_1_J_vector} will be filled */
			fVarpi_temp_matrix.Multx(u_vec,fGrad_1_J_vector, -1.0/J);
			
			
			/* {fGrad_theta_vector} will be filled */
			fShapeFluidGrad.Multx(press_vec, fGrad_theta_vector);
		
		
			/* {fGrad_phi_f_vector} will be filled */
			fGrad_phi_f_vector.SetToScaled(-1* fMaterial_Params[kPhi_s0],fGrad_1_J_vector);
		
			
			/* {fGrad_Omega_vector} will be filled */
			fTemp_nsd_vector.SetToScaled(theta/J,fGrad_phi_f_vector) ; 
			fGrad_Omega_vector = fTemp_nsd_vector;
			fTemp_nsd_vector.SetToScaled(phi_f/J,fGrad_theta_vector) ;
			fGrad_Omega_vector += fTemp_nsd_vector;
			fTemp_nsd_vector.SetToScaled(phi_f * theta,fGrad_1_J_vector) ;
			fGrad_Omega_vector += fTemp_nsd_vector;
			
			
			/* {fgrad_Omega_vector} will be formed */
			fDeformation_Gradient_Inverse_Transpose.Multx(fGrad_Omega_vector,fgrad_Omega_vector);
			
			
			/* [fJmath_temp_matrix] will be formed */
			Form_Jmath_temp_matrix(); 
			
			/* [fWp_temp_matrix] will be formed */
			Form_Wp_temp_matrix(); 
			
			
			/* [fK_thetad_H3_1_matrix] will be formed */
			const1 = fMaterial_Params[kg]*phi_f;
			if (fabs(const1) > 1e-16) 
			    scale = (-1*theta*phi_s/(fMaterial_Params[kg]*phi_f))*
				(phi_s/phi_f +1)*alpha_delta_t*scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
			fTemp_matrix_nen_press_x_ndof_se.MultABC(fTemp_matrix_nen_press_x_nsd,fChi_temp_column_matrix,fPi_temp_row_matrix);
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H3_1_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H3_2_matrix] will be formed */
			const1 = fMaterial_Params[kg]*phi_f ;
			if (fabs(const1) > 1e-16)
			    scale = 1*(theta/fMaterial_Params[kg])*
				(phi_s/phi_f -1)*alpha_delta_t*scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_ndof_se.MultAB(fTemp_matrix_nen_press_x_nsd,fVarpi_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H3_2_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H3_3_matrix] will be formed */
			const1 = fMaterial_Params[kg]*phi_f;
			if (fabs(const1) > 1e-16)
			    scale = -1*J/(fMaterial_Params[kg]*phi_f)*alpha_delta_t*scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nsd.MultATB(fShapeFluidGrad,fDeformation_Gradient_Inverse);
			fTemp_matrix_nen_press_x_ndof_se.MultABCT(fTemp_matrix_nen_press_x_nsd,fJmath_temp_matrix,fIota_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H3_3_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H3_4_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultABCT(fTemp_matrix_nen_press_x_nsd,fWp_temp_matrix,fIota_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H3_4_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetatheta_H3_1_matrix] will be formed */
			scale = (1/fMaterial_Params[kg])*alpha_delta_t*scale_const;
			fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
			fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_nsd,fShapeFluidGrad);
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H3_1_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			/* [fK_thetatheta_H3_2_matrix] will be formed */
			const1 = fMaterial_Params[kg] * phi_f ;
			if (fabs(const1)> 1e-16)
			    scale =  (1/fMaterial_Params[kg])*
				(phi_s/phi_f -1)*alpha_delta_t*scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nen_press.MultABC(fTemp_matrix_nen_press_x_nsd,fChi_temp_column_matrix,fShapeFluid_row_matrix);
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H3_2_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			
			/* Creating Second tangential elasticity tensor in the Ref. coordinate [fC_matrix] */
			Form_C_matrix(J_Prim);
 
                       /* [fEffective_Second_Piola_tensor_viscous] will be formed */
			double fAlpha = fMaterial_Params[kAlpha];
			Form_fEffective_Second_Piola_tensor_viscous(fAlpha);

			/* [fEffective_Kirchhoff_tensor_viscous] will be formed */
			fEffective_Kirchhoff_tensor_viscous.MultABCT(fDeformation_Gradient,fEffective_Second_Piola_tensor_viscous,fDeformation_Gradient);


                         /* [fEffective_Kirchhoff_tensor] will be formed */
			fEffective_Kirchhoff_tensor = fEffective_Kirchhoff_tensor_inviscid;
			fEffective_Kirchhoff_tensor += fEffective_Kirchhoff_tensor_viscous;
			
			/* [fCauchy_effective_stress_tensor_current_IP] will be formed */
			fCauchy_effective_stress_tensor_current_IP = fEffective_Kirchhoff_tensor;
			fCauchy_effective_stress_tensor_current_IP *= 1/J;
			
			
			/* extract six values of stress from symmetric cauchy stress tensor */
			Extract_six_values_from_symmetric_tensor(fCauchy_effective_stress_tensor_current_IP,fTemp_six_values);
			
			
			
			/* Save Cauchy effective stress tensor of the current IP */ 
			fCauchy_effective_stress_IPs.SetRow(IP,fTemp_six_values); 
			

			/* Creating Second tangential elasticity tensor in the Current coordinate [fc_matrix]*/
			Form_c_matrix();
		
		
			/* [fIm_Prim_temp_matrix] will be formed */
			Form_Im_Prim_temp_matrix();
			
			/* [fUpsilon_temp_matrix] will be formed */
			fUpsilon_temp_matrix.MultATB(fShapeSolid,fShapeSolid);
			
			/* [fM_dd_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se = fUpsilon_temp_matrix;
			scale = fRho_0*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fM_dd_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fC_dd_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fc_matrix,fIota_temp_matrix);
			scale = fMaterial_Params[kAlpha]*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fC_dd_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* {fFd_int_G4_vector} will be formed */
			fShapeSolid.MultTx(fGravity_vector,fTemp_vector_ndof_se);
			scale = -1*fRho_0*scale_const;
			fTemp_vector_ndof_se *= scale; 
			/* accumulate */
			fFd_int_G4_vector += fTemp_vector_ndof_se;
			
			
			/* [fM_thetad_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultAB(fLambda_temp_matrix,fShapeSolid);
			scale = 1/(fMaterial_Params[kg])*J*fRho_f*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fM_thetad_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fC_thetatheta_matrix] will be formed */
			fTemp_matrix_nen_press_x_nen_press.MultATB(fShapeFluid_row_matrix,fShapeFluid_row_matrix);
			scale = phi_f*fRho_f/(fMaterial_Params[kKf])*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fC_thetatheta_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			
			/* [fC_thetad_matrix] will be formed */		
			fTemp_matrix_nen_press_x_ndof_se.MultATBC(fShapeFluid_row_matrix,fDefGradInv_column_matrix_Transpose,fShapeSolidGrad);
			scale = J*fRho_f*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fC_thetad_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* {fFtheta_int_H4_vector} will be formed */
			fLambda_temp_matrix.Multx(fGravity_vector,fTemp_vector_nen_press);
			scale = -1/(fMaterial_Params[kg])*J*fRho_f*scale_const;
			fTemp_vector_nen_press *= scale ;
			/* accumulate */
			fFtheta_int_H4_vector += fTemp_vector_nen_press;
			
			
			/* fC1, fC2 and fC3 will be formed */
			fC1 = phi_f*fRho_f/(fMaterial_Params[kKf]*J);
			fC2 = 1/J*(fRho_f*fMaterial_Params[kPhi_s0]-
				   phi_f*fRho_f*theta/fMaterial_Params[kKf]);
			fC3 = fC2 - fMaterial_Params[kRho_sR0]*fMaterial_Params[kPhi_s0]/J;
			
			
			/* {fgradv_vector} will be formed */
			Form_gradv_vector();
			
			
			/* [fXi_temp_matrix] will be formed */
			Form_Xi_temp_matrix();
			
			
			/* [fVarsigma_temp_matrix] will be formed */
			Form_Varsigma_temp_matrix();
			
			
			/* [fI_ijkl_matrix] will be formed */
			Form_I_ijkl_matrix();
			
			
			/* [fK_dd_G2_1_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fXi_temp_matrix,fIota_temp_matrix);
			scale = -1*fMaterial_Params[kAlpha]*alpha_delta_t*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G2_1_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dd_G2_2_matrix] will be formed */
			fTemp_matrix1_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fI_ijkl_matrix,fIota_temp_matrix);
			fTemp_matrix_ndof_se_x_ndof_se.MultABC(fTemp_matrix1_ndof_se_x_ndof_se,u_dot_column_matrix,fPi_temp_row_matrix);
			scale = -2*fMaterial_Params[kAlpha]*fMaterial_Params[kLambda]*alpha_delta_t*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G2_2_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			

			/* [fK_dd_G2_3_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fc_matrix,fIota_temp_matrix);
			scale = fMaterial_Params[kAlpha]*gamma_delta_t*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G2_3_matrix += fTemp_matrix_ndof_se_x_ndof_se;

			
			/* [fK_dd_G2_4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fVarsigma_temp_matrix,fIota_temp_matrix);
			scale = -1*fMaterial_Params[kAlpha]*alpha_delta_t*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G2_4_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dd_G4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultATBC(fShapeSolid,fGravity_column_matrix,fPi_temp_row_matrix);
			scale = -1*J*(fC3 + fRho)*alpha_delta_t*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G4_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dtheta_G4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_nen_press.MultATBC(fShapeSolid,fGravity_column_matrix,fShapeFluid_row_matrix);
			scale = -1*J*fC1*alpha_delta_t*scale_const;
			fTemp_matrix_ndof_se_x_nen_press *= scale;
			/* accumulate */
			fK_dtheta_G4_matrix += fTemp_matrix_ndof_se_x_nen_press;
			

			/* [fAleph_temp_matrix] will be formed */
			Form_Aleph_temp_matrix(IP);
			

			/* [fK_thetad_H2_1_matrix] will be formed */
			fTemp_matrix1_nen_press_x_ndof_se.MultAB(press_dot_column_matrix,fPi_temp_row_matrix);
			fTemp_matrix_nen_press_x_ndof_se.MultATBC(fShapeFluid_row_matrix,fShapeFluid_row_matrix,fTemp_matrix1_nen_press_x_ndof_se);
			scale = 1/(fMaterial_Params[kKf])*fC2*alpha_delta_t*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H2_1_matrix += fTemp_matrix_nen_press_x_ndof_se;
			

			/* [fK_thetad_H2_2_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultATB(fShapeFluid_row_matrix,fPi_temp_row_matrix);
			scale = gamma_delta_t*J*fRho_f*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H2_2_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H2_3_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultAB(u_dot_column_matrix,fPi_temp_row_matrix);		
			fTemp_matrix_nen_press_x_ndof_se.MultATBC(fShapeFluid_row_matrix,fPi_temp_row_matrix,fTemp_matrix_ndof_se_x_ndof_se);
			scale = J*fRho_f*alpha_delta_t*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H2_3_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H2_4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultAB(u_dot_column_matrix,fPi_temp_row_matrix);
			fTemp_matrix_nen_press_x_ndof_se.MultATBC(fShapeFluid_row_matrix,fPi_temp_row_matrix,fTemp_matrix_ndof_se_x_ndof_se);
			scale = -1*fRho_f*theta/fMaterial_Params[kKf]*alpha_delta_t*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H2_4_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H2_5_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fShapeSolidGrad_t_Transpose,fDefGradInv_Grad_grad_Transpose,fIota_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se.MultATBC(fShapeFluid_row_matrix,u_dot_column_matrix_Transpose,fTemp_matrix_ndof_se_x_ndof_se);
			scale = -1*J*fRho_f*alpha_delta_t*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H2_5_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetatheta_H2_1_matrix] will be formed */
			fTemp_matrix_nen_press_x_nen_press.MultATB(fShapeFluid_row_matrix,fShapeFluid_row_matrix);
			scale = 1/(fMaterial_Params[kKf])*
			    gamma_delta_t*phi_f*fRho_f*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H2_1_matrix += fTemp_matrix_nen_press_x_nen_press;


			/* [fK_thetatheta_H2_2_matrix] will be formed */
			fTemp_matrix_nen_press_x_1.MultATBC(fShapeFluid_row_matrix,fShapeFluid_row_matrix,press_dot_column_matrix);
			fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_1,fShapeFluid_row_matrix);
			scale = fC1/(fMaterial_Params[kKf])*alpha_delta_t*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H2_2_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			
			/* [fK_thetatheta_H2_3_matrix] will be formed */
			fTemp_matrix_nen_press_x_1.MultATBC(fShapeFluid_row_matrix,fPi_temp_row_matrix,u_dot_column_matrix);
			fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_1,fShapeFluid_row_matrix);
			scale = fRho_f/(fMaterial_Params[kKf])*alpha_delta_t*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H2_3_matrix += fTemp_matrix_nen_press_x_nen_press;

			
			/* [fImath_temp_matrix] will be formed */
			Form_Imath_temp_matrix();
			
			
			/* [fK_thetad_H4_1_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultABC(fLambda_temp_matrix,fGravity_column_matrix,fPi_temp_row_matrix);
			scale = 1/(fMaterial_Params[kg]*fMaterial_Params[kKf])*
			    fRho_f*theta*alpha_delta_t*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H4_1_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H4_2_matrix] will be formed */
			fTemp_matrix_nsd_x_ndof_se.MultABCT(fDeformation_Gradient_Inverse,fImath_temp_matrix,fIota_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se.MultATB(fShapeFluidGrad,fTemp_matrix_nsd_x_ndof_se);
			scale = 1/(fMaterial_Params[kg])*J*fRho_f*alpha_delta_t*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H4_2_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H4_3_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultABC(fLambda_temp_matrix,fGravity_column_matrix,fPi_temp_row_matrix);
			scale = -1/(fMaterial_Params[kg])*J*fRho_f*alpha_delta_t*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H4_3_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetatheta_H4_matrix] will be formed */
			fTemp_matrix_nen_press_x_nen_press.MultABC(fLambda_temp_matrix,fGravity_column_matrix,fShapeFluid_row_matrix);
			scale = -1/(fMaterial_Params[kg]*fMaterial_Params[kKf])*
			    fRho_f*alpha_delta_t*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H4_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			

/*********************************************************************************************************************************************
 *********************************************************************************************************************************************
 *********************************************************************************************************************************************/
		
		
		
			/* implementing small strain deformation of an elastic media  for solid phase just for debugging */
			
			/* [fD_matrix] will be formed */
			Form_D_matrix();
		
			/* [fB_matrix] will be formed */
			Form_B_matrix();
			
			/* [fK_dd_BTDB_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultATBC(fB_matrix,fD_matrix,fB_matrix);
			scale = scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_BTDB_matrix += fTemp_matrix_ndof_se_x_ndof_se;
		


		
		
		
		

		
		
// defining F_1_T
			/* for debugging */
			const int ip = fShapes_displ->CurrIP()+1;
//			fs_mix_out	<< endl << "IP" << ip
//					<< setw(outputFileWidth) << ", shape function matrix for solid phase: " 
//					<< setw(outputFileWidth) << fShapeSolid;
			
//			fs_mix_out	<< endl << "terms from shape function matrix for solid phase: " 
//					<< setw(outputFileWidth) << fShapeSolid(0,0) 
//					<< setw(outputFileWidth) << fShapeSolid(0,3);	
		    }		
	    
		    /* saving eulerian effective strain for each IPs of the current element */
		    fEulerian_effective_strain_Elements_IPs.SetRow(e,fEulerian_effective_strain_IPs);
	    
	    
		    /* saving cauchy effective stress for each IPs of the current element */
		    fCauchy_effective_stress_Elements_IPs.SetRow(e,fCauchy_effective_stress_IPs);

                    /* saving physical pore water pressure for each IPs of the current element */
		    fPhysical_pore_water_pressure_Elements_IPs.SetRow(e,fPhysical_pore_water_pressure_IPs);
	    
	    
		    /* saving state variables for each IPs of the current element */
		    fState_variables_Elements_IPs.SetRow(e,fState_variables_IPs);

	    
		    /* {fFd_int_C_vector} will be formed */
		    fC_dd_matrix.Multx(u_dot_vec,fFd_int_C_vector);


		    /* {fFd_int} will be formed */
		    fFd_int = fFd_int_C_vector;
		    fFd_int += fFd_int_N1_vector;
		    fFd_int += fFd_int_N2_vector; 
		    fFd_int += fFd_int_G4_vector;
		    fFd_int *= -1;
	    
/*			fK_dd_BTDB_matrix.MultTx(u_vec,fTemp_vector_ndof_se);
			fFd_int = fTemp_vector_ndof_se;
			fFd_int *= -1.0; */


		    /* [fKdd] will be formed */
		    fKdd = fK_dd_G2_1_matrix;
		    fKdd += fK_dd_G2_2_matrix;
		    fKdd += fK_dd_G2_3_matrix;
		    fKdd += fK_dd_G2_4_matrix;
		    fKdd += fK_dd_G3_1_matrix;
		    fKdd += fK_dd_G3_2_matrix;
		    fKdd += fK_dd_G3_3_matrix;
		    fKdd += fK_dd_G3_4_matrix;
		    fKdd += fK_dd_G3_5_matrix; 
		    fKdd += fK_dd_G4_matrix; 



//		    fKdd = fK_dd_BTDB_matrix;			
	    
	    
		    /* [fKdtheta] will be formed */
		    fKdtheta = fK_dtheta_G3_matrix;
		    fKdtheta += fK_dtheta_G4_matrix;
	    
	    
		    /* [fKthetad] will be formed */
		    fKthetad = fK_thetad_H2_1_matrix;
		    fKthetad += fK_thetad_H2_2_matrix;
		    fKthetad += fK_thetad_H2_3_matrix;
		    fKthetad += fK_thetad_H2_4_matrix;
		    fKthetad += fK_thetad_H2_5_matrix;
		    fKthetad += fK_thetad_H3_1_matrix;
		    fKthetad += fK_thetad_H3_2_matrix;
		    fKthetad += fK_thetad_H3_3_matrix;
		    fKthetad += fK_thetad_H3_4_matrix;
		    fKthetad += fK_thetad_H4_1_matrix;
		    fKthetad += fK_thetad_H4_2_matrix;
		    fKthetad += fK_thetad_H4_3_matrix;


		    /* [fKthetatheta] will be formed */
		    fKthetatheta = fK_thetatheta_H2_1_matrix;
		    fKthetatheta += fK_thetatheta_H2_2_matrix;
		    fKthetatheta += fK_thetatheta_H2_3_matrix;
		    fKthetatheta += fK_thetatheta_H3_1_matrix;
		    fKthetatheta += fK_thetatheta_H3_2_matrix;
		    fKthetatheta += fK_thetatheta_H4_matrix;	


		    /* {fFtheta_int_C1_vector} will be formed */
		    fC_thetatheta_matrix.Multx(press_dot_vec,fFtheta_int_C1_vector);


		    /* {fFtheta_int_C2_vector} will be formed */
		    fC_thetad_matrix.Multx(u_dot_vec,fFtheta_int_C2_vector);


		    /* {fFtheta_int} will be formed */
		    fFtheta_int = fFtheta_int_C1_vector;
		    fFtheta_int += fFtheta_int_C2_vector; 
		    fFtheta_int += fFtheta_int_N1_vector;
		    fFtheta_int += fFtheta_int_N2_vector;
		    fFtheta_int += fFtheta_int_H4_vector;
		    fFtheta_int *= -1;
		    

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

		else if(kAnalysisType==2)
		{
		    /* residual and tangent for displacements */
                    /* solving dynamic problem */
		    const double* Det    = fShapes_displ->IPDets();
		    const double* Weight = fShapes_displ->IPWeights();
		    fShapes_displ->TopIP();
		    fShapes_press->TopIP();
		    while (fShapes_displ->NextIP() && fShapes_press->NextIP())
		    {
			double scale_const = (*Weight++)*(*Det++);
			double beta_delta_t,beta_delta_t2,gamma_delta_t;
			beta_delta_t2 = 0.3025*delta_t*delta_t;
			beta_delta_t = 0.3025*delta_t;
			gamma_delta_t = 0.6*delta_t;

		
			const int IP = fShapes_displ->CurrIP();	
			dArrayT SolidIPCoordinate(n_sd),FluidIPCoordinate(n_sd);
			fShapes_displ->IPCoords(SolidIPCoordinate);
			fShapes_press->IPCoords(FluidIPCoordinate);


			const double* shapes_displ_X = fShapes_displ->IPShapeX();
			/* [fShapeSolid]will be formed */
			Form_solid_shape_functions(shapes_displ_X);
			
			
			fShapes_displ->GradNa(fShapeSolidGrad_temp);
			/* [fShapeSolidGrad] will be formed */
			Form_Gradient_of_solid_shape_functions(fShapeSolidGrad_temp);
			
			
			/* [fShapeSolidGrad_t] and [fShapeSolidGrad_t_Transpose] will be formed */
			Form_Gradient_t_of_solid_shape_functions(fShapeSolidGrad_temp);
			fShapeSolidGrad_t_Transpose.Transpose(fShapeSolidGrad_t);
		
		
			const double* shapes_press_X = fShapes_press->IPShapeX();
			/* {fShapeFluid} will be formed */
			Form_fluid_shape_functions(shapes_press_X);
			
			/* [fShapeFluid_row_matrix] will be formed */				
			for (int i=0; i<n_en_press ; i++)
			    fShapeFluid_row_matrix(0,i) = fShapeFluid[i];
			
			
			/* [fShapeFluidGrad] will be formed */
			fShapes_press->GradNa(fShapeFluidGrad);
			
			
			
			/* [fDeformation_Gradient] will be formed */
			Form_deformation_gradient_tensor();

                       /* [fDefGradT_9x9_matrix] will be formed */
			Form_fDefGradT_9x9_matrix();
			
			
			/* [fIdentity_matrix] will be formed */
			fIdentity_matrix = 0.0;			
			for (int i=0; i<n_sd ; i++)
			    fIdentity_matrix(i,i) =1.0;
		
		
			/* [fDeformation_Gradient_Inverse] and [fDeformation_Gradient_Transpose] and [fDeformation_Gradient_Inverse_Transpose] will be formed */
			if (fDeformation_Gradient.Det()==0)
			    fDeformation_Gradient = fIdentity_matrix; 
			fDeformation_Gradient_Inverse.Inverse(fDeformation_Gradient);
			fDeformation_Gradient_Inverse_Transpose.Transpose(fDeformation_Gradient_Inverse);
			fDeformation_Gradient_Transpose.Transpose(fDeformation_Gradient);
			
			
			/* {fDefGradInv_vector} will be formed */
			Form_deformation_gradient_inv_vector();
			
			/* [fDefGradInv_column_matrix] will be formed */
			for (int i=0; i<n_sd_x_n_sd; i++)
			    fDefGradInv_column_matrix(i,0)=fDefGradInv_vector[i];
			
			/* [fDefGradInv_column_matrix_Transpose] will be formed */		
			fDefGradInv_column_matrix_Transpose.Transpose(fDefGradInv_column_matrix);
			
			/* [fDefGradInv_Grad_grad] will be formed */
			Form_Grad_grad_transformation_matrix();
			
			
			/* [fDefGradInv_Grad_grad_Transpose] will be formed */
			fDefGradInv_Grad_grad_Transpose.Transpose(fDefGradInv_Grad_grad);
			
			/* Calculating theta */
			theta = fShapeFluid[0]*press_vec[0];
			for (int i=1; i<8; i++)
			    theta += fShapeFluid[i]*press_vec[i];
			
			
			/* Calculating Jacobian */
			double J = fDeformation_Gradient.Det();


			/* Jacobian for the current IP will be saved */
			fState_variables_IPs(IP,2)=J;

			
			/* Calculating fP_f */
			double fP_f=theta/J;

                        /* Physical pore water pressure for the current IP will be saved */
			fPhysical_pore_water_pressure_IPs(IP,0)=fP_f;

			
			/* Calculating fRho_f */
			fRho_f = fMaterial_Params[kRho_fR0]*exp((fP_f-fPf_0_matrix(CurrElementNumber(),IP))/
								fMaterial_Params[kKf]);
			
			
			
			/* Calculating phi_s and phi_f, volume fractions */ 
			phi_s = fMaterial_Params[kPhi_s0]/J;
			phi_f = 1.0 - phi_s;
			
			
			/*  Calculating fRho */
			fRho = phi_f*fRho_f+ phi_s*fMaterial_Params[kRho_sR0];
			
			/* Calculating fRho_0 */
			fRho_0 = J*fRho;
			
			
			/* [fRight_Cauchy_Green_tensor] will be formed */
			fRight_Cauchy_Green_tensor.MultATB(fDeformation_Gradient, fDeformation_Gradient);
		
		
			/* [fRight_Cauchy_Green_tensor_Inverse] will be formed */
			if (fRight_Cauchy_Green_tensor.Det()==0)
			    fRight_Cauchy_Green_tensor = fIdentity_matrix;
			fRight_Cauchy_Green_tensor_Inverse.Inverse(fRight_Cauchy_Green_tensor);
			
			
			/* [fLeft_Cauchy_Green_tensor] will be formed */
			fLeft_Cauchy_Green_tensor.MultABT(fDeformation_Gradient, fDeformation_Gradient);
			/* [fLeft_Cauchy_Green_tensor_Inverse] will be formed */
			if (fLeft_Cauchy_Green_tensor.Det()==0)
			    fLeft_Cauchy_Green_tensor = fIdentity_matrix;
			fLeft_Cauchy_Green_tensor_Inverse.Inverse(fLeft_Cauchy_Green_tensor);
			
			
			/* [fEulerian_effective_strain_tensor_current_IP] will be formed */
			fEulerian_effective_strain_tensor_current_IP = fLeft_Cauchy_Green_tensor_Inverse;
			fEulerian_effective_strain_tensor_current_IP *= -1; 
			fEulerian_effective_strain_tensor_current_IP += fIdentity_matrix;
			fEulerian_effective_strain_tensor_current_IP *= 0.5;
			
			
			/* extract six values of strain from symmetric eulerian strain tensor */
			Extract_six_values_from_symmetric_tensor(fEulerian_effective_strain_tensor_current_IP,fTemp_six_values);
			
			/* Save Euilerian effective strain tensor of the current IP */ 
			fEulerian_effective_strain_IPs.SetRow(IP,fTemp_six_values);
			
			
			/* Calculating J_Prim */
			if (fRight_Cauchy_Green_tensor.Det()==0)
			    fRight_Cauchy_Green_tensor = fIdentity_matrix; 
			double TempJ_Prim=fRight_Cauchy_Green_tensor.Det();
			double J_Prim=sqrt(fabs(TempJ_Prim));
			
			/* [fEffective_Second_Piola_tensor_inviscid] will be formed */
			fEffective_Second_Piola_tensor_inviscid.SetToScaled(fMaterial_Params[kLambda]*log(J_Prim)-fMaterial_Params[kMu],fRight_Cauchy_Green_tensor_Inverse); 
			fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fIdentity_matrix);
			fEffective_Second_Piola_tensor_inviscid += fTemp_matrix_nsd_x_nsd;
			
			
			/* [fEffective_Kirchhoff_tensor_inviscid] will be formed */
			fEffective_Kirchhoff_tensor_inviscid.MultABCT(fDeformation_Gradient,fEffective_Second_Piola_tensor_inviscid,fDeformation_Gradient);
			
			
			/* {fEffective_Kirchhoff_vector_inviscid} will be formed */
			Form_effective_kirchhoff_stress_vector_inviscid();
			
			
			/* [fIota_temp_matrix] will be formed */
			fIota_temp_matrix.MultATB(fShapeSolidGrad,fDefGradInv_Grad_grad);
			
			
			/* second derivatives of solid shape functions, [fShapeSolidGradGrad] will be formed */
			fShapes_displ->Grad_GradNa(fShapeSolidGradGrad);
			
			
			
			/* [fVarpi_temp_matrix] will be formed */
			Form_Varpi_temp_matrix();
			
			
			/* hydraulic conductivity matrix in the current coordinate, [k] will be formed */
			fK_hydraulic_conductivity_matrix.SetToScaled(fMaterial_Params[kK],fIdentity_matrix); 
			fk_hydraulic_conductivity_matrix.SetToScaled(1/J,fK_hydraulic_conductivity_matrix); 
			fTemp_matrix_nsd_x_nsd.MultABCT(fDeformation_Gradient,fk_hydraulic_conductivity_matrix,fDeformation_Gradient);
			fk_hydraulic_conductivity_matrix = fTemp_matrix_nsd_x_nsd;


//			fk_hydraulic_conductivity_matrix = fK_hydraulic_conductivity_matrix;
			
			/* [fLambda_temp_matrix] will be formed */
			fLambda_temp_matrix.MultATBC(fShapeFluidGrad,fDeformation_Gradient_Inverse,fk_hydraulic_conductivity_matrix);
			
			
			/* {fChi_temp_vector} will be formed */
			fVarpi_temp_matrix.Multx(u_vec,fChi_temp_vector);
			
			
			/* [fChi_temp_column_matrix] will be formed */
			for (int i=0; i<3 ; i++)
			    fChi_temp_column_matrix(i,0)= fChi_temp_vector[i];
			
			
			/* {fFd_int_N1_vector} will be formed */
			double scale = scale_const;
			fIota_temp_matrix.Multx(fEffective_Kirchhoff_vector_inviscid,fTemp_vector_ndof_se,scale);
			/* fFd_int_N1_vector for the current IP */
			/* accumulate */
			fFd_int_N1_vector += fTemp_vector_ndof_se;
			
			/* {fFd_int_N2_vector} will be formed */
			scale = -1.0*theta*scale_const;
			fShapeSolidGrad.MultTx(fDefGradInv_vector,fTemp_vector_ndof_se);
			fTemp_vector_ndof_se *= scale;
			/* accumulate */
			fFd_int_N2_vector += fTemp_vector_ndof_se; 
			
			
			/* state vaiables(volume fractions) for the current IP will be saved */
			fState_variables_IPs(IP,0)=phi_s;
			fState_variables_IPs(IP,1)=phi_f;
			
			
			/* {fFtheta_int_N1_vector} will be formed */
			double const1 = fMaterial_Params[kg]* phi_f;
			if (fabs(const1)>1.0e-16)
			    scale = theta*(1.0/fMaterial_Params[kg]) * ((phi_s/phi_f)-1) * scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
			fTemp_matrix_nen_press_x_nsd.Multx(fChi_temp_vector, fTemp_vector_nen_press,scale);
			/* accumulate */
			fFtheta_int_N1_vector += fTemp_vector_nen_press;
			
			
			/* {fFtheta_int_N2_vector} will be formed */
			fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_nsd,fShapeFluidGrad);
			scale = 1.0/fMaterial_Params[kg]*scale_const; 
			fTemp_matrix_nen_press_x_nen_press.Multx(press_vec, fTemp_vector_nen_press,scale);
			/* accumulate */
			fFtheta_int_N2_vector += fTemp_vector_nen_press;
			
			
			/* [fIm_temp_matrix] will be formed */
			Form_Im_temp_matrix();
			
			
			/* [fHbar_temp_matrix] will be formed */
			Form_Hbar_temp_matrix();
			
			
			/* [fEll_temp_matrix] will be formed */
			Form_Ell_temp_matrix();
			
			
			/* {fPi_temp_transpose_vector} will be formed */
			fShapeSolidGrad.MultTx(fDefGradInv_vector,fPi_temp_transpose_vector);
		
		
			/* [fPi_temp_row_matrix] will be formed */
			for (int i=0; i<n_en_displ_x_n_sd; i++)
			    fPi_temp_row_matrix(0,i) = fPi_temp_transpose_vector[i];
		
			/* [fK_dd_G3_1_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fIm_temp_matrix,fIota_temp_matrix);
			scale = -1*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_1_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fI_ij_column_matrix] will be formed */
			fI_ij_column_matrix = 0.0;
			fI_ij_column_matrix(0,0) = 1.0;
			fI_ij_column_matrix(4,0) = 1.0;
			fI_ij_column_matrix(8,0) = 1.0;
			
			
			/* [fK_dd_G3_2_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fHbar_temp_matrix,fIota_temp_matrix);
			scale = fMaterial_Params[kMu] * beta_delta_t2 * scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_2_matrix += fTemp_matrix_ndof_se_x_ndof_se;
		
		
			/* [fK_dd_G3_3_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fEll_temp_matrix,fIota_temp_matrix);
			scale = fMaterial_Params[kMu] * beta_delta_t2 * scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_3_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dd_G3_4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABC(fIota_temp_matrix,fI_ij_column_matrix,fPi_temp_row_matrix);
			scale = fMaterial_Params[kLambda] * beta_delta_t2 * scale_const; 
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_4_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dd_G3_5_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fShapeSolidGrad_t_Transpose,fDefGradInv_Grad_grad_Transpose,fIota_temp_matrix);
			scale = theta * beta_delta_t2 * scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G3_5_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dtheta_G3_matrix] will be formed */
			fTemp_matrix_ndof_se_x_nen_press.MultATB(fPi_temp_row_matrix,fShapeFluid_row_matrix); 
			scale = -1*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_nen_press *= scale;
			/* accumulate */
			fK_dtheta_G3_matrix += fTemp_matrix_ndof_se_x_nen_press;
			
			/* {fGrad_1_J_vector} will be filled */
			fVarpi_temp_matrix.Multx(u_vec,fGrad_1_J_vector, -1.0/J);
			
			
			/* {fGrad_theta_vector} will be filled */
			fShapeFluidGrad.Multx(press_vec, fGrad_theta_vector);
		
		
			/* {fGrad_phi_f_vector} will be filled */
			fGrad_phi_f_vector.SetToScaled(-1* fMaterial_Params[kPhi_s0],fGrad_1_J_vector);
		
			
			/* {fGrad_Omega_vector} will be filled */
			fTemp_nsd_vector.SetToScaled(theta/J,fGrad_phi_f_vector) ; 
			fGrad_Omega_vector = fTemp_nsd_vector;
			fTemp_nsd_vector.SetToScaled(phi_f/J,fGrad_theta_vector) ;
			fGrad_Omega_vector += fTemp_nsd_vector;
			fTemp_nsd_vector.SetToScaled(phi_f * theta,fGrad_1_J_vector) ;
			fGrad_Omega_vector += fTemp_nsd_vector;
			
			
			/* {fgrad_Omega_vector} will be formed */
			fDeformation_Gradient_Inverse_Transpose.Multx(fGrad_Omega_vector,fgrad_Omega_vector);
			
			
			/* [fJmath_temp_matrix] will be formed */
			Form_Jmath_temp_matrix(); 
			
			/* [fWp_temp_matrix] will be formed */
			Form_Wp_temp_matrix(); 
			
			
			/* [fK_thetad_H3_1_matrix] will be formed */
			const1 = fMaterial_Params[kg]*phi_f;
			if (fabs(const1) > 1e-16) 
			    scale = (-1*theta*phi_s/(fMaterial_Params[kg]*phi_f))*
				(phi_s/phi_f +1)*beta_delta_t2*scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
			fTemp_matrix_nen_press_x_ndof_se.MultABC(fTemp_matrix_nen_press_x_nsd,fChi_temp_column_matrix,fPi_temp_row_matrix);
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H3_1_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H3_2_matrix] will be formed */
			const1 = fMaterial_Params[kg]*phi_f ;
			if (fabs(const1) > 1e-16)
			    scale = 1*(theta/fMaterial_Params[kg])*
				(phi_s/phi_f -1)*beta_delta_t2*scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_ndof_se.MultAB(fTemp_matrix_nen_press_x_nsd,fVarpi_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H3_2_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H3_3_matrix] will be formed */
			const1 = fMaterial_Params[kg]*phi_f;
			if (fabs(const1) > 1e-16)
			    scale = -1*J/(fMaterial_Params[kg]*phi_f)*beta_delta_t2*scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nsd.MultATB(fShapeFluidGrad,fDeformation_Gradient_Inverse);
			fTemp_matrix_nen_press_x_ndof_se.MultABCT(fTemp_matrix_nen_press_x_nsd,fJmath_temp_matrix,fIota_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H3_3_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H3_4_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultABCT(fTemp_matrix_nen_press_x_nsd,fWp_temp_matrix,fIota_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H3_4_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetatheta_H3_1_matrix] will be formed */
			scale = (1/fMaterial_Params[kg])*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
			fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_nsd,fShapeFluidGrad);
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H3_1_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			/* [fK_thetatheta_H3_2_matrix] will be formed */
			const1 = fMaterial_Params[kg] * phi_f ;
			if (fabs(const1)> 1e-16)
			    scale =  (1/fMaterial_Params[kg])*
				(phi_s/phi_f -1)*beta_delta_t2*scale_const;
			else
			    scale = 0.0;
			fTemp_matrix_nen_press_x_nen_press.MultABC(fTemp_matrix_nen_press_x_nsd,fChi_temp_column_matrix,fShapeFluid_row_matrix);
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H3_2_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			
			/* Creating Second tangential elasticity tensor in the Ref. coordinate [fC_matrix] */
			Form_C_matrix(J_Prim);

                       /* [fEffective_Second_Piola_tensor_viscous] will be formed */
			double fAlpha = fMaterial_Params[kAlpha];
			Form_fEffective_Second_Piola_tensor_viscous(fAlpha);

			/* [fEffective_Kirchhoff_tensor_viscous] will be formed */
			fEffective_Kirchhoff_tensor_viscous.MultABCT(fDeformation_Gradient,fEffective_Second_Piola_tensor_viscous,fDeformation_Gradient);


                         /* [fEffective_Kirchhoff_tensor] will be formed */
			fEffective_Kirchhoff_tensor = fEffective_Kirchhoff_tensor_inviscid;
			fEffective_Kirchhoff_tensor += fEffective_Kirchhoff_tensor_viscous;


			/* [fCauchy_effective_stress_tensor_current_IP] will be formed */
			fCauchy_effective_stress_tensor_current_IP = fEffective_Kirchhoff_tensor;
			fCauchy_effective_stress_tensor_current_IP *= 1/J;
			
			
			/* extract six values of stress from symmetric cauchy stress tensor */
			Extract_six_values_from_symmetric_tensor(fCauchy_effective_stress_tensor_current_IP,fTemp_six_values);
			
			
			
			/* Save Cauchy effective stress tensor of the current IP */ 
			fCauchy_effective_stress_IPs.SetRow(IP,fTemp_six_values); 
			

			/* Creating Second tangential elasticity tensor in the Current coordinate [fc_matrix]*/
			Form_c_matrix();
		
		
			/* [fIm_Prim_temp_matrix] will be formed */
			Form_Im_Prim_temp_matrix();
			
			/* [fUpsilon_temp_matrix] will be formed */
			fUpsilon_temp_matrix.MultATB(fShapeSolid,fShapeSolid);
			
			/* [fM_dd_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se = fUpsilon_temp_matrix;
			scale = fRho_0*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fM_dd_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fC_dd_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fc_matrix,fIota_temp_matrix);
			scale = fMaterial_Params[kAlpha]*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fC_dd_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* {fFd_int_G4_vector} will be formed */
			fShapeSolid.MultTx(fGravity_vector,fTemp_vector_ndof_se);
			scale = -1*fRho_0*scale_const;
			fTemp_vector_ndof_se *= scale; 
			/* accumulate */
			fFd_int_G4_vector += fTemp_vector_ndof_se;
			
			
			/* [fM_thetad_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultAB(fLambda_temp_matrix,fShapeSolid);
			scale = 1/(fMaterial_Params[kg])*J*fRho_f*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fM_thetad_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fC_thetatheta_matrix] will be formed */
			fTemp_matrix_nen_press_x_nen_press.MultATB(fShapeFluid_row_matrix,fShapeFluid_row_matrix);
			scale = phi_f*fRho_f/(fMaterial_Params[kKf])*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fC_thetatheta_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			
			/* [fC_thetad_matrix] will be formed */		
			fTemp_matrix_nen_press_x_ndof_se.MultATBC(fShapeFluid_row_matrix,fDefGradInv_column_matrix_Transpose,fShapeSolidGrad);
			scale = J*fRho_f*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fC_thetad_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* {fFtheta_int_H4_vector} will be formed */
			fLambda_temp_matrix.Multx(fGravity_vector,fTemp_vector_nen_press);
			scale = -1/(fMaterial_Params[kg])*J*fRho_f*scale_const;
			fTemp_vector_nen_press *= scale ;
			/* accumulate */
			fFtheta_int_H4_vector += fTemp_vector_nen_press;
			
			
			/* fC1, fC2 and fC3 will be formed */
			fC1 = phi_f*fRho_f/(fMaterial_Params[kKf]*J);
			fC2 = 1/J*(fRho_f*fMaterial_Params[kPhi_s0]-
				   phi_f*fRho_f*theta/fMaterial_Params[kKf]);
			fC3 = fC2 - fMaterial_Params[kRho_sR0]*fMaterial_Params[kPhi_s0]/J;
			
			
			/* [fK_dd_G1_1_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se = fUpsilon_temp_matrix;
			scale = fRho_0*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G1_1_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dd_G1_2_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABC(fUpsilon_temp_matrix,u_dotdot_column_matrix,fPi_temp_row_matrix);
			scale = J*(fC3 + fRho_0)*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G1_2_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dtheta_G1_matrix] will be formed */
			fTemp_matrix_ndof_se_x_nen_press.MultABC(fUpsilon_temp_matrix,u_dotdot_column_matrix,fShapeFluid_row_matrix);
			scale = J*fC1*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_nen_press *= scale;
			/* accumulate */
			fK_dtheta_G1_matrix += fTemp_matrix_ndof_se_x_nen_press;
			
			
			/* {fgradv_vector} will be formed */
			Form_gradv_vector();
			
			
			/* [fXi_temp_matrix] will be formed */
			Form_Xi_temp_matrix();
			
			
			/* [fVarsigma_temp_matrix] will be formed */
			Form_Varsigma_temp_matrix();
			
			
			/* [fI_ijkl_matrix] will be formed */
			Form_I_ijkl_matrix();
			
			
			/* [fK_dd_G2_1_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fXi_temp_matrix,fIota_temp_matrix);
			scale = -1*fMaterial_Params[kAlpha]*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G2_1_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dd_G2_2_matrix] will be formed */
			fTemp_matrix1_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fI_ijkl_matrix,fIota_temp_matrix);
			fTemp_matrix_ndof_se_x_ndof_se.MultABC(fTemp_matrix1_ndof_se_x_ndof_se,u_dot_column_matrix,fPi_temp_row_matrix);
			scale = -2*fMaterial_Params[kAlpha]*fMaterial_Params[kLambda]*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G2_2_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			

			/* [fK_dd_G2_3_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fc_matrix,fIota_temp_matrix);
			scale = fMaterial_Params[kAlpha]*gamma_delta_t*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G2_3_matrix += fTemp_matrix_ndof_se_x_ndof_se;

			
			/* [fK_dd_G2_4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fVarsigma_temp_matrix,fIota_temp_matrix);
			scale = -1*fMaterial_Params[kAlpha]*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G2_4_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dd_G4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultATBC(fShapeSolid,fGravity_column_matrix,fPi_temp_row_matrix);
			scale = -1*J*(fC3 + fRho)*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_G4_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
			/* [fK_dtheta_G4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_nen_press.MultATBC(fShapeSolid,fGravity_column_matrix,fShapeFluid_row_matrix);
			scale = -1*J*fC1*beta_delta_t2*scale_const;
			fTemp_matrix_ndof_se_x_nen_press *= scale;
			/* accumulate */
			fK_dtheta_G4_matrix += fTemp_matrix_ndof_se_x_nen_press;
			

			/* [fAleph_temp_matrix] will be formed */
			Form_Aleph_temp_matrix(IP);
			

			/* [fK_thetad_H1_1_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultAB(fLambda_temp_matrix,fShapeSolid);
			scale = 1/fMaterial_Params[kg]*J*fRho_f*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H1_1_matrix += fTemp_matrix_nen_press_x_ndof_se;

			
			/* [fK_thetad_H1_2_matrix] will be formed */
			fTemp_matrix_nsd_x_ndof_se.MultABC(fShapeSolid,u_dotdot_column_matrix,fPi_temp_row_matrix);
			fTemp_matrix_nen_press_x_ndof_se.MultAB(fLambda_temp_matrix,fTemp_matrix_nsd_x_ndof_se);
			scale = 1/fMaterial_Params[kg]*J*fRho_f*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H1_2_matrix += fTemp_matrix_nen_press_x_ndof_se;
			

			/* [fK_thetad_H1_3_matrix] will be formed */
			fTemp_matrix_nsd_x_ndof_se.MultABCT(fDeformation_Gradient_Inverse,fAleph_temp_matrix,fIota_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se.MultATB(fShapeFluidGrad,fTemp_matrix_nsd_x_ndof_se);
			scale = -1/fMaterial_Params[kg]*J*fRho_f*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H1_3_matrix += fTemp_matrix_nen_press_x_ndof_se;

 
			/* [fK_thetad_H1_4_matrix] will be formed */
			fTemp_matrix_nsd_x_ndof_se.MultABC(fShapeSolid,u_dotdot_column_matrix,fPi_temp_row_matrix);
			fTemp_matrix_nen_press_x_ndof_se.MultAB(fLambda_temp_matrix,fTemp_matrix_nsd_x_ndof_se);
			scale = 1/(fMaterial_Params[kKf]*fMaterial_Params[kg])*
			    fRho_f*theta*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H1_4_matrix += fTemp_matrix_nen_press_x_ndof_se;


			/* [fK_thetatheta_H1_matrix] will be formed */
			fTemp_matrix_nsd_x_nen_press.MultABC(fShapeSolid,u_dotdot_column_matrix,fShapeFluid_row_matrix);
			fTemp_matrix_nen_press_x_nen_press.MultAB(fLambda_temp_matrix,fTemp_matrix_nsd_x_nen_press);
			scale = 1/(fMaterial_Params[kg]*fMaterial_Params[kKf])*
			    fRho_f*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H1_matrix += fTemp_matrix_nen_press_x_nen_press;


			/* [fK_thetad_H2_1_matrix] will be formed */
			fTemp_matrix1_nen_press_x_ndof_se.MultAB(press_dot_column_matrix,fPi_temp_row_matrix);
			fTemp_matrix_nen_press_x_ndof_se.MultATBC(fShapeFluid_row_matrix,fShapeFluid_row_matrix,fTemp_matrix1_nen_press_x_ndof_se);
			scale = 1/(fMaterial_Params[kKf])*fC2*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H2_1_matrix += fTemp_matrix_nen_press_x_ndof_se;
			

			/* [fK_thetad_H2_2_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultATB(fShapeFluid_row_matrix,fPi_temp_row_matrix);
			scale = gamma_delta_t*J*fRho_f*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H2_2_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H2_3_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultAB(u_dot_column_matrix,fPi_temp_row_matrix);		
			fTemp_matrix_nen_press_x_ndof_se.MultATBC(fShapeFluid_row_matrix,fPi_temp_row_matrix,fTemp_matrix_ndof_se_x_ndof_se);
			scale = J*fRho_f*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H2_3_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H2_4_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultAB(u_dot_column_matrix,fPi_temp_row_matrix);
			fTemp_matrix_nen_press_x_ndof_se.MultATBC(fShapeFluid_row_matrix,fPi_temp_row_matrix,fTemp_matrix_ndof_se_x_ndof_se);
			scale = -1*fRho_f*theta/fMaterial_Params[kKf]*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H2_4_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H2_5_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fShapeSolidGrad_t_Transpose,fDefGradInv_Grad_grad_Transpose,fIota_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se.MultATBC(fShapeFluid_row_matrix,u_dot_column_matrix_Transpose,fTemp_matrix_ndof_se_x_ndof_se);
			scale = -1*J*fRho_f*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H2_5_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetatheta_H2_1_matrix] will be formed */
			fTemp_matrix_nen_press_x_nen_press.MultATB(fShapeFluid_row_matrix,fShapeFluid_row_matrix);
			scale = 1/(fMaterial_Params[kKf])*
			    gamma_delta_t*phi_f*fRho_f*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H2_1_matrix += fTemp_matrix_nen_press_x_nen_press;


			/* [fK_thetatheta_H2_2_matrix] will be formed */
			fTemp_matrix_nen_press_x_1.MultATBC(fShapeFluid_row_matrix,fShapeFluid_row_matrix,press_dot_column_matrix);
			fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_1,fShapeFluid_row_matrix);
			scale = fC1/(fMaterial_Params[kKf])*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H2_2_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			
			/* [fK_thetatheta_H2_3_matrix] will be formed */
			fTemp_matrix_nen_press_x_1.MultATBC(fShapeFluid_row_matrix,fPi_temp_row_matrix,u_dot_column_matrix);
			fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_1,fShapeFluid_row_matrix);
			scale = fRho_f/(fMaterial_Params[kKf])*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H2_3_matrix += fTemp_matrix_nen_press_x_nen_press;

			
			/* [fImath_temp_matrix] will be formed */
			Form_Imath_temp_matrix();
			
			
			/* [fK_thetad_H4_1_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultABC(fLambda_temp_matrix,fGravity_column_matrix,fPi_temp_row_matrix);
			scale = 1/(fMaterial_Params[kg]*fMaterial_Params[kKf])*
			    fRho_f*theta*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H4_1_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H4_2_matrix] will be formed */
			fTemp_matrix_nsd_x_ndof_se.MultABCT(fDeformation_Gradient_Inverse,fImath_temp_matrix,fIota_temp_matrix);
			fTemp_matrix_nen_press_x_ndof_se.MultATB(fShapeFluidGrad,fTemp_matrix_nsd_x_ndof_se);
			scale = 1/(fMaterial_Params[kg])*J*fRho_f*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H4_2_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetad_H4_3_matrix] will be formed */
			fTemp_matrix_nen_press_x_ndof_se.MultABC(fLambda_temp_matrix,fGravity_column_matrix,fPi_temp_row_matrix);
			scale = -1/(fMaterial_Params[kg])*J*fRho_f*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_ndof_se *= scale;
			/* accumulate */
			fK_thetad_H4_3_matrix += fTemp_matrix_nen_press_x_ndof_se;
			
			
			/* [fK_thetatheta_H4_matrix] will be formed */
			fTemp_matrix_nen_press_x_nen_press.MultABC(fLambda_temp_matrix,fGravity_column_matrix,fShapeFluid_row_matrix);
			scale = -1/(fMaterial_Params[kg]*fMaterial_Params[kKf])*
			    fRho_f*beta_delta_t2*scale_const;
			fTemp_matrix_nen_press_x_nen_press *= scale;
			/* accumulate */
			fK_thetatheta_H4_matrix += fTemp_matrix_nen_press_x_nen_press;
			
			

/*********************************************************************************************************************************************
 *********************************************************************************************************************************************
 *********************************************************************************************************************************************/
		
		
		
			/* implementing small strain deformation of an elastic media  for solid phase just for debugging */
			
			/* [fD_matrix] will be formed */
			Form_D_matrix();
		
			/* [fB_matrix] will be formed */
			Form_B_matrix();
			
			/* [fK_dd_BTDB_matrix] will be formed */
			fTemp_matrix_ndof_se_x_ndof_se.MultATBC(fB_matrix,fD_matrix,fB_matrix);
			scale = scale_const;
			fTemp_matrix_ndof_se_x_ndof_se *= scale;
			/* accumulate */
			fK_dd_BTDB_matrix += fTemp_matrix_ndof_se_x_ndof_se;
		


		
		
		
		

		
		
// defining F_1_T
			/* for debugging */
			const int ip = fShapes_displ->CurrIP()+1;
//			fs_mix_out	<< endl << "IP" << ip
//					<< setw(outputFileWidth) << ", shape function matrix for solid phase: " 
//					<< setw(outputFileWidth) << fShapeSolid;
			
//			fs_mix_out	<< endl << "terms from shape function matrix for solid phase: " 
//					<< setw(outputFileWidth) << fShapeSolid(0,0) 
//					<< setw(outputFileWidth) << fShapeSolid(0,3);	
		    }		
	    
		    /* saving eulerian effective strain for each IPs of the current element */
		    fEulerian_effective_strain_Elements_IPs.SetRow(e,fEulerian_effective_strain_IPs);
	    
	    
		    /* saving cauchy effective stress for each IPs of the current element */
		    fCauchy_effective_stress_Elements_IPs.SetRow(e,fCauchy_effective_stress_IPs);

                    /* saving physical pore water pressure for each IPs of the current element */
		    fPhysical_pore_water_pressure_Elements_IPs.SetRow(e,fPhysical_pore_water_pressure_IPs);
	    
	    
		    /* saving state variables for each IPs of the current element */
		    fState_variables_Elements_IPs.SetRow(e,fState_variables_IPs);

	    
		    /* {fFd_int_M_vector} will be formed */	    
		    fM_dd_matrix.Multx(u_dotdot_vec,fFd_int_M_vector);


		    /* {fFd_int_C_vector} will be formed */
		    fC_dd_matrix.Multx(u_dot_vec,fFd_int_C_vector);


		    /* {fFd_int} will be formed */
		    fFd_int = fFd_int_M_vector;
		    fFd_int += fFd_int_C_vector;
		    fFd_int += fFd_int_N1_vector;
		    fFd_int += fFd_int_N2_vector; 
		    fFd_int += fFd_int_G4_vector;
		    fFd_int *= -1;
	    
/*			fK_dd_BTDB_matrix.MultTx(u_vec,fTemp_vector_ndof_se);
			fFd_int = fTemp_vector_ndof_se;
			fFd_int *= -1.0; */


		    /* [fKdd] will be formed */
		    fKdd = fK_dd_G1_1_matrix;
		    fKdd += fK_dd_G1_2_matrix;
		    fKdd += fK_dd_G2_1_matrix;
		    fKdd += fK_dd_G2_2_matrix;
		    fKdd += fK_dd_G2_3_matrix;
		    fKdd += fK_dd_G2_4_matrix;
		    fKdd += fK_dd_G3_1_matrix;
		    fKdd += fK_dd_G3_2_matrix;
		    fKdd += fK_dd_G3_3_matrix;
		    fKdd += fK_dd_G3_4_matrix;
		    fKdd += fK_dd_G3_5_matrix; 
		    fKdd += fK_dd_G4_matrix; 



//		    fKdd = fK_dd_BTDB_matrix;			
	    
	    
		    /* [fKdtheta] will be formed */
		    fKdtheta = fK_dtheta_G1_matrix;	    
		    fKdtheta += fK_dtheta_G3_matrix;
		    fKdtheta += fK_dtheta_G4_matrix;
	    
	    
		    /* [fKthetad] will be formed */
		    fKthetad = fK_thetad_H1_1_matrix;	
		    fKthetad += fK_thetad_H1_2_matrix;	 
		    fKthetad += fK_thetad_H1_3_matrix;	 
		    fKthetad += fK_thetad_H1_4_matrix;
		    fKthetad += fK_thetad_H2_1_matrix;
		    fKthetad += fK_thetad_H2_2_matrix;
		    fKthetad += fK_thetad_H2_3_matrix;
		    fKthetad += fK_thetad_H2_4_matrix;
		    fKthetad += fK_thetad_H2_5_matrix;
		    fKthetad += fK_thetad_H3_1_matrix;
		    fKthetad += fK_thetad_H3_2_matrix;
		    fKthetad += fK_thetad_H3_3_matrix;
		    fKthetad += fK_thetad_H3_4_matrix;
		    fKthetad += fK_thetad_H4_1_matrix;
		    fKthetad += fK_thetad_H4_2_matrix;
		    fKthetad += fK_thetad_H4_3_matrix;


		    /* [fKthetatheta] will be formed */
		    fKthetatheta = fK_thetatheta_H1_matrix;
		    fKthetatheta += fK_thetatheta_H2_1_matrix;
		    fKthetatheta += fK_thetatheta_H2_2_matrix;
		    fKthetatheta += fK_thetatheta_H2_3_matrix;
		    fKthetatheta += fK_thetatheta_H3_1_matrix;
		    fKthetatheta += fK_thetatheta_H3_2_matrix;
		    fKthetatheta += fK_thetatheta_H4_matrix;	


		    /* {fFtheta_int_M_vector} will be formed */
		    fM_thetad_matrix.Multx(u_dotdot_vec,fFtheta_int_M_vector);


		    /* {fFtheta_int_C1_vector} will be formed */
		    fC_thetatheta_matrix.Multx(press_dot_vec,fFtheta_int_C1_vector);


		    /* {fFtheta_int_C2_vector} will be formed */
		    fC_thetad_matrix.Multx(u_dot_vec,fFtheta_int_C2_vector);


		    /* {fFtheta_int} will be formed */
		    fFtheta_int = fFtheta_int_M_vector;
		    fFtheta_int += fFtheta_int_C1_vector;
		    fFtheta_int += fFtheta_int_C2_vector; 
		    fFtheta_int += fFtheta_int_N1_vector;
		    fFtheta_int += fFtheta_int_N2_vector;
		    fFtheta_int += fFtheta_int_H4_vector;
		    fFtheta_int *= -1;
		    

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
}



/* form global shape function derivatives */
void FSSolidFluidMixT::SetGlobalShape(void)
{
    /* fetch (initial) coordinates */
    SetLocalX(fLocInitCoords);
	
    /* compute shape function derivatives */
    fShapes_displ->SetDerivatives_DN_DDN();
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
    // type of analysis
    list.AddParameter(kAnalysisType, "type_of_analysis_1consolidation_2dynamic");

    // type of initial condition
    list.AddParameter(kInitialConditionType, "initial_condition_1geostatic_2displ_vel_press");


	
    double shearMu, sLambda, sAlpha, Rho_sR0, Rho_fR0, Phi_s0, Phi_f0, bulkK, 
	newBeta, newGamma, conductivityK,gravity_g ,gravity_g1,gravity_g2,gravity_g3;
			

    // solid elasticity
    list.AddParameter(shearMu, "mu");
    list.AddParameter(sLambda, "lambda");
	
    // fluid elasticity
    list.AddParameter(bulkK, "Kf");

    // fluid hydraulic conductivity
    list.AddParameter(conductivityK, "K");

    // gravity
    list.AddParameter(gravity_g, "g");

    // gravity in each direction (depends on the coordinate system which we have chosen for the problem)
    list.AddParameter(gravity_g1, "g1");
    list.AddParameter(gravity_g2, "g2");
    list.AddParameter(gravity_g3, "g3");
	
    // viscoelasticity
    list.AddParameter(sAlpha, "alpha");
	
    // initial real mass density
    list.AddParameter(Rho_sR0, "rho_sR0");
    list.AddParameter(Rho_fR0, "rho_fR0");
	
    // initial volume fractions
    list.AddParameter(Phi_s0, "phi_s0");
    list.AddParameter(Phi_f0, "phi_f0");
	
    // Newmark time integration parameters
//    list.AddParameter(newBeta, "beta");
//    list.AddParameter(newGamma, "gamma");

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
    kAnalysisType = list.GetParameter("type_of_analysis_1consolidation_2dynamic");
    kInitialConditionType = list.GetParameter("initial_condition_1geostatic_2displ_vel_press");
  


	
    fGeometryCode_press = fGeometryCode_displ; 
    fNumIP_press = fNumIP_displ;
    fGeometryCodeSurf_press = fGeometryCodeSurf_displ;
    fNumIPSurf_press = fNumIPSurf_displ;
	
    iConstitutiveModelType = list.GetParameter("constitutive_mod_type");
	
    fMaterial_Params.Dimension ( kNUM_FMATERIAL_TERMS );
//    fIntegration_Params.Dimension ( kNUM_FINTEGRATE_TERMS );
	
    fMaterial_Params[kMu] = list.GetParameter("mu");
    fMaterial_Params[kLambda] = list.GetParameter("lambda");
    fMaterial_Params[kKf] = list.GetParameter("Kf");
    fMaterial_Params[kK] = list.GetParameter("K");
    fMaterial_Params[kg] = list.GetParameter("g");
    fMaterial_Params[kg1] = list.GetParameter("g1");
    fMaterial_Params[kg2] = list.GetParameter("g2");
    fMaterial_Params[kg3] = list.GetParameter("g3");
    fMaterial_Params[kAlpha] = list.GetParameter("alpha");
    fMaterial_Params[kRho_sR0] = list.GetParameter("rho_sR0");
    fMaterial_Params[kRho_fR0] = list.GetParameter("rho_fR0");
    fMaterial_Params[kPhi_s0] = list.GetParameter("phi_s0");
    fMaterial_Params[kPhi_f0] = list.GetParameter("phi_f0");
	
//    fIntegration_Params[kBeta] = list.GetParameter("beta");
//    fIntegration_Params[kGamma] = list.GetParameter("gamma");
	
    Echo_Input_Data();
	
    knum_d_state = 3; // evolving volume fractions for solid and fluid phases and Jacobian are 
                      // printed as internal state variables
    knum_i_state = 0; // int's needed per ip, state variables
	
    knumstrain = 6; // number of strain outputs
    knumstress = 7; // number of stress outputs + physical pore water pressure
	
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
    fShapes_displ = new ShapeFunctionT(fGeometryCode_displ, fNumIP_displ, fCurrCoords_displ,1 );
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
    u_dot.Dimension (n_en_displ, n_sd);
    u_dot_n.Dimension (n_en_displ, n_sd);
    u_dotdot_n.Dimension (n_en_displ, n_sd);
    u_dotdot.Dimension (n_en_displ, n_sd);
    u_n.Dimension (n_en_displ, n_sd);
    del_u.Dimension (n_en_displ, n_sd);
    n_en_displ_x_n_sd = n_en_displ*n_sd;
    del_u_vec.Dimension (n_en_displ_x_n_sd);
    u_vec.Dimension (n_en_displ_x_n_sd);
    u_dot_vec.Dimension (n_en_displ_x_n_sd);
    u_dotdot_vec.Dimension (n_en_displ_x_n_sd);
    //ElementSupport().RegisterCoordinates(fInitCoords_displ);
    fDispl->RegisterLocal(u);
    fDispl->RegisterLocal(u_n);



    if (fIntegrator->Order() == 1)
    {
	fDispl->RegisterLocal(u_dot);
	fDispl->RegisterLocal(u_dot_n);
    }

    if (fIntegrator->Order() == 2)
    {
	fDispl->RegisterLocal(u_dot);
	fDispl->RegisterLocal(u_dot_n);
	fDispl->RegisterLocal(u_dotdot);
	fDispl->RegisterLocal(u_dotdot_n);
    }


    /* set local arrays for pore pressure field */
    int dum=1;
    press.Dimension (n_en_press, dum);
    press_dot.Dimension (n_en_press, dum);
    press_dot_n.Dimension (n_en_press, dum);
    press_dotdot.Dimension (n_en_press, dum);
    press_dotdot_n.Dimension (n_en_press, dum);
    press_n.Dimension (n_en_press, dum);
    del_press.Dimension (n_en_press, dum);
    del_press_vec.Dimension (n_en_press);
    press_vec.Dimension (n_en_press);
    press_dot_vec.Dimension (n_en_press);
    press_dotdot_vec.Dimension (n_en_press);
    //ElementSupport().RegisterCoordinates(fInitCoords_press);

    fPress->RegisterLocal(press);
    fPress->RegisterLocal(press_n);

    if (fIntegrator->Order() == 1)
    {
	fPress->RegisterLocal(press_dot);
	fPress->RegisterLocal(press_dot_n);
    }

    if (fIntegrator->Order() == 2)
    {
	fPress->RegisterLocal(press_dot);
	fPress->RegisterLocal(press_dot_n);
	fPress->RegisterLocal(press_dotdot);
	fPress->RegisterLocal(press_dotdot_n);
    }

	
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
	
    /* workspace matricies */
    fShapeSolid.Dimension (n_sd, n_en_displ_x_n_sd);
    fShapeFluid.Dimension (n_en_press);
    n_sd_x_n_sd = n_sd*n_sd;
    fShapeSolidGrad_temp.Dimension (n_sd, n_en_displ);
    fShapeSolidGrad.Dimension (n_sd_x_n_sd, n_en_displ_x_n_sd);
    fShapeSolidGrad_t.Dimension (n_sd_x_n_sd, n_en_displ_x_n_sd);
    fShapeSolidGradGrad.Dimension (n_sd *2 , n_en_displ);
    fShapeFluidGrad.Dimension (n_sd, n_en_press);
    fDeformation_Gradient.Dimension (n_sd,n_sd);
    fGrad_disp_vector.Dimension (n_sd_x_n_sd);
    fDeformation_Gradient_Inverse.Dimension (n_sd,n_sd);
    fDeformation_Gradient_Transpose.Dimension (n_sd,n_sd);
    fDeformation_Gradient_Inverse_Transpose.Dimension (n_sd,n_sd);
    fDefGradInv_Grad_grad.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fDefGradInv_Grad_grad_Transpose.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fDefGradT_9x9_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fDefGradInv_vector.Dimension (n_sd_x_n_sd);
    fRight_Cauchy_Green_tensor.Dimension (n_sd,n_sd);
    fRight_Cauchy_Green_tensor_Inverse.Dimension (n_sd,n_sd);
    fLeft_Cauchy_Green_tensor.Dimension (n_sd,n_sd);
    fIdentity_matrix.Dimension (n_sd,n_sd);
    fEffective_Second_Piola_tensor_inviscid.Dimension (n_sd,n_sd);
    fEffective_Second_Piola_tensor_viscous.Dimension (n_sd,n_sd);
    fTemp_matrix_nsd_x_nsd.Dimension (n_sd,n_sd);
    fTemp_matrix_nen_press_x_nsd.Dimension (n_en_press,n_sd);
    fTemp_matrix_nen_press_x_nen_press.Dimension (n_en_press,n_en_press);
    fEffective_Kirchhoff_tensor_inviscid.Dimension (n_sd,n_sd);
    fEffective_Kirchhoff_tensor_viscous.Dimension (n_sd,n_sd);
    fEffective_Kirchhoff_vector_inviscid.Dimension (n_sd_x_n_sd);
    fEffective_Second_Piola_vector_viscous.Dimension (n_sd_x_n_sd);
    fIota_temp_matrix.Dimension (n_en_displ_x_n_sd,n_sd_x_n_sd);
    fVarpi_temp_matrix.Dimension (n_sd, n_en_displ_x_n_sd);
    fk_hydraulic_conductivity_matrix.Dimension (n_sd,n_sd);
    fK_hydraulic_conductivity_matrix.Dimension (n_sd,n_sd);
    fLambda_temp_matrix.Dimension (n_en_press,n_sd);
    fChi_temp_vector.Dimension (n_sd);
    fTemp_vector_9x1.Dimension (n_sd_x_n_sd);
    fFd_int_N1_vector.Dimension (n_en_displ_x_n_sd);
    fFd_int_N2_vector.Dimension (n_en_displ_x_n_sd); 
    fTemp_vector_ndof_se.Dimension (n_en_displ_x_n_sd); 
    fFtheta_int_N1_vector.Dimension (n_en_press); 
    fFtheta_int_N2_vector.Dimension (n_en_press); 
    fTemp_vector_nen_press.Dimension (n_en_press); 
    fIm_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fHbar_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fEll_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fPi_temp_transpose_vector.Dimension (n_en_displ_x_n_sd); 
    fPi_temp_row_matrix.Dimension (1,n_en_displ_x_n_sd); 
    fK_dd_G3_1_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G3_2_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G3_3_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G3_4_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G3_5_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dtheta_G3_matrix.Dimension (n_en_displ_x_n_sd,n_en_press);
    fI_ij_column_matrix.Dimension (n_sd_x_n_sd, 1);
    fShapeSolidGrad_t_Transpose.Dimension (n_en_displ_x_n_sd, n_sd_x_n_sd);
    fShapeFluid_row_matrix.Dimension (1,n_en_press); 
    fGrad_Omega_vector.Dimension (n_sd);
    fgrad_Omega_vector.Dimension (n_sd);
    fGrad_theta_vector.Dimension (n_sd);
    fGrad_phi_f_vector.Dimension (n_sd);
    fTemp_nsd_vector.Dimension (n_sd);
    fGrad_1_J_vector.Dimension (n_sd);
    fJmath_temp_matrix.Dimension (n_sd,n_sd_x_n_sd);
    fWp_temp_matrix.Dimension (n_sd,n_sd_x_n_sd);
    fK_thetad_H3_1_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetad_H3_2_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetad_H3_3_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetad_H3_4_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fChi_temp_column_matrix.Dimension (n_sd, 1);
    fTemp_matrix_nsd_x_1.Dimension (n_sd,1); 
    fTemp_matrix_nen_press_x_ndof_se.Dimension (n_en_press,n_en_displ_x_n_sd);
    fTemp_matrix_ndof_se_x_ndof_se.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fTemp_matrix_ndof_se_x_nen_press.Dimension (n_en_displ_x_n_sd,n_en_press);
    fK_thetatheta_H3_1_matrix.Dimension (n_en_press,n_en_press);
    fK_thetatheta_H3_2_matrix.Dimension (n_en_press,n_en_press);
    fc_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fC_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fIm_Prim_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fB_matrix.Dimension (6 ,n_en_displ_x_n_sd);
    fD_matrix.Dimension (6,6);
    fK_dd_BTDB_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fFd_int_smallstrain_vector.Dimension (n_en_displ_x_n_sd);
    fEulerian_effective_strain_tensor_current_IP.Dimension (n_sd,n_sd);
    fCauchy_effective_stress_tensor_current_IP.Dimension (n_sd,n_sd);
    fEulerian_effective_strain_IPs.Dimension (fNumIP_displ,6);
    fCauchy_effective_stress_IPs.Dimension (fNumIP_displ,6);
    fPhysical_pore_water_pressure_IPs.Dimension (fNumIP_displ,1);
    fState_variables_IPs.Dimension (fNumIP_displ,3);
    fTemp_six_values.Dimension (6);
    fEulerian_effective_strain_Elements_IPs.Dimension (NumElements(),fNumIP_displ*6);
    fCauchy_effective_stress_Elements_IPs.Dimension (NumElements(),fNumIP_displ*6);
    fPhysical_pore_water_pressure_Elements_IPs.Dimension (NumElements(),fNumIP_displ);
    fState_variables_Elements_IPs.Dimension (NumElements(),fNumIP_displ*3);
    fM_dd_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fUpsilon_temp_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fC_dd_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fGravity_vector.Dimension (n_sd);
    fFd_int_G4_vector.Dimension (n_en_displ_x_n_sd);
    fFd_int_M_vector.Dimension (n_en_displ_x_n_sd); 
    fFd_int_C_vector.Dimension (n_en_displ_x_n_sd); 
    fM_thetad_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fC_thetatheta_matrix.Dimension (n_en_press,n_en_press);
    fC_thetad_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fDefGradInv_column_matrix.Dimension (n_sd_x_n_sd,1);
    fDefGradInv_column_matrix_Transpose.Dimension (1,n_sd_x_n_sd);
    fFtheta_int_H4_vector.Dimension (n_en_press); 
    fFtheta_int_M_vector.Dimension (n_en_press); 
    fFtheta_int_C1_vector.Dimension (n_en_press);
    fFtheta_int_C2_vector.Dimension (n_en_press);
    fK_dd_G1_1_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G1_2_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    u_dotdot_column_matrix.Dimension (n_en_displ_x_n_sd,1);
    fK_dtheta_G1_matrix.Dimension (n_en_displ_x_n_sd,n_en_press);
    fGradv_vector.Dimension (n_sd_x_n_sd);
    fgradv_vector.Dimension (n_sd_x_n_sd);
    fXi_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fVarsigma_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fI_ijkl_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fK_dd_G2_1_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G2_2_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G2_3_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G2_4_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    u_dot_column_matrix.Dimension (n_en_displ_x_n_sd,1);
    fTemp_matrix1_ndof_se_x_ndof_se.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G4_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fGravity_column_matrix.Dimension (n_sd, 1);
    fK_dtheta_G4_matrix.Dimension (n_en_displ_x_n_sd,n_en_press);
    fAleph_temp_matrix.Dimension (n_sd,n_sd_x_n_sd);
    fK_thetad_H1_1_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetad_H1_2_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetad_H1_3_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetad_H1_4_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fTemp_matrix1_nen_press_x_ndof_se.Dimension (n_en_press,n_en_displ_x_n_sd);
    fTemp_matrix_nsd_x_ndof_se.Dimension (n_sd,n_en_displ_x_n_sd);
    fK_thetatheta_H1_matrix.Dimension (n_en_press,n_en_press);
    fTemp_matrix_nsd_x_nen_press.Dimension (n_sd,n_en_press);
    fK_thetad_H2_1_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetad_H2_2_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetad_H2_3_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetad_H2_4_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetad_H2_5_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    press_dot_column_matrix.Dimension (n_en_press,1);
    u_dot_column_matrix_Transpose.Dimension (1, n_en_displ_x_n_sd);
    fK_thetatheta_H2_1_matrix.Dimension (n_en_press,n_en_press);
    fK_thetatheta_H2_2_matrix.Dimension (n_en_press,n_en_press);
    fK_thetatheta_H2_3_matrix.Dimension (n_en_press,n_en_press);
    fTemp_matrix_nen_press_x_1.Dimension (n_en_press,1);
    fImath_temp_matrix.Dimension (n_sd,n_sd_x_n_sd);
    fK_thetad_H4_1_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetad_H4_2_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetad_H4_3_matrix.Dimension (n_en_press,n_en_displ_x_n_sd);
    fK_thetatheta_H4_matrix.Dimension (n_en_press,n_en_press);
    fPf_0_matrix.Dimension (NumElements(),fNumIP_displ);
    fP0_temp_value.Dimension (1);

    /* streams */
    ofstreamT& out = ElementSupport().Output();

    /* storage for integration point strain, stress, and ISVs*/
    fIPVariable.Dimension (n_el, fNumIP_press*(knumstrain+knumstress+knum_d_state));
    fIPVariable = 0.0;

    /* allocate storage for nodal forces */
    //fForces_at_Node.Dimension ( n_sd );
	
    /* extract natural boundary conditions */
    TakeNaturalBC(list);
	
    /* setup output file and format */
    outputPrecision = 10;
    outputFileWidth = outputPrecision + 8;
//    fs_mix_out.open("fs_mix.info");
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
	    const ParentDomainT& surf_shape = ShapeFunctionDispl().FacetShapeFunction(facet);
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

void FSSolidFluidMixT::Form_solid_shape_functions(const double* &shapes_displ_X)
{
    fShapeSolid = 0.0;
    for (int i=0; i<27; i++)
    {
	fShapeSolid(0,i*3) = shapes_displ_X[i];
	fShapeSolid(1,1+i*3) = shapes_displ_X[i];
	fShapeSolid(2,2+i*3) = shapes_displ_X[i];
    }
}

void FSSolidFluidMixT::Form_Gradient_of_solid_shape_functions(const dMatrixT &fShapeSolidGrad_temp)
{
    fShapeSolidGrad = 0.0;
    for(int i=0; i<27; i++)
    {
	fShapeSolidGrad(0,i*3) = fShapeSolidGrad_temp(0,i);
	fShapeSolidGrad(1,1+i*3) = fShapeSolidGrad_temp(0,i);
	fShapeSolidGrad(2,2+i*3) = fShapeSolidGrad_temp(0,i);

	fShapeSolidGrad(3,i*3) = fShapeSolidGrad_temp(1,i);
	fShapeSolidGrad(4,1+i*3) = fShapeSolidGrad_temp(1,i);
	fShapeSolidGrad(5,2+i*3) = fShapeSolidGrad_temp(1,i);

	fShapeSolidGrad(6,i*3) = fShapeSolidGrad_temp(2,i);
	fShapeSolidGrad(7,1+i*3) = fShapeSolidGrad_temp(2,i);
	fShapeSolidGrad(8,2+i*3) = fShapeSolidGrad_temp(2,i);
    }

}

void FSSolidFluidMixT::	Form_fluid_shape_functions(const double* &shapes_press_X)
{
    fShapeFluid = 0.0;
    for (int i=0; i<8; i++)
	fShapeFluid[i] = shapes_press_X[i];
}

void FSSolidFluidMixT::	Form_deformation_gradient_tensor(void)
{
    fShapeSolidGrad.Multx(u_vec,fGrad_disp_vector);
    fDeformation_Gradient(0,0) = fGrad_disp_vector[0]+1.0;
    fDeformation_Gradient(0,1) = fGrad_disp_vector[3]; 
    fDeformation_Gradient(0,2) = fGrad_disp_vector[6];
    fDeformation_Gradient(1,0) = fGrad_disp_vector[1];
    fDeformation_Gradient(1,1) = fGrad_disp_vector[4]+1.0;  
    fDeformation_Gradient(1,2) = fGrad_disp_vector[7];
    fDeformation_Gradient(2,0) = fGrad_disp_vector[2];
    fDeformation_Gradient(2,1) = fGrad_disp_vector[5];
    fDeformation_Gradient(2,2) = fGrad_disp_vector[8]+1.0; 
}
void FSSolidFluidMixT::Form_Grad_grad_transformation_matrix(void)
{
    fDefGradInv_Grad_grad = 0.0;
    fDefGradInv_Grad_grad(0,0) = fDeformation_Gradient_Inverse(0,0);
    fDefGradInv_Grad_grad(0,3) = fDeformation_Gradient_Inverse(0,1);
    fDefGradInv_Grad_grad(0,6) = fDeformation_Gradient_Inverse(0,2);
    fDefGradInv_Grad_grad(1,1) = fDeformation_Gradient_Inverse(0,0);
    fDefGradInv_Grad_grad(1,4) = fDeformation_Gradient_Inverse(0,1);
    fDefGradInv_Grad_grad(1,7) = fDeformation_Gradient_Inverse(0,2);
    fDefGradInv_Grad_grad(2,2) = fDeformation_Gradient_Inverse(0,0);
    fDefGradInv_Grad_grad(2,5) = fDeformation_Gradient_Inverse(0,1);
    fDefGradInv_Grad_grad(2,8) = fDeformation_Gradient_Inverse(0,2);

    fDefGradInv_Grad_grad(3,0) = fDeformation_Gradient_Inverse(1,0);
    fDefGradInv_Grad_grad(3,3) = fDeformation_Gradient_Inverse(1,1);
    fDefGradInv_Grad_grad(3,6) = fDeformation_Gradient_Inverse(1,2);
    fDefGradInv_Grad_grad(4,1) = fDeformation_Gradient_Inverse(1,0);
    fDefGradInv_Grad_grad(4,4) = fDeformation_Gradient_Inverse(1,1);
    fDefGradInv_Grad_grad(4,7) = fDeformation_Gradient_Inverse(1,2);
    fDefGradInv_Grad_grad(5,2) = fDeformation_Gradient_Inverse(1,0);
    fDefGradInv_Grad_grad(5,5) = fDeformation_Gradient_Inverse(1,1);
    fDefGradInv_Grad_grad(5,8) = fDeformation_Gradient_Inverse(1,2);

    fDefGradInv_Grad_grad(6,0) = fDeformation_Gradient_Inverse(2,0);
    fDefGradInv_Grad_grad(6,3) = fDeformation_Gradient_Inverse(2,1);
    fDefGradInv_Grad_grad(6,6) = fDeformation_Gradient_Inverse(2,2);
    fDefGradInv_Grad_grad(7,1) = fDeformation_Gradient_Inverse(2,0);
    fDefGradInv_Grad_grad(7,4) = fDeformation_Gradient_Inverse(2,1);
    fDefGradInv_Grad_grad(7,7) = fDeformation_Gradient_Inverse(2,2);
    fDefGradInv_Grad_grad(8,2) = fDeformation_Gradient_Inverse(2,0);
    fDefGradInv_Grad_grad(8,5) = fDeformation_Gradient_Inverse(2,1);
    fDefGradInv_Grad_grad(8,8) = fDeformation_Gradient_Inverse(2,2);
}

void FSSolidFluidMixT::Form_fDefGradT_9x9_matrix(void)
{
    fDefGradT_9x9_matrix = 0.0;
    fDefGradT_9x9_matrix(0,0) = fDeformation_Gradient(0,0);
    fDefGradT_9x9_matrix(0,1) = fDeformation_Gradient(1,0);
    fDefGradT_9x9_matrix(0,2) = fDeformation_Gradient(2,0);
    fDefGradT_9x9_matrix(1,3) = fDeformation_Gradient(0,0);
    fDefGradT_9x9_matrix(1,4) = fDeformation_Gradient(1,0);
    fDefGradT_9x9_matrix(1,5) = fDeformation_Gradient(2,0);
    fDefGradT_9x9_matrix(2,6) = fDeformation_Gradient(0,0);
    fDefGradT_9x9_matrix(2,7) = fDeformation_Gradient(1,0);
    fDefGradT_9x9_matrix(2,8) = fDeformation_Gradient(2,0);

    fDefGradT_9x9_matrix(3,0) = fDeformation_Gradient(0,1);
    fDefGradT_9x9_matrix(3,1) = fDeformation_Gradient(1,1);
    fDefGradT_9x9_matrix(3,2) = fDeformation_Gradient(2,1);
    fDefGradT_9x9_matrix(4,3) = fDeformation_Gradient(0,1);
    fDefGradT_9x9_matrix(4,4) = fDeformation_Gradient(1,1);
    fDefGradT_9x9_matrix(4,5) = fDeformation_Gradient(2,1);
    fDefGradT_9x9_matrix(5,6) = fDeformation_Gradient(0,1);
    fDefGradT_9x9_matrix(5,7) = fDeformation_Gradient(1,1);
    fDefGradT_9x9_matrix(5,8) = fDeformation_Gradient(2,1);

    fDefGradT_9x9_matrix(6,0) = fDeformation_Gradient(0,2);
    fDefGradT_9x9_matrix(6,1) = fDeformation_Gradient(1,2);
    fDefGradT_9x9_matrix(6,2) = fDeformation_Gradient(2,2);
    fDefGradT_9x9_matrix(7,3) = fDeformation_Gradient(0,2);
    fDefGradT_9x9_matrix(7,4) = fDeformation_Gradient(1,2);
    fDefGradT_9x9_matrix(7,5) = fDeformation_Gradient(2,2);
    fDefGradT_9x9_matrix(8,6) = fDeformation_Gradient(0,2);
    fDefGradT_9x9_matrix(8,7) = fDeformation_Gradient(1,2);
    fDefGradT_9x9_matrix(8,8) = fDeformation_Gradient(2,2);

}


void FSSolidFluidMixT::Form_deformation_gradient_inv_vector(void)
{
    fDefGradInv_vector[0] = fDeformation_Gradient_Inverse(0,0);
    fDefGradInv_vector[1] = fDeformation_Gradient_Inverse(0,1);
    fDefGradInv_vector[2] = fDeformation_Gradient_Inverse(0,2);
    fDefGradInv_vector[3] = fDeformation_Gradient_Inverse(1,0);
    fDefGradInv_vector[4] = fDeformation_Gradient_Inverse(1,1);
    fDefGradInv_vector[5] = fDeformation_Gradient_Inverse(1,2);
    fDefGradInv_vector[6] = fDeformation_Gradient_Inverse(2,0);
    fDefGradInv_vector[7] = fDeformation_Gradient_Inverse(2,1);
    fDefGradInv_vector[8] = fDeformation_Gradient_Inverse(2,2); 
    
}

void FSSolidFluidMixT::Form_effective_kirchhoff_stress_vector_inviscid()
{
    fEffective_Kirchhoff_vector_inviscid[0] = fEffective_Kirchhoff_tensor_inviscid(0,0);
    fEffective_Kirchhoff_vector_inviscid[1] = fEffective_Kirchhoff_tensor_inviscid(1,0);
    fEffective_Kirchhoff_vector_inviscid[2] = fEffective_Kirchhoff_tensor_inviscid(2,0);
    fEffective_Kirchhoff_vector_inviscid[3] = fEffective_Kirchhoff_tensor_inviscid(0,1);
    fEffective_Kirchhoff_vector_inviscid[4] = fEffective_Kirchhoff_tensor_inviscid(1,1);
    fEffective_Kirchhoff_vector_inviscid[5] = fEffective_Kirchhoff_tensor_inviscid(2,1);
    fEffective_Kirchhoff_vector_inviscid[6] = fEffective_Kirchhoff_tensor_inviscid(0,2);
    fEffective_Kirchhoff_vector_inviscid[7] = fEffective_Kirchhoff_tensor_inviscid(1,2);
    fEffective_Kirchhoff_vector_inviscid[8] = fEffective_Kirchhoff_tensor_inviscid(2,2);
}

void FSSolidFluidMixT::Form_fEffective_Second_Piola_tensor_viscous(const double& fAlpha)
{
    fShapeSolidGrad.Multx(u_dot_vec,fGradv_vector);
    fDefGradT_9x9_matrix.Multx(fGradv_vector,fTemp_vector_9x1);
    fC_matrix.Multx(fTemp_vector_9x1,fEffective_Second_Piola_vector_viscous);
    fEffective_Second_Piola_vector_viscous *= fAlpha;
    fEffective_Second_Piola_tensor_viscous(0,0) = fEffective_Second_Piola_vector_viscous[0];
    fEffective_Second_Piola_tensor_viscous(1,0) = fEffective_Second_Piola_vector_viscous[1];
    fEffective_Second_Piola_tensor_viscous(2,0) = fEffective_Second_Piola_vector_viscous[2];
    fEffective_Second_Piola_tensor_viscous(0,1) = fEffective_Second_Piola_vector_viscous[3];
    fEffective_Second_Piola_tensor_viscous(1,1) = fEffective_Second_Piola_vector_viscous[4];
    fEffective_Second_Piola_tensor_viscous(2,1) = fEffective_Second_Piola_vector_viscous[5];
    fEffective_Second_Piola_tensor_viscous(0,2) = fEffective_Second_Piola_vector_viscous[6];
    fEffective_Second_Piola_tensor_viscous(1,2) = fEffective_Second_Piola_vector_viscous[7];
    fEffective_Second_Piola_tensor_viscous(2,2) = fEffective_Second_Piola_vector_viscous[8];
}


void FSSolidFluidMixT::Form_Varpi_temp_matrix()
{
    double N_A_1I, N_A_2I, N_A_3I;
    int j, temp_j;
    j = 0 ;
    for (int A=1; A <= n_en_displ; A++)
    {
	temp_j = j;
	for (int i=0; i<3; i++)
	{
	    j = temp_j;
	    for (int n=1; n<=3; n++)
	    {
		switch (i+1)
		{
		case 1:
		{
		    N_A_1I = fShapeSolidGradGrad(0,A-1);
		    N_A_2I = fShapeSolidGradGrad(5,A-1);
		    N_A_3I = fShapeSolidGradGrad(4,A-1);
		}break;
		case 2:
		{
		    N_A_1I = fShapeSolidGradGrad(5,A-1);
		    N_A_2I = fShapeSolidGradGrad(1,A-1);
		    N_A_3I = fShapeSolidGradGrad(3,A-1);
		}break;
		case 3:
		{
		    N_A_1I = fShapeSolidGradGrad(4,A-1);
		    N_A_2I = fShapeSolidGradGrad(3,A-1);
		    N_A_3I = fShapeSolidGradGrad(2,A-1);
		}break;
		}
		fVarpi_temp_matrix(i,j) = N_A_1I * fDeformation_Gradient_Inverse(0, n-1)
		    + N_A_2I * fDeformation_Gradient_Inverse(1, n-1)
		    + N_A_3I * fDeformation_Gradient_Inverse(2, n-1);
		j += 1;
	    }
	}
    }  
}

void FSSolidFluidMixT::Form_Gradient_t_of_solid_shape_functions(const dMatrixT &fShapeSolidGrad_temp)
{
    fShapeSolidGrad_t = 0.0;
    for (int i=0; i<27; i++)
    {
	fShapeSolidGrad_t(0,i*3) = fShapeSolidGrad_temp(0,i);
	fShapeSolidGrad_t(1,i*3) = fShapeSolidGrad_temp(1,i);
	fShapeSolidGrad_t(2,i*3) = fShapeSolidGrad_temp(2,i);

	fShapeSolidGrad_t(3,1+i*3) = fShapeSolidGrad_temp(0,i);
	fShapeSolidGrad_t(4,1+i*3) = fShapeSolidGrad_temp(1,i);
	fShapeSolidGrad_t(5,1+i*3) = fShapeSolidGrad_temp(2,i);

	fShapeSolidGrad_t(6,2+i*3) = fShapeSolidGrad_temp(0,i);
	fShapeSolidGrad_t(7,2+i*3) = fShapeSolidGrad_temp(1,i);
	fShapeSolidGrad_t(8,2+i*3) = fShapeSolidGrad_temp(2,i);
    }
    
}

void FSSolidFluidMixT::Form_Im_temp_matrix()
{
    fIm_temp_matrix = 0.0;
    fIm_temp_matrix(0,0) = fEffective_Kirchhoff_tensor_inviscid(0,0);
    fIm_temp_matrix(1,0) = fEffective_Kirchhoff_tensor_inviscid(1,0);
    fIm_temp_matrix(2,0) = fEffective_Kirchhoff_tensor_inviscid(2,0);
    
    fIm_temp_matrix(3,1) = fEffective_Kirchhoff_tensor_inviscid(0,0);
    fIm_temp_matrix(4,1) = fEffective_Kirchhoff_tensor_inviscid(1,0);
    fIm_temp_matrix(5,1) = fEffective_Kirchhoff_tensor_inviscid(2,0);
    
    fIm_temp_matrix(6,2) = fEffective_Kirchhoff_tensor_inviscid(0,0);
    fIm_temp_matrix(7,2) = fEffective_Kirchhoff_tensor_inviscid(1,0);
    fIm_temp_matrix(8,2) = fEffective_Kirchhoff_tensor_inviscid(2,0);
    
    fIm_temp_matrix(0,3) = fEffective_Kirchhoff_tensor_inviscid(0,1);
    fIm_temp_matrix(1,3) = fEffective_Kirchhoff_tensor_inviscid(1,1);
    fIm_temp_matrix(2,3) = fEffective_Kirchhoff_tensor_inviscid(2,1);
    
    fIm_temp_matrix(3,4) = fEffective_Kirchhoff_tensor_inviscid(0,1);
    fIm_temp_matrix(4,4) = fEffective_Kirchhoff_tensor_inviscid(1,1);
    fIm_temp_matrix(5,4) = fEffective_Kirchhoff_tensor_inviscid(2,1);
    
    fIm_temp_matrix(6,5) = fEffective_Kirchhoff_tensor_inviscid(0,1);
    fIm_temp_matrix(7,5) = fEffective_Kirchhoff_tensor_inviscid(1,1);
    fIm_temp_matrix(8,5) = fEffective_Kirchhoff_tensor_inviscid(2,1);
    
    fIm_temp_matrix(0,6) = fEffective_Kirchhoff_tensor_inviscid(0,2);
    fIm_temp_matrix(1,6) = fEffective_Kirchhoff_tensor_inviscid(1,2);
    fIm_temp_matrix(2,6) = fEffective_Kirchhoff_tensor_inviscid(2,2);
    
    fIm_temp_matrix(3,7) = fEffective_Kirchhoff_tensor_inviscid(0,2);
    fIm_temp_matrix(4,7) = fEffective_Kirchhoff_tensor_inviscid(1,2);
    fIm_temp_matrix(5,7) = fEffective_Kirchhoff_tensor_inviscid(2,2);
    
    fIm_temp_matrix(6,8) = fEffective_Kirchhoff_tensor_inviscid(0,2);
    fIm_temp_matrix(7,8) = fEffective_Kirchhoff_tensor_inviscid(1,2);
    fIm_temp_matrix(8,8) = fEffective_Kirchhoff_tensor_inviscid(2,2);
}

void FSSolidFluidMixT::Form_Hbar_temp_matrix()
{
    fHbar_temp_matrix =0.0;
    fHbar_temp_matrix(0,0) = fLeft_Cauchy_Green_tensor(0,0);
    fHbar_temp_matrix(3,0) = fLeft_Cauchy_Green_tensor(0,1);
    fHbar_temp_matrix(6,0) = fLeft_Cauchy_Green_tensor(0,2);
    
    fHbar_temp_matrix(1,1) = fLeft_Cauchy_Green_tensor(0,0);
    fHbar_temp_matrix(4,1) = fLeft_Cauchy_Green_tensor(0,1);
    fHbar_temp_matrix(7,1) = fLeft_Cauchy_Green_tensor(0,2);
    
    fHbar_temp_matrix(2,2) = fLeft_Cauchy_Green_tensor(0,0);
    fHbar_temp_matrix(5,2) = fLeft_Cauchy_Green_tensor(0,1);
    fHbar_temp_matrix(8,2) = fLeft_Cauchy_Green_tensor(0,2);
    
    fHbar_temp_matrix(0,3) = fLeft_Cauchy_Green_tensor(1,0);
    fHbar_temp_matrix(3,3) = fLeft_Cauchy_Green_tensor(1,1);
    fHbar_temp_matrix(6,3) = fLeft_Cauchy_Green_tensor(1,2);
    
    fHbar_temp_matrix(1,4) = fLeft_Cauchy_Green_tensor(1,0);
    fHbar_temp_matrix(4,4) = fLeft_Cauchy_Green_tensor(1,1);
    fHbar_temp_matrix(7,4) = fLeft_Cauchy_Green_tensor(1,2);
    
    fHbar_temp_matrix(2,5) = fLeft_Cauchy_Green_tensor(1,0);
    fHbar_temp_matrix(5,5) = fLeft_Cauchy_Green_tensor(1,1);
    fHbar_temp_matrix(8,5) = fLeft_Cauchy_Green_tensor(1,2);
    
    fHbar_temp_matrix(0,6) = fLeft_Cauchy_Green_tensor(2,0);
    fHbar_temp_matrix(3,6) = fLeft_Cauchy_Green_tensor(2,1);
    fHbar_temp_matrix(6,6) = fLeft_Cauchy_Green_tensor(2,2);
    
    fHbar_temp_matrix(1,7) = fLeft_Cauchy_Green_tensor(2,0);
    fHbar_temp_matrix(4,7) = fLeft_Cauchy_Green_tensor(2,1);
    fHbar_temp_matrix(7,7) = fLeft_Cauchy_Green_tensor(2,2);
    
    fHbar_temp_matrix(2,8) = fLeft_Cauchy_Green_tensor(2,0);
    fHbar_temp_matrix(5,8) = fLeft_Cauchy_Green_tensor(2,1);
    fHbar_temp_matrix(8,8) = fLeft_Cauchy_Green_tensor(2,2);
    
}

void FSSolidFluidMixT::Form_Ell_temp_matrix()
{
    fEll_temp_matrix(0,0) = 0.0;
    fEll_temp_matrix(0,0) = fLeft_Cauchy_Green_tensor(0,0);
    fEll_temp_matrix(1,0) = fLeft_Cauchy_Green_tensor(1,0);
    fEll_temp_matrix(2,0) = fLeft_Cauchy_Green_tensor(2,0);
    
    fEll_temp_matrix(3,1) = fLeft_Cauchy_Green_tensor(0,0);
    fEll_temp_matrix(4,1) = fLeft_Cauchy_Green_tensor(1,0);
    fEll_temp_matrix(5,1) = fLeft_Cauchy_Green_tensor(2,0);
    
    fEll_temp_matrix(6,2) = fLeft_Cauchy_Green_tensor(0,0);
    fEll_temp_matrix(7,2) = fLeft_Cauchy_Green_tensor(1,0);
    fEll_temp_matrix(8,2) = fLeft_Cauchy_Green_tensor(2,0);
    
    fEll_temp_matrix(0,3) = fLeft_Cauchy_Green_tensor(0,1);
    fEll_temp_matrix(1,3) = fLeft_Cauchy_Green_tensor(1,1);
    fEll_temp_matrix(2,3) = fLeft_Cauchy_Green_tensor(2,1);

    fEll_temp_matrix(3,4) = fLeft_Cauchy_Green_tensor(0,1);
    fEll_temp_matrix(4,4) = fLeft_Cauchy_Green_tensor(1,1);
    fEll_temp_matrix(5,4) = fLeft_Cauchy_Green_tensor(2,1);
    
    fEll_temp_matrix(6,5) = fLeft_Cauchy_Green_tensor(0,1);
    fEll_temp_matrix(7,5) = fLeft_Cauchy_Green_tensor(1,1);
    fEll_temp_matrix(8,5) = fLeft_Cauchy_Green_tensor(2,1);
    
    fEll_temp_matrix(0,6) = fLeft_Cauchy_Green_tensor(0,2);
    fEll_temp_matrix(1,6) = fLeft_Cauchy_Green_tensor(1,2);
    fEll_temp_matrix(2,6) = fLeft_Cauchy_Green_tensor(2,2);
    
    fEll_temp_matrix(3,7) = fLeft_Cauchy_Green_tensor(0,2);
    fEll_temp_matrix(4,7) = fLeft_Cauchy_Green_tensor(1,2);
    fEll_temp_matrix(5,7) = fLeft_Cauchy_Green_tensor(2,2);
    
    fEll_temp_matrix(6,8) = fLeft_Cauchy_Green_tensor(0,2);
    fEll_temp_matrix(7,8) = fLeft_Cauchy_Green_tensor(1,2);
    fEll_temp_matrix(8,8) = fLeft_Cauchy_Green_tensor(2,2);
}

void FSSolidFluidMixT::Form_Jmath_temp_matrix(void)
{
    fJmath_temp_matrix = 0.0;
    int col =0;
    double sum;
    for (int i=0; i< 3; i++)
    {
	sum =0.0;
	for (int j=0; j< 3; j++)
	    sum += fk_hydraulic_conductivity_matrix(i,j)*fgrad_Omega_vector[j];
	for (int k=0; k<3; k++)
	{
	    fJmath_temp_matrix(k,col) = sum;
	    col += 1;
	}
	 
    }
}

void FSSolidFluidMixT::Form_Wp_temp_matrix(void)
{
    for (int i=0; i<3; i++)
    {
	for (int j=0; j<9; j++)
	{
	    switch (j)
	    {
	    case (0):
		fWp_temp_matrix(i,0) = fk_hydraulic_conductivity_matrix(i,0) * fgrad_Omega_vector[0]; 
		break;
	    case (1):
		fWp_temp_matrix(i,1) = fk_hydraulic_conductivity_matrix(i,0) * fgrad_Omega_vector[1]; 
		break;
	    case (2):
		fWp_temp_matrix(i,2) = fk_hydraulic_conductivity_matrix(i,0) * fgrad_Omega_vector[2]; 
		break;
	    case (3):
		fWp_temp_matrix(i,3) = fk_hydraulic_conductivity_matrix(i,1) * fgrad_Omega_vector[0]; 
		break;	    
	    case (4):
		fWp_temp_matrix(i,4) = fk_hydraulic_conductivity_matrix(i,1) * fgrad_Omega_vector[1]; 
		break;	    
	    case (5):
		fWp_temp_matrix(i,5) = fk_hydraulic_conductivity_matrix(i,1) * fgrad_Omega_vector[2]; 
		break;	    
	    case (6):
		fWp_temp_matrix(i,6) = fk_hydraulic_conductivity_matrix(i,2) * fgrad_Omega_vector[0]; 
		break;
	    case (7):
		fWp_temp_matrix(i,7) = fk_hydraulic_conductivity_matrix(i,2) * fgrad_Omega_vector[1]; 
		break;
	    case (8):
		fWp_temp_matrix(i,8) = fk_hydraulic_conductivity_matrix(i,2) * fgrad_Omega_vector[2]; 
		break;	    
	    }
	}
    }
}

void FSSolidFluidMixT::Form_C_matrix(const double& J_Prim)
{
    double C_IJKL;
    int row,col;
    for (int I=0; I<3; I++)
	for (int J=0; J<3; J++)
	    for (int K=0; K<3; K++)
		for (int L=0; L<3; L++)
		{
		    C_IJKL = fMaterial_Params[kLambda]*fRight_Cauchy_Green_tensor_Inverse(J,I)*fRight_Cauchy_Green_tensor_Inverse(L,K)+2*(fMaterial_Params[kMu]-fMaterial_Params[kLambda]*log(J_Prim))*1/2*(fRight_Cauchy_Green_tensor_Inverse(I,K)*fRight_Cauchy_Green_tensor_Inverse(J,L)+fRight_Cauchy_Green_tensor_Inverse(I,L)*fRight_Cauchy_Green_tensor_Inverse(J,K));

		    if (I==0 && J==0) 
			row=0;
		    else if (I==1 && J==0)
			row=1;
		    else if (I==2 && J==0)
			row=2;
		    else if (I==0 && J==1)
			row=3;
		    else if (I==1 && J==1)
			row=4;
		    else if (I==2 && J==1)
			row=5;
		    else if (I==0 && J==2)
			row=6;
		    else if (I==1 && J==2)
			row=7;
		    else 
			row=8;

		    if (K==0 && L==0) 
			col=0;
		    else if (K==1 && L==0)
			col=1;
		    else if (K==2 && L==0)
			col=2;
		    else if (K==0 && L==1)
			col=3;
		    else if (K==1 && L==1)
			col=4;
		    else if (K==2 && L==1)
			col=5;
		    else if (K==0 && L==2)
			col=6;
		    else if (K==1 && L==2)
			col=7;
		    else 
			col=8;

		    fC_matrix(row,col)= C_IJKL;
		    
		}
}

void FSSolidFluidMixT::Form_c_matrix()
{
    double c_ijkl;
    int row,col;
    for (int i=0; i<3; i++)
	for (int j=0; j<3; j++)
	    for (int k=0; k<3; k++)
		for (int l=0; l<3; l++)
		{
		    c_ijkl=0;		    
		    for (int I=0; I<3; I++)
			for (int J=0; J<3; J++)
			    for (int K=0; K<3; K++)
				for (int L=0; L<3; L++)
				{
				    if (I==0 && J==0) 
					row=0;
				    else if (I==1 && J==0)
					row=1;
				    else if (I==2 && J==0)
					row=2;
				    else if (I==0 && J==1)
					row=3;
				    else if (I==1 && J==1)
					row=4;
				    else if (I==2 && J==1)
					row=5;
				    else if (I==0 && J==2)
					row=6;
				    else if (I==1 && J==2)
					row=7;
				    else 
					     row=8;

				    if (K==0 && L==0) 
					col=0;
				    else if (K==1 && L==0)
					col=1;
				    else if (K==2 && L==0)
					col=2;
				    else if (K==0 && L==1)
					col=3;
				    else if (K==1 && L==1)
					col=4;
				    else if (K==2 && L==1)
					col=5;
				    else if (K==0 && L==2)
					col=6;
				    else if (K==1 && L==2)
					col=7;
				    else 
					     col=8;

				    c_ijkl += fDeformation_Gradient(i,I)*fDeformation_Gradient(j,J)*fDeformation_Gradient(k,K)*fDeformation_Gradient(l,L)*fC_matrix(row,col);
				}
		    if (i==0 && j==0) 
			row=0;
		    else if (i==1 && j==0)
			row=1;
		    else if (i==2 && j==0)
			row=2;
		    else if (i==0 && j==1)
			row=3;
		    else if (i==1 && j==1)
			row=4;
		    else if (i==2 && j==1)
			row=5;
		    else if (i==0 && j==2)
			row=6;
		    else if (i==1 && j==2)
			row=7;
		    else 
			     row=8;

		    if (k==0 && l==0) 
			col=0;
		    else if (k==1 && l==0)
			col=1;
		    else if (k==2 && l==0)
			col=2;
		    else if (k==0 && l==1)
			col=3;
		    else if (k==1 && l==1)
			col=4;
		    else if (k==2 && l==1)
			col=5;
		    else if (k==0 && l==2)
			col=6;
		    else if (k==1 && l==2)
			col=7;
		    else 
			     col=8;

		    fc_matrix(row,col)=c_ijkl;
		}

}

 
void FSSolidFluidMixT::Form_Im_Prim_temp_matrix()
{
    fIm_Prim_temp_matrix = 0.0;
    fIm_Prim_temp_matrix(0,0) = fEffective_Kirchhoff_tensor_inviscid(0,0);
    fIm_Prim_temp_matrix(3,0) = fEffective_Kirchhoff_tensor_inviscid(1,0);
    fIm_Prim_temp_matrix(6,0) = fEffective_Kirchhoff_tensor_inviscid(2,0);
    
    fIm_Prim_temp_matrix(1,1) = fEffective_Kirchhoff_tensor_inviscid(0,0);
    fIm_Prim_temp_matrix(4,1) = fEffective_Kirchhoff_tensor_inviscid(1,0);
    fIm_Prim_temp_matrix(7,1) = fEffective_Kirchhoff_tensor_inviscid(2,0);
    
    fIm_Prim_temp_matrix(2,2) = fEffective_Kirchhoff_tensor_inviscid(0,0);
    fIm_Prim_temp_matrix(5,2) = fEffective_Kirchhoff_tensor_inviscid(1,0);
    fIm_Prim_temp_matrix(8,2) = fEffective_Kirchhoff_tensor_inviscid(2,0);
    
    fIm_Prim_temp_matrix(0,3) = fEffective_Kirchhoff_tensor_inviscid(0,1);
    fIm_Prim_temp_matrix(3,3) = fEffective_Kirchhoff_tensor_inviscid(1,1);
    fIm_Prim_temp_matrix(6,3) = fEffective_Kirchhoff_tensor_inviscid(2,1);
    
    fIm_Prim_temp_matrix(1,4) = fEffective_Kirchhoff_tensor_inviscid(0,1);
    fIm_Prim_temp_matrix(4,4) = fEffective_Kirchhoff_tensor_inviscid(1,1);
    fIm_Prim_temp_matrix(7,4) = fEffective_Kirchhoff_tensor_inviscid(2,1);
    
    fIm_Prim_temp_matrix(2,5) = fEffective_Kirchhoff_tensor_inviscid(0,1);
    fIm_Prim_temp_matrix(5,5) = fEffective_Kirchhoff_tensor_inviscid(1,1);
    fIm_Prim_temp_matrix(8,5) = fEffective_Kirchhoff_tensor_inviscid(2,1);
    
    fIm_Prim_temp_matrix(0,6) = fEffective_Kirchhoff_tensor_inviscid(0,2);
    fIm_Prim_temp_matrix(3,6) = fEffective_Kirchhoff_tensor_inviscid(1,2);
    fIm_Prim_temp_matrix(6,6) = fEffective_Kirchhoff_tensor_inviscid(2,2);
    
    fIm_Prim_temp_matrix(1,7) = fEffective_Kirchhoff_tensor_inviscid(0,2);
    fIm_Prim_temp_matrix(4,7) = fEffective_Kirchhoff_tensor_inviscid(1,2);
    fIm_Prim_temp_matrix(7,7) = fEffective_Kirchhoff_tensor_inviscid(2,2);
    
    fIm_Prim_temp_matrix(2,8) = fEffective_Kirchhoff_tensor_inviscid(0,2);
    fIm_Prim_temp_matrix(5,8) = fEffective_Kirchhoff_tensor_inviscid(1,2);
    fIm_Prim_temp_matrix(8,8) = fEffective_Kirchhoff_tensor_inviscid(2,2);
}

void FSSolidFluidMixT::Form_D_matrix(void)
{

    fD_matrix = 0.0;
    fD_matrix(0,0) = 2*fMaterial_Params[kMu]+ fMaterial_Params[kLambda];
    fD_matrix(1,1) = 2*fMaterial_Params[kMu]+ fMaterial_Params[kLambda];
    fD_matrix(2,2) = 2*fMaterial_Params[kMu]+ fMaterial_Params[kLambda];
    fD_matrix(3,3) = fMaterial_Params[kMu];
    fD_matrix(4,4) = fMaterial_Params[kMu];
    fD_matrix(5,5) = fMaterial_Params[kMu];
    fD_matrix(0,1) = fMaterial_Params[kLambda];
    fD_matrix(1,0) = fMaterial_Params[kLambda];
    fD_matrix(0,2) = fMaterial_Params[kLambda];
    fD_matrix(2,0) = fMaterial_Params[kLambda];
    fD_matrix(1,2) = fMaterial_Params[kLambda];
    fD_matrix(2,1) = fMaterial_Params[kLambda];

}

void FSSolidFluidMixT::Form_B_matrix(void)
{

    fB_matrix = 0.0;
    for(int i=0; i<27; i++)
    {
	fB_matrix(0,i*3)=fShapeSolidGrad_temp(0,i);
	fB_matrix(1,i*3+1)=fShapeSolidGrad_temp(1,i);
	fB_matrix(2,i*3+2)=fShapeSolidGrad_temp(2,i);
	fB_matrix(3,i*3+1)=fShapeSolidGrad_temp(2,i);
	fB_matrix(3,i*3+2)=fShapeSolidGrad_temp(1,i);
	fB_matrix(4,i*3)=fShapeSolidGrad_temp(2,i);
	fB_matrix(4,i*3+2)=fShapeSolidGrad_temp(0,i);
	fB_matrix(5,i*3)=fShapeSolidGrad_temp(1,i);
	fB_matrix(5,i*3+1)=fShapeSolidGrad_temp(0,i);
    }
}


void FSSolidFluidMixT::Extract_six_values_from_symmetric_tensor(const dMatrixT &fTensor,dArrayT& fTemp_six_values)
{
    fTemp_six_values[0]=fTensor(0,0);
    fTemp_six_values[1]=fTensor(1,1);
    fTemp_six_values[2]=fTensor(2,2);
    fTemp_six_values[3]=fTensor(1,2);
    fTemp_six_values[4]=fTensor(2,0);
    fTemp_six_values[5]=fTensor(0,1);
}

void FSSolidFluidMixT::Put_values_In_dArrayT_vector(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT)
{

	fArrayT[0]=f2DArrayT(e,IP*6+0);
	fArrayT[1]=f2DArrayT(e,IP*6+1);
	fArrayT[2]=f2DArrayT(e,IP*6+2);
	fArrayT[3]=f2DArrayT(e,IP*6+3);
	fArrayT[4]=f2DArrayT(e,IP*6+4);
	fArrayT[5]=f2DArrayT(e,IP*6+5);

}

void FSSolidFluidMixT::Form_gradv_vector(void)
{
    fShapeSolidGrad.Multx(u_dot_vec,fGradv_vector);
    fDefGradInv_Grad_grad.MultTx(fGradv_vector,fgradv_vector);
}

void FSSolidFluidMixT::Form_Xi_temp_matrix(void)
{
    fXi_temp_matrix = 0.0;

    double c11kl_vkl, c21kl_vkl, c31kl_vkl;
    double c12kl_vkl, c22kl_vkl, c32kl_vkl;
    double c13kl_vkl, c23kl_vkl, c33kl_vkl;
    int col,index;
    c11kl_vkl=c21kl_vkl=c31kl_vkl=0.0;
    c12kl_vkl=c22kl_vkl=c32kl_vkl=0.0;
    c13kl_vkl=c23kl_vkl=c33kl_vkl=0.0;

    for (int k=0; k<3; k++)
	for (int l=0; l<3; l++)
	{
	    if (k==0 && l==0) 
		col=index=0;
	    else if (k==1 && l==0)
		col=index=1;
	    else if (k==2 && l==0)
		col=index=2;
	    else if (k==0 && l==1)
		col=index=3;
	    else if (k==1 && l==1)
		col=index=4;
	    else if (k==2 && l==1)
		col=index=5;
	    else if (k==0 && l==2)
		col=index=6;
	    else if (k==1 && l==2)
		col=index=7;
	    else 
		col=index=8;

	    c11kl_vkl+=fc_matrix(0,col)*fgradv_vector[index];
	    c21kl_vkl+=fc_matrix(1,col)*fgradv_vector[index];
	    c31kl_vkl+=fc_matrix(2,col)*fgradv_vector[index];
	    c12kl_vkl+=fc_matrix(3,col)*fgradv_vector[index];
	    c22kl_vkl+=fc_matrix(4,col)*fgradv_vector[index];
	    c32kl_vkl+=fc_matrix(5,col)*fgradv_vector[index];
	    c13kl_vkl+=fc_matrix(6,col)*fgradv_vector[index];
	    c23kl_vkl+=fc_matrix(7,col)*fgradv_vector[index];
	    c33kl_vkl+=fc_matrix(8,col)*fgradv_vector[index];
	}
    for (int i=0; i<3; i++)
    {
	fXi_temp_matrix(i*3,i)=c11kl_vkl;
	fXi_temp_matrix(i*3+1,i)=c21kl_vkl;
	fXi_temp_matrix(i*3+2,i)=c31kl_vkl;

	fXi_temp_matrix(i*3,i+3)=c12kl_vkl;
	fXi_temp_matrix(i*3+1,i+3)=c22kl_vkl;
	fXi_temp_matrix(i*3+2,i+3)=c32kl_vkl;

	fXi_temp_matrix(i*3,i+6)=c13kl_vkl;
	fXi_temp_matrix(i*3+1,i+6)=c23kl_vkl;
	fXi_temp_matrix(i*3+2,i+6)=c33kl_vkl;

    }
	 
}

void FSSolidFluidMixT::Form_Varsigma_temp_matrix(void)
{
    fVarsigma_temp_matrix = 0.0;
    int col,index;
    for (int i=0; i<9; i++)
    {
	for (int k=0; k<3; k++)
	{
	    if (k==0) 
	    {		
		col=0;
		index=0;
	    }
	    else if (k==1) 
	    {
		col=3;
		index=1;
	    }
	    if (k==2) 
	    {
		col=6;
		index=2;
	    }

	    fVarsigma_temp_matrix(i,0) += fc_matrix(i,col)*fgradv_vector[index];
	}

	for (int k=0; k<3; k++)
	{
	    if (k==0) 
	    {		
		col=0;
		index=3;
	    }
	    else if (k==1) 
	    {
		col=3;
		index=4;
	    }
	    if (k==2) 
	    {
		col=6;
		index=5;
	    }

	    fVarsigma_temp_matrix(i,1) += fc_matrix(i,col)*fgradv_vector[index];
	}

	for (int k=0; k<3; k++)
	{
	    if (k==0) 
	    {		
		col=0;
		index=6;
	    }
	    else if (k==1) 
	    {
		col=3;
		index=7;
	    }
	    if (k==2) 
	    {
		col=6;
		index=8;
	    }

	    fVarsigma_temp_matrix(i,2) += fc_matrix(i,col)*fgradv_vector[index];
	}

	for (int k=0; k<3; k++)
	{
	    if (k==0) 
	    {
		col=1;
		index=0;
	    }
	    else if (k==1) 
	    {
		col=4;
		index=1;
	    }
	    if (k==2) 
	    {
		col=7;
		index=2;
	    }

	    fVarsigma_temp_matrix(i,3) += fc_matrix(i,col)*fgradv_vector[index];
	}

	for (int k=0; k<3; k++)
	{
	    if (k==0) 
	    {		
		col=1;
		index=3;
	    }
	    else if (k==1) 
	    {
		col=4;
		index=4;
	    }
	    if (k==2) 
	    {
		col=7;
		index=5;
	    }

	    fVarsigma_temp_matrix(i,4) += fc_matrix(i,col)*fgradv_vector[index];
	}

	for (int k=0; k<3; k++)
	{
	    if (k==0) 
	    {		
		col=1;
		index=6;
	    }
	    else if (k==1) 
	    {
		col=4;
		index=7;
	    }
	    if (k==2) 
	    {
		col=7;
		index=8;
	    }

	    fVarsigma_temp_matrix(i,5) += fc_matrix(i,col)*fgradv_vector[index];
	}

	for (int k=0; k<3; k++)
	{
	    if (k==0) 
	    {
		col=2;
		index=0;
	    }
	    else if (k==1) 
	    {
		col=5;
		index=1;
	    }
	    if (k==2) 
	    {
		col=8;
		index=2;
	    }

	    fVarsigma_temp_matrix(i,6) += fc_matrix(i,col)*fgradv_vector[index];
	}

	for (int k=0; k<3; k++)
	{
	    if (k==0) 
	    {		
		col=2;
		index=3;
	    }
	    else if (k==1) 
	    {
		col=5;
		index=4;
	    }
	    if (k==2) 
	    {
		col=8;
		index=5;
	    }

	    fVarsigma_temp_matrix(i,7) += fc_matrix(i,col)*fgradv_vector[index];
	}

	for (int k=0; k<3; k++)
	{
	    if (k==0) 
	    {		
		col=2;
		index=6;
	    }
	    else if (k==1) 
	    {
		col=5;
		index=7;
	    }
	    if (k==2) 
	    {
		col=8;
		index=8;
	    }

	    fVarsigma_temp_matrix(i,8) += fc_matrix(i,col)*fgradv_vector[index];
	}
    }

}

void FSSolidFluidMixT::Form_I_ijkl_matrix(void)
{
    double delta_ik,delta_jl,delta_il,delta_jk;
    int row,col;
    for (int i=0; i<3; i++)
	for (int j=0; j<3; j++)
	    for (int k=0; k<3; k++)
		for (int l=0; l<3; l++)
		{
		    if (i==k)
			delta_ik =1;
		    else
			delta_ik=0;

		    if (j==l)
			delta_jl=1;
		    else
			delta_jl=0;

		    if (i==l)
			delta_il=1;
		    else
			delta_il=0;

		    if (j==k)
			delta_jk=1;
		    else
			delta_jk=0;

		    if (i==0 && j==0) 
			row=0;
		    else if (i==1 && j==0)
			row=1;
		    else if (i==2 && j==0)
			row=2;
		    else if (i==0 && j==1)
			row=3;
		    else if (i==1 && j==1)
			row=4;
		    else if (i==2 && j==1)
			row=5;
		    else if (i==0 && j==2)
			row=6;
		    else if (i==1 && j==2)
			row=7;
		    else 
			row=8;

		    if (k==0 && l==0) 
			col=0;
		    else if (k==1 && l==0)
			col=1;
		    else if (k==2 && l==0)
			col=2;
		    else if (k==0 && l==1)
			col=3;
		    else if (k==1 && l==1)
			col=4;
		    else if (k==2 && l==1)
			col=5;
		    else if (k==0 && l==2)
			col=6;
		    else if (k==1 && l==2)
			col=7;
		    else 
			col=8;

		    fI_ijkl_matrix(row,col)= 0.5*(delta_ik*delta_jl+delta_il*delta_jk);
		    
		}
}

void FSSolidFluidMixT::Form_Aleph_temp_matrix(const int& IP)
{
    fAleph_temp_matrix = 0.0;
    int col =0;
    double sum;
    for (int i=0; i< 3; i++)
    {
	sum =0.0;
	for (int j=0; j< 3; j++)
	    sum += fk_hydraulic_conductivity_matrix(i,j)*u_dotdot_vec[IP*3+j];
	for (int k=0; k<3; k++)
	{
	    fAleph_temp_matrix(k,col) = sum;
	    col += 1;
	}
	 
    }
}

void FSSolidFluidMixT::Form_Imath_temp_matrix(void)
{
    fImath_temp_matrix = 0.0;
    int col =0;
    double sum;
    for (int i=0; i< 3; i++)
    {
	sum =0.0;
	for (int j=0; j< 3; j++)
	    sum += fk_hydraulic_conductivity_matrix(i,j)*fGravity_vector[j];
	for (int k=0; k<3; k++)
	{
	    fImath_temp_matrix(k,col) = sum;
	    col += 1;
	}
	 
    }
}


void FSSolidFluidMixT::Compute_norm_of_array(double& norm,const LocalArrayT& B)
{
    int index = 0;
    double sum = 0;
    for (int i=0; i<n_en_displ; i++)
    {
	for (int j=0; j<n_sd; j++)
	    sum = sum + u_dotdot(i,j)*u_dotdot(i,j);
    }
    norm = sqrt(sum);
}
