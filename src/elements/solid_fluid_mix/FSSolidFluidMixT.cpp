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

    //-- Hydraulic Conductivity
    cout << "fMaterial_Params[kK] "  				<< fMaterial_Params[kK] 	<< endl;

    //-- gravity
    cout << "fMaterial_Params[kg] "  				<< fMaterial_Params[kg] 	<< endl;

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
	
    fs_mix_out	<< endl 
		<< setw(outputFileWidth) << "time_step"
		<< endl;
    step_number = ElementSupport().StepNumber();
    fs_mix_out	<< setw(outputFileWidth) << step_number
		<< endl;
    fs_mix_out	<< endl << "**********************************************************************************************";
    fs_mix_out	<< endl << "**********************************************************************************************" << endl;
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
		e = CurrElementNumber();

		const iArrayT& nodes_displ = fElementCards_displ[e].NodesU();
		const iArrayT& nodes_press = fElementCards_press[e].NodesU();

		u.SetLocal(nodes_displ);
		u_n.SetLocal(nodes_displ);
		press.SetLocal(nodes_press);
		press_n.SetLocal(nodes_press);
                /* applying solid boundary condition */
                /* predefined displacement */
		    u(2,1)=-0.1;
		    u(3,1)=-0.1;
		    u(6,1)=-0.1;
		    u(7,1)=-0.1;
		    u(10,1)=-0.1;
		    u(14,1)=-0.1;
		    u(18,1)=-0.1;
		    u(19,1)=-0.1;
		    u(23,1)=-0.1; 
                /* fixed dof */
		    u(0,2)=0;
		    u(1,2)=0;
		    u(2,2)=0;
		    u(3,2)=0;
		    u(8,2)=0;
		    u(9,2)=0;
		    u(10,2)=0;
		    u(11,2)=0;
		    u(20,2)=0;
		    u(0,0)=0;
		    u(3,0)=0;
		    u(7,0)=0;
		    u(4,0)=0;
		    u(11,0)=0;
		    u(19,0)=0;
		    u(15,0)=0;
		    u(16,0)=0;
		    u(24,0)=0;
		    u(0,1)=0;
		    u(1,1)=0;
		    u(4,1)=0;
		    u(5,1)=0;
		    u(8,1)=0;
		    u(16,1)=0;
		    u(12,1)=0;
		    u(17,1)=0;
		    u(22,1)=0;


		   
                /* applying fluid boundary condition */
		    press(0,0)=0;
		    press(1,0)=0;
		    press(2,0)=0;
		    press(3,0)=0;
		    press(4,0)=0;
		    press(5,0)=0;
		    press(6,0)=0;
		    press(7,0)=0;


               /* print solid displacement */
		fs_mix_out	<<"nodal solid displacement"<< endl ;
		fs_mix_out	<< u << endl ;

                /* print fluid displacement */
		fs_mix_out	<<"nodal fluid pressure"<< endl ;
		fs_mix_out	<< press << endl ;

                /* populate solid displacement in a vector */
		int index = 0;
		for (int i=0; i<n_en_displ; i++)
		{
		    for (int j=0; j<n_sd; j++)
		    {
			u_vec[index] = u(i,j);
			index += 1;
		    }
		}

                /* populate fluid displacement in a vector */
		for (int i=0; i<n_en_press; i++)
		    press_vec[i] = press(i,0);

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
				//out_variable=??;
		    } 
		}
		else 
		{ //-- Still Iterating
		    /* residual and tangent for displacements */
		    const double* Det    = fShapes_displ->IPDets();
		    const double* Weight = fShapes_displ->IPWeights();
		    fShapes_displ->TopIP();
		    fShapes_press->TopIP();
		    int IntegrationPointNumber =0;
		    while (fShapes_displ->NextIP() && fShapes_press->NextIP())
		    {
				IntegrationPointNumber++ ;
				dArrayT SolidIPCoordinate(n_sd),FluidIPCoordinate(n_sd);
				fShapes_displ->IPCoords(SolidIPCoordinate);
				fShapes_press->IPCoords(FluidIPCoordinate);
				fs_mix_out	<<"Solid Integration Point "<< IntegrationPointNumber <<endl ;
				fs_mix_out	<<SolidIPCoordinate<< endl ;
				fs_mix_out	<<"Solid Integration Point"<< endl ;
				fs_mix_out	<<"Fluid Integration Point "<< IntegrationPointNumber <<endl ;
				fs_mix_out	<<FluidIPCoordinate<< endl ;
				fs_mix_out	<<"Fluid Integration Point"<< endl ;
				fs_mix_out	<<"Weight & Jacobian determinent"<< endl ;
				fs_mix_out	<<"Weight = "<<(*Weight)<< endl ;
				fs_mix_out	<<"Jacobian determinent(for change of coordinate from global to natural) = "<<(*Det)<< endl ;
				fs_mix_out	<<"Weight & Jacobian determinent"<< endl ; 
						
				const double* shapes_displ_X = fShapes_displ->IPShapeX();
                                /* [fShapeSolid]will be formed */
				Form_solid_shape_functions(shapes_displ_X);
				fs_mix_out	<<"fShapeSolid"<< endl ;
				fs_mix_out	<<fShapeSolid<< endl ;
				fs_mix_out	<<"fShapeSolid"<< endl ;

				fShapes_displ->GradNa(fShapeSolidGrad_temp);
				/* [fShapeSolidGrad] will be formed */
				Form_Gradient_of_solid_shape_functions(fShapeSolidGrad_temp);
				fs_mix_out	<<"fShapeSolidGrad"<< endl ;
				fs_mix_out	<<fShapeSolidGrad<< endl ;
				fs_mix_out	<<"fShapeSolidGrad"<< endl ;


				/* [fShapeSolidGrad_t] and [fShapeSolidGrad_t_Transpose] will be formed */
				Form_Gradient_t_of_solid_shape_functions(fShapeSolidGrad_temp);
				fShapeSolidGrad_t_Transpose.Transpose(fShapeSolidGrad_t);
				fs_mix_out	<<"fShapeSolidGrad_t"<< endl ;
				fs_mix_out	<<fShapeSolidGrad_t<< endl ;
				fs_mix_out	<<"fShapeSolidGrad_t"<< endl ;
				fs_mix_out	<<"fShapeSolidGrad_t_Transpose"<< endl ;
				fs_mix_out	<<fShapeSolidGrad_t_Transpose<< endl ;
				fs_mix_out	<<"fShapeSolidGrad_t_Transpose"<< endl ;


				const double* shapes_press_X = fShapes_press->IPShapeX();
				/* {fShapeFluid} will be formed */
				Form_fluid_shape_functions(shapes_press_X);
				fs_mix_out	<<"fShapeFluid"<< endl ;
				fs_mix_out	<<fShapeFluid<< endl ;
				fs_mix_out	<<"fShapeFluid"<< endl ;


                                /* [fShapeFluid_row_matrix] will be formed */				
				for (int i=0; i<n_en_press ; i++)
				    fShapeFluid_row_matrix(0,i) = fShapeFluid[i];
				fs_mix_out	<<"fShapeFluid_row_matrix"<< endl ;
				fs_mix_out	<<fShapeFluid_row_matrix<< endl ;
				fs_mix_out	<<"fShapeFluid_row_matrix"<< endl ;


				/* [fShapeFluidGrad] will be formed */
				fShapes_press->GradNa(fShapeFluidGrad);
				fs_mix_out	<<"fShapeFluidGrad"<< endl ;
				fs_mix_out	<<fShapeFluidGrad<< endl ;
				fs_mix_out	<<"fShapeFluidGrad"<< endl ;


				/* [fDeformation_Gradient] will be formed */
				Form_deformation_gradient_tensor();
				fs_mix_out	<<"fDeformation_Gradient"<< endl ;
				fs_mix_out	<<fDeformation_Gradient<< endl ;
				fs_mix_out	<<"fDeformation_Gradient"<< endl ;


				/* [fIdentity_matrix] will be formed */
				fIdentity_matrix = 0.0;			
				for (int i=0; i<n_sd ; i++)
					 fIdentity_matrix(i,i) =1.0;
				fs_mix_out	<<"fIdentity_matrix"<< endl ;
				fs_mix_out	<<fIdentity_matrix<< endl ;
				fs_mix_out	<<"fIdentity_matrix"<< endl ;


				/* [fDeformation_Gradient_Inverse] and [fDeformation_Gradient_Transpose] and [fDeformation_Gradient_Inverse_Transpose] will be formed */
				if (fDeformation_Gradient.Det()==0)
				    fDeformation_Gradient = fIdentity_matrix; 
				fDeformation_Gradient_Inverse.Inverse(fDeformation_Gradient);
				fDeformation_Gradient_Inverse_Transpose.Transpose(fDeformation_Gradient_Inverse);
				fDeformation_Gradient_Transpose.Transpose(fDeformation_Gradient);
				fs_mix_out	<<"fDeformation_Gradient_Inverse"<< endl ;
				fs_mix_out	<<fDeformation_Gradient_Inverse<< endl ;
				fs_mix_out	<<"fDeformation_Gradient_Inverse"<< endl ;
				fs_mix_out	<<"fDeformation_Gradient_Transpose"<< endl ;
				fs_mix_out	<<fDeformation_Gradient_Transpose<< endl ;
				fs_mix_out	<<"fDeformation_Gradient_Transpose"<< endl ;
				fs_mix_out	<<"fDeformation_Gradient_Inverse_Transpose"<< endl ;
				fs_mix_out	<<fDeformation_Gradient_Inverse_Transpose<< endl ;
				fs_mix_out	<<"fDeformation_Gradient_Inverse_Transpose"<< endl ;




				/* {fDefGradInv_Vector} will be formed */
				Form_deformation_gradient_inv_vector();
				fs_mix_out	<<"fDefGradInv_Vector"<< endl ;
				fs_mix_out	<<fDefGradInv_Vector<< endl ;
				fs_mix_out	<<"fDefGradInv_Vector"<< endl ;


                                /* [fDefGradInv_GRAD_grad] will be formed */
				Form_GRAD_grad_transformation_matrix();
				fs_mix_out	<<"fDefGradInv_GRAD_grad"<< endl ;
				fs_mix_out	<<fDefGradInv_GRAD_grad<< endl ;
				fs_mix_out	<<"fDefGradInv_GRAD_grad"<< endl ;


                                /* [fDefGradInv_GRAD_grad_Transpose] will be formed */
				fDefGradInv_GRAD_grad_Transpose.Transpose(fDefGradInv_GRAD_grad);
				fs_mix_out	<<"fDefGradInv_GRAD_grad_Transpose"<< endl ;
				fs_mix_out	<<fDefGradInv_GRAD_grad_Transpose<< endl ;
				fs_mix_out	<<"fDefGradInv_GRAD_grad_Transpose"<< endl ;


                                /* calculating Jacobian */
				double J = fDeformation_Gradient.Det();
				fs_mix_out	<<"Jacobian"<< endl ;
				fs_mix_out	<<J<< endl ;
				fs_mix_out	<<"Jacobian"<< endl ;


				/* [fRight_Cauchy_Green_tensor] will be formed */
				fRight_Cauchy_Green_tensor.MultATB(fDeformation_Gradient, fDeformation_Gradient);
				fs_mix_out	<<"fRight_Cauchy_Green_tensor"<< endl ;
				fs_mix_out	<<fRight_Cauchy_Green_tensor<< endl ;
				fs_mix_out	<<"fRight_Cauchy_Green_tensor"<< endl ;


				/* [fRight_Cauchy_Green_tensor_Inverse] will be formed */
				if (fRight_Cauchy_Green_tensor.Det()==0)
				    fRight_Cauchy_Green_tensor = fIdentity_matrix;
				fRight_Cauchy_Green_tensor_Inverse.Inverse(fRight_Cauchy_Green_tensor);
				fs_mix_out	<<"fRight_Cauchy_Green_tensor_Inverse"<< endl ;
				fs_mix_out	<<fRight_Cauchy_Green_tensor_Inverse<< endl ;
				fs_mix_out	<<"fRight_Cauchy_Green_tensor_Inverse"<< endl ;


				/* [fLeft_Cauchy_Green_tensor] will be formed */
				fLeft_Cauchy_Green_tensor.Transpose(fRight_Cauchy_Green_tensor);
				fs_mix_out	<<"fLeft_Cauchy_Green_tensor"<< endl ;
				fs_mix_out	<<fLeft_Cauchy_Green_tensor<< endl ;
				fs_mix_out	<<"fLeft_Cauchy_Green_tensor"<< endl ;




                                /* [fEffective_Second_Piola_tensor] will be formed */
				fEffective_Second_Piola_tensor.SetToScaled(fMaterial_Params[kLambda]*log(J)-fMaterial_Params[kMu],fRight_Cauchy_Green_tensor_Inverse); 
				fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fIdentity_matrix);
				fEffective_Second_Piola_tensor += fTemp_matrix_nsd_x_nsd;
				fs_mix_out	<<"fEffective_Second_Piola_tensor"<< endl ;
				fs_mix_out	<<fEffective_Second_Piola_tensor<< endl ;
				fs_mix_out	<<"fEffective_Second_Piola_tensor"<< endl ;



                                /* [fEffective_Kirchhoff_tensor] will be formed */
				fEffective_Kirchhoff_tensor.MultABCT(fDeformation_Gradient,fEffective_Second_Piola_tensor,fDeformation_Gradient);
				fs_mix_out	<<"fEffective_Kirchhoff_tensor"<< endl ;
				fs_mix_out	<<fEffective_Kirchhoff_tensor<< endl ;
				fs_mix_out	<<"fEffective_Kirchhoff_tensor"<< endl ;


				/* {fEffective_Kirchhoff_vector} will be formed */
				Form_effective_kirchhoff_stress_vector();
				fs_mix_out	<<"fEffective_Kirchhoff_vector"<< endl ;
				fs_mix_out	<<fEffective_Kirchhoff_vector<< endl ;
				fs_mix_out	<<"fEffective_Kirchhoff_vector"<< endl ;



				/* [fIota_temp_matrix] will be formed */
				fIota_temp_matrix.MultATB(fShapeSolidGrad,fDefGradInv_GRAD_grad);
				fs_mix_out	<<"fIota_temp_matrix"<< endl ;
				fs_mix_out	<<fIota_temp_matrix<< endl ;
				fs_mix_out	<<"fIota_temp_matrix"<< endl ;


				/* second derivatives of solid shape functions, [fShapeSolidGradGrad] will be formed */
				fShapes_displ->Grad_GradNa(fShapeSolidGradGrad);
				fs_mix_out	<<"fShapeSolidGradGrad"<< endl ;
				fs_mix_out	<<fShapeSolidGradGrad<< endl ;
				fs_mix_out	<<"fShapeSolidGradGrad"<< endl ;


				/* [fVarpi_temp_matrix] will be formed */
				Form_Varpi_temp_matrix();
				fs_mix_out	<<"fVarpi_temp_matrix"<< endl ;
				fs_mix_out	<<fVarpi_temp_matrix<< endl ;
				fs_mix_out	<<"fVarpi_temp_matrix"<< endl ;



                                /* hydraulic conductivity matrix in the current coordinate, [k] will be formed */
				fK_hydraulic_conductivity_matrix.SetToScaled(fMaterial_Params[kK],fIdentity_matrix); 
				fk_hydraulic_conductivity_matrix.SetToScaled(1/J,fK_hydraulic_conductivity_matrix); 
				fs_mix_out	<<"1/J"<< endl ;
				fs_mix_out	<<1/J<< endl ;
				fs_mix_out	<<"1/J"<< endl ;
				fs_mix_out	<<"fk_hydraulic_conductivity_matrix"<< endl ;
				fs_mix_out	<<fk_hydraulic_conductivity_matrix<< endl ;
				fs_mix_out	<<"fk_hydraulic_conductivity_matrix"<< endl ;
				fTemp_matrix_nsd_x_nsd.MultABCT(fDeformation_Gradient,fk_hydraulic_conductivity_matrix,fDeformation_Gradient);
				fk_hydraulic_conductivity_matrix = fTemp_matrix_nsd_x_nsd;
				fs_mix_out	<<"fK_hydraulic_conductivity_matrix"<< endl ;
				fs_mix_out	<<fK_hydraulic_conductivity_matrix<< endl ;
				fs_mix_out	<<"fK_hydraulic_conductivity_matrix"<< endl ;
				fs_mix_out	<<"fk_hydraulic_conductivity_matrix"<< endl ;
				fs_mix_out	<<fk_hydraulic_conductivity_matrix<< endl ;
				fs_mix_out	<<"fk_hydraulic_conductivity_matrix"<< endl ;



				/* [fLambda_temp_matrix] will be formed */
				fLambda_temp_matrix.MultATBC(fShapeFluidGrad,fDeformation_Gradient_Inverse,fk_hydraulic_conductivity_matrix);
				fs_mix_out	<<"fLambda_temp_matrix"<< endl ;
				fs_mix_out	<<fLambda_temp_matrix<< endl ;
				fs_mix_out	<<"fLambda_temp_matrix"<< endl ;


				/* {fChi_temp_vector} will be formed */
				fVarpi_temp_matrix.Multx(u_vec,fChi_temp_vector);
				fs_mix_out	<<"fChi_temp_vector"<< endl ;
				fs_mix_out	<<fChi_temp_vector<< endl ;
				fs_mix_out	<<"fChi_temp_vector"<< endl ;


                                /* [fChi_temp_column_matrix] will be formed */
                                for (int i=0; i<3 ; i++)
				    fChi_temp_column_matrix(i,0)= fChi_temp_vector[i];
				fs_mix_out	<<"fChi_temp_column_matrix"<< endl ;
				fs_mix_out	<<fChi_temp_column_matrix<< endl ;
				fs_mix_out	<<"fChi_temp_column_matrix"<< endl ;


				/* {fFd_int_N1_vector} will be formed */
				double scale = (*Weight)*(*Det);
				fIota_temp_matrix.Multx(fEffective_Kirchhoff_vector,fTemp_vector_ndof_se,scale);
                                /* accumulate */
				fFd_int_N1_vector += fTemp_vector_ndof_se;
				fs_mix_out	<<"fFd_int_N1_vector"<< endl ;
				fs_mix_out	<<fFd_int_N1_vector<< endl ;
				fs_mix_out	<<"fFd_int_N1_vector"<< endl ;


                                /* {fFd_int_N2_vector} will be formed */
				theta = fShapeFluid[0]*press_vec[0];
				for (int i=1; i<8; i++)
				    theta += fShapeFluid[i]*press_vec[i];
				fs_mix_out	<<"theta"<< endl ;
				fs_mix_out	<<theta<< endl ;
				fs_mix_out	<<"theta"<< endl ;
				scale = -1*theta*(*Weight)*(*Det);
				fShapeSolidGrad.MultTx(fDefGradInv_Vector,fTemp_vector_ndof_se,scale);
				/* accumulate */
				fFd_int_N2_vector += fTemp_vector_ndof_se; 
				fs_mix_out	<<"fFd_int_N2_vector"<< endl ;
				fs_mix_out	<<fFd_int_N2_vector<< endl ;
				fs_mix_out	<<"fFd_int_N2_vector"<< endl ;



                                /* {fFtheta_int_N1_vector} will be formed */
				phi_s = fMaterial_Params[kPhi_s0]/J;
				phi_f = 1.0 - phi_s;
				double const1 = fMaterial_Params[kg]* phi_f;
				if (fabs(const1)>0.0e-16)
				    scale = theta*(-1.0/fMaterial_Params[kg]) * ((phi_s/phi_f)-1) * (*Weight)*(*Det);
				else
				    scale = 0.0;
				fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
				fTemp_matrix_nen_press_x_nsd.Multx(fChi_temp_vector, fTemp_vector_nen_press,scale);
                                /* accumulate */
				fFtheta_int_N1_vector += fTemp_vector_nen_press;
				fs_mix_out	<<"fFtheta_int_N1_vector"<< endl ;
				fs_mix_out	<<fFtheta_int_N1_vector<< endl ;
				fs_mix_out	<<"fFtheta_int_N1_vector"<< endl ;


                                /* {fFtheta_int_N2_vector} will be formed */
				fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_nsd,fShapeFluidGrad);
				scale = -1.0/fMaterial_Params[kg]*(*Weight)*(*Det); 
				fTemp_matrix_nen_press_x_nen_press.Multx(press_vec, fTemp_vector_nen_press,scale);
                                /* accumulate */
				fFtheta_int_N2_vector += fTemp_vector_nen_press;
				fs_mix_out	<<"fFtheta_int_N2_vector"<< endl ;
				fs_mix_out	<<fFtheta_int_N2_vector<< endl ;
				fs_mix_out	<<"fFtheta_int_N2_vector"<< endl ;


				/* [fIm_temp_matrix] will be formed */
				Form_Im_temp_matrix();
				fs_mix_out	<<"fIm_temp_matrix"<< endl ;
				fs_mix_out	<<fIm_temp_matrix<< endl ;
				fs_mix_out	<<"fIm_temp_matrix"<< endl ;


				/* [fHbar_temp_matrix] will be formed */
				Form_Hbar_temp_matrix();
				fs_mix_out	<<"fHbar_temp_matrix"<< endl ;
				fs_mix_out	<<fHbar_temp_matrix<< endl ;
				fs_mix_out	<<"fHbar_temp_matrix"<< endl ;


				/* [fEth_temp_matrix] will be formed */
				Form_Eth_temp_matrix();
				fs_mix_out	<<"fEth_temp_matrix"<< endl ;
				fs_mix_out	<<fEth_temp_matrix<< endl ;
				fs_mix_out	<<"fEth_temp_matrix"<< endl ;


				/* {fPi_temp_transpose_vector} will be formed */
				fShapeSolidGrad.MultTx(fDefGradInv_Vector,fPi_temp_transpose_vector);
				fs_mix_out	<<"fPi_temp_transpose_vector"<< endl ;
				fs_mix_out	<<fPi_temp_transpose_vector<< endl ;
				fs_mix_out	<<"fPi_temp_transpose_vector"<< endl ;


				/* [fPi_temp_row_matrix] will be formed */
				for (int i=0; i<n_en_displ_x_n_sd; i++)
								    fPi_temp_row_matrix(0,i) = fPi_temp_transpose_vector[i];
				fs_mix_out	<<"fPi_temp_row_matrix"<< endl ;
				fs_mix_out	<<fPi_temp_row_matrix<< endl ;
				fs_mix_out	<<"fPi_temp_row_matrix"<< endl ;

				/* [fK_dd_G3_1_matrix] will be formed */
				fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fIm_temp_matrix,fIota_temp_matrix);
				scale = -1*fIntegration_Params[kBeta]*(*Weight)*(*Det);
				fTemp_matrix_ndof_se_x_ndof_se *= scale;
                                /* accumulate */
				fK_dd_G3_1_matrix += fTemp_matrix_ndof_se_x_ndof_se;
				fs_mix_out	<<"fK_dd_G3_1_matrix"<< endl ;
				fs_mix_out	<<fK_dd_G3_1_matrix<< endl ;
				fs_mix_out	<<"fK_dd_G3_1_matrix"<< endl ;

                                /* [fI_ij_column_matrix] will be formed */
				fI_ij_column_matrix = 0.0;
				fI_ij_column_matrix(0,0) = 1.0;
				fI_ij_column_matrix(4,0) = 1.0;
				fI_ij_column_matrix(8,0) = 1.0;
				fs_mix_out	<<"fI_ij_column_matrix"<< endl ;
				fs_mix_out	<<fI_ij_column_matrix<< endl ;
				fs_mix_out	<<"fI_ij_column_matrix"<< endl ;


				/* [fK_dd_G3_2_matrix] will be formed */
				fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fHbar_temp_matrix,fIota_temp_matrix);
				scale = fMaterial_Params[kMu] * fIntegration_Params[kBeta] * (*Weight)*(*Det);
				fTemp_matrix_ndof_se_x_ndof_se *= scale;
                                /* accumulate */
				fK_dd_G3_2_matrix += fTemp_matrix_ndof_se_x_ndof_se; 
				fs_mix_out	<<"fK_dd_G3_2_matrix"<< endl ;
				fs_mix_out	<<fK_dd_G3_2_matrix<< endl ;
				fs_mix_out	<<"fK_dd_G3_2_matrix"<< endl ;

				/* [fK_dd_G3_3_matrix] will be formed */
				fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fEth_temp_matrix,fIota_temp_matrix);
				scale = fMaterial_Params[kMu] * fIntegration_Params[kBeta] * (*Weight)*(*Det);
				fTemp_matrix_ndof_se_x_ndof_se *= scale;
                                /* accumulate */
				fK_dd_G3_3_matrix += fTemp_matrix_ndof_se_x_ndof_se;
				fs_mix_out	<<"fK_dd_G3_3_matrix"<< endl ;
				fs_mix_out	<<fK_dd_G3_3_matrix<< endl ;
				fs_mix_out	<<"fK_dd_G3_3_matrix"<< endl ;

				/* [fK_dd_G3_4_matrix] will be formed */
				fTemp_matrix_ndof_se_x_ndof_se.MultABC(fIota_temp_matrix,fI_ij_column_matrix,fPi_temp_row_matrix);
				scale = fMaterial_Params[kLambda] * fIntegration_Params[kBeta] * (*Weight)*(*Det); 
				fTemp_matrix_ndof_se_x_ndof_se *= scale;
                                /* accumulate */
				fK_dd_G3_4_matrix += fTemp_matrix_ndof_se_x_ndof_se; 
				fs_mix_out	<<"fK_dd_G3_4_matrix"<< endl ;
				fs_mix_out	<<fK_dd_G3_4_matrix<< endl ;
				fs_mix_out	<<"fK_dd_G3_4_matrix"<< endl ;

				/* [fK_dd_G3_5_matrix] will be formed */
				fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fShapeSolidGrad_t_Transpose,fDefGradInv_GRAD_grad,fIota_temp_matrix);
				scale = theta * fIntegration_Params[kBeta] * (*Weight)*(*Det);
				fTemp_matrix_ndof_se_x_ndof_se *= scale;
                                /* accumulate */
				fK_dd_G3_5_matrix += fTemp_matrix_ndof_se_x_ndof_se;
				fs_mix_out	<<"fK_dd_G3_5_matrix"<< endl ;
				fs_mix_out	<<fK_dd_G3_5_matrix<< endl ;
				fs_mix_out	<<"fK_dd_G3_5_matrix"<< endl ;

                                /* [fK_dtheta_G3_matrix] will be formed */
				fTemp_matrix_ndof_se_x_nen_press.MultATB(fPi_temp_row_matrix,fShapeFluid_row_matrix); 
				scale = -1*fIntegration_Params[kBeta]*(*Weight)*(*Det);
				fTemp_matrix_ndof_se_x_nen_press *= scale;
                                /* accumulate */
				fK_dtheta_G3_matrix += fTemp_matrix_ndof_se_x_nen_press;
				fs_mix_out	<<"fK_dtheta_G3_matrix"<< endl ;
				fs_mix_out	<<fK_dtheta_G3_matrix<< endl ;
				fs_mix_out	<<"fK_dtheta_G3_matrix"<< endl ;

                                /* {fGrad_1_J_vector} will be filled */
				fVarpi_temp_matrix.Multx(u_vec,fGrad_1_J_vector, -1.0/J);
				fs_mix_out	<<"fGrad_1_J_vector"<< endl ;
				fs_mix_out	<<fGrad_1_J_vector<< endl ;
				fs_mix_out	<<"fGrad_1_J_vector"<< endl ;

                                /* {fGrad_theta_vector} will be filled */
				fShapeFluidGrad.Multx(press_vec, fGrad_theta_vector);
				fs_mix_out	<<"fGrad_theta_vector"<< endl ;
				fs_mix_out	<<fGrad_theta_vector<< endl ;
				fs_mix_out	<<"fGrad_theta_vector"<< endl ;

                                /* {fGrad_phi_f_vector} will be filled */
				fGrad_phi_f_vector.SetToScaled(-1* fMaterial_Params[kPhi_s0],fGrad_1_J_vector);
				fs_mix_out	<<"fGrad_phi_f_vector"<< endl ;
				fs_mix_out	<<fGrad_phi_f_vector<< endl ;
				fs_mix_out	<<"fGrad_phi_f_vector"<< endl ;

                                /* {fGrad_Omega_vector} will be filled */
				fTemp_nsd_vector.SetToScaled(theta/J,fGrad_phi_f_vector) ; 
				fGrad_Omega_vector = fTemp_nsd_vector;
				fTemp_nsd_vector.SetToScaled(phi_f/J,fGrad_theta_vector) ;
				fGrad_Omega_vector += fTemp_nsd_vector;
				fTemp_nsd_vector.SetToScaled(phi_f * theta,fGrad_1_J_vector) ;
				fGrad_Omega_vector += fTemp_nsd_vector;
				fs_mix_out	<<"fGrad_Omega_vector"<< endl ;
				fs_mix_out	<<fGrad_Omega_vector<< endl ;
				fs_mix_out	<<"fGrad_Omega_vector"<< endl ;

                                /* {fgrad_Omega_vector} will be formed */
				fDeformation_Gradient_Inverse_Transpose.Multx(fGrad_Omega_vector,fgrad_Omega_vector);
				fs_mix_out	<<"fgrad_Omega_vector"<< endl ;
				fs_mix_out	<<fgrad_Omega_vector<< endl ;
				fs_mix_out	<<"fgrad_Omega_vector"<< endl ;

				/* [fJmath_temp_matrix] will be formed */
				Form_Jmath_temp_matrix(); 
				fs_mix_out	<<"fJmath_temp_matrix"<< endl ;
				fs_mix_out	<<fJmath_temp_matrix<< endl ;
				fs_mix_out	<<"fJmath_temp_matrix"<< endl ;


                                /* [fWp_temp_matrix] will be formed */
				Form_Wp_temp_matrix(); 
				fs_mix_out	<<"fWp_temp_matrix"<< endl ;
				fs_mix_out	<<fWp_temp_matrix<< endl ;
				fs_mix_out	<<"fWp_temp_matrix"<< endl ;


                                /* [fK_thetad_H3_1_matrix] will be formed */
				const1 = fMaterial_Params[kg]*phi_f * J ;
				if (fabs(const1) > 0.0e-16) 
				    scale = theta*(fIntegration_Params[kBeta]/(fMaterial_Params[kg]*phi_f))*((fMaterial_Params[kPhi_s0]/J-phi_f)*(phi_s/phi_f)+2*fMaterial_Params[kPhi_s0]/J)*(*Weight)*(*Det);
				else
				    scale = 0.0;
				fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
				fTemp_matrix_nen_press_x_ndof_se.MultABC(fTemp_matrix_nen_press_x_nsd,fChi_temp_column_matrix,fPi_temp_row_matrix);
				fTemp_matrix_nen_press_x_ndof_se *= scale;
                                /* accumulate */
				fK_thetad_H3_1_matrix += fTemp_matrix_nen_press_x_ndof_se;
				fs_mix_out	<<"fK_thetad_H3_1_matrix"<< endl ;
				fs_mix_out	<<fK_thetad_H3_1_matrix<< endl ;
				fs_mix_out	<<"fK_thetad_H3_1_matrix"<< endl ;


                                /* [fK_thetad_H3_2_matrix] will be formed */
				const1 = fMaterial_Params[kg]*J*phi_f ;
				if (fabs(const1) > 0.0e-16)
				    scale = -1*(fIntegration_Params[kBeta]/fMaterial_Params[kg])*theta*(fMaterial_Params[kPhi_s0]/(J*phi_f)-1)*(*Weight)*(*Det);
				else
				    scale = 0.0;
				fTemp_matrix_nen_press_x_ndof_se.MultAB(fTemp_matrix_nen_press_x_nsd,fVarpi_temp_matrix);
				fTemp_matrix_nen_press_x_ndof_se *= scale;
                                /* accumulate */
				fK_thetad_H3_2_matrix += fTemp_matrix_nen_press_x_ndof_se;
				fs_mix_out	<<"fK_thetad_H3_2_matrix"<< endl ;
				fs_mix_out	<<fK_thetad_H3_2_matrix<< endl ;
				fs_mix_out	<<"fK_thetad_H3_2_matrix"<< endl ;


                                /* [fK_thetad_H3_3_matrix] will be formed */
				const1 = fMaterial_Params[kg]*phi_f;
				if (fabs(const1) > 0.0e-16)
				    scale = fIntegration_Params[kBeta]*J/(fMaterial_Params[kg]*phi_f)*(*Weight)*(*Det);
				else
				    scale = 0.0;
				fTemp_matrix_nen_press_x_nsd.MultATB(fShapeFluidGrad,fDeformation_Gradient_Inverse);
				fTemp_matrix_nen_press_x_ndof_se.MultABCT(fTemp_matrix_nen_press_x_nsd,fJmath_temp_matrix,fIota_temp_matrix);
				fTemp_matrix_nen_press_x_ndof_se *= scale;
                                /* accumulate */
				fK_thetad_H3_3_matrix += fTemp_matrix_nen_press_x_ndof_se;
				fs_mix_out	<<"fK_thetad_H3_3_matrix"<< endl ;
				fs_mix_out	<<fK_thetad_H3_3_matrix<< endl ;
				fs_mix_out	<<"fK_thetad_H3_3_matrix"<< endl ;


                                /* [fK_thetad_H3_4_matrix] will be formed */
				fTemp_matrix_nen_press_x_ndof_se.MultABCT(fTemp_matrix_nen_press_x_nsd,fWp_temp_matrix,fIota_temp_matrix);
				fTemp_matrix_nen_press_x_ndof_se *= scale;
                                /* accumulate */
				fK_thetad_H3_4_matrix += fTemp_matrix_nen_press_x_ndof_se;
				fs_mix_out	<<"fK_thetad_H3_4_matrix"<< endl ;
				fs_mix_out	<<fK_thetad_H3_4_matrix<< endl ;
				fs_mix_out	<<"fK_thetad_H3_4_matrix"<< endl ;


                                /* [fK_thetatheta_H3_1_matrix] will be formed */
				scale = -1*(fIntegration_Params[kBeta]/fMaterial_Params[kg])*(*Weight)*(*Det);
				fTemp_matrix_nen_press_x_nsd.MultAB(fLambda_temp_matrix,fDeformation_Gradient_Inverse_Transpose);
				fTemp_matrix_nen_press_x_nen_press.MultAB(fTemp_matrix_nen_press_x_nsd,fShapeFluidGrad);
				fTemp_matrix_nen_press_x_nen_press *= scale;
                                /* accumulate */
				fK_thetatheta_H3_1_matrix += fTemp_matrix_nen_press_x_nen_press;
				fs_mix_out	<<"fK_thetatheta_H3_1_matrix"<< endl ;
				fs_mix_out	<<fK_thetatheta_H3_1_matrix<< endl ;
				fs_mix_out	<<"fK_thetatheta_H3_1_matrix"<< endl ;


                                /* [fK_thetatheta_H3_2_matrix] will be formed */
				const1 = fMaterial_Params[kg] * J*phi_f ;
				if (fabs(const1)> 0.0e-16)
				    scale =  -1*(fIntegration_Params[kBeta]/fMaterial_Params[kg])*(fMaterial_Params[kPhi_s0]/(J*phi_f)-1)*(*Weight)*(*Det);
				else
				    scale = 0.0;
				fTemp_matrix_nen_press_x_nen_press.MultABC(fTemp_matrix_nen_press_x_nsd,fChi_temp_column_matrix,fShapeFluid_row_matrix);
				fTemp_matrix_nen_press_x_nen_press *= scale;
                                /* accumulate */
				fK_thetatheta_H3_2_matrix += fTemp_matrix_nen_press_x_nen_press;
				fs_mix_out	<<"fK_thetatheta_H3_2_matrix"<< endl ;
				fs_mix_out	<<fK_thetatheta_H3_2_matrix<< endl ;
				fs_mix_out	<<"fK_thetatheta_H3_2_matrix"<< endl ;
                                
                                /* advancing weight and jacobian determinent pointers */
				Weight++;
				Det++;

/* Test part
				Test_vector_A[0] =3;
				Test_vector_A[1]=7;
				Test_vector_A[2]=5;
				Test_vector_B.SetToScaled(5,Test_vector_A) ; 

				fTest_matrix_C.MultAB(fTest_matrix_A,fTest_matrix_B);
				double *pfTest_matrix_A	= fTest_matrix_A.Pointer();
				fs_mix_out	<<"LLLLLLLLLLLLLLLLLLLLLLL"<< endl ;			
				fs_mix_out	<<pfTest_matrix_A[0]<< endl ;
				fs_mix_out	<<pfTest_matrix_A[1]<< endl ;
				fs_mix_out	<<pfTest_matrix_A[2]<< endl ;
				fs_mix_out	<<pfTest_matrix_A[3]<< endl ;
				fs_mix_out	<<"LLLLLLLLLLLLLLLLLLLLLLLLLL"<< endl ; 
				fTest_matrix_A.SetToScaled(3,fTest_matrix_A); 
				fs_mix_out	<< "TestTestTestTestTest" << endl;
				fs_mix_out	<< Test_vector_A << endl;
				fs_mix_out	<< Test_vector_B << endl;
				fs_mix_out	<< "TestTestTestTestTest" << endl;
				Test_vector_A += Test_vector_B;
				fs_mix_out	<< "TestTestTestTestTest" << endl;
				fs_mix_out	<< Test_vector_A << endl;
				fs_mix_out	<< Test_vector_B << endl;
				fs_mix_out	<< "TestTestTestTestTest" << endl;
				fs_mix_out	<< endl << endl;                                        */







// defining F_1_T
				/* for debugging */
				const int ip = fShapes_displ->CurrIP()+1;
				fs_mix_out	<< endl << "IP" << ip
						<< setw(outputFileWidth) << ", shape function matrix for solid phase: " 
						<< setw(outputFileWidth) << fShapeSolid;
				
				fs_mix_out	<< endl << "terms from shape function matrix for solid phase: " 
						<< setw(outputFileWidth) << fShapeSolid(0,0) 
						<< setw(outputFileWidth) << fShapeSolid(0,3);			
		    }


		    fKdd = fK_dd_G3_1_matrix;
		    fKdd += fK_dd_G3_2_matrix;
		    fKdd += fK_dd_G3_3_matrix;
		    fKdd += fK_dd_G3_4_matrix;
		    fKdd += fK_dd_G3_5_matrix;

		    fKdtheta= fK_dtheta_G3_matrix;

		    fFd_int = fFd_int_N1_vector;
		    fFd_int += fFd_int_N2_vector;
				
		    fKthetad = fK_thetad_H3_1_matrix;
		    fKthetad += fK_thetad_H3_2_matrix;
		    fKthetad += fK_thetad_H3_3_matrix;
		    fKthetad += fK_thetad_H3_4_matrix;

		    fKthetatheta = fK_thetatheta_H3_1_matrix;
		    fKthetatheta += fK_thetatheta_H3_2_matrix;

		    fFtheta_int = fFtheta_int_N1_vector;
		    fFtheta_int += fFtheta_int_N2_vector;
				fs_mix_out	<<"kdd"<< endl ;
				fs_mix_out	<<fKdd<< endl ;
				fs_mix_out	<<"kdtheta"<< endl ;
				fs_mix_out	<<fKdtheta<< endl ;
				fs_mix_out	<<"Fd_int"<< endl ;
				fs_mix_out	<<fFd_int<< endl ;
				fs_mix_out	<<"Kthetad"<< endl ;
				fs_mix_out	<<fKthetad<< endl ;
				fs_mix_out	<<"Kthetatheta"<< endl ;
				fs_mix_out	<< fKthetatheta<< endl ;
				fs_mix_out	<<"Ftheta_int"<< endl ;
				fs_mix_out	<<fFtheta_int<< endl ;

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
	
    double shearMu, sLambda, sAlpha, Rho_sR0, Rho_fR0, Phi_s0, Phi_f0, bulkK, 
	newBeta, newGamma, conductivityK,gravity_g ;
			
    // solid elasticity
    list.AddParameter(shearMu, "mu");
    list.AddParameter(sLambda, "lambda");
	
    // fluid elasticity
    list.AddParameter(bulkK, "Kf");

    // fluid hydraulic conductivity
    list.AddParameter(conductivityK, "K");

    // gravity
    list.AddParameter(gravity_g, "g");
	
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
    fMaterial_Params[kK] = list.GetParameter("K");
    fMaterial_Params[kg] = list.GetParameter("g");
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
    u_n.Dimension (n_en_displ, n_sd);
    del_u.Dimension (n_en_displ, n_sd);
    n_en_displ_x_n_sd = n_en_displ*n_sd;
    del_u_vec.Dimension (n_en_displ_x_n_sd);
    u_vec.Dimension (n_en_displ_x_n_sd);
    //ElementSupport().RegisterCoordinates(fInitCoords_displ);
    fDispl->RegisterLocal(u);
    fDispl->RegisterLocal(u_n);

    /* set local arrays for pore pressure field */
    int dum=1;
    press.Dimension (n_en_press, dum);
    press_n.Dimension (n_en_press, dum);
    del_press.Dimension (n_en_press, dum);
    del_press_vec.Dimension (n_en_press);
    press_vec.Dimension (n_en_press);
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
    fGRAD_disp_vector.Dimension (n_sd_x_n_sd);
    fDeformation_Gradient_Inverse.Dimension (n_sd,n_sd);
    fDeformation_Gradient_Transpose.Dimension (n_sd,n_sd);
    fDeformation_Gradient_Inverse_Transpose.Dimension (n_sd,n_sd);
    fDefGradInv_GRAD_grad.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fDefGradInv_GRAD_grad_Transpose.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fDefGradInv_Vector.Dimension (n_sd_x_n_sd);
    fRight_Cauchy_Green_tensor.Dimension (n_sd,n_sd);
    fRight_Cauchy_Green_tensor_Inverse.Dimension (n_sd,n_sd);
    fLeft_Cauchy_Green_tensor.Dimension (n_sd,n_sd);
    fTest_matrix_A.Dimension (2,2);
    fTest_matrix_B.Dimension (2,2);
    fTest_matrix_C.Dimension (2,2);
    fIdentity_matrix.Dimension (n_sd,n_sd);
    fEffective_Second_Piola_tensor.Dimension (n_sd,n_sd);
    fTemp_matrix_nsd_x_nsd.Dimension (n_sd,n_sd);
    fTemp_matrix_nen_press_x_nsd.Dimension (n_en_press,n_sd);
    fTemp_matrix_nen_press_x_nen_press.Dimension (n_en_press,n_en_press);
    fEffective_Kirchhoff_tensor.Dimension (n_sd,n_sd);
    fEffective_Kirchhoff_vector.Dimension (n_sd_x_n_sd);
    fIota_temp_matrix.Dimension (n_en_displ_x_n_sd,n_sd_x_n_sd);
    fVarpi_temp_matrix.Dimension (n_sd, n_en_displ_x_n_sd);
    fk_hydraulic_conductivity_matrix.Dimension (n_sd,n_sd);
    fK_hydraulic_conductivity_matrix.Dimension (n_sd,n_sd);
    fLambda_temp_matrix.Dimension (n_en_press,n_sd);
    fChi_temp_vector.Dimension (n_sd);
    fFd_int_N1_vector.Dimension (n_en_displ_x_n_sd);
    fFd_int_N2_vector.Dimension (n_en_displ_x_n_sd); 
    fTemp_vector_ndof_se.Dimension (n_en_displ_x_n_sd); 
    fFtheta_int_N1_vector.Dimension (n_en_press); 
    fFtheta_int_N2_vector.Dimension (n_en_press); 
    fTemp_vector_nen_press.Dimension (n_en_press); 
    fIm_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fHbar_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fEth_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
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

    Test_vector_A.Dimension(3);
    Test_vector_B.Dimension(3);

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
    fs_mix_out.open("fs_mix.info");
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

void FSSolidFluidMixT::Form_solid_shape_functions(const double* &shapes_displ_X)
{
    fShapeSolid = 0.0;
    fShapeSolid(0,0) = shapes_displ_X[0];
    fShapeSolid(0,3) = shapes_displ_X[1];
    fShapeSolid(0,6) = shapes_displ_X[2];
    fShapeSolid(0,9) = shapes_displ_X[3];
    fShapeSolid(0,12) = shapes_displ_X[4];
    fShapeSolid(0,15) = shapes_displ_X[5];
    fShapeSolid(0,18) = shapes_displ_X[6];
    fShapeSolid(0,21) = shapes_displ_X[7];
    fShapeSolid(0,24) = shapes_displ_X[8];
    fShapeSolid(0,27) = shapes_displ_X[9];
    fShapeSolid(0,30) = shapes_displ_X[10];
    fShapeSolid(0,33) = shapes_displ_X[11];
    fShapeSolid(0,36) = shapes_displ_X[12];
    fShapeSolid(0,39) = shapes_displ_X[13];
    fShapeSolid(0,42) = shapes_displ_X[14];
    fShapeSolid(0,45) = shapes_displ_X[15];
    fShapeSolid(0,48) = shapes_displ_X[16];
    fShapeSolid(0,51) = shapes_displ_X[17];
    fShapeSolid(0,54) = shapes_displ_X[18];
    fShapeSolid(0,57) = shapes_displ_X[19];
    fShapeSolid(0,60) = shapes_displ_X[20];
    fShapeSolid(0,63) = shapes_displ_X[21];
    fShapeSolid(0,66) = shapes_displ_X[22];
    fShapeSolid(0,69) = shapes_displ_X[23];
    fShapeSolid(0,72) = shapes_displ_X[24];
    fShapeSolid(0,75) = shapes_displ_X[25];
    fShapeSolid(0,78) = shapes_displ_X[26];
    
    fShapeSolid(1,1) = shapes_displ_X[0];
    fShapeSolid(1,4) = shapes_displ_X[1];
    fShapeSolid(1,7) = shapes_displ_X[2];
    fShapeSolid(1,10) = shapes_displ_X[3];
    fShapeSolid(1,13) = shapes_displ_X[4];
    fShapeSolid(1,16) = shapes_displ_X[5];
    fShapeSolid(1,19) = shapes_displ_X[6];
    fShapeSolid(1,22) = shapes_displ_X[7];
    fShapeSolid(1,25) = shapes_displ_X[8];
    fShapeSolid(1,28) = shapes_displ_X[9];
    fShapeSolid(1,31) = shapes_displ_X[10];
    fShapeSolid(1,34) = shapes_displ_X[11];
    fShapeSolid(1,37) = shapes_displ_X[12];
    fShapeSolid(1,40) = shapes_displ_X[13];
    fShapeSolid(1,43) = shapes_displ_X[14];
    fShapeSolid(1,46) = shapes_displ_X[15];
    fShapeSolid(1,49) = shapes_displ_X[16];
    fShapeSolid(1,52) = shapes_displ_X[17];
    fShapeSolid(1,55) = shapes_displ_X[18];
    fShapeSolid(1,58) = shapes_displ_X[19];
    fShapeSolid(1,61) = shapes_displ_X[20];
    fShapeSolid(1,64) = shapes_displ_X[21];
    fShapeSolid(1,67) = shapes_displ_X[22];
    fShapeSolid(1,70) = shapes_displ_X[23];
    fShapeSolid(1,73) = shapes_displ_X[24];
    fShapeSolid(1,76) = shapes_displ_X[25];
    fShapeSolid(1,79) = shapes_displ_X[26];
    
    fShapeSolid(2,2) = shapes_displ_X[0];
    fShapeSolid(2,5) = shapes_displ_X[1];
    fShapeSolid(2,8) = shapes_displ_X[2];
    fShapeSolid(2,11) = shapes_displ_X[3];
    fShapeSolid(2,14) = shapes_displ_X[4];
    fShapeSolid(2,17) = shapes_displ_X[5];
    fShapeSolid(2,20) = shapes_displ_X[6];
    fShapeSolid(2,23) = shapes_displ_X[7];
    fShapeSolid(2,26) = shapes_displ_X[8];
    fShapeSolid(2,29) = shapes_displ_X[9];
    fShapeSolid(2,32) = shapes_displ_X[10];
    fShapeSolid(2,35) = shapes_displ_X[11];
    fShapeSolid(2,38) = shapes_displ_X[12];
    fShapeSolid(2,41) = shapes_displ_X[13];
    fShapeSolid(2,44) = shapes_displ_X[14];
    fShapeSolid(2,47) = shapes_displ_X[15];
    fShapeSolid(2,50) = shapes_displ_X[16];
    fShapeSolid(2,53) = shapes_displ_X[17];
    fShapeSolid(2,56) = shapes_displ_X[18];
    fShapeSolid(2,59) = shapes_displ_X[19];
    fShapeSolid(2,62) = shapes_displ_X[20];
    fShapeSolid(2,65) = shapes_displ_X[21];
    fShapeSolid(2,68) = shapes_displ_X[22];
    fShapeSolid(2,71) = shapes_displ_X[23];
    fShapeSolid(2,74) = shapes_displ_X[24];
    fShapeSolid(2,77) = shapes_displ_X[25];
    fShapeSolid(2,80) = shapes_displ_X[26]; 
}

void FSSolidFluidMixT::Form_Gradient_of_solid_shape_functions(const dMatrixT &fShapeSolidGrad_temp)
{
    fShapeSolidGrad = 0.0;
    fShapeSolidGrad(0,0) = fShapeSolidGrad_temp(0,0);
    fShapeSolidGrad(0,3) = fShapeSolidGrad_temp(0,1);
    fShapeSolidGrad(0,6) = fShapeSolidGrad_temp(0,2);
    fShapeSolidGrad(0,9) = fShapeSolidGrad_temp(0,3);
    fShapeSolidGrad(0,12) = fShapeSolidGrad_temp(0,4);
    fShapeSolidGrad(0,15) = fShapeSolidGrad_temp(0,5);
    fShapeSolidGrad(0,18) = fShapeSolidGrad_temp(0,6);
    fShapeSolidGrad(0,21) = fShapeSolidGrad_temp(0,7);
    fShapeSolidGrad(0,24) = fShapeSolidGrad_temp(0,8);
    fShapeSolidGrad(0,27) = fShapeSolidGrad_temp(0,9);
    fShapeSolidGrad(0,30) = fShapeSolidGrad_temp(0,10);
    fShapeSolidGrad(0,33) = fShapeSolidGrad_temp(0,11);
    fShapeSolidGrad(0,36) = fShapeSolidGrad_temp(0,12);
    fShapeSolidGrad(0,39) = fShapeSolidGrad_temp(0,13);
    fShapeSolidGrad(0,42) = fShapeSolidGrad_temp(0,14);
    fShapeSolidGrad(0,45) = fShapeSolidGrad_temp(0,15);
    fShapeSolidGrad(0,48) = fShapeSolidGrad_temp(0,16);
    fShapeSolidGrad(0,51) = fShapeSolidGrad_temp(0,17);
    fShapeSolidGrad(0,54) = fShapeSolidGrad_temp(0,18);
    fShapeSolidGrad(0,57) = fShapeSolidGrad_temp(0,19);
    fShapeSolidGrad(0,60) = fShapeSolidGrad_temp(0,20);
    fShapeSolidGrad(0,63) = fShapeSolidGrad_temp(0,21);
    fShapeSolidGrad(0,66) = fShapeSolidGrad_temp(0,22);
    fShapeSolidGrad(0,69) = fShapeSolidGrad_temp(0,23);
    fShapeSolidGrad(0,72) = fShapeSolidGrad_temp(0,24);
    fShapeSolidGrad(0,75) = fShapeSolidGrad_temp(0,25);
    fShapeSolidGrad(0,78) = fShapeSolidGrad_temp(0,26);
    
    fShapeSolidGrad(1,1) = fShapeSolidGrad_temp(0,0);
    fShapeSolidGrad(1,4) = fShapeSolidGrad_temp(0,1);
    fShapeSolidGrad(1,7) = fShapeSolidGrad_temp(0,2);
    fShapeSolidGrad(1,10) = fShapeSolidGrad_temp(0,3);
    fShapeSolidGrad(1,13) = fShapeSolidGrad_temp(0,4);
    fShapeSolidGrad(1,16) = fShapeSolidGrad_temp(0,5);
    fShapeSolidGrad(1,19) = fShapeSolidGrad_temp(0,6);
    fShapeSolidGrad(1,22) = fShapeSolidGrad_temp(0,7);
    fShapeSolidGrad(1,25) = fShapeSolidGrad_temp(0,8);
    fShapeSolidGrad(1,28) = fShapeSolidGrad_temp(0,9);
    fShapeSolidGrad(1,31) = fShapeSolidGrad_temp(0,10);
    fShapeSolidGrad(1,34) = fShapeSolidGrad_temp(0,11);
    fShapeSolidGrad(1,37) = fShapeSolidGrad_temp(0,12);
    fShapeSolidGrad(1,40) = fShapeSolidGrad_temp(0,13);
    fShapeSolidGrad(1,43) = fShapeSolidGrad_temp(0,14);
    fShapeSolidGrad(1,46) = fShapeSolidGrad_temp(0,15);
    fShapeSolidGrad(1,49) = fShapeSolidGrad_temp(0,16);
    fShapeSolidGrad(1,52) = fShapeSolidGrad_temp(0,17);
    fShapeSolidGrad(1,55) = fShapeSolidGrad_temp(0,18);
    fShapeSolidGrad(1,58) = fShapeSolidGrad_temp(0,19);
    fShapeSolidGrad(1,61) = fShapeSolidGrad_temp(0,20);
    fShapeSolidGrad(1,64) = fShapeSolidGrad_temp(0,21);
    fShapeSolidGrad(1,67) = fShapeSolidGrad_temp(0,22);
    fShapeSolidGrad(1,70) = fShapeSolidGrad_temp(0,23);
    fShapeSolidGrad(1,73) = fShapeSolidGrad_temp(0,24);
    fShapeSolidGrad(1,76) = fShapeSolidGrad_temp(0,25);
    fShapeSolidGrad(1,79) = fShapeSolidGrad_temp(0,26);
    
    
    fShapeSolidGrad(2,2) = fShapeSolidGrad_temp(0,0);
    fShapeSolidGrad(2,5) = fShapeSolidGrad_temp(0,1);
    fShapeSolidGrad(2,8) = fShapeSolidGrad_temp(0,2);
    fShapeSolidGrad(2,11) = fShapeSolidGrad_temp(0,3);
    fShapeSolidGrad(2,14) = fShapeSolidGrad_temp(0,4);
    fShapeSolidGrad(2,17) = fShapeSolidGrad_temp(0,5);
    fShapeSolidGrad(2,20) = fShapeSolidGrad_temp(0,6);
    fShapeSolidGrad(2,23) = fShapeSolidGrad_temp(0,7);
    fShapeSolidGrad(2,26) = fShapeSolidGrad_temp(0,8);
    fShapeSolidGrad(2,29) = fShapeSolidGrad_temp(0,9);
    fShapeSolidGrad(2,32) = fShapeSolidGrad_temp(0,10);
    fShapeSolidGrad(2,35) = fShapeSolidGrad_temp(0,11);
    fShapeSolidGrad(2,38) = fShapeSolidGrad_temp(0,12);
    fShapeSolidGrad(2,41) = fShapeSolidGrad_temp(0,13);
    fShapeSolidGrad(2,44) = fShapeSolidGrad_temp(0,14);
    fShapeSolidGrad(2,47) = fShapeSolidGrad_temp(0,15);
    fShapeSolidGrad(2,50) = fShapeSolidGrad_temp(0,16);
    fShapeSolidGrad(2,53) = fShapeSolidGrad_temp(0,17);
    fShapeSolidGrad(2,56) = fShapeSolidGrad_temp(0,18);
    fShapeSolidGrad(2,59) = fShapeSolidGrad_temp(0,19);
    fShapeSolidGrad(2,62) = fShapeSolidGrad_temp(0,20);
    fShapeSolidGrad(2,65) = fShapeSolidGrad_temp(0,21);
    fShapeSolidGrad(2,68) = fShapeSolidGrad_temp(0,22);
    fShapeSolidGrad(2,71) = fShapeSolidGrad_temp(0,23);
    fShapeSolidGrad(2,74) = fShapeSolidGrad_temp(0,24);
    fShapeSolidGrad(2,77) = fShapeSolidGrad_temp(0,25);
    fShapeSolidGrad(2,80) = fShapeSolidGrad_temp(0,26);
    
    fShapeSolidGrad(3,0) = fShapeSolidGrad_temp(1,0);
    fShapeSolidGrad(3,3) = fShapeSolidGrad_temp(1,1);
    fShapeSolidGrad(3,6) = fShapeSolidGrad_temp(1,2);
    fShapeSolidGrad(3,9) = fShapeSolidGrad_temp(1,3);
    fShapeSolidGrad(3,12) = fShapeSolidGrad_temp(1,4);
    fShapeSolidGrad(3,15) = fShapeSolidGrad_temp(1,5);
    fShapeSolidGrad(3,18) = fShapeSolidGrad_temp(1,6);
    fShapeSolidGrad(3,21) = fShapeSolidGrad_temp(1,7);
    fShapeSolidGrad(3,24) = fShapeSolidGrad_temp(1,8);
    fShapeSolidGrad(3,27) = fShapeSolidGrad_temp(1,9);
    fShapeSolidGrad(3,30) = fShapeSolidGrad_temp(1,10);
    fShapeSolidGrad(3,33) = fShapeSolidGrad_temp(1,11);
    fShapeSolidGrad(3,36) = fShapeSolidGrad_temp(1,12);
    fShapeSolidGrad(3,39) = fShapeSolidGrad_temp(1,13);
    fShapeSolidGrad(3,42) = fShapeSolidGrad_temp(1,14);
    fShapeSolidGrad(3,45) = fShapeSolidGrad_temp(1,15);
    fShapeSolidGrad(3,48) = fShapeSolidGrad_temp(1,16);
    fShapeSolidGrad(3,51) = fShapeSolidGrad_temp(1,17);
    fShapeSolidGrad(3,54) = fShapeSolidGrad_temp(1,18);
    fShapeSolidGrad(3,57) = fShapeSolidGrad_temp(1,19);
    fShapeSolidGrad(3,60) = fShapeSolidGrad_temp(1,20);
    fShapeSolidGrad(3,63) = fShapeSolidGrad_temp(1,21);
    fShapeSolidGrad(3,66) = fShapeSolidGrad_temp(1,22);
    fShapeSolidGrad(3,69) = fShapeSolidGrad_temp(1,23);
    fShapeSolidGrad(3,72) = fShapeSolidGrad_temp(1,24);
    fShapeSolidGrad(3,75) = fShapeSolidGrad_temp(1,25);
    fShapeSolidGrad(3,78) = fShapeSolidGrad_temp(1,26);
    
    fShapeSolidGrad(4,1) = fShapeSolidGrad_temp(1,0);
    fShapeSolidGrad(4,4) = fShapeSolidGrad_temp(1,1);
    fShapeSolidGrad(4,7) = fShapeSolidGrad_temp(1,2);
    fShapeSolidGrad(4,10) = fShapeSolidGrad_temp(1,3);
    fShapeSolidGrad(4,13) = fShapeSolidGrad_temp(1,4);
    fShapeSolidGrad(4,16) = fShapeSolidGrad_temp(1,5);
    fShapeSolidGrad(4,19) = fShapeSolidGrad_temp(1,6);
    fShapeSolidGrad(4,22) = fShapeSolidGrad_temp(1,7);
    fShapeSolidGrad(4,25) = fShapeSolidGrad_temp(1,8);
    fShapeSolidGrad(4,28) = fShapeSolidGrad_temp(1,9);
    fShapeSolidGrad(4,31) = fShapeSolidGrad_temp(1,10);
    fShapeSolidGrad(4,34) = fShapeSolidGrad_temp(1,11);
    fShapeSolidGrad(4,37) = fShapeSolidGrad_temp(1,12);
    fShapeSolidGrad(4,40) = fShapeSolidGrad_temp(1,13);
    fShapeSolidGrad(4,43) = fShapeSolidGrad_temp(1,14);
    fShapeSolidGrad(4,46) = fShapeSolidGrad_temp(1,15);
    fShapeSolidGrad(4,49) = fShapeSolidGrad_temp(1,16);
    fShapeSolidGrad(4,52) = fShapeSolidGrad_temp(1,17);
    fShapeSolidGrad(4,55) = fShapeSolidGrad_temp(1,18);
    fShapeSolidGrad(4,58) = fShapeSolidGrad_temp(1,19);
    fShapeSolidGrad(4,61) = fShapeSolidGrad_temp(1,20);
    fShapeSolidGrad(4,64) = fShapeSolidGrad_temp(1,21);
    fShapeSolidGrad(4,67) = fShapeSolidGrad_temp(1,22);
    fShapeSolidGrad(4,70) = fShapeSolidGrad_temp(1,23);
    fShapeSolidGrad(4,73) = fShapeSolidGrad_temp(1,24);
    fShapeSolidGrad(4,76) = fShapeSolidGrad_temp(1,25);
    fShapeSolidGrad(4,79) = fShapeSolidGrad_temp(1,26);
    
    fShapeSolidGrad(5,2) = fShapeSolidGrad_temp(1,0);
    fShapeSolidGrad(5,5) = fShapeSolidGrad_temp(1,1);
    fShapeSolidGrad(5,8) = fShapeSolidGrad_temp(1,2);
    fShapeSolidGrad(5,11) = fShapeSolidGrad_temp(1,3);
    fShapeSolidGrad(5,14) = fShapeSolidGrad_temp(1,4);
    fShapeSolidGrad(5,17) = fShapeSolidGrad_temp(1,5);
    fShapeSolidGrad(5,20) = fShapeSolidGrad_temp(1,6);
    fShapeSolidGrad(5,23) = fShapeSolidGrad_temp(1,7);
    fShapeSolidGrad(5,26) = fShapeSolidGrad_temp(1,8);
    fShapeSolidGrad(5,29) = fShapeSolidGrad_temp(1,9);
    fShapeSolidGrad(5,32) = fShapeSolidGrad_temp(1,10);
    fShapeSolidGrad(5,35) = fShapeSolidGrad_temp(1,11);
    fShapeSolidGrad(5,38) = fShapeSolidGrad_temp(1,12);
    fShapeSolidGrad(5,41) = fShapeSolidGrad_temp(1,13);
    fShapeSolidGrad(5,44) = fShapeSolidGrad_temp(1,14);
    fShapeSolidGrad(5,47) = fShapeSolidGrad_temp(1,15);
    fShapeSolidGrad(5,50) = fShapeSolidGrad_temp(1,16);
    fShapeSolidGrad(5,53) = fShapeSolidGrad_temp(1,17);
    fShapeSolidGrad(5,56) = fShapeSolidGrad_temp(1,18);
    fShapeSolidGrad(5,59) = fShapeSolidGrad_temp(1,19);
    fShapeSolidGrad(5,62) = fShapeSolidGrad_temp(1,20);
    fShapeSolidGrad(5,65) = fShapeSolidGrad_temp(1,21);
    fShapeSolidGrad(5,68) = fShapeSolidGrad_temp(1,22);
    fShapeSolidGrad(5,71) = fShapeSolidGrad_temp(1,23);
    fShapeSolidGrad(5,74) = fShapeSolidGrad_temp(1,24);
    fShapeSolidGrad(5,77) = fShapeSolidGrad_temp(1,25);
    fShapeSolidGrad(5,80) = fShapeSolidGrad_temp(1,26);
    
    fShapeSolidGrad(6,0) = fShapeSolidGrad_temp(2,0);
    fShapeSolidGrad(6,3) = fShapeSolidGrad_temp(2,1);
    fShapeSolidGrad(6,6) = fShapeSolidGrad_temp(2,2);
    fShapeSolidGrad(6,9) = fShapeSolidGrad_temp(2,3);
    fShapeSolidGrad(6,12) = fShapeSolidGrad_temp(2,4);
    fShapeSolidGrad(6,15) = fShapeSolidGrad_temp(2,5);
    fShapeSolidGrad(6,18) = fShapeSolidGrad_temp(2,6);
    fShapeSolidGrad(6,21) = fShapeSolidGrad_temp(2,7);
    fShapeSolidGrad(6,24) = fShapeSolidGrad_temp(2,8);
    fShapeSolidGrad(6,27) = fShapeSolidGrad_temp(2,9);
    fShapeSolidGrad(6,30) = fShapeSolidGrad_temp(2,10);
    fShapeSolidGrad(6,33) = fShapeSolidGrad_temp(2,11);
    fShapeSolidGrad(6,36) = fShapeSolidGrad_temp(2,12);
    fShapeSolidGrad(6,39) = fShapeSolidGrad_temp(2,13);
    fShapeSolidGrad(6,42) = fShapeSolidGrad_temp(2,14);
    fShapeSolidGrad(6,45) = fShapeSolidGrad_temp(2,15);
    fShapeSolidGrad(6,48) = fShapeSolidGrad_temp(2,16);
    fShapeSolidGrad(6,51) = fShapeSolidGrad_temp(2,17);
    fShapeSolidGrad(6,54) = fShapeSolidGrad_temp(2,18);
    fShapeSolidGrad(6,57) = fShapeSolidGrad_temp(2,19);
    fShapeSolidGrad(6,60) = fShapeSolidGrad_temp(2,20);
    fShapeSolidGrad(6,63) = fShapeSolidGrad_temp(2,21);
    fShapeSolidGrad(6,66) = fShapeSolidGrad_temp(2,22);
    fShapeSolidGrad(6,69) = fShapeSolidGrad_temp(2,23);
    fShapeSolidGrad(6,72) = fShapeSolidGrad_temp(2,24);
    fShapeSolidGrad(6,75) = fShapeSolidGrad_temp(2,25);
    fShapeSolidGrad(6,78) = fShapeSolidGrad_temp(2,26);
    
    fShapeSolidGrad(7,1) = fShapeSolidGrad_temp(2,0);
    fShapeSolidGrad(7,4) = fShapeSolidGrad_temp(2,1);
    fShapeSolidGrad(7,7) = fShapeSolidGrad_temp(2,2);
    fShapeSolidGrad(7,10) = fShapeSolidGrad_temp(2,3);
    fShapeSolidGrad(7,13) = fShapeSolidGrad_temp(2,4);
    fShapeSolidGrad(7,16) = fShapeSolidGrad_temp(2,5);
    fShapeSolidGrad(7,19) = fShapeSolidGrad_temp(2,6);
    fShapeSolidGrad(7,22) = fShapeSolidGrad_temp(2,7);
    fShapeSolidGrad(7,25) = fShapeSolidGrad_temp(2,8);
    fShapeSolidGrad(7,28) = fShapeSolidGrad_temp(2,9);
    fShapeSolidGrad(7,31) = fShapeSolidGrad_temp(2,10);
    fShapeSolidGrad(7,34) = fShapeSolidGrad_temp(2,11);
    fShapeSolidGrad(7,37) = fShapeSolidGrad_temp(2,12);
    fShapeSolidGrad(7,40) = fShapeSolidGrad_temp(2,13);
    fShapeSolidGrad(7,43) = fShapeSolidGrad_temp(2,14);
    fShapeSolidGrad(7,46) = fShapeSolidGrad_temp(2,15);
    fShapeSolidGrad(7,49) = fShapeSolidGrad_temp(2,16);
    fShapeSolidGrad(7,52) = fShapeSolidGrad_temp(2,17);
    fShapeSolidGrad(7,55) = fShapeSolidGrad_temp(2,18);
    fShapeSolidGrad(7,58) = fShapeSolidGrad_temp(2,19);
    fShapeSolidGrad(7,61) = fShapeSolidGrad_temp(2,20);
    fShapeSolidGrad(7,64) = fShapeSolidGrad_temp(2,21);
    fShapeSolidGrad(7,67) = fShapeSolidGrad_temp(2,22);
    fShapeSolidGrad(7,70) = fShapeSolidGrad_temp(2,23);
    fShapeSolidGrad(7,73) = fShapeSolidGrad_temp(2,24);
    fShapeSolidGrad(7,76) = fShapeSolidGrad_temp(2,25);
    fShapeSolidGrad(7,79) = fShapeSolidGrad_temp(2,26);
				
    fShapeSolidGrad(8,2) = fShapeSolidGrad_temp(2,0);
    fShapeSolidGrad(8,5) = fShapeSolidGrad_temp(2,1);
    fShapeSolidGrad(8,8) = fShapeSolidGrad_temp(2,2);
    fShapeSolidGrad(8,11) = fShapeSolidGrad_temp(2,3);
    fShapeSolidGrad(8,14) = fShapeSolidGrad_temp(2,4);
    fShapeSolidGrad(8,17) = fShapeSolidGrad_temp(2,5);
    fShapeSolidGrad(8,20) = fShapeSolidGrad_temp(2,6);
    fShapeSolidGrad(8,23) = fShapeSolidGrad_temp(2,7);
    fShapeSolidGrad(8,26) = fShapeSolidGrad_temp(2,8);
    fShapeSolidGrad(8,29) = fShapeSolidGrad_temp(2,9);
    fShapeSolidGrad(8,32) = fShapeSolidGrad_temp(2,10);
    fShapeSolidGrad(8,35) = fShapeSolidGrad_temp(2,11);
    fShapeSolidGrad(8,38) = fShapeSolidGrad_temp(2,12);
    fShapeSolidGrad(8,41) = fShapeSolidGrad_temp(2,13);
    fShapeSolidGrad(8,44) = fShapeSolidGrad_temp(2,14);
    fShapeSolidGrad(8,47) = fShapeSolidGrad_temp(2,15);
    fShapeSolidGrad(8,50) = fShapeSolidGrad_temp(2,16);
    fShapeSolidGrad(8,53) = fShapeSolidGrad_temp(2,17);
    fShapeSolidGrad(8,56) = fShapeSolidGrad_temp(2,18);
    fShapeSolidGrad(8,59) = fShapeSolidGrad_temp(2,19);
    fShapeSolidGrad(8,62) = fShapeSolidGrad_temp(2,20);
    fShapeSolidGrad(8,65) = fShapeSolidGrad_temp(2,21);
    fShapeSolidGrad(8,68) = fShapeSolidGrad_temp(2,22);
    fShapeSolidGrad(8,71) = fShapeSolidGrad_temp(2,23);
    fShapeSolidGrad(8,74) = fShapeSolidGrad_temp(2,24);
    fShapeSolidGrad(8,77) = fShapeSolidGrad_temp(2,25);
    fShapeSolidGrad(8,80) = fShapeSolidGrad_temp(2,26); 
}

void FSSolidFluidMixT::	Form_fluid_shape_functions(const double* &shapes_press_X)
{
    fShapeFluid = 0.0;
    fShapeFluid[0] = shapes_press_X[0];
    fShapeFluid[1] = shapes_press_X[1];
    fShapeFluid[2] = shapes_press_X[2];
    fShapeFluid[3] = shapes_press_X[3];
    fShapeFluid[4] = shapes_press_X[4];	
    fShapeFluid[5] = shapes_press_X[5];	
    fShapeFluid[6] = shapes_press_X[6];	
    fShapeFluid[7] = shapes_press_X[7];	
}

void FSSolidFluidMixT::	Form_deformation_gradient_tensor(void)
{
    fShapeSolidGrad.Multx(u_vec,fGRAD_disp_vector);
    fDeformation_Gradient(0,0) = fGRAD_disp_vector[0]+1.0;
    fDeformation_Gradient(0,1) = fGRAD_disp_vector[3]; 
    fDeformation_Gradient(0,2) = fGRAD_disp_vector[6];
    fDeformation_Gradient(1,0) = fGRAD_disp_vector[1];
    fDeformation_Gradient(1,1) = fGRAD_disp_vector[4]+1.0;  
    fDeformation_Gradient(1,2) = fGRAD_disp_vector[7];
    fDeformation_Gradient(2,0) = fGRAD_disp_vector[2];
    fDeformation_Gradient(2,1) = fGRAD_disp_vector[5];
    fDeformation_Gradient(2,2) = fGRAD_disp_vector[8]+1.0; 
}
void FSSolidFluidMixT::Form_GRAD_grad_transformation_matrix(void)
{
    fDefGradInv_GRAD_grad = 0.0;
    fDefGradInv_GRAD_grad(0,0) = fDeformation_Gradient_Inverse(0,0);
    fDefGradInv_GRAD_grad(0,3) = fDeformation_Gradient_Inverse(0,1);
    fDefGradInv_GRAD_grad(0,6) = fDeformation_Gradient_Inverse(0,2);
    fDefGradInv_GRAD_grad(1,1) = fDeformation_Gradient_Inverse(0,0);
    fDefGradInv_GRAD_grad(1,4) = fDeformation_Gradient_Inverse(0,1);
    fDefGradInv_GRAD_grad(1,7) = fDeformation_Gradient_Inverse(0,2);
    fDefGradInv_GRAD_grad(2,2) = fDeformation_Gradient_Inverse(0,0);
    fDefGradInv_GRAD_grad(2,5) = fDeformation_Gradient_Inverse(0,1);
    fDefGradInv_GRAD_grad(2,8) = fDeformation_Gradient_Inverse(0,2);
    fDefGradInv_GRAD_grad(3,0) = fDeformation_Gradient_Inverse(1,0);
    fDefGradInv_GRAD_grad(3,3) = fDeformation_Gradient_Inverse(1,1);
    fDefGradInv_GRAD_grad(3,6) = fDeformation_Gradient_Inverse(1,2);
    fDefGradInv_GRAD_grad(4,1) = fDeformation_Gradient_Inverse(1,0);
    fDefGradInv_GRAD_grad(4,4) = fDeformation_Gradient_Inverse(1,1);
    fDefGradInv_GRAD_grad(4,7) = fDeformation_Gradient_Inverse(1,2);
    fDefGradInv_GRAD_grad(5,2) = fDeformation_Gradient_Inverse(1,0);
    fDefGradInv_GRAD_grad(5,5) = fDeformation_Gradient_Inverse(1,1);
    fDefGradInv_GRAD_grad(5,8) = fDeformation_Gradient_Inverse(1,2);
    fDefGradInv_GRAD_grad(6,0) = fDeformation_Gradient_Inverse(2,0);
    fDefGradInv_GRAD_grad(6,3) = fDeformation_Gradient_Inverse(2,1);
    fDefGradInv_GRAD_grad(6,6) = fDeformation_Gradient_Inverse(2,2);
    fDefGradInv_GRAD_grad(7,1) = fDeformation_Gradient_Inverse(2,0);
    fDefGradInv_GRAD_grad(7,4) = fDeformation_Gradient_Inverse(2,1);
    fDefGradInv_GRAD_grad(7,7) = fDeformation_Gradient_Inverse(2,2);
    fDefGradInv_GRAD_grad(8,2) = fDeformation_Gradient_Inverse(2,0);
    fDefGradInv_GRAD_grad(8,5) = fDeformation_Gradient_Inverse(2,1);
    fDefGradInv_GRAD_grad(8,8) = fDeformation_Gradient_Inverse(2,2);
}

void FSSolidFluidMixT::Form_deformation_gradient_inv_vector(void)
{
    fDefGradInv_Vector[0] = fDeformation_Gradient_Inverse(0,0);
    fDefGradInv_Vector[1] = fDeformation_Gradient_Inverse(0,1);
    fDefGradInv_Vector[2] = fDeformation_Gradient_Inverse(0,2);
    fDefGradInv_Vector[3] = fDeformation_Gradient_Inverse(1,0);
    fDefGradInv_Vector[4] = fDeformation_Gradient_Inverse(1,1);
    fDefGradInv_Vector[5] = fDeformation_Gradient_Inverse(1,2);
    fDefGradInv_Vector[6] = fDeformation_Gradient_Inverse(2,0);
    fDefGradInv_Vector[7] = fDeformation_Gradient_Inverse(2,1);
    fDefGradInv_Vector[8] = fDeformation_Gradient_Inverse(2,2); 
    
}

void FSSolidFluidMixT::Form_effective_kirchhoff_stress_vector()
{
    fEffective_Kirchhoff_vector[0] = fEffective_Kirchhoff_tensor(0,0);
    fEffective_Kirchhoff_vector[1] = fEffective_Kirchhoff_tensor(1,0);
    fEffective_Kirchhoff_vector[2] = fEffective_Kirchhoff_tensor(2,0);
    fEffective_Kirchhoff_vector[3] = fEffective_Kirchhoff_tensor(0,1);
    fEffective_Kirchhoff_vector[4] = fEffective_Kirchhoff_tensor(1,1);
    fEffective_Kirchhoff_vector[5] = fEffective_Kirchhoff_tensor(2,1);
    fEffective_Kirchhoff_vector[6] = fEffective_Kirchhoff_tensor(0,2);
    fEffective_Kirchhoff_vector[7] = fEffective_Kirchhoff_tensor(1,2);
    fEffective_Kirchhoff_vector[8] = fEffective_Kirchhoff_tensor(2,2);
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
    fShapeSolidGrad_t(0,0) = fShapeSolidGrad_temp(0,0);
    fShapeSolidGrad_t(0,3) = fShapeSolidGrad_temp(0,1);
    fShapeSolidGrad_t(0,6) = fShapeSolidGrad_temp(0,2);
    fShapeSolidGrad_t(0,9) = fShapeSolidGrad_temp(0,3);
    fShapeSolidGrad_t(0,12) = fShapeSolidGrad_temp(0,4);
    fShapeSolidGrad_t(0,15) = fShapeSolidGrad_temp(0,5);
    fShapeSolidGrad_t(0,18) = fShapeSolidGrad_temp(0,6);
    fShapeSolidGrad_t(0,21) = fShapeSolidGrad_temp(0,7);
    fShapeSolidGrad_t(0,24) = fShapeSolidGrad_temp(0,8);
    fShapeSolidGrad_t(0,27) = fShapeSolidGrad_temp(0,9);
    fShapeSolidGrad_t(0,30) = fShapeSolidGrad_temp(0,10);
    fShapeSolidGrad_t(0,33) = fShapeSolidGrad_temp(0,11);
    fShapeSolidGrad_t(0,36) = fShapeSolidGrad_temp(0,12);
    fShapeSolidGrad_t(0,39) = fShapeSolidGrad_temp(0,13);
    fShapeSolidGrad_t(0,42) = fShapeSolidGrad_temp(0,14);
    fShapeSolidGrad_t(0,45) = fShapeSolidGrad_temp(0,15);
    fShapeSolidGrad_t(0,48) = fShapeSolidGrad_temp(0,16);
    fShapeSolidGrad_t(0,51) = fShapeSolidGrad_temp(0,17);
    fShapeSolidGrad_t(0,54) = fShapeSolidGrad_temp(0,18);
    fShapeSolidGrad_t(0,57) = fShapeSolidGrad_temp(0,19);
    fShapeSolidGrad_t(0,60) = fShapeSolidGrad_temp(0,20);
    fShapeSolidGrad_t(0,63) = fShapeSolidGrad_temp(0,21);
    fShapeSolidGrad_t(0,66) = fShapeSolidGrad_temp(0,22);
    fShapeSolidGrad_t(0,69) = fShapeSolidGrad_temp(0,23);
    fShapeSolidGrad_t(0,72) = fShapeSolidGrad_temp(0,24);
    fShapeSolidGrad_t(0,75) = fShapeSolidGrad_temp(0,25);
    fShapeSolidGrad_t(0,78) = fShapeSolidGrad_temp(0,26);
    
    fShapeSolidGrad_t(1,0) = fShapeSolidGrad_temp(1,0);
    fShapeSolidGrad_t(1,3) = fShapeSolidGrad_temp(1,1);
    fShapeSolidGrad_t(1,6) = fShapeSolidGrad_temp(1,2);
    fShapeSolidGrad_t(1,9) = fShapeSolidGrad_temp(1,3);
    fShapeSolidGrad_t(1,12) = fShapeSolidGrad_temp(1,4);
    fShapeSolidGrad_t(1,15) = fShapeSolidGrad_temp(1,5);
    fShapeSolidGrad_t(1,18) = fShapeSolidGrad_temp(1,6);
    fShapeSolidGrad_t(1,21) = fShapeSolidGrad_temp(1,7);
    fShapeSolidGrad_t(1,24) = fShapeSolidGrad_temp(1,8);
    fShapeSolidGrad_t(1,27) = fShapeSolidGrad_temp(1,9);
    fShapeSolidGrad_t(1,30) = fShapeSolidGrad_temp(1,10);
    fShapeSolidGrad_t(1,33) = fShapeSolidGrad_temp(1,11);
    fShapeSolidGrad_t(1,36) = fShapeSolidGrad_temp(1,12);
    fShapeSolidGrad_t(1,39) = fShapeSolidGrad_temp(1,13);
    fShapeSolidGrad_t(1,42) = fShapeSolidGrad_temp(1,14);
    fShapeSolidGrad_t(1,45) = fShapeSolidGrad_temp(1,15);
    fShapeSolidGrad_t(1,48) = fShapeSolidGrad_temp(1,16);
    fShapeSolidGrad_t(1,51) = fShapeSolidGrad_temp(1,17);
    fShapeSolidGrad_t(1,54) = fShapeSolidGrad_temp(1,18);
    fShapeSolidGrad_t(1,57) = fShapeSolidGrad_temp(1,19);
    fShapeSolidGrad_t(1,60) = fShapeSolidGrad_temp(1,20);
    fShapeSolidGrad_t(1,63) = fShapeSolidGrad_temp(1,21);
    fShapeSolidGrad_t(1,66) = fShapeSolidGrad_temp(1,22);
    fShapeSolidGrad_t(1,69) = fShapeSolidGrad_temp(1,23);
    fShapeSolidGrad_t(1,72) = fShapeSolidGrad_temp(1,24);
    fShapeSolidGrad_t(1,75) = fShapeSolidGrad_temp(1,25);
    fShapeSolidGrad_t(1,78) = fShapeSolidGrad_temp(1,26);
    

    fShapeSolidGrad_t(2,0) = fShapeSolidGrad_temp(2,0);
    fShapeSolidGrad_t(2,3) = fShapeSolidGrad_temp(2,1);
    fShapeSolidGrad_t(2,6) = fShapeSolidGrad_temp(2,2);
    fShapeSolidGrad_t(2,9) = fShapeSolidGrad_temp(2,3);
    fShapeSolidGrad_t(2,12) = fShapeSolidGrad_temp(2,4);
    fShapeSolidGrad_t(2,15) = fShapeSolidGrad_temp(2,5);
    fShapeSolidGrad_t(2,18) = fShapeSolidGrad_temp(2,6);
    fShapeSolidGrad_t(2,21) = fShapeSolidGrad_temp(2,7);
    fShapeSolidGrad_t(2,24) = fShapeSolidGrad_temp(2,8);
    fShapeSolidGrad_t(2,27) = fShapeSolidGrad_temp(2,9);
    fShapeSolidGrad_t(2,30) = fShapeSolidGrad_temp(2,10);
    fShapeSolidGrad_t(2,33) = fShapeSolidGrad_temp(2,11);
    fShapeSolidGrad_t(2,36) = fShapeSolidGrad_temp(2,12);
    fShapeSolidGrad_t(2,39) = fShapeSolidGrad_temp(2,13);
    fShapeSolidGrad_t(2,42) = fShapeSolidGrad_temp(2,14);
    fShapeSolidGrad_t(2,45) = fShapeSolidGrad_temp(2,15);
    fShapeSolidGrad_t(2,48) = fShapeSolidGrad_temp(2,16);
    fShapeSolidGrad_t(2,51) = fShapeSolidGrad_temp(2,17);
    fShapeSolidGrad_t(2,54) = fShapeSolidGrad_temp(2,18);
    fShapeSolidGrad_t(2,57) = fShapeSolidGrad_temp(2,19);
    fShapeSolidGrad_t(2,60) = fShapeSolidGrad_temp(2,20);
    fShapeSolidGrad_t(2,63) = fShapeSolidGrad_temp(2,21);
    fShapeSolidGrad_t(2,66) = fShapeSolidGrad_temp(2,22);
    fShapeSolidGrad_t(2,69) = fShapeSolidGrad_temp(2,23);
    fShapeSolidGrad_t(2,72) = fShapeSolidGrad_temp(2,24);
    fShapeSolidGrad_t(2,75) = fShapeSolidGrad_temp(2,25);
    fShapeSolidGrad_t(2,78) = fShapeSolidGrad_temp(2,26);
    
    fShapeSolidGrad_t(3,1) = fShapeSolidGrad_temp(0,0);
    fShapeSolidGrad_t(3,4) = fShapeSolidGrad_temp(0,1);
    fShapeSolidGrad_t(3,7) = fShapeSolidGrad_temp(0,2);
    fShapeSolidGrad_t(3,10) = fShapeSolidGrad_temp(0,3);
    fShapeSolidGrad_t(3,13) = fShapeSolidGrad_temp(0,4);
    fShapeSolidGrad_t(3,16) = fShapeSolidGrad_temp(0,5);
    fShapeSolidGrad_t(3,19) = fShapeSolidGrad_temp(0,6);
    fShapeSolidGrad_t(3,22) = fShapeSolidGrad_temp(0,7);
    fShapeSolidGrad_t(3,25) = fShapeSolidGrad_temp(0,8);
    fShapeSolidGrad_t(3,28) = fShapeSolidGrad_temp(0,9);
    fShapeSolidGrad_t(3,31) = fShapeSolidGrad_temp(0,10);
    fShapeSolidGrad_t(3,34) = fShapeSolidGrad_temp(0,11);
    fShapeSolidGrad_t(3,37) = fShapeSolidGrad_temp(0,12);
    fShapeSolidGrad_t(3,40) = fShapeSolidGrad_temp(0,13);
    fShapeSolidGrad_t(3,43) = fShapeSolidGrad_temp(0,14);
    fShapeSolidGrad_t(3,46) = fShapeSolidGrad_temp(0,15);
    fShapeSolidGrad_t(3,49) = fShapeSolidGrad_temp(0,16);
    fShapeSolidGrad_t(3,52) = fShapeSolidGrad_temp(0,17);
    fShapeSolidGrad_t(3,55) = fShapeSolidGrad_temp(0,18);
    fShapeSolidGrad_t(3,58) = fShapeSolidGrad_temp(0,19);
    fShapeSolidGrad_t(3,61) = fShapeSolidGrad_temp(0,20);
    fShapeSolidGrad_t(3,64) = fShapeSolidGrad_temp(0,21);
    fShapeSolidGrad_t(3,67) = fShapeSolidGrad_temp(0,22);
    fShapeSolidGrad_t(3,70) = fShapeSolidGrad_temp(0,23);
    fShapeSolidGrad_t(3,73) = fShapeSolidGrad_temp(0,24);
    fShapeSolidGrad_t(3,76) = fShapeSolidGrad_temp(0,25);
    fShapeSolidGrad_t(3,79) = fShapeSolidGrad_temp(0,26);
    
    fShapeSolidGrad_t(4,1) = fShapeSolidGrad_temp(1,0);
    fShapeSolidGrad_t(4,4) = fShapeSolidGrad_temp(1,1);
    fShapeSolidGrad_t(4,7) = fShapeSolidGrad_temp(1,2);
    fShapeSolidGrad_t(4,10) = fShapeSolidGrad_temp(1,3);
    fShapeSolidGrad_t(4,13) = fShapeSolidGrad_temp(1,4);
    fShapeSolidGrad_t(4,16) = fShapeSolidGrad_temp(1,5);
    fShapeSolidGrad_t(4,19) = fShapeSolidGrad_temp(1,6);
    fShapeSolidGrad_t(4,22) = fShapeSolidGrad_temp(1,7);
    fShapeSolidGrad_t(4,25) = fShapeSolidGrad_temp(1,8);
    fShapeSolidGrad_t(4,28) = fShapeSolidGrad_temp(1,9);
    fShapeSolidGrad_t(4,31) = fShapeSolidGrad_temp(1,10);
    fShapeSolidGrad_t(4,34) = fShapeSolidGrad_temp(1,11);
    fShapeSolidGrad_t(4,37) = fShapeSolidGrad_temp(1,12);
    fShapeSolidGrad_t(4,40) = fShapeSolidGrad_temp(1,13);
    fShapeSolidGrad_t(4,43) = fShapeSolidGrad_temp(1,14);
    fShapeSolidGrad_t(4,46) = fShapeSolidGrad_temp(1,15);
    fShapeSolidGrad_t(4,49) = fShapeSolidGrad_temp(1,16);
    fShapeSolidGrad_t(4,52) = fShapeSolidGrad_temp(1,17);
    fShapeSolidGrad_t(4,55) = fShapeSolidGrad_temp(1,18);
    fShapeSolidGrad_t(4,58) = fShapeSolidGrad_temp(1,19);
    fShapeSolidGrad_t(4,61) = fShapeSolidGrad_temp(1,20);
    fShapeSolidGrad_t(4,64) = fShapeSolidGrad_temp(1,21);
    fShapeSolidGrad_t(4,67) = fShapeSolidGrad_temp(1,22);
    fShapeSolidGrad_t(4,70) = fShapeSolidGrad_temp(1,23);
    fShapeSolidGrad_t(4,73) = fShapeSolidGrad_temp(1,24);
    fShapeSolidGrad_t(4,76) = fShapeSolidGrad_temp(1,25);
    fShapeSolidGrad_t(4,79) = fShapeSolidGrad_temp(1,26);
    
    fShapeSolidGrad_t(5,1) = fShapeSolidGrad_temp(2,0);
    fShapeSolidGrad_t(5,4) = fShapeSolidGrad_temp(2,1);
    fShapeSolidGrad_t(5,7) = fShapeSolidGrad_temp(2,2);
    fShapeSolidGrad_t(5,10) = fShapeSolidGrad_temp(2,3);
    fShapeSolidGrad_t(5,13) = fShapeSolidGrad_temp(2,4);
    fShapeSolidGrad_t(5,16) = fShapeSolidGrad_temp(2,5);
    fShapeSolidGrad_t(5,19) = fShapeSolidGrad_temp(2,6);
    fShapeSolidGrad_t(5,22) = fShapeSolidGrad_temp(2,7);
    fShapeSolidGrad_t(5,25) = fShapeSolidGrad_temp(2,8);
    fShapeSolidGrad_t(5,28) = fShapeSolidGrad_temp(2,9);
    fShapeSolidGrad_t(5,31) = fShapeSolidGrad_temp(2,10);
    fShapeSolidGrad_t(5,34) = fShapeSolidGrad_temp(2,11);
    fShapeSolidGrad_t(5,37) = fShapeSolidGrad_temp(2,12);
    fShapeSolidGrad_t(5,40) = fShapeSolidGrad_temp(2,13);
    fShapeSolidGrad_t(5,43) = fShapeSolidGrad_temp(2,14);
    fShapeSolidGrad_t(5,46) = fShapeSolidGrad_temp(2,15);
    fShapeSolidGrad_t(5,49) = fShapeSolidGrad_temp(2,16);
    fShapeSolidGrad_t(5,52) = fShapeSolidGrad_temp(2,17);
    fShapeSolidGrad_t(5,55) = fShapeSolidGrad_temp(2,18);
    fShapeSolidGrad_t(5,58) = fShapeSolidGrad_temp(2,19);
    fShapeSolidGrad_t(5,61) = fShapeSolidGrad_temp(2,20);
    fShapeSolidGrad_t(5,64) = fShapeSolidGrad_temp(2,21);
    fShapeSolidGrad_t(5,67) = fShapeSolidGrad_temp(2,22);
    fShapeSolidGrad_t(5,70) = fShapeSolidGrad_temp(2,23);
    fShapeSolidGrad_t(5,73) = fShapeSolidGrad_temp(2,24);
    fShapeSolidGrad_t(5,76) = fShapeSolidGrad_temp(2,25);
    fShapeSolidGrad_t(5,79) = fShapeSolidGrad_temp(2,26);
    
    fShapeSolidGrad_t(6,2) = fShapeSolidGrad_temp(0,0);
    fShapeSolidGrad_t(6,5) = fShapeSolidGrad_temp(0,1);
    fShapeSolidGrad_t(6,8) = fShapeSolidGrad_temp(0,2);
    fShapeSolidGrad_t(6,11) = fShapeSolidGrad_temp(0,3);
    fShapeSolidGrad_t(6,14) = fShapeSolidGrad_temp(0,4);
    fShapeSolidGrad_t(6,17) = fShapeSolidGrad_temp(0,5);
    fShapeSolidGrad_t(6,20) = fShapeSolidGrad_temp(0,6);
    fShapeSolidGrad_t(6,23) = fShapeSolidGrad_temp(0,7);
    fShapeSolidGrad_t(6,26) = fShapeSolidGrad_temp(0,8);
    fShapeSolidGrad_t(6,29) = fShapeSolidGrad_temp(0,9);
    fShapeSolidGrad_t(6,32) = fShapeSolidGrad_temp(0,10);
    fShapeSolidGrad_t(6,35) = fShapeSolidGrad_temp(0,11);
    fShapeSolidGrad_t(6,38) = fShapeSolidGrad_temp(0,12);
    fShapeSolidGrad_t(6,41) = fShapeSolidGrad_temp(0,13);
    fShapeSolidGrad_t(6,44) = fShapeSolidGrad_temp(0,14);
    fShapeSolidGrad_t(6,47) = fShapeSolidGrad_temp(0,15);
    fShapeSolidGrad_t(6,50) = fShapeSolidGrad_temp(0,16);
    fShapeSolidGrad_t(6,53) = fShapeSolidGrad_temp(0,17);
    fShapeSolidGrad_t(6,56) = fShapeSolidGrad_temp(0,18);
    fShapeSolidGrad_t(6,59) = fShapeSolidGrad_temp(0,19);
    fShapeSolidGrad_t(6,62) = fShapeSolidGrad_temp(0,20);
    fShapeSolidGrad_t(6,65) = fShapeSolidGrad_temp(0,21);
    fShapeSolidGrad_t(6,68) = fShapeSolidGrad_temp(0,22);
    fShapeSolidGrad_t(6,71) = fShapeSolidGrad_temp(0,23);
    fShapeSolidGrad_t(6,74) = fShapeSolidGrad_temp(0,24);
    fShapeSolidGrad_t(6,77) = fShapeSolidGrad_temp(0,25);
    fShapeSolidGrad_t(6,80) = fShapeSolidGrad_temp(0,26);
    
    fShapeSolidGrad_t(7,2) = fShapeSolidGrad_temp(1,0);
    fShapeSolidGrad_t(7,5) = fShapeSolidGrad_temp(1,1);
    fShapeSolidGrad_t(7,8) = fShapeSolidGrad_temp(1,2);
    fShapeSolidGrad_t(7,11) = fShapeSolidGrad_temp(1,3);
    fShapeSolidGrad_t(7,14) = fShapeSolidGrad_temp(1,4);
    fShapeSolidGrad_t(7,17) = fShapeSolidGrad_temp(1,5);
    fShapeSolidGrad_t(7,20) = fShapeSolidGrad_temp(1,6);
    fShapeSolidGrad_t(7,23) = fShapeSolidGrad_temp(1,7);
    fShapeSolidGrad_t(7,26) = fShapeSolidGrad_temp(1,8);
    fShapeSolidGrad_t(7,29) = fShapeSolidGrad_temp(1,9);
    fShapeSolidGrad_t(7,32) = fShapeSolidGrad_temp(1,10);
    fShapeSolidGrad_t(7,35) = fShapeSolidGrad_temp(1,11);
    fShapeSolidGrad_t(7,38) = fShapeSolidGrad_temp(1,12);
    fShapeSolidGrad_t(7,41) = fShapeSolidGrad_temp(1,13);
    fShapeSolidGrad_t(7,44) = fShapeSolidGrad_temp(1,14);
    fShapeSolidGrad_t(7,47) = fShapeSolidGrad_temp(1,15);
    fShapeSolidGrad_t(7,50) = fShapeSolidGrad_temp(1,16);
    fShapeSolidGrad_t(7,53) = fShapeSolidGrad_temp(1,17);
    fShapeSolidGrad_t(7,56) = fShapeSolidGrad_temp(1,18);
    fShapeSolidGrad_t(7,59) = fShapeSolidGrad_temp(1,19);
    fShapeSolidGrad_t(7,62) = fShapeSolidGrad_temp(1,20);
    fShapeSolidGrad_t(7,65) = fShapeSolidGrad_temp(1,21);
    fShapeSolidGrad_t(7,68) = fShapeSolidGrad_temp(1,22);
    fShapeSolidGrad_t(7,71) = fShapeSolidGrad_temp(1,23);
    fShapeSolidGrad_t(7,74) = fShapeSolidGrad_temp(1,24);
    fShapeSolidGrad_t(7,77) = fShapeSolidGrad_temp(1,25);
    fShapeSolidGrad_t(7,80) = fShapeSolidGrad_temp(1,26);
    
    fShapeSolidGrad_t(8,2) = fShapeSolidGrad_temp(2,0);
    fShapeSolidGrad_t(8,5) = fShapeSolidGrad_temp(2,1);
    fShapeSolidGrad_t(8,8) = fShapeSolidGrad_temp(2,2);
    fShapeSolidGrad_t(8,11) = fShapeSolidGrad_temp(2,3);
    fShapeSolidGrad_t(8,14) = fShapeSolidGrad_temp(2,4);
    fShapeSolidGrad_t(8,17) = fShapeSolidGrad_temp(2,5);
    fShapeSolidGrad_t(8,20) = fShapeSolidGrad_temp(2,6);
    fShapeSolidGrad_t(8,23) = fShapeSolidGrad_temp(2,7);
    fShapeSolidGrad_t(8,26) = fShapeSolidGrad_temp(2,8);
    fShapeSolidGrad_t(8,29) = fShapeSolidGrad_temp(2,9);
    fShapeSolidGrad_t(8,32) = fShapeSolidGrad_temp(2,10);
    fShapeSolidGrad_t(8,35) = fShapeSolidGrad_temp(2,11);
    fShapeSolidGrad_t(8,38) = fShapeSolidGrad_temp(2,12);
    fShapeSolidGrad_t(8,41) = fShapeSolidGrad_temp(2,13);
    fShapeSolidGrad_t(8,44) = fShapeSolidGrad_temp(2,14);
    fShapeSolidGrad_t(8,47) = fShapeSolidGrad_temp(2,15);
    fShapeSolidGrad_t(8,50) = fShapeSolidGrad_temp(2,16);
    fShapeSolidGrad_t(8,53) = fShapeSolidGrad_temp(2,17);
    fShapeSolidGrad_t(8,56) = fShapeSolidGrad_temp(2,18);
    fShapeSolidGrad_t(8,59) = fShapeSolidGrad_temp(2,19);
    fShapeSolidGrad_t(8,62) = fShapeSolidGrad_temp(2,20);
    fShapeSolidGrad_t(8,65) = fShapeSolidGrad_temp(2,21);
    fShapeSolidGrad_t(8,68) = fShapeSolidGrad_temp(2,22);
    fShapeSolidGrad_t(8,71) = fShapeSolidGrad_temp(2,23);
    fShapeSolidGrad_t(8,74) = fShapeSolidGrad_temp(2,24);
    fShapeSolidGrad_t(8,77) = fShapeSolidGrad_temp(2,25);
    fShapeSolidGrad_t(8,80) = fShapeSolidGrad_temp(2,26); 
}

void FSSolidFluidMixT::Form_Im_temp_matrix()
{
    fIm_temp_matrix = 0.0;
    fIm_temp_matrix(0,0) = fEffective_Kirchhoff_tensor(0,0);
    fIm_temp_matrix(1,0) = fEffective_Kirchhoff_tensor(1,0);
    fIm_temp_matrix(2,0) = fEffective_Kirchhoff_tensor(2,0);
    
    fIm_temp_matrix(3,1) = fEffective_Kirchhoff_tensor(0,0);
    fIm_temp_matrix(4,1) = fEffective_Kirchhoff_tensor(1,0);
    fIm_temp_matrix(5,1) = fEffective_Kirchhoff_tensor(2,0);
    
    fIm_temp_matrix(6,2) = fEffective_Kirchhoff_tensor(0,0);
    fIm_temp_matrix(7,2) = fEffective_Kirchhoff_tensor(1,0);
    fIm_temp_matrix(8,2) = fEffective_Kirchhoff_tensor(2,0);
    
    fIm_temp_matrix(0,3) = fEffective_Kirchhoff_tensor(0,1);
    fIm_temp_matrix(1,3) = fEffective_Kirchhoff_tensor(1,1);
    fIm_temp_matrix(2,3) = fEffective_Kirchhoff_tensor(2,1);
    
    fIm_temp_matrix(3,4) = fEffective_Kirchhoff_tensor(0,1);
    fIm_temp_matrix(4,4) = fEffective_Kirchhoff_tensor(1,1);
    fIm_temp_matrix(5,4) = fEffective_Kirchhoff_tensor(2,1);
    
    fIm_temp_matrix(6,5) = fEffective_Kirchhoff_tensor(0,1);
    fIm_temp_matrix(7,5) = fEffective_Kirchhoff_tensor(1,1);
    fIm_temp_matrix(8,5) = fEffective_Kirchhoff_tensor(2,1);
    
    fIm_temp_matrix(0,6) = fEffective_Kirchhoff_tensor(0,2);
    fIm_temp_matrix(1,6) = fEffective_Kirchhoff_tensor(1,2);
    fIm_temp_matrix(2,6) = fEffective_Kirchhoff_tensor(2,2);
    
    fIm_temp_matrix(3,7) = fEffective_Kirchhoff_tensor(0,2);
    fIm_temp_matrix(4,7) = fEffective_Kirchhoff_tensor(1,2);
    fIm_temp_matrix(5,7) = fEffective_Kirchhoff_tensor(2,2);
    
    fIm_temp_matrix(6,8) = fEffective_Kirchhoff_tensor(0,2);
    fIm_temp_matrix(7,8) = fEffective_Kirchhoff_tensor(1,2);
    fIm_temp_matrix(8,8) = fEffective_Kirchhoff_tensor(2,2);
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

void FSSolidFluidMixT::Form_Eth_temp_matrix()
{
    fEth_temp_matrix(0,0) = 0.0;
    fEth_temp_matrix(0,0) = fLeft_Cauchy_Green_tensor(0,0);
    fEth_temp_matrix(1,0) = fLeft_Cauchy_Green_tensor(1,0);
    fEth_temp_matrix(2,0) = fLeft_Cauchy_Green_tensor(2,0);
    
    fEth_temp_matrix(3,1) = fLeft_Cauchy_Green_tensor(0,0);
    fEth_temp_matrix(4,1) = fLeft_Cauchy_Green_tensor(1,0);
    fEth_temp_matrix(5,1) = fLeft_Cauchy_Green_tensor(2,0);
    
    fEth_temp_matrix(6,2) = fLeft_Cauchy_Green_tensor(0,0);
    fEth_temp_matrix(7,2) = fLeft_Cauchy_Green_tensor(1,0);
    fEth_temp_matrix(8,2) = fLeft_Cauchy_Green_tensor(2,0);
    
    fEth_temp_matrix(0,3) = fLeft_Cauchy_Green_tensor(0,1);
    fEth_temp_matrix(1,3) = fLeft_Cauchy_Green_tensor(1,1);
    fEth_temp_matrix(2,3) = fLeft_Cauchy_Green_tensor(2,1);

    fEth_temp_matrix(3,4) = fLeft_Cauchy_Green_tensor(0,1);
    fEth_temp_matrix(4,4) = fLeft_Cauchy_Green_tensor(1,1);
    fEth_temp_matrix(5,4) = fLeft_Cauchy_Green_tensor(2,1);
    
    fEth_temp_matrix(6,5) = fLeft_Cauchy_Green_tensor(0,1);
    fEth_temp_matrix(7,5) = fLeft_Cauchy_Green_tensor(1,1);
    fEth_temp_matrix(8,5) = fLeft_Cauchy_Green_tensor(2,1);
    
    fEth_temp_matrix(0,6) = fLeft_Cauchy_Green_tensor(0,2);
    fEth_temp_matrix(1,6) = fLeft_Cauchy_Green_tensor(1,2);
    fEth_temp_matrix(2,6) = fLeft_Cauchy_Green_tensor(2,2);
    
    fEth_temp_matrix(3,7) = fLeft_Cauchy_Green_tensor(0,2);
    fEth_temp_matrix(4,7) = fLeft_Cauchy_Green_tensor(1,2);
    fEth_temp_matrix(5,7) = fLeft_Cauchy_Green_tensor(2,2);
    
    fEth_temp_matrix(6,8) = fLeft_Cauchy_Green_tensor(0,2);
    fEth_temp_matrix(7,8) = fLeft_Cauchy_Green_tensor(1,2);
    fEth_temp_matrix(8,8) = fLeft_Cauchy_Green_tensor(2,2);
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
	for (int k; k<3; k++)
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
