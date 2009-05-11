/* Id: FSMicromorphic3DT.cpp,v 1.6 2006/10/10 19:55:23 regueiro Exp $ */
#include "FSMicromorphic3DT.h"

#include "OutputSetT.h"
#include "ParameterContainerT.h"
#include "CommunicatorT.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
FSMicromorphic3DT::FSMicromorphic3DT(const ElementSupportT& support):
    ElementBaseT(support), //pass the solid displacement field to the base class
    u(LocalArrayT::kDisp),
    u_dot(LocalArrayT::kVel),
    u_dot_n(LocalArrayT::kLastVel),
    u_dotdot(LocalArrayT::kAcc),
    u_dotdot_n(LocalArrayT::kLastAcc),
    u_n(LocalArrayT::kLastDisp),
    Phi(LocalArrayT::kDisp),
    Phi_dot(LocalArrayT::kVel),
    Phi_dot_n(LocalArrayT::kLastVel),
    Phi_dotdot(LocalArrayT::kAcc),
    Phi_dotdot_n(LocalArrayT::kLastAcc),
    Phi_n(LocalArrayT::kLastDisp),
    fInitCoords_displ(LocalArrayT::kInitCoords),
    fCurrCoords_displ(LocalArrayT::kCurrCoords),
    fInitCoords_micro(LocalArrayT::kInitCoords),
    fCurrCoords_micro(LocalArrayT::kCurrCoords),
    fTractionBCSet(0),
    fDispl(NULL),
    fMicro(NULL),
    fShapes_displ(NULL),
    fShapes_micro(NULL),
    fKdd(ElementMatrixT::kNonSymmetric),
    fKdphi(ElementMatrixT::kNonSymmetric),
    fKphid(ElementMatrixT::kNonSymmetric),
    fKphiphi(ElementMatrixT::kNonSymmetric),
    bStep_Complete(0)
{
    SetName("micromorphic_FS_3D");
}

/* destructor */
FSMicromorphic3DT::~FSMicromorphic3DT(void) 
{  
    delete fShapes_displ;
    delete fShapes_micro;
}


void FSMicromorphic3DT::Echo_Input_Data(void) 
{

    cout << "#######################################################" << endl; 
    cout << "############### ECHO FSMicromorphic3D DATA #########################" << endl; 
    cout << "#######################################################" << endl; 

    //################## material parameters ##################
    cout << "iConstitutiveModelType " 				<< iConstitutiveModelType 	<< endl; 

    //-- Elasticity parameters for solid
    cout << "fMaterial_Params[kMu] "  				<< fMaterial_Params[kMu] 	 << endl;
    cout << "fMaterial_Params[kLambda] "  			<< fMaterial_Params[kLambda] << endl;
	
}


//---------------------------------------------------------------------

void FSMicromorphic3DT::RHSDriver(void)	
{
    int curr_group = ElementSupport().CurrentGroup();

    /* traction boundary conditions acting on displacement equations */
    if (curr_group == fDispl->Group()) 
	ApplyTractionBC();

    /* choose solution method */
    if (fDispl->Group() == fMicro->Group())
		RHSDriver_monolithic();
    else
		RHSDriver_staggered();
}
//---------------------------------------------------------------------

void FSMicromorphic3DT::Equations(AutoArrayT<const iArray2DT*>& eq_d,
				 AutoArrayT<const RaggedArray2DT<int>*>& eq_phi)
{
    /* doing monolithic solution */
    if (fDispl->Group() == fMicro->Group())
    {
		int ndof_micro = fMicro->NumDOF();
		int ndof_displ = fDispl->NumDOF();
		
		/* loop over connectivity blocks */
		fEqnos_displ.Dimension(fEqnos.Length());
		fEqnos_micro.Dimension(fEqnos.Length());
		for (int i = 0; i < fEqnos.Length(); i++)
		{
		    /* connectivities */
		    const iArray2DT& connects_displ = *(fConnectivities_displ[i]);
		    const iArray2DT& connects_micro = *(fConnectivities_micro[i]);
		    int nel = connects_displ.MajorDim();
			
		    /* dimension */ 
		    fEqnos[i].Dimension(nel, n_en_displ*ndof_displ + n_en_micro*ndof_micro);
		    iArray2DT& displ_eq = fEqnos_displ[i];
		    iArray2DT& micro_eq = fEqnos_micro[i];
		    displ_eq.Dimension(nel, n_en_displ*ndof_displ);
		    micro_eq.Dimension(nel, n_en_micro*ndof_micro);
				
		    /* get equation numbers */
		    fDispl->SetLocalEqnos(connects_displ, displ_eq);
		    fMicro->SetLocalEqnos(connects_micro, micro_eq);

		    /* write into one array */
		    fEqnos[i].BlockColumnCopyAt(displ_eq, 0);
		    fEqnos[i].BlockColumnCopyAt(micro_eq, displ_eq.MinorDim());

		    /* add to list of equation numbers */
		    eq_d.Append(&fEqnos[i]);
		}
		
		/* reset pointers to element cards */
		SetElementCards(fBlockData, fConnectivities_displ, fEqnos_displ, fElementCards_displ);
		SetElementCards(fBlockData, fConnectivities_micro, fEqnos_micro, fElementCards_micro);
    }
    else
	/* doing staggered */
    {
#pragma message("initialization for staggered solution needs to be corrected")
	
		/* ElementBaseT handles equation array for displacements */
		if (ElementSupport().CurrentGroup() == fDispl->Group())
		    ElementBaseT::Equations(eq_d, eq_phi);

		/* micro-displacement-gradient equation */
		if (ElementSupport().CurrentGroup() == fMicro->Group())
		{
		    /* collect local equation numbers */
		    //fPress.SetLocalEqnos(fConnectivities_micro, fEqnos_micro);
			
		    //eq_d.Append(&fEqnos_micro);
		}
    }
	
    /* get the equation number for the nodes on the faces */
    /*
    for (int i = 0; i < fMicroFaceEqnos.Length(); i++)
    {
		iArray2DT& faces = fMicroFaces[i];
		iArray2DT& eqnos = fMicroFaceEqnos[i];
		eqnos.Dimension(faces.MajorDim(), faces.MajorDim()*fDispl->NumDOF());
		
		fDispl->SetLocalEqnos(faces, eqnos);
    }
    */
}


//---------------------------------------------------------------------

void FSMicromorphic3DT::LHSDriver(GlobalT::SystemTypeT)
{
/** Everything done in RHSDriver for efficiency */
//cout << "############### In LHS Driver ############### \n";
}

//---------------------------------------------------------------------

void FSMicromorphic3DT::Select_Equations (const int &iBalLinChoice, const int &iBalFirstMomMomChoice )
{
    /** Choices for Linear Momentum Balance Equation */

    switch ( iBalLinChoice )	{

    default :
	cout << "FSMicromorphic3DT::Select_Equations() .. currently only one linear momentum balance for micromorphic continuum \n";
	break;
    }

    /** Choices for First Moment of Momentum Balance Equation */

    switch ( iBalFirstMomMomChoice )	{

    default :
	cout << "FSMicromorphic3DT::Select_Equations() .. currently only one first moment of momentum balance equation for micromorphic continuum \n";
	break;
    }

}

//---------------------------------------------------------------------

/* return true if the element contributes to the solution of the
 * given group. */
bool FSMicromorphic3DT::InGroup(int group) const
{
    return group == fDispl->Group() || group == fMicro->Group();
}

//---------------------------------------------------------------------


/* initialize/finalize step */
void FSMicromorphic3DT::InitStep(void)
{
	/* inherited */
	ElementBaseT::InitStep();
}


/* close current time increment */
void FSMicromorphic3DT::CloseStep(void)
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
	    out_variable_all.Alias(fNumIP_displ, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
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
	int num_node_output = fDispl->NumDOF() + fMicro->NumDOF() + knumstrain + knumstress + knum_d_state;
	dArray2DT n_values(nodes_used.Length(), num_node_output);
	
	/* collect nodal values */
	const dArray2DT& fPhi = (*fMicro)[0];
	const dArray2DT& fU = (*fDispl)[0];
	for (int i = 0; i < nodes_used.Length(); i++)
	{
	    int node = nodes_used[i];
	    double* row = n_values(i);
	    for (int j = 0; j < fPhi.MinorDim(); j++)
		*row++ = fPhi(node,j);
	    
	    for (int j = 0; j < fU.MinorDim(); j++)
		*row++ = fU(node,j);
	    
	    double* p_stress = extrap_values(i);
	    for (int j = 0; j < (knumstrain+knumstress+knum_d_state); j++)
		*row++ = p_stress[j];
	}

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, fIPVariable);
	
    }

	//will not need this for quasi-static micromorphic
    /* zero first derivative of fields which are created at time=0 during calculating geostatic equilibrium (Trapezoidal rule) */    
    /*
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==1)
    {
		FieldT* fmicro = const_cast <FieldT*> (fMicro);
		FieldT* fdispl = const_cast <FieldT*> (fDispl);

		(*fdispl)[1] = 0;
		(*fmicro)[1] = 0;
    }   
    */

	//will not need this for quasi-static micromorphic
    /* zero second derivative of fields which are created at time=0 during calculating geostatic equilibrium(Newmark method) */
    /*
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==2)
    {
		FieldT* fmicro = const_cast <FieldT*> (fMicro);
		FieldT* fdispl = const_cast <FieldT*> (fDispl);
		
		(*fdispl)[2] = 0;
		(*fmicro)[2] = 0;
    }
    */

    /* reassign initial 2nd time derivative of micro-displacement-gradient to 1st derivative */
    /*
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==2)
    {
		FieldT* fmicro = const_cast <FieldT*> (fMicro);
		(*fmicro)[1] = (*fmicro)[2];
		(*fmicro)[2] = 0;
    }
    */
   
    /* store more recently updated values */
    fdState = fdState_new;
    fiState = fiState_new;
	
//    fs_micromorph3D_out	<< endl 
//		<< setw(outputFileWidth) << "time_step"
//		<< endl;
    step_number = ElementSupport().StepNumber();
//    fs_micromorph3D_out	<< setw(outputFileWidth) << step_number
//		<< endl;
//    fs_micromorph3D_out	<< endl << "**********************************************************************************************";
//    fs_micromorph3D_out	<< endl << "**********************************************************************************************" << endl;
}


/* resets to the last converged solution */
/*
GlobalT::RelaxCodeT FSMicromorphic3DT::ResetStep(void)
{
	const char caller[] = "FSMicromorphic3DT::ResetStep";
	
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
GlobalT::RelaxCodeT FSMicromorphic3DT::RelaxSystem(void)
{
	const char caller[] = "FSMicromorphic3DT::RelaxSystem";
	
	// inherited 
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	// loop over materials 
	//needs to be implemented
#pragma message("relax step for materials not implemented")	
	//ExceptionT::GeneralFail(caller, "relax step for materials not implemented");

	return relax;
}
*/


void FSMicromorphic3DT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented
}


/* return geometry and number of nodes on each facet */
void FSMicromorphic3DT::FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry, 
	iArrayT& num_facet_nodes) const
{
	/* from integration domain */
	ShapeFunctionDispl().FacetGeometry(facet_geometry, num_facet_nodes);
}


/* form of tangent matrix */
GlobalT::SystemTypeT FSMicromorphic3DT::TangentType(void) const
{
    return GlobalT::kNonSymmetric; 
}

/*
void FSMicromorphic3DT::SetStatus(const ArrayT<ElementCardT::StatusT>& status)
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
void FSMicromorphic3DT::InitialCondition(void)
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
void FSMicromorphic3DT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
    const char caller[] = "FSMicromorphic3DT::AddNodalForce";

    /* displ, micro, or neither */
    bool is_displ = false;
    dArrayT* element_force = NULL;
    int num_force = 0;
    if (field.FieldName() == fDispl->FieldName()) 
    {
		is_displ = true;
		element_force = &fFd_int;
		num_force = fDispl->NumDOF();
    }
    else if (field.FieldName() == fMicro->FieldName()) 
    {
		is_displ = false;
		element_force = &fFphi_int;
		num_force = fMicro->NumDOF();
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
		    const iArrayT& nodes_micro = fElementCards_micro[e].NodesU();

		    u.SetLocal(nodes_displ);
		    u_n.SetLocal(nodes_displ);
		    Phi.SetLocal(nodes_micro);
		    Phi_n.SetLocal(nodes_micro);

		    del_u.DiffOf (u, u_n);
		    del_Phi.DiffOf (Phi, Phi_n);

		    // calculate derivatives based on reference coordinates
		    fInitCoords_displ.SetLocal(fElementCards_displ[e].NodesX());
		    //fCurrCoords_displ.SetToCombination (1.0, fInitCoords_displ, 1.0, u); 
		    fCurrCoords_displ=fInitCoords_displ;
		    fShapes_displ->SetDerivatives_DN_DDN(); 

		    //
		    fInitCoords_micro.SetLocal(fElementCards_micro[e].NodesX());
		    fCurrCoords_micro=fInitCoords_micro;
		    //fCurrCoords_micro.SetToCombination (1.0, fInitCoords_micro, 1.0, u); 
		    fShapes_micro->SetDerivatives();
			
		    //update state variables
		    fdstatenew_all.Alias(fNumIP_displ, knum_d_state, fdState_new(CurrElementNumber()));
		    fdstate_all.Alias(fNumIP_displ, knum_d_state, fdState(CurrElementNumber()));

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
				/* residual for micro-displacement-gradient field */ 
				// generate this vector fFphi_int
				fShapes_displ->TopIP();
				while (fShapes_displ->NextIP())
				{
				    //nothing right now
				    fFphi_int=0.0;
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
		    else /* micro-displacement-gradient nodal dof */
			{
			    for (int i = 0; i < nodes_micro.Length(); i++)
			    {
					if (nodes_micro[i] == node)
					{
					    /* components for node */
					    nodalforce.Set(num_force, element_force->Pointer(dex));

					    /* accumulate */
					    force += nodalforce;
					}
					dex += fMicro->NumDOF();
			    }		
			}
		}
    }
//	cout << "F_int = \n" << fFd_int << endl;
}

//---------------------------------------------------------------------

double FSMicromorphic3DT::InternalEnergy ( void )
{
//not implemented
    return 0.0;
}

//---------------------------------------------------------------------

/* write restart data to the output stream */
void FSMicromorphic3DT::WriteRestart(ostream& out) const
{
    /* inherited */
    ElementBaseT::WriteRestart(out);

    /* write state variable data */
    out << fdState;
}

//---------------------------------------------------------------------

/* read restart data to the output stream */
void FSMicromorphic3DT::ReadRestart(istream& in)
{
    /* inherited */
    ElementBaseT::ReadRestart(in);

    /* write state variable data */
    in >> fdState;
}

//---------------------------------------------------------------------

void FSMicromorphic3DT::RegisterOutput(void)
{
    /* collect block ID's */
    ArrayT<StringT> block_ID(fBlockData.Length());
    for (int i = 0; i < block_ID.Length(); i++)
	block_ID[i] = fBlockData[i].ID();

    /* output per element - strain, stress, and ISVs at the integration points */
    ArrayT<StringT> e_labels(fNumIP_displ*(knumstrain+knumstress+knum_d_state));

    /* over integration points */
    // enter what values you need at integration points
    // stress and strain
    const char* slabels3D[] = {"s11", "s22", "s33","s23","s13","s12","e11","e22","e33","e23","e13","e12"};
    // state variables; ?
    const char* svlabels3D[] = {"thing1","thing2","J"};
    int count = 0;
    for (int j = 0; j < fNumIP_micro; j++)
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
    int num_node_output = fDispl->NumDOF() + fMicro->NumDOF() + knumstrain + knumstress + knum_d_state;
    ArrayT<StringT> n_labels(num_node_output);
    count = 0;

    /* labels from micro-displacement-gradient */
    const ArrayT<StringT>& micro_labels = fMicro->Labels();
    for (int i = 0; i < micro_labels.Length(); i++)
	n_labels[count++] = micro_labels[i];

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
#pragma message("FSMicromorphic3DT::RegisterOutput: is this right? ")
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

void FSMicromorphic3DT::WriteOutput(void)
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
		out_variable_all.Alias(fNumIP_displ, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
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
    int num_node_output = fDispl->NumDOF() + fMicro->NumDOF() + knumstrain + knumstress + knum_d_state;
    dArray2DT n_values(nodes_used.Length(), num_node_output);

    /* collect nodal values */
    const dArray2DT& fPhi = (*fMicro)[0];
    const dArray2DT& fU = (*fDispl)[0];
    for (int i = 0; i < nodes_used.Length(); i++)
    {
		int node = nodes_used[i];
		double* row = n_values(i);
		for (int j = 0; j < fPhi.MinorDim(); j++)
		    *row++ = fPhi(node,j);

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
void FSMicromorphic3DT::RHSDriver_staggered(void)
{
	const char caller[] = "FSMicromorphic3DT::RHSDriver_staggered";
#pragma message("staggered solution not implemented")
}

//---------------------------------------------------------------------
/* form group contribution to the stiffness matrix and RHS */
void FSMicromorphic3DT::RHSDriver_monolithic(void)
{
    const char caller[] = "FSMicromorphic3DT::RHSDriver_monolithic";
    if (fDispl->Group() != fMicro->Group())
	ExceptionT::GeneralFail(caller, "displacement and micro-displacement-gradient groups must be the same: %d != %d",
				fDispl->Group(), fMicro->Group());

    int curr_group = ElementSupport().CurrentGroup();

    /* stress output work space */
    dArray2DT	out_variable_all, fdstatenew_all, fdstate_all;
    dArrayT		out_variable;

    /* time Step Increment */
    double delta_t = ElementSupport().TimeStep();
    time = ElementSupport().Time();
    step_number = ElementSupport().StepNumber();

    /* print time */
//    fs_micromorph3D_out	<<"delta_t "<<delta_t << endl ;
//    fs_micromorph3D_out	<<"time "<<time << endl ;

    /* loop over elements */
    int e,l;
    Top();

//   fs_micromorph3D_out	<<"kInitialConditionType "<<kInitialConditionType << endl ;
//   fs_micromorph3D_out	<<"kAnalysisType "<<kAnalysisType << endl ;

	//don't need now for micromorphic
    /* at time=0 when geostatic initial condition is calculated, 
       trapezoidal integrator will calculate first time derivative of fields 
       which by setting alpha_delta_t = 1 will be changed to 
       displacement and pressure which should be assigned to them, 
       note that at time=0, delta_t=0 and Trapezoidal scheme 
       which is embeded in the integrator will do nothing by itself
       (in changing previous iteration values)*/
    /*   
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==1)
    {
		FieldT* fmicro = const_cast <FieldT*> (fMicro);
		FieldT* fdispl = const_cast <FieldT*> (fDispl);
		(*fmicro)[0] = (*fmicro)[1];
		(*fdispl)[0] = (*fdispl)[1];
    }  
    */

    /* at time=0 when geostatic initial condition is calculated, 
       dynamic Newmark integrator will calculate second time derivative 
       of fields which by setting beta_delta_t2 = 1 will be changed to 
       displacement and pressure which should be assigned to them, 
       note that at time=0, delta_t=0 and Newmark scheme which is embeded in 
       dynamic integrator will do nothing by itself(in changing previous iteration value)*/
    /*   
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==2)
    {
		FieldT* fmicro = const_cast <FieldT*> (fMicro);
		FieldT* fdispl = const_cast <FieldT*> (fDispl);
		(*fmicro)[0] = (*fmicro)[2];
		(*fdispl)[0] = (*fdispl)[2];
    } 
    */
   
    
    while (NextElement())
    {
    /* initialize */
	fFd_int_N1_vector = 0.0;
	fK_dd_G3_1_matrix = 0.0;
	fK_dd_G3_2_matrix = 0.0;
	fK_dd_G3_3_matrix = 0.0;
	fK_dd_G3_4_matrix = 0.0;
	fK_dd_BTDB_matrix = 0.0;
	fFd_int_smallstrain_vector = 0.0;
	fFd_int_G4_vector = 0.0;
	fK_dd_G4_matrix = 0.0;	
		
	e = CurrElementNumber();
	const iArrayT& nodes_displ = fElementCards_displ[e].NodesU();
	const iArrayT& nodes_micro = fElementCards_micro[e].NodesU();

	u.SetLocal(nodes_displ);
	u_n.SetLocal(nodes_displ);
	if (u_dot.IsRegistered()) u_dot.SetLocal(nodes_displ);
	if (u_dot_n.IsRegistered()) u_dot_n.SetLocal(nodes_displ);
	if (u_dotdot.IsRegistered()) u_dotdot.SetLocal(nodes_displ);
	if (u_dotdot_n.IsRegistered())u_dotdot_n.SetLocal(nodes_displ);

	Phi.SetLocal(nodes_micro);
	Phi_n.SetLocal(nodes_micro);
	if (Phi_dot.IsRegistered()) Phi_dot.SetLocal(nodes_micro);
	if (Phi_dot_n.IsRegistered()) Phi_dot_n.SetLocal(nodes_micro);
	if (Phi_dotdot.IsRegistered()) Phi_dotdot.SetLocal(nodes_micro);
	if (Phi_dotdot_n.IsRegistered()) Phi_dotdot_n.SetLocal(nodes_micro);

	/* print solid displacement at current step (u)*/
/*	fs_micromorph3D_out	<<"nodal solid displacement at current step(u)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    fs_micromorph3D_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		fs_micromorph3D_out << u(i,j) << "\t";
	    fs_micromorph3D_out	<< endl ;
	}
*/

	/* print solid displacement from previous step (u_n)*/
/*	fs_micromorph3D_out	<<"nodal solid displacement from previous step(u_n)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    fs_micromorph3D_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		fs_micromorph3D_out << u_n(i,j) << "\t";
	    fs_micromorph3D_out	<< endl ;
	}
*/
	
	/* print solid velocity at current step (u_dot)*/
/*	fs_micromorph3D_out	<<"nodal solid velocity at current step(u_dot)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    fs_micromorph3D_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		fs_micromorph3D_out << u_dot(i,j) << "\t";
	    fs_micromorph3D_out	<< endl ;
	}
*/

	/* print solid velocity from previous step (u_dot_n)*/
/*	fs_micromorph3D_out	<<"nodal solid velocity from previous step(u_dot_n)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    fs_micromorph3D_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		fs_micromorph3D_out << u_dot_n(i,j) << "\t";
	    fs_micromorph3D_out	<< endl ;
	}
*/


	/* print solid acceleration at current step (u_dotdot)*/
/*	fs_micromorph3D_out	<<"nodal solid velocity at current step(u_dotdot)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    fs_micromorph3D_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		fs_micromorph3D_out << u_dotdot(i,j) << "\t";
	    fs_micromorph3D_out	<< endl ;
	}
*/

	/* print solid acceleration from previous step (u_dotdot_n)*/
/*	fs_micromorph3D_out	<<"nodal solid velocity from previous step(u_dotdot_n)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    fs_micromorph3D_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		fs_micromorph3D_out << u_dotdot_n(i,j) << "\t";
	    fs_micromorph3D_out	<< endl ;
	}
*/
		 
	/* print micro-displacement-gradient at current step (micro)*/
/*	fs_micromorph3D_out	<<"nodal micro-displacement-gradient at current step(micro)"<< endl ;
	for (int i=0; i<n_en_micro; i++)
	{
	    fs_micromorph3D_out	<< "node number " << i+1 <<" :  " ;	
	    fs_micromorph3D_out	<< Phi(i,0) << endl;
	}
*/

	/* print micro-displacement-gradient from previous step (micro_n)*/
/*	fs_micromorph3D_out	<<"nodal micro-displacement-gradient from previous step(micro_n)"<< endl ;
	for (int i=0; i<n_en_micro; i++)
	{
	    fs_micromorph3D_out	<< "node number " << i+1 <<" :  " ;	
	    fs_micromorph3D_out	<< Phi_n(i,0) << endl;
	}
*/
	
	/* print first derivative of micro-displacement-gradient at current step (micro_dot)*/
/*	fs_micromorph3D_out	<<"first derivative of nodal micro-displacement-gradient at current step(micro_dot)"<< endl ;
	for (int i=0; i<n_en_micro; i++)
	{
	    fs_micromorph3D_out	<< "node number " << i+1 <<" :  " ;	
	    fs_micromorph3D_out	<< Phi_dot(i,0) << endl;
	}
*/

	/* print first derivative of micro-displacement-gradient from previous step (micro_dot_n)*/
/*	fs_micromorph3D_out	<<"first derivative of nodal micro-displacement-gradient from previous step(micro_dot_n)"<< endl ;
	for (int i=0; i<n_en_micro; i++)
	{
	    fs_micromorph3D_out	<< "node number " << i+1 <<" :  " ;	
	    fs_micromorph3D_out	<< Phi_dot_n(i,0) << endl;
	}
*/


	/* print second derivative of micro-displacement-gradient at current step (micro_dotdot)*/
/*	fs_micromorph3D_out	<<"second derivative of nodal micro-displacement-gradient at current step(micro_dotdot)"<< endl ;
	for (int i=0; i<n_en_micro; i++)
	{
	    fs_micromorph3D_out	<< "node number " << i+1 <<" :  " ;	
	    fs_micromorph3D_out	<< Phi_dotdot(i,0) << endl;
	}
*/

	/* print second derivative of micro-displacement-gradient from previous step (micro_dotdot_n)*/
/*	fs_micromorph3D_out	<<"second derivative of nodal micro-displacement-gradient from previous step(micro_dotdot_n)"<< endl ;
	for (int i=0; i<n_en_micro; i++)
	{
	    fs_micromorph3D_out	<< "node number " << i+1 <<" :  " ;	
	    fs_micromorph3D_out	<< Phi_dotdot_n(i,0) << endl;
	}
*/
	
	
	/* populate solid displacement,solid velocity and 
	   solid accelration in vector form*/
	int index_u = 0;
	for (int i=0; i<n_en_displ; i++)
	{
	    for (int j=0; j<n_sd; j++)
	    {
			u_vec[index_u] = u(i,j);
			u_dot_vec[index_u] = u_dot(i,j);
			u_dotdot_vec[index_u] = u_dotdot(i,j);
			index_u += 1;
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
	
	/* populate micro-displacement-gradient, first and second time derivatives of 
	   micro-displacement-gradient in vector form*/
	int index_Phi = 0;
	for (int i=0; i<n_en_micro; i++)
	{
	    for (int j=0; j<ndof_per_nd_micro; j++)
	    {
			Phi_vec[index_Phi] = Phi(i,j);
			Phi_dot_vec[index_Phi] = Phi_dot(i,j);
			Phi_dotdot_vec[index_Phi] = Phi_dotdot(i,j);
			index_Phi += 1;
	    }
	}

	/* [micro_dot_column_matrix] will be formed */
	/*	
	for (int i=0; i<n_en_micro_x_ndof_per_nd_micro; i++)
	    Phi_dot_column_matrix(i,0) = Phi_dot_vec[i];
	*/

	del_u.DiffOf (u, u_n);
	del_Phi.DiffOf (Phi, Phi_n);
	
	// calculate derivatives based on reference coordinates
	fInitCoords_displ.SetLocal(fElementCards_displ[e].NodesX());
	fCurrCoords_displ=fInitCoords_displ;
	//fCurrCoords_displ.SetToCombination (1.0, fInitCoords_displ, 1.0, u); 
	fShapes_displ->SetDerivatives_DN_DDN(); 
	//
	fInitCoords_micro.SetLocal(fElementCards_micro[e].NodesX());
	fCurrCoords_micro=fInitCoords_micro;
	//fCurrCoords_micro.SetToCombination (1.0, fInitCoords_micro, 1.0, u); 
	fShapes_micro->SetDerivatives(); 
	
	//update state variables
	fdstatenew_all.Alias(fNumIP_displ, knum_d_state, fdState_new(CurrElementNumber()));
	fdstate_all.Alias(fNumIP_displ, knum_d_state, fdState(CurrElementNumber()));
	/*
	fdstatenew_all.Alias(fNumIP_micro, knum_d_state, fdState_new(CurrElementNumber()));
	fdstate_all.Alias(fNumIP_micro, knum_d_state, fdState(CurrElementNumber()));
	*/
	
	if (bStep_Complete) 
	{ 
	    //-- Store/Register data in classic tahoe manner 
	    out_variable_all.Alias(fNumIP_displ, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
	    //out_variable_all.Alias(fNumIP_micro, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
	    for (l=0; l < fNumIP_displ; l++) 
	    //for (l=0; l < fNumIP_micro; l++) 
	    {
			out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(l));
			Put_values_In_dArrayT_vector(fCauchy_stress_Elements_IPs, e,l,fTemp_six_values);
			out_variable.CopyIn(0,fTemp_six_values);
			Put_values_In_dArrayT_vector(fEulerian_strain_Elements_IPs, e,l,fTemp_six_values);
			out_variable.CopyIn(6,fTemp_six_values);
			/*
			out_variable[13]=fState_variables_Elements_IPs(e,l*3+0);
			out_variable[14]=fState_variables_Elements_IPs(e,l*3+1);
			out_variable[15]=fState_variables_Elements_IPs(e,l*3+2);
			*/
	    } 
	}
	else 
	{ //-- Still Iterating
		
		    /* residual and tangent for displacements */
            
		    const double* Det    = fShapes_displ->IPDets();
		    const double* Weight = fShapes_displ->IPWeights();
		    fShapes_displ->TopIP();
		    fShapes_micro->TopIP();
		    
		    while (fShapes_displ->NextIP() && fShapes_micro->NextIP())
		    {
				double scale_const = (*Weight++)*(*Det++);
			
				const int IP = fShapes_displ->CurrIP();	
				dArrayT DisplIPCoordinate(n_sd), MicroIPCoordinate(n_sd);
				fShapes_displ->IPCoords(DisplIPCoordinate);
				fShapes_micro->IPCoords(MicroIPCoordinate);

				const double* shapes_displ_X = fShapes_displ->IPShapeX();
				/* [fShapeDispl]will be formed */
				Form_solid_shape_functions(shapes_displ_X);
				
				fShapes_displ->GradNa(fShapeDisplGrad_temp);
				/* [fShapeDisplGrad] will be formed */
				Form_Gradient_of_solid_shape_functions(fShapeDisplGrad_temp);
				
				/* [fShapeDisplGrad_t] and [fShapeDisplGrad_t_Transpose] will be formed */
				Form_Gradient_t_of_solid_shape_functions(fShapeDisplGrad_temp);
				fShapeDisplGrad_t_Transpose.Transpose(fShapeDisplGrad_t);
			
				const double* shapes_micro_X = fShapes_micro->IPShapeX();
				/* {fShapeMicro} will be formed */
				Form_micro_shape_functions(shapes_micro_X);
				
				/* [fShapeMicro_row_matrix] will be formed */
				//need?				
				for (int i=0; i<n_en_micro ; i++)
				    fShapeMicro_row_matrix(0,i) = fShapeMicro[i];
				
				/* [fShapeMicroGrad] will be formed */
				fShapes_micro->GradNa(fShapeMicroGrad);			
				
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
				
				/* Calculating Jacobian */
				double J = fDeformation_Gradient.Det();
				/* Jacobian for the current IP will be saved */
				fState_variables_IPs(IP,2)=J;
					
				/*  fRho */
				fRho_0 = fMaterial_Params[kRho_0];
				
				/* Calculating fRho */
				fRho = fRho_0/J;
				
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

				/* [fEulerian_strain_tensor_current_IP] will be formed */
				fEulerian_strain_tensor_current_IP = fLeft_Cauchy_Green_tensor_Inverse;
				fEulerian_strain_tensor_current_IP *= -1; 
				fEulerian_strain_tensor_current_IP += fIdentity_matrix;
				fEulerian_strain_tensor_current_IP *= 0.5;
				
				/* extract six values of strain from symmetric eulerian strain tensor */
				Extract_six_values_from_symmetric_tensor(fEulerian_strain_tensor_current_IP,fTemp_six_values);
				
				/* Save Eulerian strain tensor of the current IP */ 
				fEulerian_strain_IPs.SetRow(IP,fTemp_six_values);
				
				/* Calculating J_Prim */
				if (fRight_Cauchy_Green_tensor.Det()==0)
				    fRight_Cauchy_Green_tensor = fIdentity_matrix; 
				double TempJ_Prim=fRight_Cauchy_Green_tensor.Det();
				double J_Prim=sqrt(fabs(TempJ_Prim));
				
				/* [fSecond_Piola_tensor] will be formed */
				fSecond_Piola_tensor.SetToScaled(fMaterial_Params[kLambda]*log(J_Prim)-fMaterial_Params[kMu],fRight_Cauchy_Green_tensor_Inverse); 
				fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fIdentity_matrix);
				fSecond_Piola_tensor += fTemp_matrix_nsd_x_nsd;
				
				/* [fKirchhoff_tensor] will be formed */
				fKirchhoff_tensor.MultABCT(fDeformation_Gradient,fSecond_Piola_tensor,fDeformation_Gradient);
				
				/* [fCauchy_effective_stress_tensor_current_IP] will be formed */
				fCauchy_stress_tensor_current_IP = fKirchhoff_tensor;
				fCauchy_stress_tensor_current_IP *= 1/J;
				/* extract six values of stress from symmetric cauchy stress tensor */
				Extract_six_values_from_symmetric_tensor(fCauchy_stress_tensor_current_IP,fTemp_six_values);
				/* Save Cauchy effective stress tensor of the current IP */ 
				fCauchy_stress_IPs.SetRow(IP,fTemp_six_values); 
			
				/* {fKirchhoff_vector} will be formed */
				Form_kirchhoff_stress_vector();
				
				/* [fIota_temp_matrix] will be formed */
				fIota_temp_matrix.MultATB(fShapeDisplGrad,fDefGradInv_Grad_grad);
				
				/* second derivatives of solid shape functions, [fShapeDisplGradGrad] will be formed */
				fShapes_displ->Grad_GradNa(fShapeDisplGradGrad);
				
				
				//need?
				/* [fVarpi_temp_matrix] will be formed */
				Form_Varpi_temp_matrix();
				
				/* {fChi_temp_vector} will be formed */
				fVarpi_temp_matrix.Multx(u_vec,fChi_temp_vector);
				
				/* [fChi_temp_column_matrix] will be formed */
				for (int i=0; i<3 ; i++)
				    fChi_temp_column_matrix(i,0)= fChi_temp_vector[i];
				
				
				
				/* {fFd_int_N1_vector} will be formed */
				//pg33,34 of Davoud's thesis
				double scale = scale_const;
				fIota_temp_matrix.Multx(fKirchhoff_vector,fTemp_vector_ndof_se,scale);
				/* fFd_int_N1_vector for the current IP */
				/* accumulate */
				fFd_int_N1_vector += fTemp_vector_ndof_se;
				
				/* [fIm_temp_matrix] will be formed */
				Form_Im_temp_matrix();
				
				/* [fHbar_temp_matrix] will be formed */
				Form_Hbar_temp_matrix();
				
				/* [fEll_temp_matrix] will be formed */
				Form_Ell_temp_matrix();
				
				/* {fPi_temp_transpose_vector} will be formed */
				fShapeDisplGrad.MultTx(fDefGradInv_vector,fPi_temp_transpose_vector);
			
				/* [fPi_temp_row_matrix] will be formed */
				for (int i=0; i<n_en_displ_x_n_sd; i++)
				    fPi_temp_row_matrix(0,i) = fPi_temp_transpose_vector[i];
			
				/* [fK_dd_G3_1_matrix] will be formed */
				//pg57 of Davoud's thesis
				fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fIm_temp_matrix,fIota_temp_matrix);
				scale = -1*scale_const;
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
				scale = fMaterial_Params[kMu] * scale_const;
				fTemp_matrix_ndof_se_x_ndof_se *= scale;
				/* accumulate */
				fK_dd_G3_2_matrix += fTemp_matrix_ndof_se_x_ndof_se;
			
			
				/* [fK_dd_G3_3_matrix] will be formed */
				fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fEll_temp_matrix,fIota_temp_matrix);
				scale = fMaterial_Params[kMu] * scale_const;
				fTemp_matrix_ndof_se_x_ndof_se *= scale;
				/* accumulate */
				fK_dd_G3_3_matrix += fTemp_matrix_ndof_se_x_ndof_se;
				
				
				/* [fK_dd_G3_4_matrix] will be formed */
				fTemp_matrix_ndof_se_x_ndof_se.MultABC(fIota_temp_matrix,fI_ij_column_matrix,fPi_temp_row_matrix);
				scale = fMaterial_Params[kLambda] * scale_const; 
				fTemp_matrix_ndof_se_x_ndof_se *= scale;
				/* accumulate */
				fK_dd_G3_4_matrix += fTemp_matrix_ndof_se_x_ndof_se;
				
				
				//need?
				/* {fGrad_1_J_vector} will be filled */
				fVarpi_temp_matrix.Multx(u_vec,fGrad_1_J_vector, -1.0/J);
				
				
				/* Creating Second tangential elasticity tensor in the Ref. coordinate [fC_matrix] */
				Form_C_matrix(J_Prim);
	 
				
			
				/* Creating Second tangential elasticity tensor in the Current coordinate [fc_matrix]*/
				Form_c_matrix();
			
				/* [fIm_Prim_temp_matrix] will be formed */
				Form_Im_Prim_temp_matrix();
				
				/* {fFd_int_G4_vector} will be formed */
				fShapeDispl.MultTx(fGravity_vector,fTemp_vector_ndof_se);
				scale = -1*fRho_0*scale_const;
				fTemp_vector_ndof_se *= scale; 
				/* accumulate */
				fFd_int_G4_vector += fTemp_vector_ndof_se;
				
				
				//need??
				/* {fgradv_vector} will be formed */
				Form_gradv_vector();
				
				/* [fXi_temp_matrix] will be formed */
				Form_Xi_temp_matrix();
				
				/* [fVarsigma_temp_matrix] will be formed */
				Form_Varsigma_temp_matrix();
				
				
				/* [fI_ijkl_matrix] will be formed */
				Form_I_ijkl_matrix();
				
				
				/* [fK_dd_G4_matrix] will be formed */
				//pg57 Davoud's thesis
				fTemp_matrix_ndof_se_x_ndof_se.MultATBC(fShapeDispl,fGravity_column_matrix,fPi_temp_row_matrix);
				scale = -1*J*(fRho)*scale_const;
				fTemp_matrix_ndof_se_x_ndof_se *= scale;
				/* accumulate */
				fK_dd_G4_matrix += fTemp_matrix_ndof_se_x_ndof_se;
				
				
				
				/*************************************************/
				/* implementing small strain deformation of an isotropic linear elastic media for debugging */
				
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
				
				/* end of small strain  */
				/*************************************************/
			
			
				/* for debugging */
				/*
				const int ip = fShapes_displ->CurrIP()+1;
				fs_micromorph3D_out	<< endl << "IP" << ip
						<< setw(outputFileWidth) << ", shape function matrix for solid phase: " 
						<< setw(outputFileWidth) << fShapeDispl;
				
				fs_micromorph3D_out	<< endl << "terms from shape function matrix for solid phase: " 
						<< setw(outputFileWidth) << fShapeDispl(0,0) 
						<< setw(outputFileWidth) << fShapeDispl(0,3);	
				*/
			
		    } //end Gauss integration loop
		    		    
	    
			
			/* saving eulerian strain for each IPs of the current element */
		    fEulerian_strain_Elements_IPs.SetRow(e,fEulerian_strain_IPs);
	    
		    /* saving cauchy stress for each IPs of the current element */
		    fCauchy_stress_Elements_IPs.SetRow(e,fCauchy_stress_IPs);

		    /* saving state variables for each IPs of the current element */
		    fState_variables_Elements_IPs.SetRow(e,fState_variables_IPs);

			/* {fFd_int} will be formed */
		    fFd_int = fFd_int_N1_vector;
		    fFd_int += fFd_int_G4_vector;
		    fFd_int *= -1;
	    
		    /* [fKdd] will be formed */
		    fKdd = fK_dd_G3_1_matrix;
		    fKdd += fK_dd_G3_2_matrix;
		    fKdd += fK_dd_G3_3_matrix;
		    fKdd += fK_dd_G3_4_matrix;
		    fKdd += fK_dd_G4_matrix;
		    
		    /* for small strain check */
			/*		
			fK_dd_BTDB_matrix.MultTx(u_vec,fTemp_vector_ndof_se);
			fFd_int = fTemp_vector_ndof_se;
			fFd_int *= -1.0; 
			
			fKdd = fK_dd_BTDB_matrix;	
			*/
		    
		    /* [fKdphi] will be formed */
		    //need to code
		    fKdphi = 0.0;
	    
		    /* [fKphid] will be formed */
		    //need to code
		    fKphid = 0.0;

		    /* [fKphiphi] will be formed */
		    //need to code
		    fKphiphi = 0.0;	

		    /* {fFphi_int} will be formed */
		    //need to code
		    fFphi_int = 0.0;

		    /* equations numbers */
		    const iArrayT& displ_eq = fElementCards_displ[e].Equations();
		    const iArrayT& micro_eq = fElementCards_micro[e].Equations();
		    
		    /* assemble residuals */
		    ElementSupport().AssembleRHS(curr_group, fFd_int, displ_eq);
		    ElementSupport().AssembleRHS(curr_group, fFphi_int, micro_eq);
		    
		    /* assemble components of the tangent */
		    ElementSupport().AssembleLHS(curr_group, fKdd, displ_eq);
		    ElementSupport().AssembleLHS(curr_group, fKphiphi, micro_eq);
		    ElementSupport().AssembleLHS(curr_group, fKdphi, displ_eq, micro_eq);
		    ElementSupport().AssembleLHS(curr_group, fKphid, micro_eq, displ_eq);

	}	
    }
}



/* form global shape function derivatives */
void FSMicromorphic3DT::SetGlobalShape(void)
{
    /* fetch (initial) coordinates */
    SetLocalX(fLocInitCoords);
	
    /* compute shape function derivatives */
    fShapes_displ->SetDerivatives_DN_DDN();
    fShapes_micro->SetDerivatives();
}



/* describe the parameters needed by the interface */
void FSMicromorphic3DT::DefineParameters(ParameterListT& list) const
{
    /* inherited */
    ElementBaseT::DefineParameters(list);

    /* displacement field */
    //already done in ElementBaseT
    //list.AddParameter(ParameterT::Word, "displ_field_name");
	
    /* micro-displacement-gradient field */
    list.AddParameter(ParameterT::Word, "micro_field_name");
	
    list.AddParameter(fGeometryCode_displ_int, "GeometryCode_displ");
    list.AddParameter(fNumIP_displ, "NumIP_displ");
    list.AddParameter(fGeometryCodeSurf_displ_int, "GeometryCodeSurf_displ");
    list.AddParameter(fNumIPSurf_displ, "NumIPSurf_displ");
    list.AddParameter(n_en_displ, "n_en_displ");
    list.AddParameter(n_en_micro, "n_en_micro");
    list.AddParameter(ndof_per_nd_micro, "ndof_per_nd_micro");
	
    list.AddParameter(iConstitutiveModelType, "constitutive_mod_type");

    double shearMu, sLambda, Rho_0, gravity_g, gravity_g1, gravity_g2, gravity_g3;

    // solid elasticity
    list.AddParameter(shearMu, "mu");
    list.AddParameter(sLambda, "lambda");
	
    // gravity
    list.AddParameter(gravity_g, "g");

    // gravity in each direction (depends on the coordinate system which we have chosen for the problem)
    list.AddParameter(gravity_g1, "g1");
    list.AddParameter(gravity_g2, "g2");
    list.AddParameter(gravity_g3, "g3");
	
    // reference mass density
    list.AddParameter(Rho_0, "rho_0");
	
    // Newmark time integration parameters
//    list.AddParameter(newBeta, "beta");
//    list.AddParameter(newGamma, "gamma");

}


/* accept parameter list */
void FSMicromorphic3DT::TakeParameterList(const ParameterListT& list)
{
    const char caller[] = "FSMicromorphic3DT::TakeParameterList";
	
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

    /* get micro-displacement-gradient field */
    const StringT& micro_field_name = list.GetParameter("micro_field_name");
    fMicro = ElementSupport().Field(micro_field_name);
    if (!fMicro)
	ExceptionT::GeneralFail(caller, "could not resolve \"%s\" micro_field", 
				micro_field_name.Pointer());

    fGeometryCode_displ_int = list.GetParameter("GeometryCode_displ");
    fGeometryCode_displ = GeometryT::int2CodeT(fGeometryCode_displ_int);
    fNumIP_displ = list.GetParameter("NumIP_displ");
    fGeometryCodeSurf_displ_int = list.GetParameter("GeometryCodeSurf_displ");
    fGeometryCodeSurf_displ = GeometryT::int2CodeT(fGeometryCodeSurf_displ_int);
    fNumIPSurf_displ = list.GetParameter("NumIPSurf_displ");
    n_en_displ = list.GetParameter("n_en_displ");
    n_en_micro = list.GetParameter("n_en_micro");
    ndof_per_nd_micro = list.GetParameter("ndof_per_nd_micro");
    
    //kAnalysisType = list.GetParameter("type_of_analysis_1consolidation_2dynamic");
    //kInitialConditionType = list.GetParameter("initial_condition_1geostatic_2displ_vel_press");
  
    fGeometryCode_micro = fGeometryCode_displ; 
    fNumIP_micro = fNumIP_displ;
    fGeometryCodeSurf_micro = fGeometryCodeSurf_displ;
    fNumIPSurf_micro = fNumIPSurf_displ;
	
    iConstitutiveModelType = list.GetParameter("constitutive_mod_type");
	
    fMaterial_Params.Dimension ( kNUM_FMATERIAL_TERMS );
//    fIntegration_Params.Dimension ( kNUM_FINTEGRATE_TERMS );
	
    fMaterial_Params[kMu] = list.GetParameter("mu");
    fMaterial_Params[kLambda] = list.GetParameter("lambda");
    fMaterial_Params[kg] = list.GetParameter("g");
    fMaterial_Params[kg1] = list.GetParameter("g1");
    fMaterial_Params[kg2] = list.GetParameter("g2");
    fMaterial_Params[kg3] = list.GetParameter("g3");
    fMaterial_Params[kRho_0] = list.GetParameter("rho_0");
	
//    fIntegration_Params[kBeta] = list.GetParameter("beta");
//    fIntegration_Params[kGamma] = list.GetParameter("gamma");
	
    Echo_Input_Data();
	
	//need to change for appropriate number of ISVs for micromorphic model
    knum_d_state = 3; // #? internal state variables
    knum_i_state = 0; // int's needed per ip, state variables
	
	//need to change these for non-symmetric stress, and higher order couple stress output
    knumstrain = 6; // number of strain outputs
    knumstress = 6; // number of stress outputs + higher order = ??
	
    output = "out";
	
    /* dimensions (notation as per Hughes' Book) */
    int& n_ip_displ = fNumIP_displ;
    int& n_ip_micro = fNumIP_micro;
    n_sd = NumSD();
    int nen = NumElementNodes(); /* number of nodes/element in the mesh */

    /* initialize connectivities */
    fConnectivities_displ.Alias(fConnectivities);
    fConnectivities_micro.Alias(fConnectivities);

    /* pick element interpolations based on available number of element nodes
     * and the specified number of integration points */
    // only implemented for 3D, quadratic hexs
    //if (n_sd == 2 && n_en_micro != n_en_displ && fGeometryCode_displ == GeometryT::kQuadrilateral) 
    if (n_sd == 3 && n_en_micro != n_en_displ && fGeometryCode_displ == GeometryT::kHexahedron) 
    {
		// don't expect reduced integration for both fields 
		// if (n_ip_displ == 4 && n_ip_micro == 4)
		//	ExceptionT::GeneralFail(caller, "not expecting 4 ips for both fields");
		//else if (n_ip_displ == 4 || n_ip_micro == 4) // create reduced connectivities
		//{ 
		// reduce the number of element nodes based on the number ip's
		int& nen_red = (n_ip_displ == 8) ? n_en_displ : n_en_micro;
		nen_red = 8;
		ArrayT<const iArray2DT*>& connects_red = (n_ip_displ == 8) ? 
		    fConnectivities_displ : 
		    fConnectivities_micro;
			
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
    // phi
    fInitCoords_micro.Dimension(n_en_micro, n_sd);
    ElementSupport().RegisterCoordinates(fInitCoords_micro);	
    fCurrCoords_micro.Dimension(n_en_micro, n_sd);
    fShapes_micro = new ShapeFunctionT(fGeometryCode_micro, fNumIP_micro, fCurrCoords_micro);
    //fShapes_micro = new ShapeFunctionT(fGeometryCode_micro, fNumIP_micro, fCurrCoords_displ);
    //fShapes_micro->Initialize();

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

    /* set local arrays for micro-displacement-gradient field */
    Phi.Dimension (n_en_micro, ndof_per_nd_micro);
    Phi_dot.Dimension (n_en_micro, ndof_per_nd_micro);
    Phi_dot_n.Dimension (n_en_micro, ndof_per_nd_micro);
    Phi_dotdot.Dimension (n_en_micro, ndof_per_nd_micro);
    Phi_dotdot_n.Dimension (n_en_micro, ndof_per_nd_micro);
    Phi_n.Dimension (n_en_micro, ndof_per_nd_micro);
    del_Phi.Dimension (n_en_micro, ndof_per_nd_micro);
    n_en_micro_x_ndof_per_nd_micro = n_en_micro*ndof_per_nd_micro;
    del_Phi_vec.Dimension (n_en_micro_x_ndof_per_nd_micro);
    Phi_vec.Dimension (n_en_micro_x_ndof_per_nd_micro);
    Phi_dot_vec.Dimension (n_en_micro_x_ndof_per_nd_micro);
    Phi_dotdot_vec.Dimension (n_en_micro_x_ndof_per_nd_micro);
    //ElementSupport().RegisterCoordinates(fInitCoords_micro);

    fMicro->RegisterLocal(Phi);
    fMicro->RegisterLocal(Phi_n);

    if (fIntegrator->Order() == 1)
    {
		fMicro->RegisterLocal(Phi_dot);
		fMicro->RegisterLocal(Phi_dot_n);
    }

    if (fIntegrator->Order() == 2)
    {
		fMicro->RegisterLocal(Phi_dot);
		fMicro->RegisterLocal(Phi_dot_n);
		fMicro->RegisterLocal(Phi_dotdot);
		fMicro->RegisterLocal(Phi_dotdot_n);
    }

	
    /* allocate state variable storage */
    // state variables are calculated at IPs for displacement field
    int num_ip = fNumIP_displ;
    //int num_ip = fNumIP_micro;
    fdState_new.Dimension(n_el, num_ip*knum_d_state);
    fdState.Dimension(n_el, num_ip*knum_d_state);
    fiState_new.Dimension(n_el, num_ip*knum_i_state);
    fiState.Dimension(n_el, num_ip*knum_i_state);
	
    /* initialize equations */
    fEqnos_displ.Alias(fEqnos_displ);
    fEqnos_micro.Dimension(fConnectivities_micro.Length());

    /* initialize state variables */
    fdState = 0;
    fdState_new = 0;
    fiState = 0;
    fiState_new = 0;

    /* initialize element cards */
    fElementCards_displ.Alias(fElementCards);
    fElementCards_micro.Dimension(fElementCards.Length());
	
    /* set cards to data in array - NOT NEEDED IF YOU'RE NOT
     * GOING TO USE THE ElementCardT ARRAY? */
    for (int i= 0; i < fElementCards.Length(); i++)
	fElementCards[i].Set(fiState.MinorDim(), fiState(i), fdState.MinorDim(), fdState(i));

    fKdd.Dimension 			( n_en_displ_x_n_sd, n_en_displ_x_n_sd );
    fKdphi.Dimension 		( n_en_displ_x_n_sd, n_en_micro_x_ndof_per_nd_micro );
    fKphid.Dimension 		( n_en_micro_x_ndof_per_nd_micro, n_en_displ_x_n_sd );
    fKphiphi.Dimension 		( n_en_micro_x_ndof_per_nd_micro, n_en_micro_x_ndof_per_nd_micro );

    fFd_int.Dimension 		( n_en_displ_x_n_sd );
    fFd_ext.Dimension 		( n_en_displ_x_n_sd );
    fFphi_int.Dimension 	( n_en_micro_x_ndof_per_nd_micro );
    fFphi_ext.Dimension 	( n_en_micro_x_ndof_per_nd_micro );
	
    /* workspace matricies */
    fShapeDispl.Dimension (n_sd, n_en_displ_x_n_sd);
    fShapeMicro.Dimension (ndof_per_nd_micro, n_en_micro_x_ndof_per_nd_micro);
    n_sd_x_n_sd = n_sd*n_sd;
    fShapeDisplGrad_temp.Dimension (n_sd, n_en_displ);
    fShapeDisplGrad.Dimension (n_sd_x_n_sd, n_en_displ_x_n_sd);
    fShapeDisplGrad_t.Dimension (n_sd_x_n_sd, n_en_displ_x_n_sd);
    fShapeDisplGradGrad.Dimension (n_sd*2 , n_en_displ);
    ndof_per_nd_micro_x_n_sd = ndof_per_nd_micro*n_sd;
    fShapeMicroGrad.Dimension (ndof_per_nd_micro_x_n_sd, n_en_micro_x_ndof_per_nd_micro);
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
    fSecond_Piola_tensor.Dimension (n_sd,n_sd);
    fTemp_matrix_nsd_x_nsd.Dimension (n_sd,n_sd);
    fKirchhoff_tensor.Dimension (n_sd,n_sd);
    fKirchhoff_vector.Dimension (n_sd_x_n_sd);
    fIota_temp_matrix.Dimension (n_en_displ_x_n_sd,n_sd_x_n_sd);
    fVarpi_temp_matrix.Dimension (n_sd, n_en_displ_x_n_sd);
    
    fChi_temp_vector.Dimension (n_sd);
    fTemp_vector_9x1.Dimension (n_sd_x_n_sd);
    fFd_int_N1_vector.Dimension (n_en_displ_x_n_sd);
    fTemp_vector_ndof_se.Dimension (n_en_displ_x_n_sd); 
    fTemp_vector_nen_micro.Dimension (n_en_micro); 
    fIm_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fHbar_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fEll_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fPi_temp_transpose_vector.Dimension (n_en_displ_x_n_sd); 
    fPi_temp_row_matrix.Dimension (1,n_en_displ_x_n_sd); 
    
    fK_dd_G3_1_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G3_2_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G3_3_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G3_4_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fI_ij_column_matrix.Dimension (n_sd_x_n_sd, 1);
    fShapeDisplGrad_t_Transpose.Dimension (n_en_displ_x_n_sd, n_sd_x_n_sd);
    fShapeMicro_row_matrix.Dimension (1,n_en_micro); 
    fTemp_nsd_vector.Dimension (n_sd);
    fGrad_1_J_vector.Dimension (n_sd);
    fChi_temp_column_matrix.Dimension (n_sd, 1);
    fTemp_matrix_nsd_x_1.Dimension (n_sd,1); 
    fTemp_matrix_ndof_se_x_ndof_se.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fTemp_matrix_ndof_se_x_nen_micro.Dimension (n_en_displ_x_n_sd,n_en_micro);
    fc_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fC_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fIm_Prim_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fB_matrix.Dimension (6 ,n_en_displ_x_n_sd);
    fD_matrix.Dimension (6,6);
    fK_dd_BTDB_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fFd_int_smallstrain_vector.Dimension (n_en_displ_x_n_sd);
    fEulerian_strain_tensor_current_IP.Dimension (n_sd,n_sd);
    fCauchy_stress_tensor_current_IP.Dimension (n_sd,n_sd);
    fEulerian_strain_IPs.Dimension (fNumIP_displ,knumstrain);
    fCauchy_stress_IPs.Dimension (fNumIP_displ,knumstress);
    fState_variables_IPs.Dimension (fNumIP_displ,knum_d_state);
    fTemp_six_values.Dimension (6);
    fEulerian_strain_Elements_IPs.Dimension (NumElements(),fNumIP_displ*knumstrain);
    fCauchy_stress_Elements_IPs.Dimension (NumElements(),fNumIP_displ*knumstress);
    fState_variables_Elements_IPs.Dimension (NumElements(),fNumIP_displ*knum_d_state);
    fGravity_vector.Dimension (n_sd);
    fFd_int_G4_vector.Dimension (n_en_displ_x_n_sd);
    fDefGradInv_column_matrix.Dimension (n_sd_x_n_sd,1);
    fDefGradInv_column_matrix_Transpose.Dimension (1,n_sd_x_n_sd);
    u_dotdot_column_matrix.Dimension (n_en_displ_x_n_sd,1);
    fGradv_vector.Dimension (n_sd_x_n_sd);
    fgradv_vector.Dimension (n_sd_x_n_sd);
    fXi_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fVarsigma_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fI_ijkl_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    u_dot_column_matrix.Dimension (n_en_displ_x_n_sd,1);
    fTemp_matrix1_ndof_se_x_ndof_se.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G4_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fGravity_column_matrix.Dimension (n_sd, 1);
    fTemp_matrix_nsd_x_ndof_se.Dimension (n_sd,n_en_displ_x_n_sd);
    fTemp_matrix_nsd_x_nen_micro.Dimension (n_sd,n_en_micro);
    micro_dot_column_matrix.Dimension (n_en_micro,1);
    u_dot_column_matrix_Transpose.Dimension (1, n_en_displ_x_n_sd);
    fImath_temp_matrix.Dimension (n_sd,n_sd_x_n_sd);

    /* streams */
    ofstreamT& out = ElementSupport().Output();

    /* storage for integration point strain, stress, and ISVs*/
    fIPVariable.Dimension (n_el, fNumIP_displ*(knumstrain+knumstress+knum_d_state));
    //fIPVariable.Dimension (n_el, fNumIP_micro*(knumstrain+knumstress+knum_d_state));
    fIPVariable = 0.0;

    /* allocate storage for nodal forces */
    //fForces_at_Node.Dimension ( n_sd );
	
    /* extract natural boundary conditions */
    TakeNaturalBC(list);
	
    /* setup output file and format */
    outputPrecision = 10;
    outputFileWidth = outputPrecision + 8;
//    fs_micromorph3D_out.open("fs_micromorph3D.info");
}


/* information about subordinate parameter lists */
void FSMicromorphic3DT::DefineSubs(SubListT& sub_list) const
{
    /* inherited */
    ElementBaseT::DefineSubs(sub_list);

    /* element blocks */
    sub_list.AddSub("micromorphic_FS_3D_element_block");
	
    /* tractions */
    sub_list.AddSub("micromorphic_FS_3D_natural_bc", ParameterListT::Any);
}



/* return the description of the given inline subordinate parameter list */
void FSMicromorphic3DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
				       SubListT& sub_lists) const
{
    ElementBaseT::DefineInlineSub(name, order, sub_lists);
}



/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSMicromorphic3DT::NewSub(const StringT& name) const
{
    /* create non-const this */
    FSMicromorphic3DT* non_const_this = const_cast<FSMicromorphic3DT*>(this);

    if (name == "micromorphic_FS_3D_natural_bc") /* traction bc */
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
    else if (name == "micromorphic_FS_3D_element_block")
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
void FSMicromorphic3DT::SetTractionBC(void)
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
void FSMicromorphic3DT::TakeNaturalBC(const ParameterListT& list)
{
    const char caller[] = "FSMicromorphic3DT::TakeTractionBC";

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
void FSMicromorphic3DT::ApplyTractionBC(void)
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

void FSMicromorphic3DT::Form_solid_shape_functions(const double* &shapes_displ_X)
{
    fShapeDispl = 0.0;
    for (int i=0; i<27; i++)
    {
		fShapeDispl(0,i*3) = shapes_displ_X[i];
		fShapeDispl(1,1+i*3) = shapes_displ_X[i];
		fShapeDispl(2,2+i*3) = shapes_displ_X[i];
    }
}

void FSMicromorphic3DT::Form_Gradient_of_solid_shape_functions(const dMatrixT &fShapeDisplGrad_temp)
{
    fShapeDisplGrad = 0.0;
    for(int i=0; i<27; i++)
    {
		fShapeDisplGrad(0,i*3) = fShapeDisplGrad_temp(0,i);
		fShapeDisplGrad(1,1+i*3) = fShapeDisplGrad_temp(0,i);
		fShapeDisplGrad(2,2+i*3) = fShapeDisplGrad_temp(0,i);

		fShapeDisplGrad(3,i*3) = fShapeDisplGrad_temp(1,i);
		fShapeDisplGrad(4,1+i*3) = fShapeDisplGrad_temp(1,i);
		fShapeDisplGrad(5,2+i*3) = fShapeDisplGrad_temp(1,i);

		fShapeDisplGrad(6,i*3) = fShapeDisplGrad_temp(2,i);
		fShapeDisplGrad(7,1+i*3) = fShapeDisplGrad_temp(2,i);
		fShapeDisplGrad(8,2+i*3) = fShapeDisplGrad_temp(2,i);
    }

}

void FSMicromorphic3DT::Form_micro_shape_functions(const double* &shapes_micro_X)
{
    fShapeMicro = 0.0;
    //hard coded for n_en_micro=8; can change
    for (int i=0; i<8; i++)
		fShapeMicro[i] = shapes_micro_X[i];
}

void FSMicromorphic3DT::Form_deformation_gradient_tensor(void)
{
    fShapeDisplGrad.Multx(u_vec,fGrad_disp_vector);
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
void FSMicromorphic3DT::Form_Grad_grad_transformation_matrix(void)
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

void FSMicromorphic3DT::Form_fDefGradT_9x9_matrix(void)
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


void FSMicromorphic3DT::Form_deformation_gradient_inv_vector(void)
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

void FSMicromorphic3DT::Form_kirchhoff_stress_vector()
{
    fKirchhoff_vector[0] = fKirchhoff_tensor(0,0);
    fKirchhoff_vector[1] = fKirchhoff_tensor(1,0);
    fKirchhoff_vector[2] = fKirchhoff_tensor(2,0);
    fKirchhoff_vector[3] = fKirchhoff_tensor(0,1);
    fKirchhoff_vector[4] = fKirchhoff_tensor(1,1);
    fKirchhoff_vector[5] = fKirchhoff_tensor(2,1);
    fKirchhoff_vector[6] = fKirchhoff_tensor(0,2);
    fKirchhoff_vector[7] = fKirchhoff_tensor(1,2);
    fKirchhoff_vector[8] = fKirchhoff_tensor(2,2);
}

void FSMicromorphic3DT::Form_Varpi_temp_matrix()
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
		    N_A_1I = fShapeDisplGradGrad(0,A-1);
		    N_A_2I = fShapeDisplGradGrad(5,A-1);
		    N_A_3I = fShapeDisplGradGrad(4,A-1);
		}break;
		case 2:
		{
		    N_A_1I = fShapeDisplGradGrad(5,A-1);
		    N_A_2I = fShapeDisplGradGrad(1,A-1);
		    N_A_3I = fShapeDisplGradGrad(3,A-1);
		}break;
		case 3:
		{
		    N_A_1I = fShapeDisplGradGrad(4,A-1);
		    N_A_2I = fShapeDisplGradGrad(3,A-1);
		    N_A_3I = fShapeDisplGradGrad(2,A-1);
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

void FSMicromorphic3DT::Form_Gradient_t_of_solid_shape_functions(const dMatrixT &fShapeDisplGrad_temp)
{
    fShapeDisplGrad_t = 0.0;
    for (int i=0; i<27; i++)
    {
	fShapeDisplGrad_t(0,i*3) = fShapeDisplGrad_temp(0,i);
	fShapeDisplGrad_t(1,i*3) = fShapeDisplGrad_temp(1,i);
	fShapeDisplGrad_t(2,i*3) = fShapeDisplGrad_temp(2,i);

	fShapeDisplGrad_t(3,1+i*3) = fShapeDisplGrad_temp(0,i);
	fShapeDisplGrad_t(4,1+i*3) = fShapeDisplGrad_temp(1,i);
	fShapeDisplGrad_t(5,1+i*3) = fShapeDisplGrad_temp(2,i);

	fShapeDisplGrad_t(6,2+i*3) = fShapeDisplGrad_temp(0,i);
	fShapeDisplGrad_t(7,2+i*3) = fShapeDisplGrad_temp(1,i);
	fShapeDisplGrad_t(8,2+i*3) = fShapeDisplGrad_temp(2,i);
    }
    
}

void FSMicromorphic3DT::Form_Im_temp_matrix()
{
    fIm_temp_matrix = 0.0;
    fIm_temp_matrix(0,0) = fKirchhoff_tensor(0,0);
    fIm_temp_matrix(1,0) = fKirchhoff_tensor(1,0);
    fIm_temp_matrix(2,0) = fKirchhoff_tensor(2,0);
    
    fIm_temp_matrix(3,1) = fKirchhoff_tensor(0,0);
    fIm_temp_matrix(4,1) = fKirchhoff_tensor(1,0);
    fIm_temp_matrix(5,1) = fKirchhoff_tensor(2,0);
    
    fIm_temp_matrix(6,2) = fKirchhoff_tensor(0,0);
    fIm_temp_matrix(7,2) = fKirchhoff_tensor(1,0);
    fIm_temp_matrix(8,2) = fKirchhoff_tensor(2,0);
    
    fIm_temp_matrix(0,3) = fKirchhoff_tensor(0,1);
    fIm_temp_matrix(1,3) = fKirchhoff_tensor(1,1);
    fIm_temp_matrix(2,3) = fKirchhoff_tensor(2,1);
    
    fIm_temp_matrix(3,4) = fKirchhoff_tensor(0,1);
    fIm_temp_matrix(4,4) = fKirchhoff_tensor(1,1);
    fIm_temp_matrix(5,4) = fKirchhoff_tensor(2,1);
    
    fIm_temp_matrix(6,5) = fKirchhoff_tensor(0,1);
    fIm_temp_matrix(7,5) = fKirchhoff_tensor(1,1);
    fIm_temp_matrix(8,5) = fKirchhoff_tensor(2,1);
    
    fIm_temp_matrix(0,6) = fKirchhoff_tensor(0,2);
    fIm_temp_matrix(1,6) = fKirchhoff_tensor(1,2);
    fIm_temp_matrix(2,6) = fKirchhoff_tensor(2,2);
    
    fIm_temp_matrix(3,7) = fKirchhoff_tensor(0,2);
    fIm_temp_matrix(4,7) = fKirchhoff_tensor(1,2);
    fIm_temp_matrix(5,7) = fKirchhoff_tensor(2,2);
    
    fIm_temp_matrix(6,8) = fKirchhoff_tensor(0,2);
    fIm_temp_matrix(7,8) = fKirchhoff_tensor(1,2);
    fIm_temp_matrix(8,8) = fKirchhoff_tensor(2,2);
}

void FSMicromorphic3DT::Form_Hbar_temp_matrix()
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

void FSMicromorphic3DT::Form_Ell_temp_matrix()
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


void FSMicromorphic3DT::Form_C_matrix(const double& J_Prim)
{
    double C_IJKL;
    int row,col;
    for (int I=0; I<3; I++)
	for (int J=0; J<3; J++)
	    for (int K=0; K<3; K++)
		for (int L=0; L<3; L++)
		{
		    C_IJKL = fMaterial_Params[kLambda]*fRight_Cauchy_Green_tensor_Inverse(J,I)
		    *fRight_Cauchy_Green_tensor_Inverse(L,K)+2*(fMaterial_Params[kMu]
		    -fMaterial_Params[kLambda]*log(J_Prim))*1/2*(fRight_Cauchy_Green_tensor_Inverse(I,K)
		    *fRight_Cauchy_Green_tensor_Inverse(J,L)+fRight_Cauchy_Green_tensor_Inverse(I,L)
		    *fRight_Cauchy_Green_tensor_Inverse(J,K));

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

void FSMicromorphic3DT::Form_c_matrix()
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

				    c_ijkl += fDeformation_Gradient(i,I)*fDeformation_Gradient(j,J)
				    *fDeformation_Gradient(k,K)*fDeformation_Gradient(l,L)*fC_matrix(row,col);
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

 
void FSMicromorphic3DT::Form_Im_Prim_temp_matrix()
{
    fIm_Prim_temp_matrix = 0.0;
    fIm_Prim_temp_matrix(0,0) = fKirchhoff_tensor(0,0);
    fIm_Prim_temp_matrix(3,0) = fKirchhoff_tensor(1,0);
    fIm_Prim_temp_matrix(6,0) = fKirchhoff_tensor(2,0);
    
    fIm_Prim_temp_matrix(1,1) = fKirchhoff_tensor(0,0);
    fIm_Prim_temp_matrix(4,1) = fKirchhoff_tensor(1,0);
    fIm_Prim_temp_matrix(7,1) = fKirchhoff_tensor(2,0);
    
    fIm_Prim_temp_matrix(2,2) = fKirchhoff_tensor(0,0);
    fIm_Prim_temp_matrix(5,2) = fKirchhoff_tensor(1,0);
    fIm_Prim_temp_matrix(8,2) = fKirchhoff_tensor(2,0);
    
    fIm_Prim_temp_matrix(0,3) = fKirchhoff_tensor(0,1);
    fIm_Prim_temp_matrix(3,3) = fKirchhoff_tensor(1,1);
    fIm_Prim_temp_matrix(6,3) = fKirchhoff_tensor(2,1);
    
    fIm_Prim_temp_matrix(1,4) = fKirchhoff_tensor(0,1);
    fIm_Prim_temp_matrix(4,4) = fKirchhoff_tensor(1,1);
    fIm_Prim_temp_matrix(7,4) = fKirchhoff_tensor(2,1);
    
    fIm_Prim_temp_matrix(2,5) = fKirchhoff_tensor(0,1);
    fIm_Prim_temp_matrix(5,5) = fKirchhoff_tensor(1,1);
    fIm_Prim_temp_matrix(8,5) = fKirchhoff_tensor(2,1);
    
    fIm_Prim_temp_matrix(0,6) = fKirchhoff_tensor(0,2);
    fIm_Prim_temp_matrix(3,6) = fKirchhoff_tensor(1,2);
    fIm_Prim_temp_matrix(6,6) = fKirchhoff_tensor(2,2);
    
    fIm_Prim_temp_matrix(1,7) = fKirchhoff_tensor(0,2);
    fIm_Prim_temp_matrix(4,7) = fKirchhoff_tensor(1,2);
    fIm_Prim_temp_matrix(7,7) = fKirchhoff_tensor(2,2);
    
    fIm_Prim_temp_matrix(2,8) = fKirchhoff_tensor(0,2);
    fIm_Prim_temp_matrix(5,8) = fKirchhoff_tensor(1,2);
    fIm_Prim_temp_matrix(8,8) = fKirchhoff_tensor(2,2);
}

void FSMicromorphic3DT::Form_D_matrix(void)
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

void FSMicromorphic3DT::Form_B_matrix(void)
{
    fB_matrix = 0.0;
    for(int i=0; i<27; i++)
    {
		fB_matrix(0,i*3)=fShapeDisplGrad_temp(0,i);
		fB_matrix(1,i*3+1)=fShapeDisplGrad_temp(1,i);
		fB_matrix(2,i*3+2)=fShapeDisplGrad_temp(2,i);
		fB_matrix(3,i*3+1)=fShapeDisplGrad_temp(2,i);
		fB_matrix(3,i*3+2)=fShapeDisplGrad_temp(1,i);
		fB_matrix(4,i*3)=fShapeDisplGrad_temp(2,i);
		fB_matrix(4,i*3+2)=fShapeDisplGrad_temp(0,i);
		fB_matrix(5,i*3)=fShapeDisplGrad_temp(1,i);
		fB_matrix(5,i*3+1)=fShapeDisplGrad_temp(0,i);
    }
}

void FSMicromorphic3DT::Extract_six_values_from_symmetric_tensor(const dMatrixT &fTensor,dArrayT& fTemp_six_values)
{
    fTemp_six_values[0]=fTensor(0,0);
    fTemp_six_values[1]=fTensor(1,1);
    fTemp_six_values[2]=fTensor(2,2);
    fTemp_six_values[3]=fTensor(1,2);
    fTemp_six_values[4]=fTensor(2,0);
    fTemp_six_values[5]=fTensor(0,1);
}

void FSMicromorphic3DT::Put_values_In_dArrayT_vector(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT)
{
	fArrayT[0]=f2DArrayT(e,IP*6+0);
	fArrayT[1]=f2DArrayT(e,IP*6+1);
	fArrayT[2]=f2DArrayT(e,IP*6+2);
	fArrayT[3]=f2DArrayT(e,IP*6+3);
	fArrayT[4]=f2DArrayT(e,IP*6+4);
	fArrayT[5]=f2DArrayT(e,IP*6+5);
}

void FSMicromorphic3DT::Form_gradv_vector(void)
{
    fShapeDisplGrad.Multx(u_dot_vec,fGradv_vector);
    fDefGradInv_Grad_grad.MultTx(fGradv_vector,fgradv_vector);
}

void FSMicromorphic3DT::Form_Xi_temp_matrix(void)
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

void FSMicromorphic3DT::Form_Varsigma_temp_matrix(void)
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

void FSMicromorphic3DT::Form_I_ijkl_matrix(void)
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




void FSMicromorphic3DT::Compute_norm_of_array(double& norm,const LocalArrayT& B)
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
