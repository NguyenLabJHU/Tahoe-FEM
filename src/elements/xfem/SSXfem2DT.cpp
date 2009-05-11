/* Id: SSXfem2DT.cpp,v 1.6 2006/10/10 19:55:23 regueiro Exp $ */
#include "SSXfem2DT.h"

#include "OutputSetT.h"
#include "ParameterContainerT.h"
#include "CommunicatorT.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
SSXfem2DT::SSXfem2DT(const ElementSupportT& support):
    ElementBaseT(support), //pass the solid displacement field to the base class
    u(LocalArrayT::kDisp),
    u_dot(LocalArrayT::kVel),
    u_dot_n(LocalArrayT::kLastVel),
    u_dotdot(LocalArrayT::kAcc),
    u_dotdot_n(LocalArrayT::kLastAcc),
    u_n(LocalArrayT::kLastDisp),
    ujump(LocalArrayT::kDisp),
    ujump_dot(LocalArrayT::kVel),
    ujump_dot_n(LocalArrayT::kLastVel),
    ujump_dotdot(LocalArrayT::kAcc),
    ujump_dotdot_n(LocalArrayT::kLastAcc),
    ujump_n(LocalArrayT::kLastDisp),
    fInitCoords_displ(LocalArrayT::kInitCoords),
    fCurrCoords_displ(LocalArrayT::kCurrCoords),
    fInitCoords_enhan(LocalArrayT::kInitCoords),
    fCurrCoords_enhan(LocalArrayT::kCurrCoords),
    fTractionBCSet(0),
    fDispl(NULL),
    fEnhan(NULL),
    fShapes_displ(NULL),
    fShapes_enhan(NULL),
    fKdd(ElementMatrixT::kNonSymmetric),
    fKdq(ElementMatrixT::kNonSymmetric),
    fKqd(ElementMatrixT::kNonSymmetric),
    fKqq(ElementMatrixT::kNonSymmetric),
    bStep_Complete(0)
{
    SetName("xfem_SS_2D");
}

/* destructor */
SSXfem2DT::~SSXfem2DT(void) 
{  
    delete fShapes_displ;
    delete fShapes_enhan;
}


void SSXfem2DT::Echo_Input_Data(void) 
{

    cout << "#######################################################" << endl; 
    cout << "############### ECHO SSXfem2D DATA #########################" << endl; 
    cout << "#######################################################" << endl; 

    //################## material parameters ##################
    cout << "iConstitutiveModelType " 				<< iConstitutiveModelType 	<< endl; 

    //-- Elasticity parameters for solid
    cout << "fMaterial_Params[kMu] "  				<< fMaterial_Params[kMu] 	 << endl;
    cout << "fMaterial_Params[kLambda] "  			<< fMaterial_Params[kLambda] << endl;
	
}


//---------------------------------------------------------------------

void SSXfem2DT::RHSDriver(void)	
{
    int curr_group = ElementSupport().CurrentGroup();

    /* traction boundary conditions acting on displacement equations */
    if (curr_group == fDispl->Group()) 
	ApplyTractionBC();

    /* choose solution method */
    if (fDispl->Group() == fEnhan->Group())
		RHSDriver_monolithic();
    else
		RHSDriver_staggered();
}
//---------------------------------------------------------------------

void SSXfem2DT::Equations(AutoArrayT<const iArray2DT*>& eq_d,
				 AutoArrayT<const RaggedArray2DT<int>*>& eq_q)
{
    /* doing monolithic solution */
    if (fDispl->Group() == fEnhan->Group())
    {
		int ndof_enhan = fEnhan->NumDOF();
		int ndof_displ = fDispl->NumDOF();
		
		/* loop over connectivity blocks */
		fEqnos_displ.Dimension(fEqnos.Length());
		fEqnos_enhan.Dimension(fEqnos.Length());
		for (int i = 0; i < fEqnos.Length(); i++)
		{
		    /* connectivities */
		    const iArray2DT& connects_displ = *(fConnectivities_displ[i]);
		    const iArray2DT& connects_enhan = *(fConnectivities_enhan[i]);
		    int nel = connects_displ.MajorDim();
			
		    /* dimension */ 
		    fEqnos[i].Dimension(nel, n_en_displ*ndof_displ + n_en_enhan*ndof_enhan);
		    iArray2DT& displ_eq = fEqnos_displ[i];
		    iArray2DT& enhan_eq = fEqnos_enhan[i];
		    displ_eq.Dimension(nel, n_en_displ*ndof_displ);
		    enhan_eq.Dimension(nel, n_en_enhan*ndof_enhan);
				
		    /* get equation numbers */
		    fDispl->SetLocalEqnos(connects_displ, displ_eq);
		    fEnhan->SetLocalEqnos(connects_enhan, enhan_eq);

		    /* write into one array */
		    fEqnos[i].BlockColumnCopyAt(displ_eq, 0);
		    fEqnos[i].BlockColumnCopyAt(enhan_eq, displ_eq.MinorDim());

		    /* add to list of equation numbers */
		    eq_d.Append(&fEqnos[i]);
		}
		
		/* reset pointers to element cards */
		SetElementCards(fBlockData, fConnectivities_displ, fEqnos_displ, fElementCards_displ);
		SetElementCards(fBlockData, fConnectivities_enhan, fEqnos_enhan, fElementCards_enhan);
    }
    else
	/* doing staggered */
    {
#pragma message("initialization for staggered solution needs to be corrected")
	
		/* ElementBaseT handles equation array for displacements */
		if (ElementSupport().CurrentGroup() == fDispl->Group())
		    ElementBaseT::Equations(eq_d, eq_q);

		/* enhan-displacement-gradient equation */
		if (ElementSupport().CurrentGroup() == fEnhan->Group())
		{
		    /* collect local equation numbers */
		    //fPress.SetLocalEqnos(fConnectivities_enhan, fEqnos_enhan);
			
		    //eq_d.Append(&fEqnos_enhan);
		}
    }
	
    /* get the equation number for the nodes on the faces */
    /*
    for (int i = 0; i < fEnhanFaceEqnos.Length(); i++)
    {
		iArray2DT& faces = fEnhanFaces[i];
		iArray2DT& eqnos = fEnhanFaceEqnos[i];
		eqnos.Dimension(faces.MajorDim(), faces.MajorDim()*fDispl->NumDOF());
		
		fDispl->SetLocalEqnos(faces, eqnos);
    }
    */
}


//---------------------------------------------------------------------

void SSXfem2DT::LHSDriver(GlobalT::SystemTypeT)
{
/** Everything done in RHSDriver for efficiency */
//cout << "############### In LHS Driver ############### \n";
}

//---------------------------------------------------------------------

void SSXfem2DT::Select_Equations (const int &iBalLinChoice, const int &iBalLinEnhChoice )
{
    /** Choices for Linear Momentum Balance Equation */

    switch ( iBalLinChoice )	{

    default :
	cout << "SSXfem2DT::Select_Equations() .. currently only one linear momentum balance for enhanced xfem continuum \n";
	break;
    }

    /** Choices for Linear Momentum Balance Enhanced Equation */

    switch ( iBalLinEnhChoice )	{

    default :
	cout << "SSXfem2DT::Select_Equations() .. currently only one linear momentum balance equation for enhanced xfem continuum \n";
	break;
    }

}

//---------------------------------------------------------------------

/* return true if the element contributes to the solution of the
 * given group. */
bool SSXfem2DT::InGroup(int group) const
{
    return group == fDispl->Group() || group == fEnhan->Group();
}

//---------------------------------------------------------------------


/* initialize/finalize step */
void SSXfem2DT::InitStep(void)
{
	/* inherited */
	ElementBaseT::InitStep();
}


/* close current time increment */
void SSXfem2DT::CloseStep(void)
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
	int num_node_output = fDispl->NumDOF() + fEnhan->NumDOF() + knumstrain + knumstress + knum_d_state;
	dArray2DT n_values(nodes_used.Length(), num_node_output);
	
	/* collect nodal values */
	const dArray2DT& fUjump = (*fEnhan)[0];
	const dArray2DT& fU = (*fDispl)[0];
	for (int i = 0; i < nodes_used.Length(); i++)
	{
	    int node = nodes_used[i];
	    double* row = n_values(i);
	    for (int j = 0; j < fUjump.MinorDim(); j++)
		*row++ = fUjump(node,j);
	    
	    for (int j = 0; j < fU.MinorDim(); j++)
		*row++ = fU(node,j);
	    
	    double* p_stress = extrap_values(i);
	    for (int j = 0; j < (knumstrain+knumstress+knum_d_state); j++)
		*row++ = p_stress[j];
	}

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, fIPVariable);
	
    }

	//will not need this for quasi-static xfem
    /* zero first derivative of fields which are created at time=0 during calculating geostatic equilibrium (Trapezoidal rule) */    
    /*
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==1)
    {
		FieldT* fenhan = const_cast <FieldT*> (fEnhan);
		FieldT* fdispl = const_cast <FieldT*> (fDispl);

		(*fdispl)[1] = 0;
		(*fenhan)[1] = 0;
    }   
    */

	//will not need this for quasi-static xfem
    /* zero second derivative of fields which are created at time=0 during calculating geostatic equilibrium(Newmark method) */
    /*
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==2)
    {
		FieldT* fenhan = const_cast <FieldT*> (fEnhan);
		FieldT* fdispl = const_cast <FieldT*> (fDispl);
		
		(*fdispl)[2] = 0;
		(*fenhan)[2] = 0;
    }
    */

    /* reassign initial 2nd time derivative of enhan-displacement-gradient to 1st derivative */
    /*
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==2)
    {
		FieldT* fenhan = const_cast <FieldT*> (fEnhan);
		(*fenhan)[1] = (*fenhan)[2];
		(*fenhan)[2] = 0;
    }
    */
   
    /* store more recently updated values */
    fdState = fdState_new;
    fiState = fiState_new;
	
//    ss_xfem2D_out	<< endl 
//		<< setw(outputFileWidth) << "time_step"
//		<< endl;
    step_number = ElementSupport().StepNumber();
//    ss_xfem2D_out	<< setw(outputFileWidth) << step_number
//		<< endl;
//    ss_xfem2D_out	<< endl << "**********************************************************************************************";
//    ss_xfem2D_out	<< endl << "**********************************************************************************************" << endl;
}


/* resets to the last converged solution */
/*
GlobalT::RelaxCodeT SSXfem2DT::ResetStep(void)
{
	const char caller[] = "SSXfem2DT::ResetStep";
	
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
GlobalT::RelaxCodeT SSXfem2DT::RelaxSystem(void)
{
	const char caller[] = "SSXfem2DT::RelaxSystem";
	
	// inherited 
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	// loop over materials 
	//needs to be implemented
#pragma message("relax step for materials not implemented")	
	//ExceptionT::GeneralFail(caller, "relax step for materials not implemented");

	return relax;
}
*/


void SSXfem2DT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented
}


/* return geometry and number of nodes on each facet */
void SSXfem2DT::FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry, 
	iArrayT& num_facet_nodes) const
{
	/* from integration domain */
	ShapeFunctionDispl().FacetGeometry(facet_geometry, num_facet_nodes);
}


/* form of tangent matrix */
GlobalT::SystemTypeT SSXfem2DT::TangentType(void) const
{
    return GlobalT::kNonSymmetric; 
}

/*
void SSXfem2DT::SetStatus(const ArrayT<ElementCardT::StatusT>& status)
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
void SSXfem2DT::InitialCondition(void)
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
void SSXfem2DT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
    const char caller[] = "SSXfem2DT::AddNodalForce";

    /* displ, enhan, or neither */
    bool is_displ = false;
    dArrayT* element_force = NULL;
    int num_force = 0;
    if (field.FieldName() == fDispl->FieldName()) 
    {
		is_displ = true;
		element_force = &fFd_int;
		num_force = fDispl->NumDOF();
    }
    else if (field.FieldName() == fEnhan->FieldName()) 
    {
		is_displ = false;
		element_force = &fFq_int;
		num_force = fEnhan->NumDOF();
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
		    const iArrayT& nodes_enhan = fElementCards_enhan[e].NodesU();

		    u.SetLocal(nodes_displ);
		    u_n.SetLocal(nodes_displ);
		    ujump.SetLocal(nodes_enhan);
		    ujump_n.SetLocal(nodes_enhan);

		    del_u.DiffOf (u, u_n);
		    del_ujump.DiffOf (ujump, ujump_n);

		    // calculate derivatives based on reference coordinates
		    fInitCoords_displ.SetLocal(fElementCards_displ[e].NodesX());
		    //fCurrCoords_displ.SetToCombination (1.0, fInitCoords_displ, 1.0, u); 
		    fCurrCoords_displ=fInitCoords_displ;
		    fShapes_displ->SetDerivatives_DN_DDN(); 

		    //
		    fInitCoords_enhan.SetLocal(fElementCards_enhan[e].NodesX());
		    fCurrCoords_enhan=fInitCoords_enhan;
		    //fCurrCoords_enhan.SetToCombination (1.0, fInitCoords_enhan, 1.0, u); 
		    fShapes_enhan->SetDerivatives();
			
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
				/* residual for enhan-displacement-gradient field */ 
				// generate this vector fFq_int
				fShapes_displ->TopIP();
				while (fShapes_displ->NextIP())
				{
				    //nothing right now
				    fFq_int=0.0;
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
		    else /* enhan-displacement-gradient nodal dof */
			{
			    for (int i = 0; i < nodes_enhan.Length(); i++)
			    {
					if (nodes_enhan[i] == node)
					{
					    /* components for node */
					    nodalforce.Set(num_force, element_force->Pointer(dex));

					    /* accumulate */
					    force += nodalforce;
					}
					dex += fEnhan->NumDOF();
			    }		
			}
		}
    }
//	cout << "F_int = \n" << fFd_int << endl;
}

//---------------------------------------------------------------------

double SSXfem2DT::InternalEnergy ( void )
{
//not implemented
    return 0.0;
}

//---------------------------------------------------------------------

/* write restart data to the output stream */
void SSXfem2DT::WriteRestart(ostream& out) const
{
    /* inherited */
    ElementBaseT::WriteRestart(out);

    /* write state variable data */
    out << fdState;
}

//---------------------------------------------------------------------

/* read restart data to the output stream */
void SSXfem2DT::ReadRestart(istream& in)
{
    /* inherited */
    ElementBaseT::ReadRestart(in);

    /* write state variable data */
    in >> fdState;
}

//---------------------------------------------------------------------

void SSXfem2DT::RegisterOutput(void)
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
    const char* slabels3D[] = {"s11", "s22","s12","e11","e22","e12"};
    // state variables; ?
    const char* svlabels3D[] = {"thing1","thing2"};
    int count = 0;
    for (int j = 0; j < fNumIP_enhan; j++)
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
    int num_node_output = fDispl->NumDOF() + fEnhan->NumDOF() + knumstrain + knumstress + knum_d_state;
    ArrayT<StringT> n_labels(num_node_output);
    count = 0;

    /* labels from enhan-displacement-gradient */
    const ArrayT<StringT>& enhan_labels = fEnhan->Labels();
    for (int i = 0; i < enhan_labels.Length(); i++)
	n_labels[count++] = enhan_labels[i];

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
#pragma message("SSXfem2DT::RegisterOutput: is this right? ")
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

void SSXfem2DT::WriteOutput(void)
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
    int num_node_output = fDispl->NumDOF() + fEnhan->NumDOF() + knumstrain + knumstress + knum_d_state;
    dArray2DT n_values(nodes_used.Length(), num_node_output);

    /* collect nodal values */
    const dArray2DT& fUjump = (*fEnhan)[0];
    const dArray2DT& fU = (*fDispl)[0];
    for (int i = 0; i < nodes_used.Length(); i++)
    {
		int node = nodes_used[i];
		double* row = n_values(i);
		for (int j = 0; j < fUjump.MinorDim(); j++)
		    *row++ = fUjump(node,j);

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
void SSXfem2DT::RHSDriver_staggered(void)
{
	const char caller[] = "SSXfem2DT::RHSDriver_staggered";
#pragma message("staggered solution not implemented")
}

//---------------------------------------------------------------------
/* form group contribution to the stiffness matrix and RHS */
void SSXfem2DT::RHSDriver_monolithic(void)
{
    const char caller[] = "SSXfem2DT::RHSDriver_monolithic";
    if (fDispl->Group() != fEnhan->Group())
	ExceptionT::GeneralFail(caller, "displacement and enhan-displacement-gradient groups must be the same: %d != %d",
				fDispl->Group(), fEnhan->Group());

    int curr_group = ElementSupport().CurrentGroup();

    /* stress output work space */
    dArray2DT	out_variable_all, fdstatenew_all, fdstate_all;
    dArrayT		out_variable;

    /* time Step Increment */
    double delta_t = ElementSupport().TimeStep();
    time = ElementSupport().Time();
    step_number = ElementSupport().StepNumber();

    /* print time */
//    ss_xfem2D_out	<<"delta_t "<<delta_t << endl ;
//    ss_xfem2D_out	<<"time "<<time << endl ;

    /* loop over elements */
    int e,l;
    Top();

//   ss_xfem2D_out	<<"kInitialConditionType "<<kInitialConditionType << endl ;
//   ss_xfem2D_out	<<"kAnalysisType "<<kAnalysisType << endl ;

	//don't need now for xfem
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
		FieldT* fenhan = const_cast <FieldT*> (fEnhan);
		FieldT* fdispl = const_cast <FieldT*> (fDispl);
		(*fenhan)[0] = (*fenhan)[1];
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
		FieldT* fenhan = const_cast <FieldT*> (fEnhan);
		FieldT* fdispl = const_cast <FieldT*> (fDispl);
		(*fenhan)[0] = (*fenhan)[2];
		(*fdispl)[0] = (*fdispl)[2];
    } 
    */
   
    
    while (NextElement())
    {
    /* initialize */
	fK_dd_BTDB_matrix = 0.0;
		
	e = CurrElementNumber();
	const iArrayT& nodes_displ = fElementCards_displ[e].NodesU();
	const iArrayT& nodes_enhan = fElementCards_enhan[e].NodesU();

	u.SetLocal(nodes_displ);
	u_n.SetLocal(nodes_displ);
	if (u_dot.IsRegistered()) u_dot.SetLocal(nodes_displ);
	if (u_dot_n.IsRegistered()) u_dot_n.SetLocal(nodes_displ);
	if (u_dotdot.IsRegistered()) u_dotdot.SetLocal(nodes_displ);
	if (u_dotdot_n.IsRegistered())u_dotdot_n.SetLocal(nodes_displ);

	ujump.SetLocal(nodes_enhan);
	ujump_n.SetLocal(nodes_enhan);
	if (ujump_dot.IsRegistered()) ujump_dot.SetLocal(nodes_enhan);
	if (ujump_dot_n.IsRegistered()) ujump_dot_n.SetLocal(nodes_enhan);
	if (ujump_dotdot.IsRegistered()) ujump_dotdot.SetLocal(nodes_enhan);
	if (ujump_dotdot_n.IsRegistered()) ujump_dotdot_n.SetLocal(nodes_enhan);

	/* print solid displacement at current step (u)*/
/*	ss_xfem2D_out	<<"nodal solid displacement at current step(u)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    ss_xfem2D_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		ss_xfem2D_out << u(i,j) << "\t";
	    ss_xfem2D_out	<< endl ;
	}
*/

	/* print solid displacement from previous step (u_n)*/
/*	ss_xfem2D_out	<<"nodal solid displacement from previous step(u_n)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    ss_xfem2D_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		ss_xfem2D_out << u_n(i,j) << "\t";
	    ss_xfem2D_out	<< endl ;
	}
*/
	
	/* print solid velocity at current step (u_dot)*/
/*	ss_xfem2D_out	<<"nodal solid velocity at current step(u_dot)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    ss_xfem2D_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		ss_xfem2D_out << u_dot(i,j) << "\t";
	    ss_xfem2D_out	<< endl ;
	}
*/

	/* print solid velocity from previous step (u_dot_n)*/
/*	ss_xfem2D_out	<<"nodal solid velocity from previous step(u_dot_n)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    ss_xfem2D_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		ss_xfem2D_out << u_dot_n(i,j) << "\t";
	    ss_xfem2D_out	<< endl ;
	}
*/


	/* print solid acceleration at current step (u_dotdot)*/
/*	ss_xfem2D_out	<<"nodal solid velocity at current step(u_dotdot)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    ss_xfem2D_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		ss_xfem2D_out << u_dotdot(i,j) << "\t";
	    ss_xfem2D_out	<< endl ;
	}
*/

	/* print solid acceleration from previous step (u_dotdot_n)*/
/*	ss_xfem2D_out	<<"nodal solid velocity from previous step(u_dotdot_n)"<< endl ;
	for (int i=0; i<n_en_displ; i++)
	{
	    ss_xfem2D_out	<< "node number " << i+1 <<" :  " ;		
	    for (int j=0; j<n_sd; j++)
		ss_xfem2D_out << u_dotdot_n(i,j) << "\t";
	    ss_xfem2D_out	<< endl ;
	}
*/
		 
	/* print enhan-displacement-gradient at current step (enhan)*/
/*	ss_xfem2D_out	<<"nodal enhan-displacement-gradient at current step(enhan)"<< endl ;
	for (int i=0; i<n_en_enhan; i++)
	{
	    ss_xfem2D_out	<< "node number " << i+1 <<" :  " ;	
	    ss_xfem2D_out	<< ujump(i,0) << endl;
	}
*/

	/* print enhan-displacement-gradient from previous step (enhan_n)*/
/*	ss_xfem2D_out	<<"nodal enhan-displacement-gradient from previous step(enhan_n)"<< endl ;
	for (int i=0; i<n_en_enhan; i++)
	{
	    ss_xfem2D_out	<< "node number " << i+1 <<" :  " ;	
	    ss_xfem2D_out	<< ujump_n(i,0) << endl;
	}
*/
	
	/* print first derivative of enhan-displacement-gradient at current step (enhan_dot)*/
/*	ss_xfem2D_out	<<"first derivative of nodal enhan-displacement-gradient at current step(enhan_dot)"<< endl ;
	for (int i=0; i<n_en_enhan; i++)
	{
	    ss_xfem2D_out	<< "node number " << i+1 <<" :  " ;	
	    ss_xfem2D_out	<< ujump_dot(i,0) << endl;
	}
*/

	/* print first derivative of enhan-displacement-gradient from previous step (enhan_dot_n)*/
/*	ss_xfem2D_out	<<"first derivative of nodal enhan-displacement-gradient from previous step(enhan_dot_n)"<< endl ;
	for (int i=0; i<n_en_enhan; i++)
	{
	    ss_xfem2D_out	<< "node number " << i+1 <<" :  " ;	
	    ss_xfem2D_out	<< ujump_dot_n(i,0) << endl;
	}
*/


	/* print second derivative of enhan-displacement-gradient at current step (enhan_dotdot)*/
/*	ss_xfem2D_out	<<"second derivative of nodal enhan-displacement-gradient at current step(enhan_dotdot)"<< endl ;
	for (int i=0; i<n_en_enhan; i++)
	{
	    ss_xfem2D_out	<< "node number " << i+1 <<" :  " ;	
	    ss_xfem2D_out	<< ujump_dotdot(i,0) << endl;
	}
*/

	/* print second derivative of enhan-displacement-gradient from previous step (enhan_dotdot_n)*/
/*	ss_xfem2D_out	<<"second derivative of nodal enhan-displacement-gradient from previous step(enhan_dotdot_n)"<< endl ;
	for (int i=0; i<n_en_enhan; i++)
	{
	    ss_xfem2D_out	<< "node number " << i+1 <<" :  " ;	
	    ss_xfem2D_out	<< ujump_dotdot_n(i,0) << endl;
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
			//u_dot_vec[index_u] = u_dot(i,j);
			//u_dotdot_vec[index_u] = u_dotdot(i,j);
			index_u += 1;
	    }
	}

	/* populate enhan-displacement-gradient, first and second time derivatives of 
	   enhan-displacement-gradient in vector form*/
	int index_ujump = 0;
	for (int i=0; i<n_en_enhan; i++)
	{
	    for (int j=0; j<ndof_per_nd_enhan; j++)
	    {
			ujump_vec[index_ujump] = ujump(i,j);
			//ujump_dot_vec[index_ujump] = ujump_dot(i,j);
			//ujump_dotdot_vec[index_ujump] = ujump_dotdot(i,j);
			index_ujump += 1;
	    }
	}

	del_u.DiffOf (u, u_n);
	del_ujump.DiffOf (ujump, ujump_n);
	
	// calculate derivatives based on reference coordinates
	fInitCoords_displ.SetLocal(fElementCards_displ[e].NodesX());
	fCurrCoords_displ=fInitCoords_displ;
	//fCurrCoords_displ.SetToCombination (1.0, fInitCoords_displ, 1.0, u); 
	//fShapes_displ->SetDerivatives_DN_DDN(); 
	fShapes_displ->SetDerivatives();
	//
	fInitCoords_enhan.SetLocal(fElementCards_enhan[e].NodesX());
	fCurrCoords_enhan=fInitCoords_enhan;
	//fCurrCoords_enhan.SetToCombination (1.0, fInitCoords_enhan, 1.0, u); 
	fShapes_enhan->SetDerivatives(); 
	
	//update state variables
	fdstatenew_all.Alias(fNumIP_displ, knum_d_state, fdState_new(CurrElementNumber()));
	fdstate_all.Alias(fNumIP_displ, knum_d_state, fdState(CurrElementNumber()));
	
	if (bStep_Complete) 
	{ 
	    //-- Store/Register data in classic tahoe manner 
	    out_variable_all.Alias(fNumIP_displ, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
	    for (l=0; l < fNumIP_displ; l++) 
	    {
			out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(l));
			out_variable[0]=fStress_Elements_IPs(e,l*3+0);
			out_variable[1]=fStress_Elements_IPs(e,l*3+1);
			out_variable[2]=fStress_Elements_IPs(e,l*3+2);
			out_variable[3]=fStrain_Elements_IPs(e,l*3+0);
			out_variable[4]=fStrain_Elements_IPs(e,l*3+1);
			out_variable[5]=fStrain_Elements_IPs(e,l*3+2);
			/* 
			// if 3 state variables
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
		    fShapes_enhan->TopIP();
		    
		    while (fShapes_displ->NextIP() && fShapes_enhan->NextIP())
		    {
				double scale_const = (*Weight++)*(*Det++);
			
				const int IP = fShapes_displ->CurrIP();	
				dArrayT DisplIPCoordinate(n_sd), EnhanIPCoordinate(n_sd);
				fShapes_displ->IPCoords(DisplIPCoordinate);
				fShapes_enhan->IPCoords(EnhanIPCoordinate);

				const double* shapes_displ_X = fShapes_displ->IPShapeX();
				/* [fShapeDispl]will be formed */
				Form_solid_shape_functions(shapes_displ_X);
				
				fShapes_displ->GradNa(fShapeDisplGrad);
				/* [fShapeDisplGrad] will be formed */
				//Form_Gradient_of_solid_shape_functions(fShapeDisplGrad);
				
				const double* shapes_enhan_X = fShapes_enhan->IPShapeX();
				/* {fShapeEnhan} will be formed */
				Form_enhan_shape_functions(shapes_enhan_X);
				
				fShapes_enhan->GradNa(fShapeEnhanGrad);		
				/* [fShapeEnhanGrad] will be formed */
				//Form_Gradient_of_enhance_shape_functions(fShapeEnhanGrad);
				
				/* [fD_matrix] will be formed */
				Form_D_matrix();
			
				/* [fB_matrix] will be formed */
				Form_B_matrix();
				
				/* [fK_dd_BTDB_matrix] will be formed */
				fTemp_matrix_ndof_se_x_ndof_se.MultATBC(fB_matrix,fD_matrix,fB_matrix);
				double scale = scale_const;
				fTemp_matrix_ndof_se_x_ndof_se *= scale;
				/* accumulate */
				fK_dd_BTDB_matrix += fTemp_matrix_ndof_se_x_ndof_se;
				
				/* [fStrain_vector_current_IP] will be formed */
				fB_matrix.Multx(u_vec,fStrain_vector_current_IP);
				/* Save strain tensor of the current IP */ 
				fStrain_IPs.SetRow(IP,fStrain_vector_current_IP);
				
				/* [fStress_vector_current_IP] will be formed */
				fD_matrix.Multx(fStrain_vector_current_IP,fStress_vector_current_IP);
				/* Save stress vector of the current IP */ 
				fStress_IPs.SetRow(IP,fStress_vector_current_IP); 
			

				/* for debugging */
				/*
				const int ip = fShapes_displ->CurrIP()+1;
				ss_xfem2D_out	<< endl << "IP" << ip
						<< setw(outputFileWidth) << ", shape function matrix for solid phase: " 
						<< setw(outputFileWidth) << fShapeDispl;
				
				ss_xfem2D_out	<< endl << "terms from shape function matrix for solid phase: " 
						<< setw(outputFileWidth) << fShapeDispl(0,0) 
						<< setw(outputFileWidth) << fShapeDispl(0,1);	
				*/
			
		    } //end Gauss integration loop
		    		    
			/* saving strain for each IPs of the current element */
		    fStrain_Elements_IPs.SetRow(e,fStrain_IPs);
	    
		    /* saving stress for each IPs of the current element */
		    fStress_Elements_IPs.SetRow(e,fStress_IPs);

		    /* saving state variables for each IPs of the current element */
		    fState_variables_Elements_IPs.SetRow(e,fState_variables_IPs);

			//stiffness matrix
			fKdd = fK_dd_BTDB_matrix;
		    fK_dd_BTDB_matrix.MultTx(u_vec,fTemp_vector_ndof_se);
		    //internal force vector
			fFd_int = fTemp_vector_ndof_se;
			fFd_int *= -1.0; 
			
		    /* [fKdq] will be formed */
		    //need to code
		    fKdq = 0.0;
	    
		    /* [fKqd] will be formed */
		    //need to code
		    fKqd = 0.0;

		    /* [fKqq] will be formed */
		    //need to code
		    fKqq = 0.0;	

		    /* {fFq_int} will be formed */
		    //need to code
		    fFq_int = 0.0;

		    /* equations numbers */
		    const iArrayT& displ_eq = fElementCards_displ[e].Equations();
		    const iArrayT& enhan_eq = fElementCards_enhan[e].Equations();
		    
		    /* assemble residuals */
		    ElementSupport().AssembleRHS(curr_group, fFd_int, displ_eq);
		    ElementSupport().AssembleRHS(curr_group, fFq_int, enhan_eq);
		    
		    /* assemble components of the tangent */
		    ElementSupport().AssembleLHS(curr_group, fKdd, displ_eq);
		    ElementSupport().AssembleLHS(curr_group, fKqq, enhan_eq);
		    ElementSupport().AssembleLHS(curr_group, fKdq, displ_eq, enhan_eq);
		    ElementSupport().AssembleLHS(curr_group, fKqd, enhan_eq, displ_eq);

	}	
    }
}



/* form global shape function derivatives */
void SSXfem2DT::SetGlobalShape(void)
{
    /* fetch (initial) coordinates */
    SetLocalX(fLocInitCoords);
	
    /* compute shape function derivatives */
    fShapes_displ->SetDerivatives_DN_DDN();
    fShapes_enhan->SetDerivatives();
}


/* describe the parameters needed by the interface */
void SSXfem2DT::DefineParameters(ParameterListT& list) const
{
    /* inherited */
    ElementBaseT::DefineParameters(list);

    /* displacement field */
    //already done in ElementBaseT
    //list.AddParameter(ParameterT::Word, "displ_field_name");
	
    /* enhan-displacement-gradient field */
    list.AddParameter(ParameterT::Word, "enhan_field_name");
	
    list.AddParameter(fGeometryCode_displ_int, "GeometryCode_displ");
    list.AddParameter(fNumIP_displ, "NumIP_displ");
    list.AddParameter(fGeometryCodeSurf_displ_int, "GeometryCodeSurf_displ");
    list.AddParameter(fNumIPSurf_displ, "NumIPSurf_displ");
    list.AddParameter(n_en_displ, "n_en_displ");
    list.AddParameter(n_en_enhan, "n_en_enhan");
    list.AddParameter(ndof_per_nd_enhan, "ndof_per_nd_enhan");
	
    list.AddParameter(iConstitutiveModelType, "constitutive_mod_type");

    double shearMu, sLambda;

    // solid elasticity
    list.AddParameter(shearMu, "mu");
    list.AddParameter(sLambda, "lambda");
	
    // Newmark time integration parameters
//    list.AddParameter(newBeta, "beta");
//    list.AddParameter(newGamma, "gamma");

}


/* accept parameter list */
void SSXfem2DT::TakeParameterList(const ParameterListT& list)
{
    const char caller[] = "SSXfem2DT::TakeParameterList";
	
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

    /* get enhan-displacement-gradient field */
    const StringT& enhan_field_name = list.GetParameter("enhan_field_name");
    fEnhan = ElementSupport().Field(enhan_field_name);
    if (!fEnhan)
	ExceptionT::GeneralFail(caller, "could not resolve \"%s\" enhan_field", 
				enhan_field_name.Pointer());

    fGeometryCode_displ_int = list.GetParameter("GeometryCode_displ");
    fGeometryCode_displ = GeometryT::int2CodeT(fGeometryCode_displ_int);
    fNumIP_displ = list.GetParameter("NumIP_displ");
    fGeometryCodeSurf_displ_int = list.GetParameter("GeometryCodeSurf_displ");
    fGeometryCodeSurf_displ = GeometryT::int2CodeT(fGeometryCodeSurf_displ_int);
    fNumIPSurf_displ = list.GetParameter("NumIPSurf_displ");
    n_en_displ = list.GetParameter("n_en_displ");
    n_en_enhan = list.GetParameter("n_en_enhan");
    ndof_per_nd_enhan = list.GetParameter("ndof_per_nd_enhan");
    
    fGeometryCode_enhan = fGeometryCode_displ; 
    fNumIP_enhan = fNumIP_displ;
    fGeometryCodeSurf_enhan = fGeometryCodeSurf_displ;
    fNumIPSurf_enhan = fNumIPSurf_displ;
	
    iConstitutiveModelType = list.GetParameter("constitutive_mod_type");
	
    fMaterial_Params.Dimension ( kNUM_FMATERIAL_TERMS );
//    fIntegration_Params.Dimension ( kNUM_FINTEGRATE_TERMS );
	
    fMaterial_Params[kMu] = list.GetParameter("mu");
    fMaterial_Params[kLambda] = list.GetParameter("lambda");
	
//    fIntegration_Params[kBeta] = list.GetParameter("beta");
//    fIntegration_Params[kGamma] = list.GetParameter("gamma");
	
    Echo_Input_Data();
	
	//need to change for appropriate number of ISVs for xfem model
    knum_d_state = 2; // #? internal state variables
    knum_i_state = 0; // int's needed per ip, state variables
	
	//need to change these for non-symmetric stress, and higher order couple stress output
    knumstrain = 3; // number of strain outputs
    knumstress = 3; // number of stress outputs
	
    output = "out";
	
    /* dimensions (notation as per Hughes' Book) */
    int& n_ip_displ = fNumIP_displ;
    int& n_ip_enhan = fNumIP_enhan;
    n_sd = NumSD();
    int nen = NumElementNodes(); /* number of nodes/element in the mesh */

    /* initialize connectivities */
    fConnectivities_displ.Alias(fConnectivities);
    fConnectivities_enhan.Alias(fConnectivities);

    /* pick element interpolations based on available number of element nodes
     * and the specified number of integration points */
    // only implemented for 2D, quadrilaterals
    if (n_sd == 2 && n_en_enhan != n_en_displ && fGeometryCode_displ == GeometryT::kQuadrilateral) 
    {
		// don't expect reduced integration for both fields 
		// if (n_ip_displ == 4 && n_ip_enhan == 4)
		//	ExceptionT::GeneralFail(caller, "not expecting 4 ips for both fields");
		//else if (n_ip_displ == 4 || n_ip_enhan == 4) // create reduced connectivities
		//{ 
		// reduce the number of element nodes based on the number ip's
		int& nen_red = (n_ip_displ == 4) ? n_en_displ : n_en_enhan;
		nen_red = 4;
		ArrayT<const iArray2DT*>& connects_red = (n_ip_displ == 4) ? 
		    fConnectivities_displ : 
		    fConnectivities_enhan;
			
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
    
    // q
    fInitCoords_enhan.Dimension(n_en_enhan, n_sd);
    ElementSupport().RegisterCoordinates(fInitCoords_enhan);	
    fCurrCoords_enhan.Dimension(n_en_enhan, n_sd);
    fShapes_enhan = new ShapeFunctionT(fGeometryCode_enhan, fNumIP_enhan, fCurrCoords_enhan);
    //fShapes_enhan = new ShapeFunctionT(fGeometryCode_enhan, fNumIP_enhan, fCurrCoords_displ);
    //fShapes_enhan->Initialize();

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

    /* set local arrays for enhan-displacement-gradient field */
    ujump.Dimension (n_en_enhan, ndof_per_nd_enhan);
    ujump_dot.Dimension (n_en_enhan, ndof_per_nd_enhan);
    ujump_dot_n.Dimension (n_en_enhan, ndof_per_nd_enhan);
    ujump_dotdot.Dimension (n_en_enhan, ndof_per_nd_enhan);
    ujump_dotdot_n.Dimension (n_en_enhan, ndof_per_nd_enhan);
    ujump_n.Dimension (n_en_enhan, ndof_per_nd_enhan);
    del_ujump.Dimension (n_en_enhan, ndof_per_nd_enhan);
    n_en_enhan_x_ndof_per_nd_enhan = n_en_enhan*ndof_per_nd_enhan;
    del_ujump_vec.Dimension (n_en_enhan_x_ndof_per_nd_enhan);
    ujump_vec.Dimension (n_en_enhan_x_ndof_per_nd_enhan);
    ujump_dot_vec.Dimension (n_en_enhan_x_ndof_per_nd_enhan);
    ujump_dotdot_vec.Dimension (n_en_enhan_x_ndof_per_nd_enhan);
    //ElementSupport().RegisterCoordinates(fInitCoords_enhan);

    fEnhan->RegisterLocal(ujump);
    fEnhan->RegisterLocal(ujump_n);

    if (fIntegrator->Order() == 1)
    {
		fEnhan->RegisterLocal(ujump_dot);
		fEnhan->RegisterLocal(ujump_dot_n);
    }

    if (fIntegrator->Order() == 2)
    {
		fEnhan->RegisterLocal(ujump_dot);
		fEnhan->RegisterLocal(ujump_dot_n);
		fEnhan->RegisterLocal(ujump_dotdot);
		fEnhan->RegisterLocal(ujump_dotdot_n);
    }

    /* allocate state variable storage */
    // state variables are calculated at IPs for displacement field
    int num_ip = fNumIP_displ;
    fdState_new.Dimension(n_el, num_ip*knum_d_state);
    fdState.Dimension(n_el, num_ip*knum_d_state);
    fiState_new.Dimension(n_el, num_ip*knum_i_state);
    fiState.Dimension(n_el, num_ip*knum_i_state);
	
    /* initialize equations */
    fEqnos_displ.Alias(fEqnos_displ);
    fEqnos_enhan.Dimension(fConnectivities_enhan.Length());

    /* initialize state variables */
    fdState = 0;
    fdState_new = 0;
    fiState = 0;
    fiState_new = 0;

    /* initialize element cards */
    fElementCards_displ.Alias(fElementCards);
    fElementCards_enhan.Dimension(fElementCards.Length());
	
    /* set cards to data in array - NOT NEEDED IF YOU'RE NOT
     * GOING TO USE THE ElementCardT ARRAY? */
    for (int i= 0; i < fElementCards.Length(); i++)
	fElementCards[i].Set(fiState.MinorDim(), fiState(i), fdState.MinorDim(), fdState(i));

    fKdd.Dimension 		( n_en_displ_x_n_sd, n_en_displ_x_n_sd );
    fKdq.Dimension 		( n_en_displ_x_n_sd, n_en_enhan_x_ndof_per_nd_enhan );
    fKqd.Dimension 		( n_en_enhan_x_ndof_per_nd_enhan, n_en_displ_x_n_sd );
    fKqq.Dimension 		( n_en_enhan_x_ndof_per_nd_enhan, n_en_enhan_x_ndof_per_nd_enhan );

    fFd_int.Dimension 	( n_en_displ_x_n_sd );
    fFd_ext.Dimension 	( n_en_displ_x_n_sd );
    fFq_int.Dimension 	( n_en_enhan_x_ndof_per_nd_enhan );
    fFq_ext.Dimension 	( n_en_enhan_x_ndof_per_nd_enhan );
	
    /* workspace matricies */
    fShapeDispl.Dimension (n_sd, n_en_displ_x_n_sd);
    fShapeEnhan.Dimension (ndof_per_nd_enhan, n_en_enhan_x_ndof_per_nd_enhan);
    
    //n_sd_x_n_sd = n_sd*n_sd;
    fShapeDisplGrad.Dimension (n_sd, n_en_displ);    
    fShapeEnhanGrad.Dimension (n_sd, n_en_enhan);

    fB_matrix.Dimension (3,n_en_displ_x_n_sd);
    fD_matrix.Dimension (3,3);
    fK_dd_BTDB_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fTemp_matrix_ndof_se_x_ndof_se.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fTemp_vector_ndof_se.Dimension (n_en_displ_x_n_sd);

    fStrain_vector_current_IP.Dimension (knumstrain);
    fStress_vector_current_IP.Dimension (knumstress);
    
    fStrain_IPs.Dimension (fNumIP_displ,knumstrain);
    fStress_IPs.Dimension (fNumIP_displ,knumstress);
    fState_variables_IPs.Dimension (fNumIP_displ,knum_d_state);
    
    fStrain_Elements_IPs.Dimension (NumElements(),fNumIP_displ*knumstrain);
    fStress_Elements_IPs.Dimension (NumElements(),fNumIP_displ*knumstress);
    fState_variables_Elements_IPs.Dimension (NumElements(),fNumIP_displ*knum_d_state);
    
    /* streams */
    ofstreamT& out = ElementSupport().Output();

    /* storage for integration point strain, stress, and ISVs*/
    fIPVariable.Dimension (n_el, fNumIP_displ*(knumstrain+knumstress+knum_d_state));
    fIPVariable = 0.0;

    /* allocate storage for nodal forces */
    //fForces_at_Node.Dimension ( n_sd );
	
    /* extract natural boundary conditions */
    TakeNaturalBC(list);
	
    /* setup output file and format */
    outputPrecision = 10;
    outputFileWidth = outputPrecision + 8;
//    ss_xfem2D_out.open("ss_xfem2D.info");
}


/* information about subordinate parameter lists */
void SSXfem2DT::DefineSubs(SubListT& sub_list) const
{
    /* inherited */
    ElementBaseT::DefineSubs(sub_list);

    /* element blocks */
    sub_list.AddSub("xfem_SS_2D_element_block");
	
    /* tractions */
    sub_list.AddSub("xfem_SS_2D_natural_bc", ParameterListT::Any);
}



/* return the description of the given inline subordinate parameter list */
void SSXfem2DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
				       SubListT& sub_lists) const
{
    ElementBaseT::DefineInlineSub(name, order, sub_lists);
}



/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSXfem2DT::NewSub(const StringT& name) const
{
    /* create non-const this */
    SSXfem2DT* non_const_this = const_cast<SSXfem2DT*>(this);

    if (name == "xfem_SS_2D_natural_bc") /* traction bc */
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
    else if (name == "xfem_SS_2D_element_block")
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
void SSXfem2DT::SetTractionBC(void)
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
void SSXfem2DT::TakeNaturalBC(const ParameterListT& list)
{
    const char caller[] = "SSXfem2DT::TakeTractionBC";

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
void SSXfem2DT::ApplyTractionBC(void)
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

void SSXfem2DT::Form_solid_shape_functions(const double* &shapes_displ_X)
{
    fShapeDispl = 0.0;
    for (int i=0; i<n_en_displ; i++)
    {
    	//for 2D
		fShapeDispl(0,i*n_sd) = shapes_displ_X[i];
		fShapeDispl(1,1+i*n_sd) = shapes_displ_X[i];
    }
}


void SSXfem2DT::Form_enhan_shape_functions(const double* &shapes_enhan_X)
{
    fShapeEnhan = 0.0;
	for (int i=0; i<n_en_enhan; i++)
    {
    	//for 2D
		fShapeEnhan(0,i*n_sd) = shapes_enhan_X[i];
		fShapeEnhan(1,1+i*n_sd) = shapes_enhan_X[i];
    }	
}





void SSXfem2DT::Form_D_matrix(void)
{
    fD_matrix = 0.0;
    fD_matrix(0,0) = 2*fMaterial_Params[kMu]+ fMaterial_Params[kLambda];
    fD_matrix(1,1) = 2*fMaterial_Params[kMu]+ fMaterial_Params[kLambda];
    fD_matrix(2,2) = fMaterial_Params[kMu];
    fD_matrix(0,1) = fMaterial_Params[kLambda];
    fD_matrix(1,0) = fMaterial_Params[kLambda];
}

void SSXfem2DT::Form_B_matrix(void)
{
    fB_matrix = 0.0;
    for(int i=0; i<n_en_displ; i++)
    {
		fB_matrix(0,i*n_sd)=fShapeDisplGrad(0,i);
		fB_matrix(1,i*n_sd+1)=fShapeDisplGrad(1,i);
		fB_matrix(2,i*n_sd)=fShapeDisplGrad(1,i);
		fB_matrix(2,i*n_sd+1)=fShapeDisplGrad(0,i);
    }
}


