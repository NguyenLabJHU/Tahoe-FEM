/* $Id: BridgingScaleManagerT.cpp,v 1.12 2005-07-21 21:57:55 d-farrell2 Exp $ */

#include "BridgingScaleManagerT.h"

#if defined(BRIDGING_ELEMENT) && defined(BRIDGING_ELEMENT_DEV)

#include "ifstreamT.h"
#include "SolverT.h"
#include "FEManagerT_THK.h"
#include "NodeManagerT.h"
#include "TimeManagerT.h"
#include "ParticlePairT.h"
#include "FieldT.h"
#include "dSPMatrixT.h"
#include "CommunicatorT.h"
#include "CommManagerT.h"

using namespace Tahoe;

// constructor
BridgingScaleManagerT::BridgingScaleManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
	const ArrayT<StringT>& argv, TaskT task):
	MultiManagerT(input_file, output, comm, argv, task),
	fFine_THK(NULL)
{
	SetName("tahoe_bridging_scale");
	
	// split communicator for the coarse scale into N single processor communicators
	fCoarseComm = new CommunicatorT(fComm, fComm.Rank());
}

BridgingScaleManagerT::~BridgingScaleManagerT(void)
{
	// clean up communicators
	if (fCoarseComm != &fComm) delete fCoarseComm;
	if (fFineComm != &fComm) delete fFineComm;
}
	
// (re-)set system to initial conditions
ExceptionT::CodeT BridgingScaleManagerT::InitialCondition(void)
{
	// shouldn't enter here. BridgingScaleManagerT::Solve overrides all
	ExceptionT::GeneralFail("BridgingScaleManagerT::InitialCondition", "should not be called");
	return ExceptionT::kNoError;
}

// solve all the time sequences
void BridgingScaleManagerT::Solve(void)
{
	const char caller[] = "BridgingScaleManagerT::Solve";
	
	// Set up the required parameters, etc, and find initial condition
	Initialize();
	
	if (fignore == false)	// BSM simulation solver
	{
		SolveBSM();
	}
	else if	(fignore == true)	// MD/THK simulation solver
	{
		SolveMDTHK();
	}
}

/* accept parameter list */
void BridgingScaleManagerT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "BridgingScaleManagerT::TakeParameterList";

	/* inherited */
	MultiManagerT::TakeParameterList(list);

	/* resolve atomistic solver */
	fFine_THK = TB_DYNAMIC_CAST(FEManagerT_THK*, fFine);
	if (!fFine_THK || fFine_THK->Name() != "tahoe_THK")
		ExceptionT::GeneralFail(caller, "expecting \"tahoe_THK\" fine scale not \"%s\"",
			fFine_THK->Name().Pointer());
}

/**********************************************************************
 * Private
 **********************************************************************/

// Determine basic solution information
void BridgingScaleManagerT::Initialize(void)
{
	const char caller[] = "BridgingScaleManagerT::Initialize";
	
	// get size and dimension information from node manager
	fNSD = (*(fFine_THK->NodeManager())).NumSD();
	fNND = (*(fFine_THK->NodeManager())).NumNodes();
	
#pragma message ("FIX ME!! - finish implementation")	
	// Set up THK, etc
	if (fignore == true) // MDTHK
	{
		if (fNSD == 2 || fNSD == 3)
			fFine_THK->InitializeMDTHK();
		else
			ExceptionT::GeneralFail(caller, "%d dimensions not supported", fNSD);
	}
	else	// BSM by default
	{
		if (fNSD == 2)
			fFine_THK->Initialize2D();
		else if (fNSD == 3)
			fFine_THK->Initialize3D();
		else
			ExceptionT::GeneralFail(caller, "%d dimensions not supported", fNSD);
	}
	
	/* figure out boundary atoms for use with THK boundary conditions, 
	   ghost atoms for usage with MD force calculations */
	if (fNSD == 2)
		fBoundaryghostatoms = fFine_THK->InterpolationNodes2D();
	else if (fNSD == 3)
		fBoundaryghostatoms = fFine_THK->InterpolationNodes3D();
	else
		ExceptionT::GeneralFail(caller, "%d dimensions not implemented", fNSD);
	
	fNumgatoms = (fFine_THK->GhostNodes()).Length();	// total number of ghost atoms
	fNumbatoms = fBoundaryghostatoms.Length() - fNumgatoms;	// total number of boundary atoms
	
	// Dimension ghost atom kinematic info
	fGadisp.Dimension(fNumgatoms,fNSD);
	fGavel.Dimension(fNumgatoms,fNSD);
	fGaacc.Dimension(fNumgatoms,fNSD);
	
	// Dimension boundary atom kinematic info 
	fBadisp.Dimension(fNumbatoms,fNSD);
	fBavel.Dimension(fNumbatoms,fNSD);
	fBaacc.Dimension(fNumbatoms,fNSD);
	
	// Dimension some arrays of atom set numbers
	fAllatoms.Dimension(fBoundaryghostatoms.Length());
	fGatoms.Dimension(fNumgatoms);
	fBatoms.Dimension(fNumbatoms);
	fBoundatoms.Dimension(fNumbatoms);
	
	// Populate the atom set numbers
	fAllatoms.SetValueToPosition();
	fBatoms.CopyPart(0, fAllatoms, fNumgatoms, fNumbatoms);
	fGatoms.CopyPart(0, fAllatoms, 0, fNumgatoms);      
	fBoundatoms.CopyPart(0, fBoundaryghostatoms, fNumgatoms, fNumbatoms);
	
	
	// now initialize BSM or MD/THK stuff as needed
	if (fignore == true) // don't use continuum, i.e. MD/THK
	{
		InitMDTHK();
	}
	else // use continuum, i.e. BSM by default
	{
		InitBSM();
	}
}

// Determine basic BSM solution information, initial conditions
void BridgingScaleManagerT::InitBSM(void)
{
	const char caller[] = "BridgingScaleManagerT::InitBSM";
	
	const StringT& bridging_field = fFineField->FieldName();
	
	// internal force vector - communicated within atomistic side
	fRHS_2D_true.Dimension(fNND,fNSD);
	fRHS_2D_true = 0.0e0;
	
	fFubig.Dimension(fNND,fNSD);
	
	fFine_comm_manager = fFine_THK->CommManager();
	if (!fFine_comm_manager) ExceptionT::GeneralFail(caller, "could not resolve fine scale comm manager");
	fFubig_ID = fFine_comm_manager->Init_AllGather(fFubig);
	
	// reset interpolation data set in MultiManagerT::TakeParameterList
	fCoarse->InitInterpolation(bridging_field, fBoundaryghostatoms, (*(fFine_THK->NodeManager())).InitialCoordinates());
	
	// reset interpolation data set in MultiManagerT::TakeParameterList with makeinactive = false
	bool makeinactive = false;	
	fCoarse->InitProjection(bridging_field, *(fFine_THK->CommManager()), fFine_THK->NonGhostNodes(), 
		*(fFine_THK->NodeManager()), makeinactive, true);
	
	// define time managers
	fFine_time_manager = fFine_THK->TimeManager();
	fCoarse_time_manager = fCoarse->TimeManager();
	
	// compute initial conditions
	fFine->SetComputeInitialCondition(true);
	fCoarse->SetComputeInitialCondition(true);
	
	// set to initial condition
	fFine_THK->InitialCondition();
	
	// calculate fine scale part of MD displacement and total displacement u
	fCoarse->InitialProject(bridging_field, *(fFine_THK->NodeManager()), fProjectedu, 0, 0);
	
	// solve for initial FEM force f(u) as function of fine scale + FEM
	// use projected totalu instead of totalu for initial FEM displacements
	const int promap_dim = (fFine_THK->PropertiesMap(0)).Rows();	// assumes square property matrix,  element group for particles = 0
	
	// now dimension the ghost on/off property mappings
	fGhostonmap.Dimension(promap_dim);
	fGhostoffmap.Dimension(promap_dim);
	
	// now define the ghost on/off property mappings
	fGhostoffmap = 0;
	fGhostoffmap = fFine_THK->GetGhostMap();
	
	fGhostonmap = (fFine_THK->PropertiesMap(0));	// copy original properties map
	(fFine_THK->PropertiesMap(0)) = fGhostoffmap;	// turn ghost atoms off for f(u) calculation
	
	fFubig = TotalForce(bridging_field, fProjectedu, *fFine_THK, fRHS_2D_true);
	
	fFine_comm_manager->AllGather(fFubig_ID, fFubig);	
	(fFine_THK->PropertiesMap(0)) = fGhostonmap;   // turn ghost atoms back on for MD force calculations

	// calculate global interpolation matrix NtF
	fCoarse->Ntf(fNtF, fFine_THK->NonGhostNodes(), fActiveFENodes);
	
	// compute FEM RHS force as product of ntf and fu
	fFx.Dimension(fNtF.Rows());
	fTempx.Dimension(fNtF.Cols());
	fNtfproduct.Dimension(fNtF.Rows(), fNSD);
	fFu.Dimension((fFine_THK->NonGhostNodes()).Length(), fNSD);
	fFu.RowCollect(fFine_THK->NonGhostNodes(), fFubig); 
	
	for (int i = 0; i < fNSD; i++)
	{
		fFu.ColumnCopy(i, fTempx);
		fNtF.Multx(fTempx, fFx);
		fNtfproduct.SetColumn(i, fFx);
	}
	
	// Add FEM RHS force to RHS using SetExternalForce
	fCoarse->SetExternalForce(fFineField->FieldName(), fNtfproduct, fActiveFENodes);
	
	// now d0, v0 and a0 are known after InitialCondition
	fCoarse->InitialCondition();
		
	// interpolate FEM values to MD ghost nodes which will act as MD boundary conditions (0 - disp, 1-vel, 2-acc)
	fCoarse->InterpolateField(bridging_field, 0, fBoundghostdisp);
	fCoarse->InterpolateField(bridging_field, 1, fBoundghostvel);
	fCoarse->InterpolateField(bridging_field, 2, fBoundghostacc);
	
	// sort boundary + ghost atom info into separate arrays
	fGadisp.RowCollect(fGatoms, fBoundghostdisp);
	fGavel.RowCollect(fGatoms, fBoundghostvel);
	fGaacc.RowCollect(fGatoms, fBoundghostacc);
	fBadisp.RowCollect(fBatoms, fBoundghostdisp);
	fBavel.RowCollect(fBatoms, fBoundghostvel);
	fBaacc.RowCollect(fBatoms, fBoundghostacc);
			
	if (fNSD == 2)
	{
		// store initial MD boundary displacement histories
		fTHKforce = fFine_THK->THKForce2D(bridging_field, fBadisp);
		fFine_THK->SetExternalForce(bridging_field, fTHKforce, fBoundatoms);  // sets pointer to thkforce 
	}
	else
	{
		// thkdisp = fine scale part of ghost atom displacement
		fTHKforce = fFine_THK->THKForce3D(bridging_field, fBadisp);
		fFine_THK->SetExternalForce(bridging_field, fTHKforce, fBoundatoms);
	}
}

// Determine basic MD/THK solution information, initial conditions
void BridgingScaleManagerT::InitMDTHK(void)
{
	const StringT& bridging_field = fFineField->FieldName();
	
	// time manager
	fFine_time_manager = fFine_THK->TimeManager();
	
	// compute initial conditions
	fFine->SetComputeInitialCondition(true);
	
	// set to initial condition
	fFine_THK->InitialCondition();

	// Dimension the arrays to hold boundary atom information
	fBoundghostdisp.Dimension(fBoundaryghostatoms.Length(),fNSD);
	fBoundghostvel.Dimension(fBoundaryghostatoms.Length(),fNSD);
	fBoundghostacc.Dimension(fBoundaryghostatoms.Length(),fNSD);
	
	// set the boundary information to zero displacement, vel, accel.
	fBoundghostdisp = 0.0;
	fBoundghostvel = 0.0;
	fBoundghostacc = 0.0;
	
	// sort boundary + ghost atom info into separate arrays
	fGadisp.RowCollect(fGatoms, fBoundghostdisp);
	fGavel.RowCollect(fGatoms, fBoundghostvel);
	fGaacc.RowCollect(fGatoms, fBoundghostacc);
	fBadisp.RowCollect(fBatoms, fBoundghostdisp);
	fBavel.RowCollect(fBatoms, fBoundghostvel);
	fBaacc.RowCollect(fBatoms, fBoundghostacc);
			
	if (fNSD == 2)
	{
		// Determine the ghost atom displacement from displacement formulation
		fGadisp = fFine_THK->THKDisp(bridging_field, fBadisp);		

		// Write interpolated FEM values at MD ghost nodes into MD field - displacements only
		fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), 0, fGadisp); 
	}
	else
	{
		// Determine the ghost atom displacement from displacement formulation
		fGadisp = fFine_THK->THKDisp(bridging_field, fBadisp);		

		// Write interpolated FEM values at MD ghost nodes into MD field - displacements only
		fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), 0, fGadisp);
	}
}

// Determine BSM solution
void BridgingScaleManagerT::SolveBSM(void)
{

	const StringT& bridging_field = fFineField->FieldName();
	
	// figure out timestep ratio between fem and md simulations
	int nfesteps = fCoarse_time_manager->NumberOfSteps();
	double mddt = fFine_time_manager->TimeStep();
	double fedt = fCoarse_time_manager->TimeStep();
	double d_ratio = fedt/mddt;		
	int ratio = int((2.0*d_ratio + 1.0)/2.0);
	
	// running status flag
	ExceptionT::CodeT error = ExceptionT::kNoError;	

	for (int i = 0; i < nfesteps; i++)	
	{
		for (int j = 0; j < ratio; j++)	// MD update first
		{
			fFine_time_manager->Step();	
								
			// initialize step
			if (1 || error == ExceptionT::kNoError) error = fFine_THK->InitStep();
				
			/* update FEM solution interpolated at boundary atoms and ghost atoms assuming 
			constant acceleration - because of constant acceleration assumption, predictor and 
			corrector are combined into one function */
			fFine_THK->BAPredictAndCorrect(mddt, fBadisp, fBavel, fBaacc);
			fFine_THK->BAPredictAndCorrect(mddt, fGadisp, fGavel, fGaacc);

			if (fNSD == 2)
			{
				// Write interpolated FEM values at MD ghost nodes into MD field - displacements only
				fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), 0, fGadisp);	
				
				// calculate THK force on boundary atoms, update displacement histories
				fTHKforce = fFine_THK->THKForce2D(bridging_field, fBadisp);  // SetExternalForce set via pointer
			}
			else
			{
				// Write interpolated FEM values at MD ghost nodes into MD field - displacements only
				fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), 0, fGadisp);
					
				// calculate thkforces
				fTHKforce = fFine_THK->THKForce3D(bridging_field, fBadisp);
			}
				
			// solve MD equations of motion
			if (1 || error == ExceptionT::kNoError) {
					fFine_THK->ResetCumulativeUpdate(0);
					error = fFine_THK->SolveStep();
			}

			// close  md step
			if (1 || error == ExceptionT::kNoError) error = fFine_THK->CloseStep();    
		}

		fCoarse_time_manager->Step();

		// initialize step
		if (1 || error == ExceptionT::kNoError) error = fCoarse->InitStep();
            
		// calculate total displacement u = FE + fine scale here using updated FEM displacement
		fCoarse->BridgingFields(bridging_field, *(fFine_THK->NodeManager()), *(fCoarse->NodeManager()), fTotalu);
		
		// calculate FE internal force as function of total displacement u here
		(fFine_THK->PropertiesMap(0)) = fGhostoffmap;  // turn off ghost atoms for f(u) calculations
		
		fFubig = TotalForce(bridging_field, fTotalu, *fFine_THK, fRHS_2D_true);
		
		fFine_comm_manager->AllGather(fFubig_ID, fFubig);	
		(fFine_THK->PropertiesMap(0)) = fGhostonmap;   // turn on ghost atoms for MD force calculations
		fFu.RowCollect(fFine_THK->NonGhostNodes(), fFubig); 

		for (int i = 0; i < fNSD; i++)
		{
			fFu.ColumnCopy(i, fTempx);
			fNtF.Multx(fTempx, fFx);
			fNtfproduct.SetColumn(i, fFx);
		}
			
		// solve FE equation of motion using internal force just calculated
		if (1 || error == ExceptionT::kNoError)
		{
			fCoarse->ResetCumulativeUpdate(0);
			error = fCoarse->SolveStep();
		}

		// Interpolate FEM values to MD ghost nodes which will act as MD boundary conditions (0 - disp, 1-vel, 2-acc)
		fCoarse->InterpolateField(bridging_field, 0, fBoundghostdisp);
		fCoarse->InterpolateField(bridging_field, 1, fBoundghostvel);
		fCoarse->InterpolateField(bridging_field, 2, fBoundghostacc);
			
		// sort boundary + ghost atom info into separate arrays
		fGadisp.RowCollect(fGatoms, fBoundghostdisp);
		fGavel.RowCollect(fGatoms, fBoundghostvel);
		fGaacc.RowCollect(fGatoms, fBoundghostacc);
		fBadisp.RowCollect(fBatoms, fBoundghostdisp);
		fBavel.RowCollect(fBatoms, fBoundghostvel);
		fBaacc.RowCollect(fBatoms, fBoundghostacc);
			
		// close fe step
		if (1 || error == ExceptionT::kNoError) error = fCoarse->CloseStep();                   
	}
}

// Determine MD/THK solution
void BridgingScaleManagerT::SolveMDTHK(void)
{
	// Note: we use the displacement formulation for MD/THK problems
	
	const StringT& bridging_field = fFineField->FieldName();
	
	// figure out timestep ratio between fem and md simulations
	int nfesteps = fFine_time_manager->NumberOfSteps();
	double mddt = fFine_time_manager->TimeStep();
	double d_ratio = 1.0;		
	int ratio = 1;
	dArray2DT gadisp_coarse = fGadisp;		
	
	// running status flag
	ExceptionT::CodeT error = ExceptionT::kNoError;	

	for (int i = 0; i < nfesteps; i++)	
	{
		fFine_time_manager->Step();	
							
		// initialize step
		if (1 || error == ExceptionT::kNoError) error = fFine_THK->InitStep();

		/* update MD solution at boundary atoms assuming constant acceleration -
		because of constant acceleration assumption, predictor and 
		corrector are combined into one function */
		fFine_THK->BAPredictAndCorrect(mddt, fBadisp, fBavel, fBaacc);
		
		// Determine coarse scale displacement of ghost atoms (zero for now)
		gadisp_coarse = 0.0;
		fGadisp = gadisp_coarse;
		
		// Determine coarse scale displacement of boundary atoms (zero for now)
		fBadisp = 0.0;
		
		// Determine the ghost atom displacement from displacement formulation (fluctuation due to THK)
		// then add to coarse ghost atom displacement
		fGadisp += fFine_THK->THKDisp(bridging_field, fBadisp);		

		// Write displacement at MD ghost nodes into MD field - displacements only
		fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), 0, fGadisp);	
			
		// solve MD equations of motion
		if (1 || error == ExceptionT::kNoError) {
				fFine_THK->ResetCumulativeUpdate(0);
				error = fFine_THK->SolveStep();
		}

		// close  md step
		if (1 || error == ExceptionT::kNoError) error = fFine_THK->CloseStep(); 
	}
}
 
// calculate total finescale force for the given displacement u
const dArray2DT& BridgingScaleManagerT::TotalForce(const StringT& field_name, const dArray2DT& field_values, 
	FEManagerT_bridging& bridging, dArray2DT& rhs_2D) const
{
// DEBUG
//cout << "BridgingScaleManagerT::InternalForce : Start" << endl;

	NodeManagerT& fine_node_manager = *(fFine_THK->NodeManager());
	int fNSD = fine_node_manager.NumSD(); // get the number of spatial dims

	/* first obtain the MD displacement field */
	FieldT* field = bridging.NodeManager()->Field(field_name);
	if (!field) 
		ExceptionT::GeneralFail("BridgingScaleManagerT::TotalForce", "field \"%s\" not found",
			field_name.Pointer());
	
	dArray2DT disp_0 = (*field)[0];	// temporarily store current displacements
	int group = 0; // assume particle group number = 0
	
	/* obtain atom node list - can calculate once and store... */
	int fNND = field_values.MajorDim();
	iArrayT nodes(fNND);
	nodes.SetValueToPosition();
	
	// get the ghost nodes
	iArrayT& ghostnodes = (iArrayT&) bridging.GhostNodes();
		
	/* now write total bridging scale displacement u into field */
	int order = 0;	// write displacement only
	bridging.SetFieldValues(field_name, nodes, order, field_values);
	
	/* compute RHS - ParticlePairT fForce calculated by this call */
	bridging.FormRHS(0);
	
	// lets try to get the external force vector information from FieldT
	// fext* will be NULL if no external force (i.e. length = 0)
	ArrayT<FieldT*> finefields;
	fine_node_manager.CollectFields(0, finefields);
	FieldT* fields_ref;
	for (int ifield = 0; ifield < finefields.Length(); ifield++)
	{
		fields_ref = finefields[ifield]; // not sure what to do with this when multiple fields...
	}
	const dArrayT& fextvals = fields_ref->FieldT::GetfFBCValues();
	const iArrayT& fexteqns = fields_ref->FieldT::GetfFBCEqnos(); 	
	
	/* write actual MD displacements back into field */
	bridging.SetFieldValues(field_name, nodes, order, disp_0);

	// get the internal force contribution associated with the last call to FormRHS
	dArray2DT& internalforce = (dArray2DT&) bridging.InternalForce(0);
	
	// get the equation numbers from the MD displacement field (1 is first equation # (corresponds to first entry in RHS)
	iArray2DT& eq_nos = field->Equations();
	
	//get the internal and external (interatomic & BC) force into the projection
	dArrayT& rhs = (dArrayT&) bridging.RHS(0);
	dArrayT rhs_temp = rhs;
	
	// insert the internal force into the total force
	for (int i = 0; i < fNND; i++)
	{
		for (int j = 0; j < fNSD; j++)
		{
			rhs_2D(i,j) = internalforce(i,j);
		}
	}
	
	// now add in the external force contribution -> roll this into above once it works (maybe make it faster..)
	// since eqn numbering in 2D force array is row order (xcomp # < ycomp # < zcomp#)
	
	// first translate Fext array into 2D with dimensions (fNND,fNSD) -> mapping in general???
	dArray2DT external_force_vals(fNND,fNSD);
	external_force_vals = 0.0e0;
	dArray2DT& external_force_ref = external_force_vals; 
	int extmarker = -1;
	
	
	if (fextvals.Length() != 0)
	{	
		for (int i = 0; i < fNND; i++)
		{
			for (int j = 0; j < fNSD; j++)
			{
				// this searches the external force equations... needs better implementation.
				for (int k = 0; k < fextvals.Length(); k++) // assumes length(fextvals) = length(fexteqns)
				{
					if (fexteqns[k] == eq_nos(i,j)) // if the equation number matches up
					{
						// record k, break out of loop (try to save some time...)
						extmarker = k;
						break;
					}
				}				
				if (extmarker > -1)
				{
					rhs_2D(i,j) += fextvals[extmarker];
				}
				else
				{
					cout << "BridgingScaleManagerT::TotalForce, extmarker = -1" << endl;
				}
			}
		}
	}
	
	return rhs_2D;
}

#endif /* BRIDGING_ELEMENT && BRIDGING_ELEMENT_DEV */
