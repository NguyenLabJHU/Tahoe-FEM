/* $Id: BridgingScaleManagerT.cpp,v 1.8 2005-04-09 18:27:33 d-farrell2 Exp $ */
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

/* constructor */
BridgingScaleManagerT::BridgingScaleManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
	const ArrayT<StringT>& argv, TaskT task):
	MultiManagerT(input_file, output, comm, argv, task),
	fFine_THK(NULL)
{
	SetName("tahoe_bridging_scale");
	
	/* split communicator for the coarse scale into N single processor communicators */
	fCoarseComm = new CommunicatorT(fComm, fComm.Rank());
}

BridgingScaleManagerT::~BridgingScaleManagerT(void)
{
	/* clean up communicators */
	if (fCoarseComm != &fComm) delete fCoarseComm;
	if (fFineComm != &fComm) delete fFineComm;
}
	
/* (re-)set system to initial conditions */
ExceptionT::CodeT BridgingScaleManagerT::InitialCondition(void)
{
	/* shouldn't enter here. BridgingScaleManagerT::Solve overrides all */
	ExceptionT::GeneralFail("BridgingScaleManagerT::InitialCondition", "should not be called");
	return ExceptionT::kNoError;
}

/* (re-)set the equation number for the given group */
void BridgingScaleManagerT::Solve(void)
{
	const char caller[] = "FEExecutionManagerT::RunDynamicBridging";

	/* configure ghost nodes */
	NodeManagerT& fine_node_manager = *(fFine_THK->NodeManager());
	int nsd = fine_node_manager.NumSD();
	int nnd = fine_node_manager.NumNodes(); // DEF added
	int group = 0;
	int order1 = 0;	// For InterpolateField, 3 calls to obtain displacement/velocity/acceleration
	int order2 = 1;
	int order3 = 2;
	double dissipation = 0.0;
	
	/* internal force vector - communicated within atomistic side */
	dArray2DT rhs_2D_true(nnd,nsd);
	rhs_2D_true = 0.0e0;
	dArray2DT& rhs_2D = rhs_2D_true; 
	dArray2DT fubig(nnd,nsd);

	CommManagerT* fine_comm_manager = fFine_THK->CommManager();
	if (!fine_comm_manager) ExceptionT::GeneralFail(caller, "could not resolve fine scale comm manager");
	int fubig_ID = fine_comm_manager->Init_AllGather(fubig);
	
	/* other work space */
	dArray2DT field_at_ghosts, totalu, fu, projectedu, boundghostdisp, boundghostvel, boundghostacc;
	dArray2DT thkforce, totaldisp, elecdens, embforce, wavedisp;
	dSPMatrixT ntf;
	iArrayT activefenodes, boundaryghostatoms;
	const StringT& bridging_field = fFineField->FieldName();

	/* figure out boundary atoms for use with THK boundary conditions, 
	   ghost atoms for usage with MD force calculations */
	if (nsd == 2)
		boundaryghostatoms = fFine_THK->InterpolationNodes2D();
	else if (nsd == 3)
		boundaryghostatoms = fFine_THK->InterpolationNodes3D();
	else
		ExceptionT::GeneralFail(caller, "%d dimensions not implemented", nsd);
		
	int numgatoms = (fFine_THK->GhostNodes()).Length();	// total number of ghost atoms
	int numbatoms = boundaryghostatoms.Length() - numgatoms;	// total number of boundary atoms
	dArray2DT gadisp(numgatoms,nsd), gavel(numgatoms,nsd), gaacc(numgatoms,nsd);
	dArray2DT badisp(numbatoms,nsd), bavel(numbatoms,nsd), baacc(numbatoms,nsd);
	iArrayT allatoms(boundaryghostatoms.Length()), gatoms(numgatoms), batoms(numbatoms), boundatoms(numbatoms);
	allatoms.SetValueToPosition();
	batoms.CopyPart(0, allatoms, numgatoms, numbatoms);
	gatoms.CopyPart(0, allatoms, 0, numgatoms);      
	boundatoms.CopyPart(0, boundaryghostatoms, numgatoms, numbatoms);
	//elecdens.Dimension(gatoms.Length(), 1);
	//embforce.Dimension(gatoms.Length(), 1);
	TimeManagerT* atom_time = NULL;		// Declare the time managers, which are defined later.
	TimeManagerT* continuum_time = NULL;
	dArrayT fx, tempx;				// Declare some more arrays, TBD later
	dArray2DT ntfproduct;				
	
	if (fignore == false)	// if we are to use the continuum
	{
		/* reset interpolation data set in MultiManagerT::TakeParameterList */
		fCoarse->InitInterpolation(bridging_field, boundaryghostatoms, fine_node_manager.InitialCoordinates());
		
		/* reset interpolation data set in MultiManagerT::TakeParameterList with makeinactive = false */
		bool makeinactive = false;	
		fCoarse->InitProjection(bridging_field, *(fFine_THK->CommManager()), fFine_THK->NonGhostNodes(), 
			fine_node_manager, makeinactive);
		
		if (nsd == 3)
		{
			/* set pointers to embedding force/electron density in FEManagerT_bridging atoms */
			//fFine_THK->SetExternalElecDensity(elecdens, fFine_THK->GhostNodes());
			//fFine_THK->SetExternalEmbedForce(embforce, fFine_THK->GhostNodes());
			//fCoarse->ElecDensity(gatoms.Length(), elecdens, embforce);
		}	
		
		/* time managers */
		atom_time = fFine_THK->TimeManager();
		continuum_time = fCoarse->TimeManager();
		
		/* compute initial conditions */
		fFine->SetComputeInitialCondition(true);
		fCoarse->SetComputeInitialCondition(true);
		
		/* set to initial condition */
		fFine_THK->InitialCondition();
		
		/* calculate fine scale part of MD displacement and total displacement u */
		fCoarse->InitialProject(bridging_field, fine_node_manager, projectedu, order1);
	}
	else if	(fignore == true)	// if we are to ignore the continuum
	{
		if (nsd == 3)
		{
			/* set pointers to embedding force/electron density in FEManagerT_bridging atoms */
			//fFine_THK->SetExternalElecDensity(elecdens, fFine_THK->GhostNodes());
			//fFine_THK->SetExternalEmbedForce(embforce, fFine_THK->GhostNodes());
			//fCoarse->ElecDensity(gatoms.Length(), elecdens, embforce);
		}
		
		/* time managers */
		atom_time = fFine_THK->TimeManager();
		
		/* compute initial conditions */
		fFine->SetComputeInitialCondition(true);
		
		/* set to initial condition */
		fFine_THK->InitialCondition();
#pragma message("not exactly sure what to do here")		
//		/* calculate fine scale part of MD displacement and total displacement u */
//		fCoarse->InitialProject(bridging_field, fine_node_manager, projectedu, order1);
	}
	
	/* solve for initial FEM force f(u) as function of fine scale + FEM */
	/* use projected totalu instead of totalu for initial FEM displacements */
	nMatrixT<int>& promap = fFine_THK->PropertiesMap(0);   // element group for particles = 0
	const int promap_dim = promap.Rows(); // assumes square property matrix
	
	// now defne the mappings
	nMatrixT<int> ghostonmap(promap_dim),ghostoffmap(promap_dim);

	if (fignore == false)	// if we are to use the continuum
	{
		ghostoffmap = 0;
		
		ghostoffmap = fFine_THK->GetGhostMap();
		
		ghostonmap = promap; // copy original properties map
		promap = ghostoffmap;  // turn ghost atoms off for f(u) calculation
		
		fubig = TotalForce(bridging_field, projectedu, *fFine_THK, rhs_2D) ;
		
		fine_comm_manager->AllGather(fubig_ID, fubig);	
		promap = ghostonmap;   // turn ghost atoms back on for MD force calculations
		
		/* calculate global interpolation matrix ntf */
		fCoarse->Ntf(ntf, fFine_THK->NonGhostNodes(), activefenodes);
		
		/* compute FEM RHS force as product of ntf and fu */
		fx.Dimension(ntf.Rows());
		tempx.Dimension(ntf.Cols());
		ntfproduct.Dimension(ntf.Rows(), nsd);
		fu.Dimension((fFine_THK->NonGhostNodes()).Length(), nsd);
		fu.RowCollect(fFine_THK->NonGhostNodes(), fubig); 
		
		for (int i = 0; i < nsd; i++)
		{
			fu.ColumnCopy(i, tempx);
			ntf.Multx(tempx, fx);
			ntfproduct.SetColumn(i, fx);
		}
	
		/* Add FEM RHS force to RHS using SetExternalForce */
		fCoarse->SetExternalForce(bridging_field, ntfproduct, activefenodes);
		
		/* now d0, v0 and a0 are known after InitialCondition */
		fCoarse->InitialCondition();
			
		if (nsd == 3)
		{
			/* Calculate EAM electron density/embedding terms for ghost atoms using continuum information */
			//fCoarse->ElecDensity(gatoms.Length(), elecdens, embforce);
		}
			
		/* interpolate FEM values to MD ghost nodes which will act as MD boundary conditions */
		fCoarse->InterpolateField(bridging_field, order1, boundghostdisp);
		fCoarse->InterpolateField(bridging_field, order2, boundghostvel);
		fCoarse->InterpolateField(bridging_field, order3, boundghostacc);
	}
	else if	(fignore == true)	// if we are to ignore the continuum
	{
		// Dimension the arrays to hold boundary atom information
		boundghostdisp.Dimension(boundaryghostatoms.Length(),nsd);
		boundghostvel.Dimension(boundaryghostatoms.Length(),nsd);
		boundghostacc.Dimension(boundaryghostatoms.Length(),nsd);
		
		// set the boundary information to zero displacement, vel, accel.
		boundghostdisp = 0.0;
		boundghostvel = 0.0;
		boundghostacc = 0.0;
	}
	
		
	/* sort boundary + ghost atom info into separate arrays */
	gadisp.RowCollect(gatoms, boundghostdisp);
	gavel.RowCollect(gatoms, boundghostvel);
	gaacc.RowCollect(gatoms, boundghostacc);
	badisp.RowCollect(batoms, boundghostdisp);
	bavel.RowCollect(batoms, boundghostvel);
	baacc.RowCollect(batoms, boundghostacc);
			
	if (nsd == 2)
	{
		/* store initial MD boundary displacement histories */
		thkforce = fFine_THK->THKForce2D(bridging_field, badisp);
		fFine_THK->SetExternalForce(bridging_field, thkforce, boundatoms);  // sets pointer to thkforce 
	}
	else
	{
		/* thkdisp = fine scale part of ghost atom displacement */
		thkforce = fFine_THK->THKForce3D(bridging_field, badisp);
		fFine_THK->SetExternalForce(bridging_field, thkforce, boundatoms);
	}
#pragma message("roll up redundancy here when working, perhaps put into separate method?")	
	if (fignore == false)	// if we are to use the continuum
	{
		/* figure out timestep ratio between fem and md simulations */
		int nfesteps = continuum_time->NumberOfSteps();
		double mddt = atom_time->TimeStep();
		double fedt = continuum_time->TimeStep();
		double d_ratio = fedt/mddt;		
		int ratio = int((2.0*d_ratio + 1.0)/2.0);
		
		/* running status flag */
		ExceptionT::CodeT error = ExceptionT::kNoError;	

		for (int i = 0; i < nfesteps; i++)	
		{
			for (int j = 0; j < ratio; j++)	// MD update first
			{
				atom_time->Step();	
									
				/* initialize step */
				if (1 || error == ExceptionT::kNoError) error = fFine_THK->InitStep();
					
				/* update FEM solution interpolated at boundary atoms and ghost atoms assuming 
				constant acceleration - because of constant acceleration assumption, predictor and 
				corrector are combined into one function */
				fFine_THK->BAPredictAndCorrect(mddt, badisp, bavel, baacc);
				fFine_THK->BAPredictAndCorrect(mddt, gadisp, gavel, gaacc);

				if (nsd == 2)
				{
					/* Write interpolated FEM values at MD ghost nodes into MD field - displacements only */
					fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), order1, gadisp);	
					
					/* calculate THK force on boundary atoms, update displacement histories */
					thkforce = fFine_THK->THKForce2D(bridging_field, badisp);  // SetExternalForce set via pointer
				}
				else
				{
					/* Write interpolated FEM values at MD ghost nodes into MD field - displacements only */
					fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), order1, gadisp);
						
					/* calculate thkforces */
					thkforce = fFine_THK->THKForce3D(bridging_field, badisp);
				}
					
				/* solve MD equations of motion */
				if (1 || error == ExceptionT::kNoError) {
						fFine_THK->ResetCumulativeUpdate(group);
						error = fFine_THK->SolveStep();
				}

				/* close  md step */
				if (1 || error == ExceptionT::kNoError) error = fFine_THK->CloseStep();    
			}

			continuum_time->Step();

			/* initialize step */
			if (1 || error == ExceptionT::kNoError) error = fCoarse->InitStep();
                
			/* calculate total displacement u = FE + fine scale here using updated FEM displacement */
			fCoarse->BridgingFields(bridging_field, fine_node_manager, *(fCoarse->NodeManager()), totalu);
			
			/* calculate FE internal force as function of total displacement u here */
			promap = ghostoffmap;  // turn off ghost atoms for f(u) calculations
			
			fubig = TotalForce(bridging_field, totalu, *fFine_THK, rhs_2D);
			
			fine_comm_manager->AllGather(fubig_ID, fubig);	
			promap = ghostonmap;   // turn on ghost atoms for MD force calculations
			fu.RowCollect(fFine_THK->NonGhostNodes(), fubig); 

			for (int i = 0; i < nsd; i++)
			{
				fu.ColumnCopy(i, tempx);
				ntf.Multx(tempx, fx);
				ntfproduct.SetColumn(i, fx);
			}
				
			/* solve FE equation of motion using internal force just calculated */
			if (1 || error == ExceptionT::kNoError)
			{
				fCoarse->ResetCumulativeUpdate(group);
				error = fCoarse->SolveStep();
			}

			/* Interpolate FEM values to MD ghost nodes which will act as MD boundary conditions */
			fCoarse->InterpolateField(bridging_field, order1, boundghostdisp);
			fCoarse->InterpolateField(bridging_field, order2, boundghostvel);
			fCoarse->InterpolateField(bridging_field, order3, boundghostacc);
				
			/* sort boundary + ghost atom info into separate arrays */
			gadisp.RowCollect(gatoms, boundghostdisp);
			gavel.RowCollect(gatoms, boundghostvel);
			gaacc.RowCollect(gatoms, boundghostacc);
			badisp.RowCollect(batoms, boundghostdisp);
			bavel.RowCollect(batoms, boundghostvel);
			baacc.RowCollect(batoms, boundghostacc);
				
			if (nsd == 3)
			{
				/* Calculate EAM electron density/embedding terms for ghost atoms using updated continuum information */
				//fCoarse->ElecDensity(gatoms.Length(), elecdens, embforce);
			}
				
			/* close fe step */
			if (1 || error == ExceptionT::kNoError) error = fCoarse->CloseStep();                   
		}
	}
	else if	(fignore == true)	// if we are to ignore the continuum
	{
		/* figure out timestep ratio between fem and md simulations */
		int nfesteps = atom_time->NumberOfSteps();
		double mddt = atom_time->TimeStep();
		double d_ratio = 1.0;		
		int ratio = 1;		
		
		/* running status flag */
		ExceptionT::CodeT error = ExceptionT::kNoError;	

		for (int i = 0; i < nfesteps; i++)	
		{
			atom_time->Step();	
								
			/* initialize step */
			if (1 || error == ExceptionT::kNoError) error = fFine_THK->InitStep();
				
			/* update MD solution at boundary atoms and ghost atoms assuming 
			constant acceleration - because of constant acceleration assumption, predictor and 
			corrector are combined into one function */
			fFine_THK->BAPredictAndCorrect(mddt, badisp, bavel, baacc);
			fFine_THK->BAPredictAndCorrect(mddt, gadisp, gavel, gaacc);

			if (nsd == 2)
			{
				/* Write interpolated FEM values at MD ghost nodes into MD field - displacements only */
				fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), order1, gadisp);	
				
				/* calculate THK force on boundary atoms, update displacement histories */
				thkforce = fFine_THK->THKForce2D(bridging_field, badisp);  // SetExternalForce set via pointer
			}
			else
			{
				/* Write interpolated FEM values at MD ghost nodes into MD field - displacements only */
				fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), order1, gadisp);
					
				/* calculate thkforces */
				thkforce = fFine_THK->THKForce3D(bridging_field, badisp);
			}
				
			/* solve MD equations of motion */
			if (1 || error == ExceptionT::kNoError) {
					fFine_THK->ResetCumulativeUpdate(group);
					error = fFine_THK->SolveStep();
			}

			/* close  md step */
			if (1 || error == ExceptionT::kNoError) error = fFine_THK->CloseStep(); 
		}
	} // end ignore continuum if loop

	/* check for error */
//	if (0) ExceptionT::GeneralFail(caller, "hit error %d", error);
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
 * Protected
 **********************************************************************/
// calculate total finescale force for the given displacement u, DEF added
const dArray2DT& BridgingScaleManagerT::TotalForce(const StringT& field_name, const dArray2DT& field_values, 
	FEManagerT_bridging& bridging, dArray2DT& rhs_2D) const
{
// DEBUG
//cout << "BridgingScaleManagerT::InternalForce : Start" << endl;

	NodeManagerT& fine_node_manager = *(fFine_THK->NodeManager());
	int nsd = fine_node_manager.NumSD(); // get the number of spatial dims

	/* first obtain the MD displacement field */
	FieldT* field = bridging.NodeManager()->Field(field_name);
	if (!field) 
		ExceptionT::GeneralFail("BridgingScaleManagerT::TotalForce", "field \"%s\" not found",
			field_name.Pointer());
	
	dArray2DT disp_0 = (*field)[0];	// temporarily store current displacements
	int group = 0; // assume particle group number = 0
	
	/* obtain atom node list - can calculate once and store... */
	int nnd = field_values.MajorDim();
	iArrayT nodes(nnd);
	nodes.SetValueToPosition();
	
	// get the ghost nodes
	iArrayT& ghostnodes = (iArrayT&) bridging.GhostNodes();
		
	/* now write total bridging scale displacement u into field */
	int order = 0;	// write displacement only
	bridging.SetFieldValues(field_name, nodes, order, field_values);
	
	/* compute RHS - ParticlePairT fForce calculated by this call */
	bridging.FormRHS(group);
	
	// lets try to get the external force vector information from FieldT
	// fext* will be NULL if no external force (i.e. length = 0)
	ArrayT<FieldT*> finefields;
	fine_node_manager.CollectFields(group, finefields);
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
	dArray2DT& internalforce = (dArray2DT&) bridging.InternalForce(group);
	
	// get the equation numbers from the MD displacement field (1 is first equation # (corresponds to first entry in RHS)
	iArray2DT& eq_nos = field->Equations();
	
	//get the internal and external (interatomic & BC) force into the projection
	dArrayT& rhs = (dArrayT&) bridging.RHS(group);
	dArrayT rhs_temp = rhs;
	
	// insert the internal force into the total force
	for (int i = 0; i < nnd; i++)
	{
		for (int j = 0; j < nsd; j++)
		{
			rhs_2D(i,j) = internalforce(i,j);
		}
	}
	
	// now add in the external force contribution -> roll this into above once it works (maybe make it faster..)
	// since eqn numbering in 2D force array is row order (xcomp # < ycomp # < zcomp#)
	
	// first translate Fext array into 2D with dimensions (nnd,nsd) -> mapping in general???
	dArray2DT external_force_vals(nnd,nsd);
	external_force_vals = 0.0e0;
	dArray2DT& external_force_ref = external_force_vals; 
	int extmarker = -1;
	
	
	if (fextvals.Length() != 0)
	{	
		for (int i = 0; i < nnd; i++)
		{
			for (int j = 0; j < nsd; j++)
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

/* calculate internal force for the given displacement u */
const dArray2DT& BridgingScaleManagerT::InternalForce(const StringT& field_name, const dArray2DT& field_values, 
	FEManagerT_bridging& bridging) const
{
	/* first obtain the MD displacement field */
	FieldT* field = bridging.NodeManager()->Field(field_name);
	if (!field) 
		ExceptionT::GeneralFail("BridgingScaleManagerT::InternalForce", "field \"%s\" not found",
			field_name.Pointer());
	
	dArray2DT disp_0 = (*field)[0];	// temporarily store current displacements
	int group = 0; // assume particle group number = 0
	
	/* obtain atom node list - can calculate once and store... */
	int nnd = field_values.MajorDim();
	iArrayT nodes(nnd);
	nodes.SetValueToPosition();
		
	/* now write total bridging scale displacement u into field */
	int order = 0;	// write displacement only
	bridging.SetFieldValues(field_name, nodes, order, field_values);
		
	/* compute RHS - ParticlePairT fForce calculated by this call */
	bridging.FormRHS(group);
	
	/* write actual MD displacements back into field */
	bridging.SetFieldValues(field_name, nodes, order, disp_0);

	/* get the internal force contribution associated with the last call to FormRHS */
	return bridging.InternalForce(group);
}

#endif /* BRIDGING_ELEMENT && BRIDGING_ELEMENT_DEV */
